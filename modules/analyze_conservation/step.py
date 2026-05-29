"""Orchestration for analyze_conservation: per-track conservation analysis plus
the preflight (FASTA inspection) and postflight (summary) hooks."""

import datetime
import json
from collections import Counter
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from rich import box
from rich.panel import Panel
from rich.table import Table
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

import config
from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin, is_interactive_session
from utils.csv_write import write_user_facing_csv
from utils.naming import (
    COLUMN_ALLELES_UNITED,
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_PEPTIDE,
    STAR_MARKER,
    get_step_filename,
)

from .charts import write_conservation_dual_panel_png
from .core import (
    compute_epitope_conservation,
    compute_mutation_records,
    _LENGTH_TOLERANCE,
    _build_iedb_summary_rows,
)
from .io import (
    load_fasta_sequences,
    write_conservation_summary_xlsx,
    write_mutations_xlsx,
)
from .prompts import _ask_for_local_fasta_overrides, prompt_analysis_threshold, prompt_length_filter
from .render import _render_track_status_table, print_conservation_rich_table

def _inspect_track_fasta_status(
    project_name: str,
    track_id: str,
) -> dict:
    """
    Return a status dict describing the variant FASTA situation for one track.

    Used by `preflight()` to render the consolidated status table before the
    user is asked whether they want to provide local FASTA files. The dict
    has the shape:

        {
            "track_id":      str,
            "ref_length":    int | None,
            "fasta_path":    Path | None,
            "fasta_exists":  bool,
            "n_records_raw": int,
            "n_records_kept": int,
            "status":        str,   # "ok" | "missing_fasta" | "no_reference" | "length_mismatch"
        }
    """
    project_root = Path("projects") / project_name
    track_dir = project_root / "data" / "intermediate" / track_id
    input_track_dir = project_root / "data" / "input" / track_id

    reference_fasta_path = input_track_dir / get_step_filename("SEQUENCES", track_id, ext="fasta")
    reference_length: Optional[int] = None
    if reference_fasta_path.exists():
        reference_records = list(SeqIO.parse(str(reference_fasta_path), "fasta"))
        if reference_records:
            reference_length = len(reference_records[0].seq)

    variants_fasta_path = track_dir / "variants" / get_step_filename("VARIANTS", track_id, ext="fasta")
    fasta_exists = variants_fasta_path.exists() and variants_fasta_path.stat().st_size > 0

    n_records_raw = 0
    n_records_kept = 0
    if fasta_exists:
        kept_records, excluded_count = load_fasta_sequences(variants_fasta_path, reference_length or 0)
        n_records_kept = len(kept_records)
        n_records_raw = n_records_kept + excluded_count

    if not reference_length:
        status_label = "no_reference"
    elif not fasta_exists:
        status_label = "missing_fasta"
    elif n_records_kept == 0 and n_records_raw > 0:
        status_label = "length_mismatch"
    else:
        status_label = "ok"

    return {
        "track_id":       track_id,
        "ref_length":     reference_length,
        "fasta_path":     variants_fasta_path if fasta_exists else None,
        "fasta_exists":   fasta_exists,
        "n_records_raw":  n_records_raw,
        "n_records_kept": n_records_kept,
        "status":         status_label,
    }


def _identify_problematic_tracks(
    project_name: str,
    track_outcomes: dict,
    minimum_acceptable_variants: int = 5,
) -> list[dict]:
    """
    Return tracks whose conservation analysis did not produce a meaningful
    result. A track is problematic when (a) it errored, (b) it ran but had
    fewer than `minimum_acceptable_variants` usable sequences, or (c) it ran
    with no variants at all (conservation_unknown).
    """
    problematic_track_rows: list[dict] = []

    for track_id, outcome_dict in track_outcomes.items():
        outcome_status = outcome_dict.get("status")

        if outcome_status == "error":
            problematic_track_rows.append({
                "track_id": track_id,
                "problem":  outcome_dict.get("error_message", "unknown error"),
                "n_used":   None,
            })
            continue

        audit_path = (
            Path("projects") / project_name / "data" / "intermediate" / track_id /
            "conservation" / get_step_filename("CONSERVATION_AUDIT", track_id, ext="json")
        )
        if not audit_path.exists():
            continue

        try:
            audit_payload = json.loads(audit_path.read_text(encoding="utf-8"))
        except Exception:
            continue

        n_variants_used = int(audit_payload.get("n_variants_used", 0))
        if n_variants_used < minimum_acceptable_variants:
            problematic_track_rows.append({
                "track_id": track_id,
                "problem":  f"only {n_variants_used} variants used"
                            f" (label_counts={audit_payload.get('label_counts', {})})",
                "n_used":   n_variants_used,
            })

    return problematic_track_rows


# ── Step ──────────────────────────────────────────────────────────────────────

class AnalyzeConservationStep(BaseTrackStep):
    step_name   = "analyze_conservation"
    description = (
        "Slides every ★ epitope across each variant sequence, records best-match "
        "identity, and labels each peptide as perfect / high / moderate / low. "
        "Also emits per-variant mutation verdicts (BLOSUM62 + MHC-I anchors P2 "
        "and PΩ). Qualitative: never removes epitopes."
    )
    long_description = (
        "Tells you whether your selected epitopes will still work against "
        "real-world strain variation. For each ★ epitope and each variant "
        "sequence, slides a window of the epitope's length across the variant "
        "and records the best-matching identity.\n\n"
        "Results are [bold]qualitative annotations[/bold]; this step never "
        "removes epitopes. It tells downstream decision-makers which epitopes "
        "are 'safe across strains' vs. which are 'present only in the reference'."
    )
    methodology = (
        "1. Loads ★ representatives from select_representatives.\n"
        "2. Loads variant sequences from search_variants (or a user-supplied "
        "FASTA; the preflight asks).\n"
        "3. Filters variants by length tolerance (default ± 20% of reference).\n"
        "4. For each (epitope, variant) pair: sliding-window identity at the "
        "epitope's length; records max identity.\n"
        "5. Aggregates per epitope: mean of max-identities, min/max, tier "
        "fractions (≥ 100% / 90% / 80% / threshold).\n"
        "6. Labels: [bold]perfect[/bold] (mean=100%), [bold]high[/bold] (≥90%), "
        "[bold]moderate[/bold] (≥80%), [bold]low[/bold] (<80%), and "
        "[bold]conservation_unknown[/bold] when no variant FASTA is available.\n"
        "7. Mutation verdicts per (epitope, variant) using BLOSUM62 substitution "
        "scores at MHC-I anchor positions (P2 + PΩ)."
    )
    references = [
        {
            'authors': 'Bui HH, Sidney J, Li W, Fusseder N, Sette A.',
            'title':   'Development of an epitope conservancy analysis tool to facilitate the design of epitope-based diagnostics and vaccines',
            'journal': 'BMC Bioinformatics',
            'year':    2007,
            'doi':     '10.1186/1471-2105-8-361',
        },
        {
            'authors': 'Henikoff S, Henikoff JG.',
            'title':   'Amino acid substitution matrices from protein blocks',
            'journal': 'PNAS',
            'year':    1992,
            'doi':     '10.1073/pnas.89.22.10915',
        },
    ]
    data_format = (
        "Inputs are automatic:\n"
        "  • CLUSTER_REPR_{track_id}.csv from select_representatives (only ★ rows).\n"
        "  • VARIANTS_{track_id}.fasta from search_variants, or, if missing or "
        "you have a curated set, the preflight prompt asks for a path.\n\n"
        "You will be asked once for the [bold]identity threshold[/bold] "
        "(default 0.80 = 80%)."
    )
    outputs_overview = (
        "[bold]CONSERVATION_VIEW_{track_id}.csv[/bold]      slim per-step view (peptide + min/max/avg identity + label).\n"
        "[bold]CONSERVATION_{track_id}.csv/.xlsx[/bold]      IEDB-style summary, one row per ★ representative.\n"
        "[bold]CONSERVATION_MUTATIONS_{track_id}.xlsx[/bold] per (epitope, variant) breakdown with anchor flags + MHC verdict.\n"
        "[bold]CONSERVATION_HEATMAP_{track_id}.png[/bold]    dual-panel heatmap visualisation.\n"
        "[bold]CONSERVATION_AUDIT_{track_id}.json[/bold]     threshold, FASTA source, label counts, verdict counts."
    )
    tips = [
        "This step never removes epitopes; it only annotates. Use the labels downstream to prioritise.",
        "Low-conservation epitopes are not necessarily bad; they may be species-specific by design.",
        "If search_variants returned few or zero variants, supply a curated FASTA at the preflight prompt.",
        "The ±20% length filter avoids comparing partial/incomplete variant sequences.",
    ]

    @classmethod
    def preflight(cls, project_name: str, project_config: dict, track_ids: list[str],
                  is_rerun: bool = False) -> Optional[dict]:
        if not track_ids:
            return None

        track_status_rows = [
            _inspect_track_fasta_status(project_name, track_id)
            for track_id in track_ids
        ]

        status_panel = Panel(
            _render_track_status_table(track_status_rows),
            title="[bold]Variant FASTA status across all tracks[/bold]",
            border_style="cyan", box=box.ROUNDED,
        )

        if not is_interactive_session():
            console.print(status_panel)
            return None

        # Ask the local-FASTA question up front (before the full status table) so the
        # option isn't buried; the table then follows as context for picking track ids.
        try:
            response = input(
                "  Do you want to attach a local FASTA to any track? [y/N]: "
            ).strip().lower()
        except EOFError:
            response = ""

        console.print(status_panel)

        if response not in {"y", "yes"}:
            return None

        console.print(
            "[dim]  Type the track id (must match the table above), then the FASTA path. "
            "Leave the track id blank to stop.[/dim]"
        )
        fasta_overrides = _ask_for_local_fasta_overrides(track_status_rows)
        return fasta_overrides if fasta_overrides else None

    @classmethod
    def postflight(cls, project_name: str, project_config: dict, track_outcomes: dict) -> None:
        if not track_outcomes or not is_interactive_session():
            return

        problematic_track_rows = _identify_problematic_tracks(project_name, track_outcomes)
        if not problematic_track_rows:
            return

        recovery_table = Table(box=box.SIMPLE, header_style="bold yellow")
        recovery_table.add_column("Track", style="cyan")
        recovery_table.add_column("Variants used", justify="right")
        recovery_table.add_column("Problem")
        for problem_row in problematic_track_rows:
            recovery_table.add_row(
                problem_row["track_id"],
                "-" if problem_row["n_used"] is None else str(problem_row["n_used"]),
                problem_row["problem"],
            )

        console.print(Panel(
            recovery_table,
            title="[bold]Tracks with weak conservation results[/bold]",
            border_style="yellow", box=box.ROUNDED,
        ))

        console.print(
            "  Re-run weak tracks?  [cyan]f[/cyan] = with a local FASTA   "
            "[cyan]l[/cyan] = with the length filter OFF   [cyan]N[/cyan] = no"
        )
        try:
            response = input("  > ").strip().lower()
        except EOFError:
            response = ""
        if response not in {"f", "l"}:
            return

        from utils.pipeline_state import reset_track_step

        problematic_track_id_set = {row["track_id"] for row in problematic_track_rows}
        while True:
            try:
                target_track = input("  Track id to re-run (or blank to finish): ").strip()
            except EOFError:
                target_track = ""
            if not target_track:
                break
            if target_track not in problematic_track_id_set:
                console.print(f"  [yellow]Not in the list above: '{target_track}'.[/yellow]")
                continue

            reset_track_step(project_name, target_track, cls.step_name)
            rerun_step_instance = cls(
                project_name=project_name,
                project_config=project_config,
                track_id=target_track,
            )

            if response == "f":
                try:
                    fasta_path_input = input(f"  Local FASTA path for {target_track}: ").strip()
                except EOFError:
                    fasta_path_input = ""
                override_path = Path(fasta_path_input).expanduser() if fasta_path_input else None
                if override_path is None or not override_path.exists() or override_path.stat().st_size == 0:
                    console.print("  [red]File not found or empty. Skipping.[/red]")
                    continue
                rerun_step_instance.preflight_config = {target_track: str(override_path)}
                rerun_step_instance.execute(force_rerun=False)
            else:  # response == "l" — re-run with the length filter disabled
                rerun_step_instance._force_length_filter_off = True
                rerun_step_instance.execute(reconfigure=True)

    def describe_outputs(self) -> dict[Path, str]:
        conservation_dir = self.track_dir / "conservation"
        return {
            conservation_dir / get_step_filename("CONSERVATION_VIEW", self.track_id):
                "Slim per-step view: peptide + length + min/max/avg identity + conservation_label only.",
            conservation_dir / get_step_filename("CONSERVATION", self.track_id):
                "IEDB-style summary: per ★ rep with tier %, fractions, min/max/avg identity, label.",
            conservation_dir / get_step_filename("CONSERVATION", self.track_id, ext="xlsx"):
                "Same table with header-only styling; conservation_label cell coloured.",
            conservation_dir / get_step_filename("CONSERVATION_HEATMAP", self.track_id, ext="png"):
                "Dual-panel heatmap: position conservation + identity tiers (filtered to ≤2-mut variants).",
            conservation_dir / get_step_filename("CONSERVATION_MUTATIONS", self.track_id, ext="xlsx"):
                "Per (epitope, variant) breakdown with anchor flag, BLOSUM62 score, MHC verdict.",
            conservation_dir / get_step_filename("CONSERVATION_AUDIT", self.track_id, ext="json"):
                "Run audit: threshold, FASTA source, variants used, label counts, verdict counts.",
        }

    def run(self, input_data=None):
        analysis_threshold = prompt_analysis_threshold(
            self.project_name, self.project_config, is_rerun=self.is_rerun
        )
        # postflight recovery can force the filter off deterministically (no re-prompt)
        if getattr(self, "_force_length_filter_off", False):
            apply_length_filter = False
            console.print("[dim]→ Length filter OFF (postflight recovery).[/dim]")
        else:
            apply_length_filter = prompt_length_filter(self.is_rerun)

        # ── Load ★ representatives ────────────────────────────────────────────
        clusters_dir        = self.track_dir / "clusters"
        representatives_csv = clusters_dir / get_step_filename("CLUSTER_REPR", self.track_id)
        if not representatives_csv.exists():
            raise FileNotFoundError(
                f"select_representatives output not found: {representatives_csv}\n"
                "Run 'select_representatives' before 'analyze_conservation'."
            )

        df_repr  = pd.read_csv(representatives_csv)
        df_stars = df_repr[df_repr[COLUMN_BEST_REPRESENTATIVE] == STAR_MARKER].copy()
        if df_stars.empty:
            raise ValueError(
                f"No ★ representatives found in {representatives_csv.name}. "
                "Run 'select_representatives' first."
            )

        # ── Reference length (for variant length filtering) ───────────────────
        input_dir  = self.track_dir.parent.parent / "input" / self.track_id
        ref_fasta  = input_dir / get_step_filename("SEQUENCES", self.track_id, ext="fasta")
        ref_length = 0
        if ref_fasta.exists():
            ref_recs = list(SeqIO.parse(str(ref_fasta), "fasta"))
            if ref_recs:
                ref_length = len(ref_recs[0].seq)
                filter_note = (
                    f"±{int(_LENGTH_TOLERANCE*100)}% filter applied"
                    if apply_length_filter else "length filter OFF"
                )
                console.print(
                    f"[dim]→ Reference length: {ref_length} aa ({filter_note})[/dim]"
                )

        # ── Resolve FASTA source ──────────────────────────────────────────────
        variants_dir = self.track_dir / "variants"
        fasta_path   = variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta")
        fasta_source = "none"
        records: list = []
        n_excl = 0

        override_fasta_path: Optional[Path] = None
        if isinstance(self.preflight_config, dict):
            override_value = self.preflight_config.get(self.track_id)
            if override_value:
                override_fasta_path = Path(override_value).expanduser()

        if override_fasta_path is not None and override_fasta_path.exists() and override_fasta_path.stat().st_size > 0:
            records, n_excl = load_fasta_sequences(override_fasta_path, ref_length, apply_length_filter)
            fasta_source = "user_provided"
            msg = f"[dim]→ Using user-provided FASTA {override_fasta_path.name} ({len(records)} sequences"
            if n_excl:
                msg += f", {n_excl} excluded by length filter"
            console.print(msg + ")[/dim]")
        elif fasta_path.exists() and fasta_path.stat().st_size > 0:
            records, n_excl = load_fasta_sequences(fasta_path, ref_length, apply_length_filter)
            fasta_source    = "variants_cache"
            msg = f"[dim]→ Using {fasta_path.name} ({len(records)} sequences"
            if n_excl:
                msg += f", {n_excl} excluded by length filter"
            console.print(msg + ")[/dim]")
        else:
            console.print(
                f"[yellow]⚠ No variant FASTA available for {self.track_id}: "
                "epitopes will be labelled 'conservation_unknown'. Provide a "
                "local FASTA after the loop completes if needed.[/yellow]"
            )

        # Loud warning when the length filter wiped out every variant — common for
        # divergent proteins (e.g. membrane proteins) whose variants differ in length.
        if apply_length_filter and not records and n_excl > 0:
            console.print(Panel(
                f"[bold]Length filter removed all {n_excl} variant(s)[/bold] for "
                f"{self.track_id}.\nThe protein's variants differ too much in length "
                f"(±{int(_LENGTH_TOLERANCE*100)}% of {ref_length} aa). Conservation will be "
                f"'unknown'.\n[dim]Fix: redo this step ([cyan]r[/cyan] in the menu) and turn "
                f"the length filter OFF.[/dim]",
                title="[yellow]⚠ Divergent protein: length filter wiped variants[/yellow]",
                border_style="yellow", box=box.ROUNDED,
            ))

        # ── Analyse each peptide ──────────────────────────────────────────────
        metrics_rows:        list[dict] = []
        alignment_data:      list[dict] = []
        all_mutation_records: list[dict] = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(
                style="yellow",
                complete_style="green",
                finished_style="green",
                pulse_style="bold yellow",
            ),
            MofNCompleteColumn(),
            console=console,
            transient=True,
        ) as conservation_progress_bar:
            conservation_task_id = conservation_progress_bar.add_task(
                "  Analyzing per-epitope conservation",
                total=len(df_stars),
            )
            for _, repr_row in df_stars.iterrows():
                peptide = repr_row[COLUMN_PEPTIDE]
                alleles_united_value = repr_row.get(COLUMN_ALLELES_UNITED, "")
                alleles_united_str   = str(alleles_united_value) if pd.notna(alleles_united_value) else ""

                metrics, alignment_tuples = compute_epitope_conservation(
                    peptide, records, analysis_threshold
                )

                peptide_mutation_records = compute_mutation_records(
                    peptide, alignment_tuples, alleles_united_str
                )
                all_mutation_records.extend(peptide_mutation_records)

                variants_exact_labels = [
                    label for label, identity, _window in alignment_tuples
                    if identity == 1.0
                ]
                tolerable_records = [
                    r for r in peptide_mutation_records
                    if r["mhc_verdict"] in {"excellent_match", "tolerated"}
                ]
                n_excellent = sum(1 for r in peptide_mutation_records if r["mhc_verdict"] == "excellent_match")
                n_tolerated = sum(1 for r in peptide_mutation_records if r["mhc_verdict"] == "tolerated")

                metrics["n_excellent_match"]    = n_excellent
                metrics["n_tolerated"]          = n_tolerated
                metrics["variants_exact_match"] = "; ".join(variants_exact_labels)
                metrics["variants_tolerable"]   = "; ".join(
                    f"{r['variant_accession']}[{r['mutations']}]" for r in tolerable_records
                )

                metrics_rows.append(metrics)

                alignment_data.append({
                    "peptide":            peptide,
                    "alignment_tuples":   alignment_tuples,
                    "conservation_label": metrics["conservation_label"],
                    "n_exact_match":      metrics["n_exact_match"],
                    "n_identity_90":      metrics["n_identity_90"],
                    "n_identity_80":      metrics["n_identity_80"],
                    "n_passed_threshold": metrics["n_passed_threshold"],
                })

                conservation_progress_bar.advance(conservation_task_id)

        flush_stdin()

        # ── Build full result DataFrame (metrics + repr columns) ──────────────
        result_df = pd.concat(
            [df_stars.reset_index(drop=True), pd.DataFrame(metrics_rows)],
            axis=1,
        )

        # ── Save outputs ──────────────────────────────────────────────────────
        conservation_dir = self.track_dir / "conservation"
        conservation_dir.mkdir(parents=True, exist_ok=True)

        output_csv             = conservation_dir / get_step_filename("CONSERVATION",            self.track_id)
        output_xlsx            = conservation_dir / get_step_filename("CONSERVATION",            self.track_id, ext="xlsx")
        output_heatmap_png     = conservation_dir / get_step_filename("CONSERVATION_HEATMAP",    self.track_id, ext="png")
        output_mutations_xlsx  = conservation_dir / get_step_filename("CONSERVATION_MUTATIONS",  self.track_id, ext="xlsx")
        audit_path             = conservation_dir / get_step_filename("CONSERVATION_AUDIT",      self.track_id, ext="json")

        # Drop deprecated files if present from previous runs
        for deprecated_name in ("CONSERVATION_VISUAL", "CONSERVATION_POSITIONS"):
            deprecated_path = conservation_dir / get_step_filename(deprecated_name, self.track_id, ext="xlsx")
            if deprecated_path.exists():
                deprecated_path.unlink()
            deprecated_png_path = conservation_dir / get_step_filename(deprecated_name, self.track_id, ext="png")
            if deprecated_png_path.exists():
                deprecated_png_path.unlink()

        iedb_summary_df = _build_iedb_summary_rows(result_df)
        write_user_facing_csv(iedb_summary_df, output_csv)
        write_conservation_summary_xlsx(iedb_summary_df, output_xlsx)

        view_columns = ["peptide", "length", "min_identity", "max_identity", "avg_identity", "conservation_label"]
        view_csv = conservation_dir / get_step_filename("CONSERVATION_VIEW", self.track_id)
        write_user_facing_csv(iedb_summary_df[view_columns], view_csv)
        write_mutations_xlsx(all_mutation_records, output_mutations_xlsx)
        write_conservation_dual_panel_png(
            alignment_data, self.track_id, len(records), analysis_threshold, output_heatmap_png,
        )

        # ── Console summary ───────────────────────────────────────────────────
        print_conservation_rich_table(result_df, self.track_id, analysis_threshold)

        n_analyzed   = len(result_df)
        label_counts = result_df["conservation_label"].value_counts().to_dict()
        mean_all     = round(float(result_df["mean_max_identity"].mean()), 4) if n_analyzed else 0.0

        verdict_counts = dict(Counter(r["mhc_verdict"] for r in all_mutation_records))
        n_epitopes_with_tolerant = len({r["peptide"] for r in all_mutation_records})

        console.print(
            f"\n[bold green]RESULT: {n_analyzed} representatives analysed "
            f"({fasta_source}, threshold={int(round(analysis_threshold * 100))}%).[/bold green]\n"
        )

        # ── Audit JSON ────────────────────────────────────────────────────────
        audit = {
            "timestamp":                         datetime.datetime.now().isoformat(),
            "track_id":                          self.track_id,
            "analysis_threshold":                analysis_threshold,
            "fasta_source":                      fasta_source,
            "n_variants_used":                   len(records),
            "ref_length_aa":                     ref_length,
            "length_filter_tolerance":           _LENGTH_TOLERANCE,
            "length_filter_applied":             apply_length_filter,
            "n_representatives_analyzed":        n_analyzed,
            "label_counts":                      label_counts,
            "mean_identity_across_all_epitopes": mean_all,
            "n_pairs_within_2mut":               len(all_mutation_records),
            "verdict_counts":                    verdict_counts,
            "n_epitopes_with_tolerant_variants": n_epitopes_with_tolerant,
            "substitution_criterion":            "BLOSUM62 >= 0",
            "anchor_positions":                  "P2 + PΩ (universal HLA-I)",
            "reference_thresholds": {
                "high":     config.CONSERVATION_HIGH_THRESHOLD,
                "moderate": config.CONSERVATION_MODERATE_THRESHOLD,
            },
            "outputs": {
                "csv":            str(output_csv),
                "summary_xlsx":   str(output_xlsx),
                "heatmap_png":    str(output_heatmap_png),
                "mutations_xlsx": str(output_mutations_xlsx),
                "audit":          str(audit_path),
            },
        }
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False))

        return {
            "output_csv":            str(output_csv),
            "output_xlsx":           str(output_xlsx),
            "output_heatmap_png":    str(output_heatmap_png),
            "output_mutations_xlsx": str(output_mutations_xlsx),
            "n_analyzed":            n_analyzed,
            "fasta_source":          fasta_source,
            "label_counts":          label_counts,
        }
