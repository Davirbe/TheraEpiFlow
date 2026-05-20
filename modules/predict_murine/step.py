"""Orchestration for predict_murine: runs NetMHCpan + MHCFlurry against H-2
alleles for the ★ epitopes and writes the murine prediction tables."""

import datetime
import json
import time
from pathlib import Path

import pandas as pd
from rich import box
from rich.panel import Panel
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn
from rich.text import Text

import config
from modules.base_step import BaseTrackStep
from modules.predict_binding import (
    _run_mhcflurry_with_progress,
    _run_netmhcpan_iedb_silent,
)
from utils.console import console, flush_stdin
from utils.naming import COLUMN_PEPTIDE, get_step_filename
from utils.project_manager import save_project_config
from utils.step_summary import print_step_summary

from .core import (
    _LABEL_BORDERLINE,
    _LABEL_GOOD,
    _LABEL_NON_BINDER,
    _LABEL_OPTIMAL,
    _aggregate_per_peptide,
    _assign_binder_label,
    _build_synthetic_seqrecords,
    _combine_per_tool_long_tables,
    _load_star_peptides,
)
from .prompts import _ask_murine_strain

# ── Step class ────────────────────────────────────────────────────────────────

class PredictMurineStep(BaseTrackStep):
    step_name = 'predict_murine'
    description = (
        "Runs NetMHCpan EL + MHCFlurry on the ★ representatives against "
        "H-2 (murine) alleles to verify that the human-selected epitopes "
        "would also bind in a mouse-model strain. Qualitative — no "
        "peptides are removed; every ★ gets per-allele records plus an "
        "aggregated row with best percentile, alleles bound (best first), "
        "and the four-tier binder label."
    )
    long_description = (
        "In vivo vaccine validation typically starts in mice — but a peptide "
        "that binds human HLA-I doesn't automatically bind murine H-2. This "
        "step re-runs the two binding predictors with H-2 alleles instead of "
        "HLA-I and labels each ★ epitope with how well it translates to the "
        "murine model.\n\n"
        "Results are [bold]qualitative annotations[/bold] — never removes "
        "epitopes. Used by curate_murine downstream to flag which candidates "
        "are immediately testable in vivo."
    )
    methodology = (
        "1. Loads ★ representatives from select_representatives.\n"
        "2. For each ★ epitope: peptide-direct mode (no protein context "
        "re-enumeration), feeds the peptide list to NetMHCpan and MHCFlurry "
        "against the configured H-2 allele set.\n"
        "3. Strain groups available: [bold]C57BL/6[/bold] (H-2Kb, H-2Db), "
        "[bold]BALB/c[/bold] (H-2Kd, H-2Dd, H-2Ld), or [bold]complete[/bold] "
        "(all H-2 alleles in both predictors).\n"
        "4. Binder tier per peptide: [bold]optimal[/bold] (best %rank ≤ 0.5), "
        "[bold]strong[/bold] (≤ 2), [bold]borderline[/bold] (≤ 2.5), "
        "[bold]none[/bold] (above).\n"
        "5. Promiscuous criterion: ≥ 2 H-2 alleles bound at strong-level."
    )
    references = [
        {
            'authors': 'Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M.',
            'title':   'NetMHCpan-4.1 and NetMHCIIpan-4.0',
            'journal': 'Nucleic Acids Research',
            'year':    2020,
            'doi':     '10.1093/nar/gkaa379',
        },
        {
            'authors': "O'Donnell TJ, Rubinsteyn A, Laserson U.",
            'title':   'MHCflurry 2.0',
            'journal': 'Cell Systems',
            'year':    2020,
            'doi':     '10.1016/j.cels.2020.06.010',
        },
    ]
    data_format = (
        "Input is automatic — ★ representatives from select_representatives. "
        "You will be asked once for the [bold]strain group[/bold] (C57BL/6 / "
        "BALB/c / complete) — choose based on which mouse strain your lab "
        "uses for immunisation."
    )
    outputs_overview = (
        "[bold]MURINE_{track_id}.csv[/bold]      — long format: one row per (peptide, allele, tool).\n"
        "[bold]MURINE_AGG_{track_id}.csv[/bold]  — one row per ★ peptide with best percentile, H-2 alleles bound, count, tier.\n"
        "[bold]MURINE_VIEW_{track_id}.csv[/bold] — slim view: peptide + tier + best percentile + num alleles bound.\n"
        "[bold]MURINE_AUDIT_{track_id}.json[/bold] — strain group, allele set, totals per tier."
    )
    tips = [
        "Pick the strain group your wet lab actually uses — running 'complete' is heavier and rarely needed.",
        "Peptides labelled 'optimal' or 'strong' in murine are immediate candidates for in vivo testing.",
        "The MHCFlurry first-run delay applies here too — same TF model load as predict_binding.",
        "Inputs are peptide-direct, not from FASTA — no length re-enumeration happens at this step.",
    ]

    def describe_outputs(self) -> dict:
        murine_dir = self.track_dir / "murine"
        return {
            murine_dir / get_step_filename("MURINE", self.track_id):
                "Long-format table — one row per (peptide, allele, tool).",
            murine_dir / get_step_filename("MURINE_AGG", self.track_id):
                "Aggregated table — one row per ★ peptide with best percentile, "
                "alleles bound (best first), and binder label.",
            murine_dir / get_step_filename("MURINE_VIEW", self.track_id):
                "Slim view — peptide + tier + best percentile + num alleles bound.",
            murine_dir / get_step_filename("MURINE_AUDIT", self.track_id, ext='json'):
                "Run audit — strain, alleles, lengths, label counts.",
        }

    def run(self, input_data=None):
        run_start_time = time.time()

        strain_name, h2_alleles = _ask_murine_strain(self.project_name, self.project_config)
        peptide_lengths = list(self.project_config.get('peptide_lengths') or [9])

        star_peptides_df = _load_star_peptides(self.track_dir, self.track_id)
        star_peptide_sequences = star_peptides_df[COLUMN_PEPTIDE].astype(str).tolist()
        synthetic_seqrecords = _build_synthetic_seqrecords(star_peptide_sequences)

        console.print(Panel(
            Text.from_markup(
                f"[bold]Track:[/bold] {self.track_id}\n"
                f"[bold]Strain:[/bold] {strain_name}  "
                f"[dim]({', '.join(h2_alleles)})[/dim]\n"
                f"[bold]★ peptides:[/bold] {len(star_peptide_sequences)}\n"
                f"[bold]Peptide lengths:[/bold] {peptide_lengths}"
            ),
            title="Murine prediction",
            border_style="cyan",
            box=box.ROUNDED,
        ))

        console.print()
        console.print(Panel(
            Text.from_markup(
                "[bold yellow]⚠ MHCFlurry first-time load: ~30-60 seconds.[/bold yellow]\n"
                "TensorFlow models are being read from disk and warmed up. The bar may "
                "briefly pause during native imports — that's expected.\n"
                "[dim]Please don't press any keys until the load finishes.[/dim]"
            ),
            border_style="yellow", box=box.ROUNDED, padding=(0, 1),
        ))

        # ── Run both predictors with a shared Progress UI ─────────────────────
        with Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.fields[label]:<10}[/bold blue]"),
            BarColumn(
                bar_width=30,
                style="yellow",
                complete_style="green",
                finished_style="green",
                pulse_style="bold yellow",
            ),
            MofNCompleteColumn(),
            TextColumn("{task.description}"),
            console=console,
            transient=False,
        ) as progress_bar:
            netmhcpan_task_id = progress_bar.add_task(
                "sending request to IEDB",
                total=1,
                label="NetMHCpan",
            )
            mhcflurry_task_id = progress_bar.add_task(
                "preparing",
                total=None,
                label="MHCFlurry",
            )

            netmhcpan_raw_df = _run_netmhcpan_iedb_silent(
                synthetic_seqrecords, h2_alleles, peptide_lengths,
            )
            progress_bar.update(
                netmhcpan_task_id,
                completed=1,
                description=f"[green]done — {len(netmhcpan_raw_df):,} rows[/green]",
            )

            mhcflurry_raw_df, _mhcflurry_capture_log = _run_mhcflurry_with_progress(
                sequence_records=synthetic_seqrecords,
                hla_alleles=h2_alleles,
                peptide_lengths=peptide_lengths,
                progress_handle=progress_bar,
                progress_task_id=mhcflurry_task_id,
            )
            progress_bar.update(
                mhcflurry_task_id,
                description=f"[green]done — {len(mhcflurry_raw_df):,} rows[/green]",
            )

        flush_stdin()

        # ── Combine per (peptide, allele), then aggregate per peptide ─────────
        best_per_pair_df, full_long_df = _combine_per_tool_long_tables(
            netmhcpan_raw_df, mhcflurry_raw_df,
        )
        best_per_pair_df['binder_label'] = best_per_pair_df['percentile'].apply(_assign_binder_label)
        full_long_df['binder_label']     = full_long_df['percentile'].apply(_assign_binder_label)

        aggregated_per_peptide_df = _aggregate_per_peptide(
            best_per_pair_df, star_peptide_sequences,
        )

        # ── Persist outputs ───────────────────────────────────────────────────
        murine_dir = self.track_dir / "murine"
        murine_dir.mkdir(parents=True, exist_ok=True)

        long_csv_path       = murine_dir / get_step_filename("MURINE",       self.track_id)
        aggregated_csv_path = murine_dir / get_step_filename("MURINE_AGG",   self.track_id)
        audit_json_path     = murine_dir / get_step_filename("MURINE_AUDIT", self.track_id, ext='json')
        view_csv_path       = murine_dir / get_step_filename("MURINE_VIEW", self.track_id)

        if 'percentile' in full_long_df.columns:
            full_long_df['percentile'] = pd.to_numeric(
                full_long_df['percentile'], errors='coerce'
            ).round(2)
        full_long_df.to_csv(long_csv_path, index=False)
        aggregated_per_peptide_df.to_csv(aggregated_csv_path, index=False)

        # Slim VIEW: only the columns produced here — drops alleles_bound string for terminal browsing.
        view_columns = ['peptide', 'best_percentile_label', 'best_percentile_value', 'num_murine_alleles_bound']
        present_view_columns = [c for c in view_columns if c in aggregated_per_peptide_df.columns]
        aggregated_per_peptide_df[present_view_columns].to_csv(view_csv_path, index=False)

        label_counts = aggregated_per_peptide_df['best_percentile_label'].value_counts().to_dict()
        for label_name in (_LABEL_OPTIMAL, _LABEL_GOOD, _LABEL_BORDERLINE, _LABEL_NON_BINDER):
            label_counts.setdefault(label_name, 0)

        num_with_any_binder = int(
            (aggregated_per_peptide_df['num_murine_alleles_bound'] > 0).sum()
        )

        audit_payload = {
            "timestamp":          datetime.datetime.now().isoformat(),
            "track_id":           self.track_id,
            "strain":             strain_name,
            "h2_alleles":         h2_alleles,
            "peptide_lengths":    peptide_lengths,
            "n_star_peptides":    len(star_peptide_sequences),
            "label_counts": {
                _LABEL_OPTIMAL:    int(label_counts[_LABEL_OPTIMAL]),
                _LABEL_GOOD:       int(label_counts[_LABEL_GOOD]),
                _LABEL_BORDERLINE: int(label_counts[_LABEL_BORDERLINE]),
                _LABEL_NON_BINDER: int(label_counts[_LABEL_NON_BINDER]),
            },
            "num_with_any_binder": num_with_any_binder,
            "netmhcpan_rows":      int(len(netmhcpan_raw_df)),
            "mhcflurry_rows":      int(len(mhcflurry_raw_df)),
            "output_long_csv":     str(long_csv_path),
            "output_agg_csv":      str(aggregated_csv_path),
        }
        audit_json_path.write_text(json.dumps(audit_payload, indent=2, ensure_ascii=False))

        # ── Step summary ──────────────────────────────────────────────────────
        elapsed_seconds = time.time() - run_start_time
        print_step_summary(
            step_title="Murine prediction complete",
            elapsed_seconds=elapsed_seconds,
            narrative_lines=[
                f"Strain: [bold]{strain_name}[/bold] "
                f"({len(h2_alleles)} H-2 alleles).",
                f"Scored [bold]{len(star_peptide_sequences)}[/bold] ★ peptides "
                f"with NetMHCpan EL + MHCFlurry.",
                f"Labels: "
                f"optimal={label_counts[_LABEL_OPTIMAL]}, "
                f"good={label_counts[_LABEL_GOOD]}, "
                f"borderline={label_counts[_LABEL_BORDERLINE]}, "
                f"non_binder={label_counts[_LABEL_NON_BINDER]}.",
                f"[bold]{num_with_any_binder}[/bold] of "
                f"{len(star_peptide_sequences)} peptides bind at least one H-2 allele.",
            ],
            output_files=[long_csv_path, aggregated_csv_path, audit_json_path],
        )

        return {
            "output_long_csv":    str(long_csv_path),
            "output_agg_csv":     str(aggregated_csv_path),
            "output_audit_json":  str(audit_json_path),
            "n_star_peptides":    len(star_peptide_sequences),
            "n_alleles":          len(h2_alleles),
            "label_counts":       audit_payload['label_counts'],
        }
