"""
Orchestration for the population_coverage step: reads the ★ representatives,
drives the coverage math, writes every output and records the audit JSON.
"""

import datetime
import json
from pathlib import Path
from typing import Optional

import pandas as pd
from rich import box
from rich.panel import Panel
from rich.table import Table

from modules.base_step import BaseTrackStep
from utils.console import console, is_interactive_session
from utils.naming import (
    COLUMN_ALLELES_UNITED,
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_NUM_ALLELES_UNITED,
    COLUMN_PEPTIDE,
    STAR_MARKER,
    get_step_filename,
)

from .charts import write_coverage_matrix_png, write_hit_chart_png
from .core import (
    _MHC_CLASS,
    _build_population_freq_map,
    _compute_epitope_coverage,
    _list_available_populations,
    _load_population_database,
    _parse_alleles_united,
    _safe_filename_token,
)
from .io import write_detail_csv, write_summary_csv, write_summary_xlsx, write_view_csv
from .prompts import _prompt_coverage_cutoff, _prompt_populations
from .render import _print_console_summary


class PopulationCoverageStep(BaseTrackStep):
    step_name   = "population_coverage"
    description = (
        "For each ★ epitope, computes the fraction of one or more human "
        "populations that carries at least one of the epitope's HLA alleles, "
        "using the vendored IEDB allele-frequency database (diploid coverage "
        "model). Qualitative — never removes epitopes."
    )
    long_description = (
        "Tells you how many people in a given population are expected to "
        "respond to each ★ epitope. Uses the IEDB Population Coverage Tool "
        "methodology (Bui 2006), shipped as a vendored allele-frequency "
        "database — no internet calls.\n\n"
        "Results are [bold]qualitative annotations[/bold] — coverage is "
        "appended to each ★ epitope, not used to filter. The decision of "
        "which populations to target belongs to you (and your collaborators)."
    )
    methodology = (
        "1. Reads ★ representatives from select_representatives (uses the "
        "alleles_united union column generated there).\n"
        "2. For each (population, epitope) pair: looks up the diploid "
        "frequency of every HLA allele the epitope binds, in the vendored "
        "IEDB pickle.\n"
        "3. Coverage formula (Bui 2006): probability that a random individual "
        "in the population carries ≥ 1 of the epitope's alleles, computed "
        "from per-allele frequencies under Hardy-Weinberg equilibrium.\n"
        "4. Per-population CSV + a per-(epitope, population) PNG hit chart + "
        "a global matrix PNG (when ≥ 2 populations are selected)."
    )
    references = [
        {
            'authors': 'Bui HH, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A.',
            'title':   'Predicting population coverage of T-cell epitope-based diagnostics and vaccines',
            'journal': 'BMC Bioinformatics',
            'year':    2006,
            'doi':     '10.1186/1471-2105-7-153',
        },
    ]
    data_format = (
        "Input is automatic — uses CLUSTER_REPR_{track_id}.csv from "
        "select_representatives (★ rows only). You will be asked once for the "
        "populations to evaluate (e.g. 'World', 'Europe', 'East Asia', "
        "'South America') — multi-select."
    )
    outputs_overview = (
        "[bold]COVERAGE_{track_id}.csv[/bold]                — long format: one row per (peptide, population) with coverage % + metadata.\n"
        "[bold]COVERAGE_{track_id}.xlsx[/bold]               — same data, formatted spreadsheet.\n"
        "[bold]COVERAGE_VIEW_{track_id}.csv[/bold]           — slim view: peptide + population + coverage_pct only.\n"
        "[bold]COVERAGE_DETAIL_{population}_{track_id}.csv[/bold] — per-population IEDB-style detail.\n"
        "[bold]COVERAGE_HIT_CHART_{population}_{track_id}.png[/bold] — bar chart per population.\n"
        "[bold]COVERAGE_MATRIX_{track_id}.png[/bold]         — heatmap (when ≥ 2 populations).\n"
        "[bold]COVERAGE_AUDIT_{track_id}.json[/bold]         — populations queried, allele DB version, totals."
    )
    tips = [
        "Coverage is the EXPECTED fraction of responders — not a guarantee per individual.",
        "Pick 'World' for global vaccine candidates; pick region-specific populations to optimise local deployment.",
        "Alleles missing from the IEDB database default to 0 frequency in that population — they don't error out.",
        "This step runs fully offline using the vendored pickle in modules/population_coverage/.",
    ]

    @classmethod
    def preflight(
        cls,
        project_name: str,
        project_config: dict,
        track_ids: list[str],
        is_rerun: bool = False,
    ) -> Optional[dict]:
        if not track_ids:
            return None
        # Load the DB once at preflight to validate population names against it.
        population_db = _load_population_database()
        _prompt_populations(project_name, project_config, population_db, is_rerun=is_rerun)
        _prompt_coverage_cutoff(project_name, project_config, is_rerun=is_rerun)
        return None

    @classmethod
    def postflight(
        cls,
        project_name: str,
        project_config: dict,
        track_outcomes: dict,
    ) -> None:
        cutoff = project_config.get("coverage_minimum_pct")
        if cutoff is None or not track_outcomes or not is_interactive_session():
            return

        low_coverage_rows: list[tuple[str, str, float]] = []
        for track_id, outcome in track_outcomes.items():
            if outcome.get("status") != "completed":
                continue
            mean_coverage_per_population = outcome.get("mean_coverage_per_population", {})
            for population, mean_value in mean_coverage_per_population.items():
                if float(mean_value) < float(cutoff):
                    low_coverage_rows.append((track_id, population, float(mean_value)))

        if not low_coverage_rows:
            return

        low_table = Table(box=box.SIMPLE, header_style="bold yellow")
        low_table.add_column("Track", style="cyan")
        low_table.add_column("Population")
        low_table.add_column("Mean coverage", justify="right")
        for track_id, population, mean_value in low_coverage_rows:
            low_table.add_row(track_id, population, f"{mean_value:.2f}%")
        console.print(Panel(
            low_table,
            title=f"[bold]Below {cutoff:g}% cutoff (informational)[/bold]",
            border_style="yellow", box=box.ROUNDED,
        ))

    def describe_outputs(self) -> dict[Path, str]:
        coverage_dir = self.track_dir / "coverage"
        described: dict[Path, str] = {
            coverage_dir / get_step_filename("COVERAGE", self.track_id):
                "Long-format coverage: one row per (peptide × population) + metadata columns.",
            coverage_dir / get_step_filename("COVERAGE", self.track_id, ext="xlsx"):
                "Same long-format table, coloured by coverage band.",
            coverage_dir / get_step_filename("COVERAGE_VIEW", self.track_id):
                "Slim view — peptide + population + coverage_pct only (no metadata).",
            coverage_dir / get_step_filename("COVERAGE_AUDIT", self.track_id, ext="json"):
                "Run metadata, populations selected, allele-match stats.",
        }
        populations = self.project_config.get("coverage_populations", ["World"])
        for population in populations:
            safe_name = _safe_filename_token(population)
            described[coverage_dir / f"COVERAGE_DETAIL_{safe_name}_{self.track_id}.csv"] = (
                f"IEDB-style per-epitope-per-HLA table for {population}."
            )
            described[coverage_dir / f"COVERAGE_HIT_CHART_{safe_name}_{self.track_id}.png"] = (
                f"IEDB-style hit chart PNG for {population}."
            )
        if len(populations) >= 2:
            described[coverage_dir / f"COVERAGE_MATRIX_{self.track_id}.png"] = (
                "Comparative heatmap across selected populations."
            )
        return described

    def run(self, input_data=None):
        coverage_dir = self.track_dir / "coverage"
        coverage_dir.mkdir(parents=True, exist_ok=True)

        # ── Load inputs ──────────────────────────────────────────────────────
        cluster_repr_csv = (
            self.track_dir / "clusters" /
            get_step_filename("CLUSTER_REPR", self.track_id)
        )
        if not cluster_repr_csv.exists():
            raise FileNotFoundError(
                f"select_representatives output not found: {cluster_repr_csv}\n"
                "Run 'select_representatives' before 'population_coverage'."
            )

        repr_df = pd.read_csv(cluster_repr_csv)
        star_df = repr_df[repr_df[COLUMN_BEST_REPRESENTATIVE] == STAR_MARKER].copy()
        if star_df.empty:
            raise ValueError(
                f"No ★ representatives found in {cluster_repr_csv.name}. "
                "Run 'select_representatives' first."
            )

        # ── Resolve runtime configuration ────────────────────────────────────
        population_db = _load_population_database()
        populations: list[str] = list(
            self.project_config.get("coverage_populations") or ["World"]
        )
        cutoff: Optional[float] = self.project_config.get("coverage_minimum_pct")

        valid_populations = set(_list_available_populations(population_db))
        unknown_populations = [p for p in populations if p not in valid_populations]
        if unknown_populations:
            raise ValueError(
                "Unknown populations in project_config: "
                f"{unknown_populations}. Edit project_config.json or rerun preflight."
            )

        console.print(
            f"[dim]→ Computing coverage for {len(star_df)} ★ epitopes "
            f"× {len(populations)} population(s): "
            f"{', '.join(populations)}[/dim]"
        )

        # ── Compute coverage for every (peptide, population) ─────────────────
        summary_rows:           list[dict]              = []
        per_population_details: dict[str, list[dict]]   = {pop: [] for pop in populations}
        per_population_matched: dict[str, int]          = {pop: 0  for pop in populations}
        per_population_missed:  dict[str, int]          = {pop: 0  for pop in populations}

        for population in populations:
            freq_map = _build_population_freq_map(population_db, population)

            for _, peptide_row in star_df.iterrows():
                peptide          = peptide_row[COLUMN_PEPTIDE]
                epitope_alleles  = _parse_alleles_united(
                    peptide_row.get(COLUMN_ALLELES_UNITED, "")
                )
                coverage_pct, matched_alleles = _compute_epitope_coverage(
                    epitope_alleles, freq_map
                )

                summary_rows.append({
                    "peptide":        peptide,
                    "population":     population,
                    "mhc_class":      _MHC_CLASS,
                    "coverage_pct":   coverage_pct,
                    "n_hlas_used":    int(peptide_row.get(COLUMN_NUM_ALLELES_UNITED, len(epitope_alleles)))
                                       if pd.notna(peptide_row.get(COLUMN_NUM_ALLELES_UNITED, None))
                                       else len(epitope_alleles),
                    "n_hlas_in_db":   len(matched_alleles),
                    "alleles_united": ";".join(sorted(epitope_alleles)),
                })

                per_population_details[population].append({
                    "peptide":          peptide,
                    "coverage_pct":     coverage_pct,
                    "all_alleles":      epitope_alleles,
                    "matched_alleles":  matched_alleles,
                })

                per_population_matched[population] += len(matched_alleles)
                per_population_missed[population]  += len(epitope_alleles - matched_alleles)

        # Sort summary by peptide, then by configured population order.
        population_order_index = {p: i for i, p in enumerate(populations)}
        summary_rows.sort(key=lambda r: (r["peptide"], population_order_index[r["population"]]))

        # ── Write outputs ────────────────────────────────────────────────────
        output_csv  = coverage_dir / get_step_filename("COVERAGE", self.track_id)
        output_xlsx = coverage_dir / get_step_filename("COVERAGE", self.track_id, ext="xlsx")
        view_csv    = coverage_dir / get_step_filename("COVERAGE_VIEW", self.track_id)
        audit_path  = coverage_dir / get_step_filename("COVERAGE_AUDIT", self.track_id, ext="json")

        write_summary_csv(summary_rows, output_csv)
        write_summary_xlsx(summary_rows, output_xlsx, cutoff)
        write_view_csv(summary_rows, view_csv)

        detail_csv_paths:    dict[str, Path] = {}
        hit_chart_png_paths: dict[str, Path] = {}

        for population in populations:
            safe_name     = _safe_filename_token(population)
            detail_csv    = coverage_dir / f"COVERAGE_DETAIL_{safe_name}_{self.track_id}.csv"
            hit_chart_png = coverage_dir / f"COVERAGE_HIT_CHART_{safe_name}_{self.track_id}.png"

            # HLAs appearing in any ★ epitope AND present in the population DB,
            # sorted by per-locus-normalized frequency descending (most frequent
            # alleles appear leftmost in the hit chart — matches the prototype).
            freq_map = _build_population_freq_map(population_db, population)
            seen_alleles: set[str] = set()
            for result in per_population_details[population]:
                seen_alleles.update(result["all_alleles"])
            hlas_with_freq = sorted(
                [(allele, round(freq_map[allele] * 100, 2))
                 for allele in seen_alleles if allele in freq_map],
                key=lambda pair: pair[1], reverse=True,
            )

            write_detail_csv(
                population, per_population_details[population], hlas_with_freq, detail_csv,
            )
            write_hit_chart_png(
                population, self.track_id,
                per_population_details[population], hlas_with_freq, cutoff, hit_chart_png,
            )
            detail_csv_paths[population]    = detail_csv
            hit_chart_png_paths[population] = hit_chart_png

        matrix_png_path: Optional[Path] = None
        if len(populations) >= 2:
            matrix_png_path = coverage_dir / f"COVERAGE_MATRIX_{self.track_id}.png"
            write_coverage_matrix_png(
                self.track_id, summary_rows, populations, matrix_png_path,
            )

        # ── Console summary ──────────────────────────────────────────────────
        _print_console_summary(
            self.track_id, populations, summary_rows, len(star_df), cutoff,
        )

        # ── Audit JSON ───────────────────────────────────────────────────────
        summary_df = pd.DataFrame(summary_rows)
        mean_coverage_per_population: dict[str, float] = {}
        for population in populations:
            rows_for_pop = summary_df[summary_df["population"] == population]
            mean_coverage_per_population[population] = round(
                float(rows_for_pop["coverage_pct"].mean()) if len(rows_for_pop) else 0.0, 4
            )

        audit_payload = {
            "timestamp":                       datetime.datetime.now().isoformat(),
            "track_id":                        self.track_id,
            "mhc_class":                       _MHC_CLASS,
            "populations":                     populations,
            "coverage_minimum_pct":            cutoff,
            "n_star_epitopes":                 int(len(star_df)),
            "matched_alleles_per_population":  per_population_matched,
            "missed_alleles_per_population":   per_population_missed,
            "mean_coverage_per_population":    mean_coverage_per_population,
            "outputs": {
                "summary_csv":   str(output_csv),
                "summary_xlsx":  str(output_xlsx),
                "detail_csvs":   {pop: str(path) for pop, path in detail_csv_paths.items()},
                "hit_charts":    {pop: str(path) for pop, path in hit_chart_png_paths.items()},
                "matrix_png":    str(matrix_png_path) if matrix_png_path else None,
                "audit":         str(audit_path),
            },
        }
        audit_path.write_text(json.dumps(audit_payload, indent=2, ensure_ascii=False))

        return {
            "output_csv":                    str(output_csv),
            "output_xlsx":                   str(output_xlsx),
            "populations":                   populations,
            "n_star_epitopes":               int(len(star_df)),
            "mean_coverage_per_population":  mean_coverage_per_population,
        }
