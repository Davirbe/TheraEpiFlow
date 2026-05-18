"""curate_murine step.

Re-joins the qualitative annotations onto the human ★ table from
select_representatives. The cumulative-column pattern (each step appending
columns) stopped at select_representatives because conservation/coverage/
murine wrote side files. This step resumes it.

No filtering or ranking — every ★ peptide kept; missing annotations stay empty.

Inputs (track-relative):
    clusters/CLUSTER_REPR_{track_id}.csv       — ★ table (required)
    conservation/CONSERVATION_{track_id}.csv   — required
    coverage/COVERAGE_{track_id}.csv           — required (pivoted long→wide)
    murine/MURINE_AGG_{track_id}.csv           — optional

Outputs (murine/):
    CURATE_MURINE_{track_id}.csv       — full per-★-peptide master table
    CURATE_MURINE_VIEW_{track_id}.csv  — slim view
    CURATE_MURINE_AUDIT_{track_id}.json
"""

import datetime
import json
import time
from pathlib import Path

import pandas as pd
from rich import box
from rich.panel import Panel
from rich.text import Text

from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.naming import (
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_PEPTIDE,
    STAR_MARKER,
    get_step_filename,
)
from utils.step_summary import print_step_summary


_MURINE_COLUMN_RENAME = {
    'best_percentile_label': 'murine_label',
    'best_percentile_value': 'murine_best_percentile',
}

_MURINE_OUTPUT_COLUMNS = [
    'murine_label',
    'murine_best_percentile',
    'murine_alleles_bound',
    'num_murine_alleles_bound',
]

_VIEW_COLUMNS = [
    COLUMN_PEPTIDE,
    'alleles_united',
    'final_score',
    'conservation_label',
    'murine_label',
    'murine_best_percentile',
    'num_murine_alleles_bound',
]

# Conservation duplicates these columns from CLUSTER_REPR; drop on join.
_CONSERVATION_DROP_ON_JOIN = ['#', 'length', 'alleles_united', 'num_alleles_united']


def _load_star_human_table(track_dir: Path, track_id: str) -> pd.DataFrame:
    """Loads CLUSTER_REPR_{track_id}.csv and keeps only the ★ rows."""
    representatives_csv = track_dir / "clusters" / get_step_filename("CLUSTER_REPR", track_id)
    if not representatives_csv.exists():
        raise FileNotFoundError(
            f"select_representatives output not found: {representatives_csv}\n"
            "Run 'select_representatives' before 'curate_murine'."
        )

    representatives_df = pd.read_csv(representatives_csv)
    if COLUMN_BEST_REPRESENTATIVE not in representatives_df.columns:
        raise ValueError(
            f"Column '{COLUMN_BEST_REPRESENTATIVE}' not found in {representatives_csv.name}."
        )

    star_df = representatives_df[
        representatives_df[COLUMN_BEST_REPRESENTATIVE] == STAR_MARKER
    ].copy()
    if star_df.empty:
        raise ValueError(
            "No ★ representatives found in select_representatives output."
        )
    return star_df.reset_index(drop=True)


def _load_murine_aggregate_if_present(track_dir: Path, track_id: str) -> pd.DataFrame | None:
    """Loads MURINE_AGG_{track_id}.csv, or returns None if predict_murine was skipped.
    predict_murine is optional — when absent the master table just lacks the four murine columns."""
    murine_agg_csv = track_dir / "murine" / get_step_filename("MURINE_AGG", track_id)
    if not murine_agg_csv.exists():
        return None
    return pd.read_csv(murine_agg_csv)


def _load_conservation(track_dir: Path, track_id: str) -> pd.DataFrame:
    """Loads CONSERVATION_{track_id}.csv ready for left-join; raises if missing.
    Drops columns that would collide with CLUSTER_REPR (#, length, alleles_united, num_alleles_united)."""
    conservation_csv = track_dir / "conservation" / get_step_filename("CONSERVATION", track_id)
    if not conservation_csv.exists():
        raise FileNotFoundError(
            f"analyze_conservation output not found: {conservation_csv}\n"
            "Run 'analyze_conservation' before 'curate_murine'."
        )
    conservation_df = pd.read_csv(conservation_csv)
    if COLUMN_PEPTIDE not in conservation_df.columns:
        raise ValueError(
            f"Column '{COLUMN_PEPTIDE}' not found in {conservation_csv.name}."
        )
    drop_columns = [c for c in _CONSERVATION_DROP_ON_JOIN if c in conservation_df.columns]
    return conservation_df.drop(columns=drop_columns)


def _load_coverage_pivoted(track_dir: Path, track_id: str) -> pd.DataFrame:
    """Loads COVERAGE_{track_id}.csv (long) and pivots to wide; raises if missing.
    Returns one row per peptide with `coverage_{population}` columns holding coverage_pct."""
    coverage_csv = track_dir / "coverage" / get_step_filename("COVERAGE", track_id)
    if not coverage_csv.exists():
        raise FileNotFoundError(
            f"population_coverage output not found: {coverage_csv}\n"
            "Run 'population_coverage' before 'curate_murine'."
        )
    long_df = pd.read_csv(coverage_csv)
    required_columns = {COLUMN_PEPTIDE, 'population', 'coverage_pct'}
    if not required_columns.issubset(long_df.columns):
        raise ValueError(
            f"{coverage_csv.name} is missing required columns "
            f"{required_columns - set(long_df.columns)}."
        )
    wide_df = long_df.pivot_table(
        index=COLUMN_PEPTIDE,
        columns='population',
        values='coverage_pct',
        aggfunc='first',
    ).reset_index()
    wide_df.columns.name = None
    rename_map = {c: f'coverage_{c}' for c in wide_df.columns if c != COLUMN_PEPTIDE}
    return wide_df.rename(columns=rename_map)


def _build_full_master_table(
    human_star_df: pd.DataFrame,
    conservation_df: pd.DataFrame,
    coverage_pivoted_df: pd.DataFrame,
    murine_agg_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """Sequentially left-joins conservation, coverage, and (optional) murine onto the human ★ table.
    Conservation/coverage are required and always joined; murine is optional."""
    master_df = human_star_df.copy()
    master_df = master_df.merge(conservation_df, on=COLUMN_PEPTIDE, how='left')
    master_df = master_df.merge(coverage_pivoted_df, on=COLUMN_PEPTIDE, how='left')

    if murine_agg_df is not None:
        murine_renamed_df = murine_agg_df.rename(columns=_MURINE_COLUMN_RENAME)
        murine_keep_columns = [COLUMN_PEPTIDE] + _MURINE_OUTPUT_COLUMNS
        murine_renamed_df = murine_renamed_df[
            [c for c in murine_keep_columns if c in murine_renamed_df.columns]
        ]
        master_df = master_df.merge(murine_renamed_df, on=COLUMN_PEPTIDE, how='left')

    return master_df


def _build_label_histogram(joined_df: pd.DataFrame) -> dict:
    """Returns count per murine_label, with explicit 0 for missing tiers."""
    raw_counts = joined_df['murine_label'].fillna('null').value_counts().to_dict()
    histogram = {label: int(raw_counts.get(label, 0))
                 for label in ('optimal', 'good', 'borderline', 'non_binder', 'null')}
    return histogram


class CurateMurineStep(BaseTrackStep):
    step_name = 'curate_murine'
    description = (
        "Reassembles the per-track master table by left-joining conservation, "
        "population-coverage (pivoted long→wide), and murine annotations onto "
        "the human ★ table from select_representatives. No filtering — every "
        "★ peptide is kept; missing annotations stay as empty cells."
    )
    long_description = (
        "The cumulative-column pattern (each step appending columns to the "
        "previous CSV) ran until select_representatives. After that, "
        "analyze_conservation, population_coverage, and predict_murine wrote "
        "side files instead of extending the main table.\n\n"
        "This step resumes the pattern: it takes the ★ rows from "
        "select_representatives and left-joins the qualitative annotations "
        "back in, so a single per-track table carries all the per-epitope "
        "evidence (human binding + cluster scores + conservation + coverage "
        "per population + murine).\n\n"
        "No ranking decision is made here — that depends on the user's "
        "selection priorities and happens in the report. curate_murine just "
        "guarantees data completeness for integrate_data and generate_report."
    )
    methodology = (
        "1. Reads `clusters/CLUSTER_REPR_{track_id}.csv` and keeps only ★ rows.\n"
        "2. Optionally reads `conservation/CONSERVATION_{track_id}.csv` and "
        "left-joins it (drops `#`, `length`, and the alleles columns that "
        "would collide with the ★ table).\n"
        "3. Optionally reads `coverage/COVERAGE_{track_id}.csv` (long format) "
        "and pivots [italic]population → coverage_{population}[/italic] before "
        "left-joining.\n"
        "4. Reads `murine/MURINE_AGG_{track_id}.csv`, renames "
        "[italic]best_percentile_label → murine_label[/italic] and "
        "[italic]best_percentile_value → murine_best_percentile[/italic], "
        "then left-joins.\n"
        "5. Steps that the user did not run (conservation/coverage) simply "
        "produce no columns — peptides are never dropped.\n"
        "6. Writes the full master table + a slim VIEW + a JSON audit."
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
        "All inputs are picked up automatically — no prompts. Required: "
        "[bold]CLUSTER_REPR_{track_id}.csv[/bold] (★ table), "
        "[bold]CONSERVATION_{track_id}.csv[/bold] (per-epitope conservation), "
        "and [bold]COVERAGE_{track_id}.csv[/bold] (per-population coverage). "
        "Optional: [bold]MURINE_AGG_{track_id}.csv[/bold] — if you did not "
        "run predict_murine (mouse model not relevant for your project), the "
        "four murine columns are skipped, no error."
    )
    outputs_overview = (
        "[bold]CURATE_MURINE_{track_id}.csv[/bold] — one row per ★ peptide "
        "with: CLUSTER_REPR columns + conservation columns + "
        "coverage_{population} columns + 4 murine columns (if predict_murine "
        "was run).\n"
        "[bold]CURATE_MURINE_VIEW_{track_id}.csv[/bold] — slim view: peptide, "
        "alleles_united, final_score, conservation_label, murine_label, "
        "murine_best_percentile, num_murine_alleles_bound (murine fields are "
        "still listed but empty if predict_murine was skipped).\n"
        "[bold]CURATE_MURINE_AUDIT_{track_id}.json[/bold] — total column "
        "count, list of coverage populations, murine_present flag, label "
        "histogram (when murine is present)."
    )
    tips = [
        "predict_murine is optional — skip it if you don't have a mouse model in mind. The other 3 inputs are required.",
        "Coverage columns are named `coverage_<population>` — e.g. `coverage_World`, `coverage_Brazil`.",
        "Empty murine cells (when predict_murine WAS run) mean the peptide has no H-2 binder above threshold.",
        "This step is non-destructive: rerunning it just overwrites the join; upstream files are untouched.",
        "Murine `optimal`/`good`/`borderline` bands match `predict_murine`'s cutoffs (0.5 / 2.0 / 2.5).",
    ]

    def describe_outputs(self) -> dict:
        murine_dir = self.track_dir / "murine"
        return {
            murine_dir / get_step_filename("CURATE_MURINE", self.track_id):
                "Full per-★-peptide master table — CLUSTER_REPR columns + conservation + "
                "coverage_{population} + 4 murine columns. One row per ★ peptide.",
            murine_dir / get_step_filename("CURATE_MURINE_VIEW", self.track_id):
                "Slim view — peptide, alleles_united, final_score, conservation_label, "
                "and the three headline murine columns.",
            murine_dir / get_step_filename("CURATE_MURINE_AUDIT", self.track_id, ext='json'):
                "Run audit — conservation/coverage presence flags, total column count, "
                "murine label histogram, murine match counts.",
        }

    def run(self, input_data=None):
        run_start_time = time.time()

        human_star_df       = _load_star_human_table(self.track_dir, self.track_id)
        conservation_df     = _load_conservation(self.track_dir, self.track_id)
        coverage_pivoted_df = _load_coverage_pivoted(self.track_dir, self.track_id)
        murine_agg_df       = _load_murine_aggregate_if_present(self.track_dir, self.track_id)

        murine_present       = murine_agg_df is not None
        coverage_populations = [
            c.removeprefix('coverage_') for c in coverage_pivoted_df.columns
            if c != COLUMN_PEPTIDE
        ]

        console.print(Panel(
            Text.from_markup(
                f"[bold]Track:[/bold] {self.track_id}\n"
                f"[bold]★ peptides (human):[/bold] {len(human_star_df)}\n"
                f"[bold]Conservation:[/bold] attached "
                f"[dim]({len(conservation_df)} rows)[/dim]\n"
                f"[bold]Coverage:[/bold] attached "
                f"[dim]({len(coverage_populations)} populations)[/dim]\n"
                f"[bold]Murine:[/bold] "
                f"{f'attached ({len(murine_agg_df)} rows)' if murine_present else '[yellow]skipped — predict_murine not run[/yellow]'}"
            ),
            title="Curate murine",
            border_style="cyan",
            box=box.ROUNDED,
        ))

        master_df = _build_full_master_table(
            human_star_df       = human_star_df,
            conservation_df     = conservation_df,
            coverage_pivoted_df = coverage_pivoted_df,
            murine_agg_df       = murine_agg_df,
        )

        if murine_present:
            # Fill explicit empty defaults for peptides that have no murine row at all.
            master_df['murine_label'] = master_df['murine_label'].fillna('non_binder')
            master_df['murine_alleles_bound'] = master_df['murine_alleles_bound'].fillna('')
            master_df['num_murine_alleles_bound'] = (
                master_df['num_murine_alleles_bound'].fillna(0).astype(int)
            )

        murine_dir = self.track_dir / "murine"
        murine_dir.mkdir(parents=True, exist_ok=True)

        full_csv_path   = murine_dir / get_step_filename("CURATE_MURINE",       self.track_id)
        view_csv_path   = murine_dir / get_step_filename("CURATE_MURINE_VIEW",  self.track_id)
        audit_json_path = murine_dir / get_step_filename("CURATE_MURINE_AUDIT", self.track_id, ext='json')

        master_df.to_csv(full_csv_path, index=False)

        view_present_columns = [c for c in _VIEW_COLUMNS if c in master_df.columns]
        master_df[view_present_columns].to_csv(view_csv_path, index=False)

        if murine_present:
            label_histogram = _build_label_histogram(master_df)
            peptides_with_murine_match    = int((master_df['num_murine_alleles_bound'] > 0).sum())
            peptides_without_murine_match = len(master_df) - peptides_with_murine_match
        else:
            label_histogram               = None
            peptides_with_murine_match    = None
            peptides_without_murine_match = None

        audit_payload = {
            'timestamp':            datetime.datetime.now().isoformat(),
            'track_id':             self.track_id,
            'n_star_peptides':      len(master_df),
            'n_total_columns':      int(master_df.shape[1]),
            'coverage_populations': coverage_populations,
            'murine_present':       murine_present,
            'outputs': {
                'full_csv':   str(full_csv_path),
                'view_csv':   str(view_csv_path),
                'audit_json': str(audit_json_path),
            },
        }
        if murine_present:
            audit_payload['peptides_with_murine_match']    = peptides_with_murine_match
            audit_payload['peptides_without_murine_match'] = peptides_without_murine_match
            audit_payload['murine_label_histogram']        = label_histogram
        with open(audit_json_path, 'w', encoding='utf-8') as audit_file:
            json.dump(audit_payload, audit_file, indent=2, ensure_ascii=False)

        flush_stdin()

        elapsed_seconds = time.time() - run_start_time

        narrative_lines = [
            f"[bold]★ peptides joined:[/bold] {len(master_df):,}  "
            f"[dim]({master_df.shape[1]} columns total)[/dim]",
            f"  conservation columns: attached",
            f"  coverage columns:     attached ({len(coverage_populations)} populations)",
            f"  murine columns:       "
            f"{'attached' if murine_present else '[yellow]skipped — predict_murine not run[/yellow]'}",
        ]
        if murine_present:
            narrative_lines.append(
                f"[bold]Murine match:[/bold] {peptides_with_murine_match:,} with / "
                f"{peptides_without_murine_match:,} without"
            )
            narrative_lines.append(
                f"[bold]Murine label histogram:[/bold] "
                f"optimal={label_histogram['optimal']:,}  "
                f"good={label_histogram['good']:,}  "
                f"borderline={label_histogram['borderline']:,}  "
                f"non_binder={label_histogram['non_binder']:,}"
            )

        print_step_summary(
            step_title=f"Curate murine complete for {self.track_id}",
            elapsed_seconds=elapsed_seconds,
            narrative_lines=narrative_lines,
            output_files=[full_csv_path, view_csv_path, audit_json_path],
        )

        return {
            'full_csv':   str(full_csv_path),
            'view_csv':   str(view_csv_path),
            'audit_json': str(audit_json_path),
        }
