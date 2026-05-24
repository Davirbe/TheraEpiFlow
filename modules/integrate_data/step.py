"""integrate_data — global step orchestration.

Stacks every track's CURATE_MURINE table into a project-wide master table,
asks the user (once) which columns the VIEW projection should contain, and
writes FULL/VIEW XLSX + a VIEW CSV + an audit JSON to data/output/.

The customization choice is persisted to
`project_config.step_overrides.integrate_data.view_columns` so reruns reuse
it without re-prompting; pass --reconfigure to re-open the prompt.
"""

from __future__ import annotations

import datetime
import json
import time
from typing import ClassVar

from rich import box
from rich.panel import Panel
from rich.text import Text

from modules.base_step import BaseGlobalStep
from utils.console import console
from utils.project_manager import save_project_config
from utils.step_summary import print_step_summary

from .core import (
    attach_anchor_aggregates,
    build_column_catalog,
    coverage_populations_in,
    default_column_selection,
    project_view,
    stack_tracks,
    view_headers,
)
from .io import write_full_xlsx, write_view_csv, write_view_xlsx
from .prompts import prompt_view_customization


_ANCHOR_COLUMNS = {'n_mutations_safe', 'n_mutations_risky'}


def _filename(stem: str, project_name: str, ext: str) -> str:
    """Project-scoped output filename. Mirrors get_step_filename's style but
    keyed on project_name since global outputs aren't per-track."""
    return f"{stem}_{project_name}.{ext}"


class IntegrateDataStep(BaseGlobalStep):
    step_name = 'integrate_data'

    description = (
        "Stacks every track's curated master table into one project-wide table, "
        "then projects a user-configurable VIEW slim enough to feed the report's "
        "calculator. Identity collapses from track_id to organism+protein."
    )
    long_description: ClassVar[str] = (
        "By the time `curate_murine` finishes, every track has a single CSV that "
        "carries the full per-★-peptide evidence (binding, conservation, "
        "coverage, murine). This step joins those into one project-wide table "
        "and asks one question: which columns should the report's calculator "
        "see?\n\n"
        "Two files are always written: a FULL XLSX with every column for "
        "auditing/research, and a VIEW XLSX/CSV with only the columns the "
        "researcher will use to pick epitopes for the construct. The VIEW is "
        "the input for `generate_report`."
    )
    methodology: ClassVar[str] = (
        "1. Reads `CURATE_MURINE_{track_id}.csv` from each track in the "
        "project. Tracks without that file are recorded in the audit JSON "
        "and skipped (no error).\n"
        "2. Prepends `track_id`, `organism`, `protein` (the latter two pulled "
        "from `project_config` so labels are canonical, even for organisms "
        "with dashes in their name).\n"
        "3. Concatenates all tracks (outer-style) into a single dataframe.\n"
        "4. Builds the column catalog: a curated list grouped by topic "
        "(identity / binding / conservation / coverage / murine), with each "
        "entry marked default-on or optional.\n"
        "5. On first run, asks the user to accept the default selection or "
        "open a checkbox UI to toggle individual columns. The choice is "
        "saved under `project_config.step_overrides.integrate_data.view_columns` "
        "and honored on subsequent runs (use `--reconfigure` to re-ask).\n"
        "6. If anchor-mutation columns are selected, aggregates per peptide "
        "from `CONSERVATION_MUTATIONS_*.xlsx` and joins them in.\n"
        "7. Writes FULL XLSX (all columns) + VIEW XLSX/CSV (selected columns "
        "with the canonical palette) + an audit JSON."
    )
    references: ClassVar[list] = []
    data_format: ClassVar[str] = (
        "Inputs are picked up automatically — one "
        "[bold]CURATE_MURINE_{track_id}.csv[/bold] per track under "
        "`data/intermediate/{track_id}/murine/`. Tracks missing that file are "
        "skipped with a warning, not an error.\n\n"
        "VIEW headers are in English and renamed for readability (e.g. "
        "`best_combined_percentile` → 'Best percentile (consensus)'). The "
        "internal column names are preserved in the FULL XLSX and the CSV "
        "sidecar."
    )
    outputs_overview: ClassVar[str] = (
        "[bold]MASTER_TABLE_FULL_{project}.xlsx[/bold] — every column from "
        "every track, plus `track_id`/`organism`/`protein` for traceability. "
        "No cell coloring — this is the audit-grade dump.\n"
        "[bold]MASTER_TABLE_VIEW_{project}.xlsx[/bold] — selected columns "
        "with display headers and the canonical palette (conservation/murine "
        "labels colored per band, coverage cells colored by threshold, "
        "percentile cols orange, HLA-count cols pink, tooltip-only cols "
        "rendered narrow). No `track_id` — identity is organism + protein.\n"
        "[bold]MASTER_TABLE_VIEW_{project}.csv[/bold] — same VIEW with "
        "display headers, in CSV.\n"
        "[bold]MASTER_TABLE_AUDIT_{project}.json[/bold] — track counts, "
        "skipped tracks, chosen VIEW columns, generation timestamp."
    )
    tips: ClassVar[list] = [
        "Run `curate_murine` on every track first — tracks without it are skipped.",
        "Pass `--reconfigure` to re-open the column-customization prompt.",
        "The VIEW is the input to `generate_report`; choose columns to match what users will need in the calculator.",
        "FULL keeps every column so you can always re-derive a custom VIEW later.",
    ]

    # ── Reconfigure plumbing ──────────────────────────────────────────────────
    #
    # BaseGlobalStep doesn't expose `is_rerun` like BaseTrackStep does, so we
    # override execute() just to capture the reconfigure flag and stash it for
    # run() to inspect.

    _force_reprompt: bool = False

    def execute(self, input_data=None, force_rerun: bool = False, reconfigure: bool = False) -> dict:
        self._force_reprompt = bool(force_rerun or reconfigure)
        return super().execute(input_data, force_rerun=force_rerun, reconfigure=reconfigure)

    # ── Main entry ────────────────────────────────────────────────────────────

    def run(self, input_data=None):
        run_start_time = time.time()

        full_df, present_tracks, skipped_tracks = stack_tracks(
            self.intermediate_dir, self.project_config,
        )

        if full_df.empty:
            raise FileNotFoundError(
                "No CURATE_MURINE_*.csv found in any track. "
                "Run 'curate_murine' before 'integrate_data'."
            )

        coverage_populations = coverage_populations_in(full_df)
        catalog              = build_column_catalog(coverage_populations)

        selected_columns = self._resolve_selected_columns(catalog)

        if any(name in selected_columns for name in _ANCHOR_COLUMNS):
            full_df = attach_anchor_aggregates(full_df, self.intermediate_dir, present_tracks)

        view_df = project_view(full_df, selected_columns, catalog)
        headers = view_headers(list(view_df.columns), catalog)

        self.output_dir.mkdir(parents=True, exist_ok=True)
        full_xlsx_path  = self.output_dir / _filename('MASTER_TABLE_FULL',  self.project_name, 'xlsx')
        view_xlsx_path  = self.output_dir / _filename('MASTER_TABLE_VIEW',  self.project_name, 'xlsx')
        view_csv_path   = self.output_dir / _filename('MASTER_TABLE_VIEW',  self.project_name, 'csv')
        audit_json_path = self.output_dir / _filename('MASTER_TABLE_AUDIT', self.project_name, 'json')

        console.print(Panel(
            Text.from_markup(
                f"[bold]Project:[/bold] {self.project_name}\n"
                f"[bold]Tracks integrated:[/bold] {len(present_tracks)} "
                f"[dim]({', '.join(present_tracks) or '—'})[/dim]\n"
                f"[bold]Tracks skipped:[/bold] "
                f"{len(skipped_tracks)}"
                + (f" [dim]({', '.join(skipped_tracks)})[/dim]" if skipped_tracks else '')
                + f"\n"
                f"[bold]Coverage populations:[/bold] "
                f"{', '.join(coverage_populations) or '—'}\n"
                f"[bold]VIEW columns:[/bold] {len(view_df.columns)} "
                f"[dim](FULL has {full_df.shape[1]})[/dim]"
            ),
            title="Integrate data",
            border_style="cyan",
            box=box.ROUNDED,
        ))

        write_full_xlsx(full_df, full_xlsx_path)
        write_view_xlsx(view_df, view_xlsx_path, headers=headers)
        write_view_csv (view_df, view_csv_path,  headers=headers)

        audit_payload = {
            'timestamp':              datetime.datetime.now().isoformat(),
            'project_name':           self.project_name,
            'tracks_integrated':      present_tracks,
            'tracks_skipped':         skipped_tracks,
            'coverage_populations':   coverage_populations,
            'n_rows_full':            int(full_df.shape[0]),
            'n_columns_full':         int(full_df.shape[1]),
            'n_columns_view':         int(view_df.shape[1]),
            'view_columns_selected':  list(view_df.columns),
            'view_columns_persisted': self._persisted_setting_repr(),
            'outputs': {
                'full_xlsx':  str(full_xlsx_path),
                'view_xlsx':  str(view_xlsx_path),
                'view_csv':   str(view_csv_path),
                'audit_json': str(audit_json_path),
            },
        }
        with open(audit_json_path, 'w', encoding='utf-8') as audit_file:
            json.dump(audit_payload, audit_file, indent=2, ensure_ascii=False)

        elapsed_seconds = time.time() - run_start_time

        narrative_lines = [
            f"[bold]Tracks integrated:[/bold] {len(present_tracks):,}  "
            f"[dim](skipped {len(skipped_tracks)})[/dim]",
            f"[bold]★ peptides total:[/bold] {full_df.shape[0]:,}",
            f"[bold]FULL columns:[/bold] {full_df.shape[1]}",
            f"[bold]VIEW columns:[/bold] {view_df.shape[1]} "
            f"[dim]({'default' if self._persisted_setting_repr() == 'default' else 'custom'})[/dim]",
            f"[bold]Coverage populations:[/bold] "
            f"{', '.join(coverage_populations) if coverage_populations else 'none'}",
        ]
        if skipped_tracks:
            narrative_lines.append(
                f"[yellow]Skipped (no CURATE_MURINE_*.csv):[/yellow] {', '.join(skipped_tracks)}"
            )

        print_step_summary(
            step_title      = f"Master tables built for {self.project_name}",
            elapsed_seconds = elapsed_seconds,
            narrative_lines = narrative_lines,
            output_files    = [full_xlsx_path, view_xlsx_path, view_csv_path, audit_json_path],
        )

        return {
            'full_xlsx':  str(full_xlsx_path),
            'view_xlsx':  str(view_xlsx_path),
            'view_csv':   str(view_csv_path),
            'audit_json': str(audit_json_path),
            'n_tracks':   len(present_tracks),
            'n_rows':     int(full_df.shape[0]),
        }

    # ── Column-selection helpers ──────────────────────────────────────────────

    def _resolve_selected_columns(self, catalog: list[dict]) -> list[str]:
        """Returns the final VIEW column list. Reads persisted choice unless
        the user explicitly asked to reconfigure; otherwise prompts."""
        persisted = (
            self.project_config
                .get('step_overrides', {})
                .get('integrate_data', {})
                .get('view_columns')
        )

        if persisted is not None and not self._force_reprompt:
            if persisted == 'default':
                return default_column_selection(catalog)
            if isinstance(persisted, list):
                # Filter out anything the catalog no longer knows about, but
                # always include required identity columns.
                valid_names = {entry['name'] for entry in catalog}
                required_names = [entry['name'] for entry in catalog if entry['required']]
                kept = [c for c in persisted if c in valid_names]
                for required_name in required_names:
                    if required_name not in kept:
                        kept.append(required_name)
                return kept

        chosen = prompt_view_customization(catalog)
        self._persist_selection(chosen, catalog)
        return chosen

    def _persist_selection(self, chosen: list[str], catalog: list[dict]) -> None:
        """Writes the chosen list (or 'default' sentinel) into project_config."""
        is_default = set(chosen) == set(default_column_selection(catalog))
        overrides = self.project_config.setdefault('step_overrides', {})
        step_overrides = overrides.setdefault('integrate_data', {})
        step_overrides['view_columns'] = 'default' if is_default else chosen
        save_project_config(self.project_name, self.project_config)

    def _persisted_setting_repr(self) -> object:
        """Returns the value currently saved in project_config for audit purposes."""
        return (
            self.project_config
                .get('step_overrides', {})
                .get('integrate_data', {})
                .get('view_columns', 'unknown')
        )
