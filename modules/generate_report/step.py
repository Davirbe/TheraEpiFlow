"""generate_report — global step orchestration.

Reads the three master-table artifacts written by `integrate_data` and the
IEDB allele-frequency pickle, builds the per-peptide JSON + project-meta +
coverage-DB payloads, and renders a single self-contained HTML calculator
to `data/output/REPORT_{project}.html`.

The rendered HTML is fully offline (JSZip vendored inline; Google Fonts is
the only remote dependency — it degrades gracefully to system fonts).
"""

from __future__ import annotations

import time
from typing import ClassVar

from rich import box
from rich.panel import Panel
from rich.text import Text

from modules.base_step import BaseGlobalStep
from utils.console import console
from utils.step_summary import print_step_summary

from .core import (
    assert_master_tables_exist,
    build_coverage_payload,
    build_epitopes_payload,
    build_project_meta,
    load_audit,
    load_full_dataframe,
)
from .io import render_report


class GenerateReportStep(BaseGlobalStep):
    step_name = 'generate_report'

    description = (
        "Renders an offline interactive HTML calculator from the master "
        "tables — researchers filter/sort epitopes, configure the vaccine "
        "construct (linker, TAG, adjuvant) and export a ZIP bundle with "
        "FASTA, stats and coverage."
    )
    long_description: ClassVar[str] = (
        "The pipeline's final deliverable. Consumes the three artifacts "
        "produced by `integrate_data` (VIEW CSV, FULL XLSX, AUDIT JSON) "
        "and the vendored IEDB allele-frequency pickle. Renders a single "
        "self-contained HTML file using Jinja2: every dataset is inlined "
        "as JSON, JSZip is embedded for the ZIP export, so the report "
        "works offline on any browser.\n\n"
        "Inside the HTML, the user filters by organism/protein/conservation, "
        "sorts every numeric column, selects a set of epitopes, then opens "
        "the 'Finalize construct' modal — which renders an organism × "
        "protein heatmap, cumulative population coverage (computed in JS "
        "via the standard `1 - Π(1 - f_i)` formula over the union of "
        "alleles), and construct stats (length, MW, mean conservation). A "
        "single 'Download ZIP' button packages the FASTA, the full-detail "
        "CSV of selected peptides, a stats JSON, the heatmap PNG and a "
        "human-readable summary."
    )
    methodology: ClassVar[str] = (
        "1. Reads `MASTER_TABLE_FULL_{project}.xlsx` (46 raw columns) and "
        "`MASTER_TABLE_AUDIT_{project}.json` (to know which optional "
        "columns the user opted into during integrate_data).\n"
        "2. Builds the per-peptide JSON list. Each entry carries the "
        "organism/protein/peptide identity, the consensus best percentile, "
        "the bound HLA allele list (parsed into an array), the two "
        "conservation fractions + label, the per-population coverage, and "
        "the four murine fields.\n"
        "3. Reduces the IEDB pickle to a small "
        "`{population: {allele: frequency}}` dict, restricted to alleles "
        "actually present in the project (drops 3 MB to ~5 KB).\n"
        "4. Renders `calculator.html.j2` via Jinja2, inlining every dataset "
        "and the vendored JSZip script. Output: `REPORT_{project}.html`."
    )
    references: ClassVar[list] = []
    data_format: ClassVar[str] = (
        "Inputs are picked up automatically — the three "
        "[bold]MASTER_TABLE_*_{project}[/bold] files under `data/output/`. "
        "If any is missing, the step raises a clear FileNotFoundError "
        "pointing back to `integrate_data`. No prompts."
    )
    outputs_overview: ClassVar[str] = (
        "[bold]REPORT_{project}.html[/bold] — a single self-contained HTML "
        "file (typically 300–500 KB depending on n peptides). Opens "
        "directly in any browser, works offline. Inside the report the "
        "user assembles a construct and downloads a ZIP bundle (FASTA + "
        "selected_epitopes.csv + construction_stats.json + "
        "coverage_heatmap.png + selection_summary.txt) — that ZIP is the "
        "scientific deliverable the researcher hands off."
    )
    tips: ClassVar[list] = [
        "[bold yellow]This step is the pipeline's final deliverable.[/bold yellow] The "
        "HTML is a snapshot in time — every dataset (epitopes, allele frequencies, "
        "selected populations) is baked into the file at render time. To add a new "
        "population (e.g. 'Brazil'), edit `project_config.coverage_populations`, then "
        "re-run population_coverage → curate_murine → integrate_data → generate_report.",
        "Run `integrate_data` first — generate_report depends on its output.",
        "The report is fully offline once rendered; copy the single HTML to share with collaborators.",
        "Click any column header to sort; click again to reverse.",
        "Hover over Conservation / HLA / Murine cells to see the underlying numbers and lists.",
        "Use the Finalize construct button only after picking at least one epitope.",
    ]

    def describe_outputs(self) -> dict:
        return {
            self.output_dir / f'REPORT_{self.project_name}.html':
                "Self-contained interactive HTML calculator — opens in any browser, "
                "works offline. The pipeline's final deliverable. Press [o] to open in your browser.",
        }

    def run(self, input_data=None):
        run_start_time = time.time()

        paths = assert_master_tables_exist(self.output_dir, self.project_name)

        full_df = load_full_dataframe(paths['full_xlsx'])
        audit   = load_audit(paths['audit_json'])

        selected_view_columns = audit.get('view_columns_selected', [])

        project_meta = build_project_meta(
            project_name   = self.project_name,
            project_config = self.project_config,
            full_df        = full_df,
            audit          = audit,
        )

        epitopes    = build_epitopes_payload(full_df, selected_view_columns)
        coverage_db = build_coverage_payload(full_df, project_meta['populations'])

        report_path = self.output_dir / f'REPORT_{self.project_name}.html'

        console.print(Panel(
            Text.from_markup(
                f"[bold]Project:[/bold] {self.project_name}\n"
                f"[bold]★ peptides:[/bold] {project_meta['n_peptides_total']}\n"
                f"[bold]Organisms × proteins:[/bold] "
                f"{project_meta['n_organisms']} × {project_meta['n_proteins']}\n"
                f"[bold]Populations:[/bold] {', '.join(project_meta['populations']) or '—'}\n"
                f"[bold]COVERAGE_DB alleles inlined:[/bold] "
                f"{sum(len(v) for v in coverage_db.values())} "
                f"[dim](across {len(coverage_db)} population(s))[/dim]"
            ),
            title="Generate report",
            border_style="cyan",
            box=box.ROUNDED,
        ))

        render_report(
            output_html_path = report_path,
            project_meta     = project_meta,
            epitopes         = epitopes,
            coverage_db      = coverage_db,
            full_df          = full_df,
        )

        elapsed_seconds = time.time() - run_start_time

        narrative_lines = [
            f"[bold]Peptides exposed:[/bold] {len(epitopes):,}",
            f"[bold]Organisms in report:[/bold] {project_meta['n_organisms']}",
            f"[bold]Proteins in report:[/bold] {project_meta['n_proteins']}",
            f"[bold]Populations with coverage:[/bold] {len(coverage_db)}",
            f"[bold]Alleles inlined for coverage math:[/bold] "
            f"{sum(len(v) for v in coverage_db.values())}",
            f"[bold]HTML size:[/bold] "
            f"{report_path.stat().st_size / 1024:.1f} KB "
            f"[dim](includes JSZip + inline data)[/dim]",
        ]

        print_step_summary(
            step_title      = f"Interactive report built for {self.project_name}",
            elapsed_seconds = elapsed_seconds,
            narrative_lines = narrative_lines,
            output_files    = [report_path],
        )

        return {
            'report_html':  str(report_path),
            'n_peptides':   len(epitopes),
            'n_organisms':  project_meta['n_organisms'],
            'n_proteins':   project_meta['n_proteins'],
            'n_populations': len(coverage_db),
        }
