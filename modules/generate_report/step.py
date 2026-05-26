"""generate_report — global step orchestration.

Reads the three master-table artifacts written by `integrate_data` and the
IEDB allele-frequency pickle, builds the per-peptide JSON + project-meta +
coverage-DB payloads, and renders a single self-contained HTML calculator
to `data/output/REPORT_{project}.html`.

The rendered HTML is fully offline. Google Fonts is the only remote
dependency and degrades gracefully to system fonts. Downloads (TSV, FASTA,
matrix) are produced inline in JS via Blob + URL.createObjectURL — no
external JSZip dependency.
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
        "tables — researchers filter/sort epitopes, build a vaccine "
        "construct (linker + 6×His tag) and download TSV / FASTA / "
        "matrix inline."
    )
    long_description: ClassVar[str] = (
        "The pipeline's final deliverable. Consumes the three artifacts "
        "produced by `integrate_data` (VIEW CSV, FULL XLSX, AUDIT JSON) "
        "and the vendored IEDB allele-frequency pickle. Renders a single "
        "self-contained HTML file using Jinja2: every dataset is inlined "
        "as JSON so the report works offline on any browser.\n\n"
        "Inside the HTML the user filters by Protein / Organism / Flag "
        "buttons, sorts every numeric column, and picks epitopes to "
        "assemble a construct (linker + 6×His). The right pane shows the "
        "selected × populations coverage heatmap (with a Cumulative row "
        "computed via the diploid IEDB model, locus-aware) plus the "
        "organism × protein distribution heatmap. Three buttons in the "
        "footer download the selection as TSV, the construct as FASTA, "
        "and the distribution as CSV — no modal, no JSZip dependency."
    )
    methodology: ClassVar[str] = (
        "1. Reads `MASTER_TABLE_FULL_{project}.xlsx` and "
        "`MASTER_TABLE_AUDIT_{project}.json` (the second tells us which "
        "optional columns the user opted into during `integrate_data`).\n"
        "2. Builds the per-peptide JSON list. Each entry carries the "
        "organism / protein / peptide identity, the consensus best "
        "percentile, the bound HLA allele list (parsed into an array), "
        "the two conservation fractions + label, the per-population "
        "coverage, and the four murine fields.\n"
        "3. Reduces the IEDB pickle to a per-locus map "
        "`{population: {locus: {allele: frequency}}}`, restricted to "
        "alleles actually present in the project (drops 3 MB to a few "
        "KB). The locus grouping is required by the diploid coverage "
        "formula that runs in JS — it mirrors "
        "`modules/population_coverage/core.py:_compute_epitope_coverage`.\n"
        "4. Renders `calculator.html.j2` via Jinja2, inlining every "
        "dataset. Output: `REPORT_{project}.html`."
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
        "file (typically 100–300 KB). Opens directly in any browser, works "
        "offline. Inside the report the user assembles a construct and "
        "downloads three artifacts inline: a TSV of selected epitopes, a "
        "FASTA of the assembled construct sequence, and a CSV of the "
        "organism × protein distribution matrix."
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
        "Hover over MP / HLA / MHC 🐭 cells to see the per-allele percentile breakdown.",
        "Coverage math in JS uses the diploid IEDB model (locus-aware) — values match the population_coverage step exactly.",
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

        n_alleles_inlined = sum(
            len(allele_map)
            for per_locus in coverage_db.values()
            for allele_map in per_locus.values()
        )
        console.print(Panel(
            Text.from_markup(
                f"[bold]Project:[/bold] {self.project_name}\n"
                f"[bold]★ peptides:[/bold] {project_meta['n_peptides_total']}\n"
                f"[bold]Organisms × proteins:[/bold] "
                f"{project_meta['n_organisms']} × {project_meta['n_proteins']}\n"
                f"[bold]Populations:[/bold] {', '.join(project_meta['populations']) or '—'}\n"
                f"[bold]COVERAGE_DB alleles inlined:[/bold] "
                f"{n_alleles_inlined} "
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
            f"[bold]Alleles inlined for coverage math:[/bold] {n_alleles_inlined}",
            f"[bold]HTML size:[/bold] "
            f"{report_path.stat().st_size / 1024:.1f} KB "
            f"[dim](self-contained, no external dependencies)[/dim]",
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
