"""generate_report — Jinja2 HTML render.

Single render function that pulls the calculator template and inlines the
three datasets the JS expects on `window.*`: per-peptide records, project
metadata, and the per-locus coverage frequency map.
"""

from __future__ import annotations

import json
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape


_TEMPLATES_DIR = Path(__file__).parent / 'templates'


def render_report(
    output_html_path: Path,
    project_meta:     dict,
    epitopes:         list[dict],
    coverage_db:      dict,
    full_df=None,
) -> None:
    """Renders calculator.html.j2 to `output_html_path` with every dataset inlined.

    `full_df` is accepted for backwards compatibility with the step's call site
    but no longer used — the in-browser downloads (TSV / FASTA / matrix) derive
    everything they need from the per-peptide payload.
    """
    jinja_env = Environment(
        loader     = FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape = select_autoescape(disabled_extensions=('j2',)),
    )
    template = jinja_env.get_template('calculator.html.j2')

    rendered_html = template.render(
        project_meta      = project_meta,
        project_meta_json = json.dumps(project_meta, ensure_ascii=False),
        epitopes_json     = json.dumps(epitopes,     ensure_ascii=False),
        coverage_db_json  = json.dumps(coverage_db,  ensure_ascii=False),
    )

    output_html_path.parent.mkdir(parents=True, exist_ok=True)
    output_html_path.write_text(rendered_html, encoding='utf-8')
