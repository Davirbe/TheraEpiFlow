"""generate_report — Jinja2 HTML render.

Single render function that pulls the calculator template + the vendored
JSZip (via the `{% include %}` directive) and inlines every dataset that
the JS expects on `window.*`.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape


_TEMPLATES_DIR = Path(__file__).parent / 'templates'


def _coalesce(value):
    """JSON-safe scalar coercion (NaN → None, numpy scalar → Python)."""
    if pd.isna(value):
        return None
    if hasattr(value, 'item'):
        return value.item()
    return value


def _full_rows_for_export(full_df: pd.DataFrame) -> list[dict]:
    """Materializes the FULL dataframe as a list of dicts the JS can filter
    by `_id` (1-based, aligned with the epitopes payload). Used to build
    `selected_epitopes.csv` inside the ZIP, so the user gets EVERY column
    (including BEST_REPRESENTATIVE, ToxinPred3, cluster, etc.) for just the
    peptides they selected."""
    records: list[dict] = []
    for one_based_index, (_, row) in enumerate(full_df.iterrows(), start=1):
        record = {'_id': one_based_index}
        for column_name in full_df.columns:
            record[column_name] = _coalesce(row[column_name])
        records.append(record)
    return records


def render_report(
    output_html_path: Path,
    project_meta:     dict,
    epitopes:         list[dict],
    coverage_db:      dict,
    full_df:          pd.DataFrame,
) -> None:
    """Renders calculator.html.j2 to `output_html_path` with every dataset inlined."""
    jinja_env = Environment(
        loader     = FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape = select_autoescape(disabled_extensions=('j2',)),
    )
    template = jinja_env.get_template('calculator.html.j2')

    rendered_html = template.render(
        project_meta      = project_meta,
        project_meta_json = json.dumps(project_meta,    ensure_ascii=False),
        epitopes_json     = json.dumps(epitopes,        ensure_ascii=False),
        coverage_db_json  = json.dumps(coverage_db,     ensure_ascii=False),
        full_rows_json    = json.dumps(_full_rows_for_export(full_df), ensure_ascii=False),
    )

    output_html_path.parent.mkdir(parents=True, exist_ok=True)
    output_html_path.write_text(rendered_html, encoding='utf-8')
