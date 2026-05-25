"""generate_report — VIEW/FULL/AUDIT readers + per-peptide JSON builder.

Pure data prep. No Rich, no Jinja, no openpyxl writes — just pandas/json/csv
reads and dict construction. Outputs the three JSON blobs the calculator
template needs:

  - epitopes_json  — one object per ★ peptide for the table/calculator
  - project_meta_json — header/filter constants (organisms, proteins, n_*, etc.)
  - coverage_db_json  — {population: {allele: freq}} (via coverage_db.py)
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path

import pandas as pd

from .coverage_db import build_coverage_db


# Internal master-table column names (NOT the English display headers in the CSV).
# We read the XLSX so we get raw column names instead of headers.
_FULL_PEPTIDE_COL  = 'peptide'
_FULL_ORGANISM_COL = 'organism'
_FULL_PROTEIN_COL  = 'protein'
_FULL_ALLELES_COL  = 'alleles_united'

# Internal → JSON key map for the per-peptide payload.
_BASE_JSON_KEYS = {
    'peptide':                   'peptide',
    'organism':                  'organism',
    'protein':                   'protein',
    'best_combined_percentile':  'best_pct',
    'num_alleles_united':        'n_hla',
    'alleles_united':            'hla_list',   # split to array
    'pct_identity_100':          'cons_100',
    'pct_identity_threshold':    'cons_thr',
    'conservation_label':        'cons_label',
    'murine_label':              'murine_label',
    'murine_best_percentile':    'murine_pct',
    'num_murine_alleles_bound':  'murine_n',
}


def _master_table_paths(output_dir: Path, project_name: str) -> dict[str, Path]:
    """Returns the four expected MASTER_TABLE_* paths for this project."""
    return {
        'view_csv':   output_dir / f'MASTER_TABLE_VIEW_{project_name}.csv',
        'full_xlsx':  output_dir / f'MASTER_TABLE_FULL_{project_name}.xlsx',
        'audit_json': output_dir / f'MASTER_TABLE_AUDIT_{project_name}.json',
    }


def assert_master_tables_exist(output_dir: Path, project_name: str) -> dict[str, Path]:
    """Raises FileNotFoundError with a clear next-step message if any input is missing."""
    paths = _master_table_paths(output_dir, project_name)
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "generate_report needs the master tables produced by integrate_data.\n"
            f"Missing: {', '.join(missing)}\n"
            f"Run: python main.py --project {project_name} --step integrate_data"
        )
    return paths


def load_full_dataframe(full_xlsx: Path) -> pd.DataFrame:
    """Reads MASTER_TABLE_FULL_*.xlsx (raw internal column names)."""
    return pd.read_excel(full_xlsx, sheet_name='Master Full')


def load_audit(audit_json: Path) -> dict:
    """Reads MASTER_TABLE_AUDIT_*.json so we know which optional columns the user opted into."""
    with open(audit_json, encoding='utf-8') as audit_file:
        return json.load(audit_file)


def _split_alleles(allele_blob: object) -> list[str]:
    """Turns the semicolon-separated string into a list; tolerates NaN/empty."""
    if pd.isna(allele_blob):
        return []
    return [a.strip() for a in str(allele_blob).split(';') if a.strip()]


def _split_floats(blob: object) -> list[float | None]:
    """Turns a semicolon-separated number blob into a list of floats; tolerates
    NaN/empty (whole column missing → []; individual non-numeric token → None)."""
    if pd.isna(blob):
        return []
    out: list[float | None] = []
    for token in str(blob).split(';'):
        token = token.strip()
        if not token:
            continue
        try:
            out.append(float(token))
        except ValueError:
            out.append(None)
    return out


def _coalesce_scalar(value):
    """JSON-safe scalar: NaN → None; numpy scalars → Python."""
    if pd.isna(value):
        return None
    if hasattr(value, 'item'):
        return value.item()
    return value


def build_epitopes_payload(
    full_df: pd.DataFrame,
    selected_view_columns: list[str],
) -> list[dict]:
    """Turns each FULL row into the per-peptide JSON dict the template consumes.

    `selected_view_columns` are the internal column names the user opted into
    via `integrate_data`'s prompt. Anything outside the BASE_JSON_KEYS map is
    surfaced under the `optional` key so the HTML can render extra columns
    automatically without code changes.
    """
    payload: list[dict] = []
    base_internal_cols  = set(_BASE_JSON_KEYS.keys())

    # Coverage columns are dynamic per project: coverage_World, coverage_Brazil, ...
    coverage_columns = [c for c in full_df.columns if c.startswith('coverage_')]

    optional_columns = [
        c for c in selected_view_columns
        if c not in base_internal_cols and not c.startswith('coverage_')
        and c in full_df.columns
    ]

    for one_based_index, (_, row) in enumerate(full_df.iterrows(), start=1):
        epitope_record: dict = {'id': one_based_index}

        for internal_col, json_key in _BASE_JSON_KEYS.items():
            if internal_col not in full_df.columns:
                continue
            value = row[internal_col]
            if json_key == 'hla_list':
                epitope_record[json_key] = _split_alleles(value)
            else:
                epitope_record[json_key] = _coalesce_scalar(value)

        epitope_record['coverage'] = {
            col.removeprefix('coverage_'): _coalesce_scalar(row[col])
            for col in coverage_columns
        }

        # Per-method per-allele percentiles for the Best %ile tooltip (Round 2).
        # Arrays are positionally aligned with `hla_list` (same order in
        # alleles_united / netmhcpan_el_percentiles_all / mhcflurry_presentation_percentiles_all).
        if 'netmhcpan_el_percentiles_all' in full_df.columns:
            epitope_record['net_pct_per_allele'] = _split_floats(row['netmhcpan_el_percentiles_all'])
        if 'mhcflurry_presentation_percentiles_all' in full_df.columns:
            epitope_record['flurry_pct_per_allele'] = _split_floats(row['mhcflurry_presentation_percentiles_all'])

        if optional_columns:
            epitope_record['optional'] = {
                c: _coalesce_scalar(row[c]) for c in optional_columns
            }

        payload.append(epitope_record)

    return payload


def build_project_meta(
    project_name: str,
    project_config: dict,
    full_df: pd.DataFrame,
    audit: dict,
) -> dict:
    """Header/filter constants the JS reads to populate the UI dynamically."""
    organisms = sorted(full_df[_FULL_ORGANISM_COL].dropna().unique().tolist())
    proteins  = sorted(full_df[_FULL_PROTEIN_COL].dropna().unique().tolist())
    populations = audit.get('coverage_populations') or list(
        project_config.get('coverage_populations', [])
    )

    conservation_threshold = project_config.get('conservation_threshold')

    return {
        'project_name':            project_name,
        'project_description':     project_config.get('description', ''),
        'generated_at':            datetime.datetime.now().isoformat(timespec='seconds'),
        'organisms':               organisms,
        'proteins':                proteins,
        'populations':             populations,
        'n_organisms':             len(organisms),
        'n_proteins':              len(proteins),
        'n_peptides_total':        int(len(full_df)),
        'conservation_threshold':  conservation_threshold,
        'optional_columns':        [
            c for c in audit.get('view_columns_selected', [])
            if c not in set(_BASE_JSON_KEYS.keys()) and not c.startswith('coverage_')
        ],
    }


def collect_alleles_in_project(full_df: pd.DataFrame) -> set[str]:
    """Returns the union of every allele appearing across all peptides."""
    alleles_in_project: set[str] = set()
    for blob in full_df[_FULL_ALLELES_COL].dropna():
        alleles_in_project.update(_split_alleles(blob))
    return alleles_in_project


def build_coverage_payload(full_df: pd.DataFrame, populations: list[str]) -> dict:
    """{population: {allele: freq}} — wrapper for coverage_db.build_coverage_db."""
    alleles_in_project = collect_alleles_in_project(full_df)
    return build_coverage_db(populations, alleles_in_project)
