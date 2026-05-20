"""
Domain logic for predict_murine: ★ peptide loading, synthetic SeqRecords, the
binder-tier label, the per-tool collapse and the per-peptide aggregation. Pure
data functions — no Rich, no prompts, no tool invocation.
"""

from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import config
from utils.naming import (
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_PEPTIDE,
    STAR_MARKER,
    get_step_filename,
)

_LABEL_OPTIMAL     = 'optimal'
_LABEL_GOOD        = 'good'
_LABEL_BORDERLINE  = 'borderline'
_LABEL_NON_BINDER  = 'non_binder'

_BOUND_LABELS = {_LABEL_OPTIMAL, _LABEL_GOOD, _LABEL_BORDERLINE}
# ── Star peptide loader ───────────────────────────────────────────────────────

def _load_star_peptides(track_dir: Path, track_id: str) -> pd.DataFrame:
    """Returns a DataFrame with the ★ representatives. Raises if input missing."""
    representatives_csv = track_dir / "clusters" / get_step_filename("CLUSTER_REPR", track_id)
    if not representatives_csv.exists():
        raise FileNotFoundError(
            f"select_representatives output not found: {representatives_csv}\n"
            "Run 'select_representatives' before 'predict_murine'."
        )

    representatives_df = pd.read_csv(representatives_csv)
    if COLUMN_BEST_REPRESENTATIVE not in representatives_df.columns:
        raise ValueError(
            f"Column '{COLUMN_BEST_REPRESENTATIVE}' not found in "
            f"{representatives_csv.name}."
        )
    if COLUMN_PEPTIDE not in representatives_df.columns:
        raise ValueError(
            f"Column '{COLUMN_PEPTIDE}' not found in {representatives_csv.name}."
        )

    star_representatives_df = representatives_df[
        representatives_df[COLUMN_BEST_REPRESENTATIVE] == STAR_MARKER
    ].copy()
    if star_representatives_df.empty:
        raise ValueError(
            "No ★ representatives found in select_representatives output. "
            "Re-run select_representatives or check that the marker column is correct."
        )

    return star_representatives_df.reset_index(drop=True)


# ── Synthetic SeqRecord builder ───────────────────────────────────────────────

def _build_synthetic_seqrecords(peptide_sequences: list[str]) -> list[SeqRecord]:
    """One SeqRecord per peptide — len(seq) == len(peptide), so no k-mer regeneration."""
    return [
        SeqRecord(seq=Seq(peptide_string), id=peptide_string, description="")
        for peptide_string in peptide_sequences
    ]


# ── Tier label ────────────────────────────────────────────────────────────────

def _assign_binder_label(percentile_value) -> str:
    if percentile_value is None or pd.isna(percentile_value):
        return _LABEL_NON_BINDER
    if percentile_value <= config.MURINE_OPTIMAL_BINDER_RANK_MAX:
        return _LABEL_OPTIMAL
    if percentile_value <= config.MURINE_STRONG_BINDER_EL_RANK:
        return _LABEL_GOOD
    if percentile_value <= config.MURINE_BORDERLINE_BINDER_RANK_MAX:
        return _LABEL_BORDERLINE
    return _LABEL_NON_BINDER


# ── Per-(peptide, allele) collapse between tools ──────────────────────────────

def _combine_per_tool_long_tables(
    netmhcpan_raw_df: pd.DataFrame,
    mhcflurry_raw_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (best_per_pair_df, full_long_df).
    full_long_df keeps every (peptide, allele, tool) for audit; best_per_pair_df
    collapses to one row per (peptide, allele) with the smaller percentile."""
    netmhcpan_long_df = netmhcpan_raw_df.rename(columns={
        'netmhcpan_el_percentile': 'percentile',
    })[['peptide', 'allele', 'percentile']].copy()
    netmhcpan_long_df['tool'] = 'netmhcpan_el'

    mhcflurry_long_df = mhcflurry_raw_df.rename(columns={
        'mhcflurry_presentation_percentile': 'percentile',
    })[['peptide', 'allele', 'percentile']].copy()
    mhcflurry_long_df['tool'] = 'mhcflurry'

    full_long_df = pd.concat([netmhcpan_long_df, mhcflurry_long_df], ignore_index=True)

    best_per_pair_df = (
        full_long_df
        .sort_values('percentile', ascending=True)
        .groupby(['peptide', 'allele'], as_index=False)
        .first()
        .drop(columns=['tool'])
    )
    return best_per_pair_df, full_long_df


# ── Aggregation per peptide ───────────────────────────────────────────────────

def _aggregate_per_peptide(
    best_per_pair_df: pd.DataFrame,
    all_star_peptides: list[str],
) -> pd.DataFrame:
    """One row per ★ peptide with best_percentile_* + murine_alleles_bound."""
    aggregated_rows: list[dict] = []

    for peptide_string in all_star_peptides:
        peptide_rows = best_per_pair_df[best_per_pair_df['peptide'] == peptide_string]
        bound_rows = peptide_rows[peptide_rows['binder_label'].isin(_BOUND_LABELS)]
        bound_rows = bound_rows.sort_values('percentile', ascending=True)

        if bound_rows.empty:
            aggregated_rows.append({
                'peptide':                  peptide_string,
                'best_percentile_label':    _LABEL_NON_BINDER,
                'best_percentile_value':    None,
                'murine_alleles_bound':     '',
                'num_murine_alleles_bound': 0,
            })
            continue

        alleles_sorted_by_percentile = bound_rows['allele'].tolist()
        best_percentile_value = float(bound_rows.iloc[0]['percentile'])
        best_percentile_label = bound_rows.iloc[0]['binder_label']

        aggregated_rows.append({
            'peptide':                  peptide_string,
            'best_percentile_label':    best_percentile_label,
            'best_percentile_value':    round(best_percentile_value, 2),
            'murine_alleles_bound':     ';'.join(alleles_sorted_by_percentile),
            'num_murine_alleles_bound': len(alleles_sorted_by_percentile),
        })

    return pd.DataFrame(aggregated_rows)


