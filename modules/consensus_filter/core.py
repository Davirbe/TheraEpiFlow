"""
Filtering logic for consensus_filter: column resolution, the Stage-1
presentation filter (load / consolidate / intersect / finalize) and the
Stage-2 Calis 2013 immunogenicity scoring. Heavy compute, no Rich panels.
"""

import contextlib
import io
from pathlib import Path

import pandas as pd
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

from utils.console import console, flush_stdin

_NET_PERCENTILE_SUBSTRINGS    = ['netmhcpan_el', 'percentile']
_FLURRY_PERCENTILE_SUBSTRINGS = ['mhcflurry', 'presentation', 'percentile']
# ── DataFrame helpers ──────────────────────────────────────────────────────────

def _find_col(dataframe: pd.DataFrame, substrings: list) -> str | None:
    """Returns the first column name that contains all given substrings."""
    for column_name in dataframe.columns:
        if all(substring in column_name for substring in substrings):
            return column_name
    return None


def _normalize_col_name(raw_name: str) -> str:
    """Lowercase + underscores. 'seq #' → 'sequence_number'."""
    cleaned = raw_name.strip()
    if cleaned.lower() == 'seq #':
        return 'sequence_number'
    return cleaned.replace(' ', '_').replace('#', 'num').lower()


# ── Stage 1 — Presentation filter ─────────────────────────────────────────────

def _load_and_filter(csv_path: Path, percentile_substrings: list, threshold: float) -> dict:
    """Loads a prediction CSV and applies sub-phase 0a (drop NaNs) and 0b (percentile ≤ threshold).
    Returns both intermediate DataFrames plus audit counts."""
    dataframe = pd.read_csv(csv_path, sep=',', dtype=str, keep_default_na=False,
                            na_values=['', 'NA', 'N/A'], engine='python')
    if dataframe.shape[1] == 1:
        dataframe = pd.read_csv(csv_path, sep=';', dtype=str, keep_default_na=False,
                                na_values=['', 'NA', 'N/A'], engine='python')

    dataframe.columns = [_normalize_col_name(col) for col in dataframe.columns]
    if 'peptide' in dataframe.columns:
        dataframe['peptide'] = dataframe['peptide'].str.strip()

    numeric_columns = [col for col in dataframe.columns
                       if any(keyword in col for keyword in ('percentile', 'score', 'ic50', 'rank', 'affinity'))]
    for numeric_col in numeric_columns:
        dataframe[numeric_col] = pd.to_numeric(
            dataframe[numeric_col].str.replace(',', '.', regex=False), errors='coerce'
        )

    percentile_column = _find_col(dataframe, percentile_substrings)
    if percentile_column is None:
        raise ValueError(
            f'Percentile column not found in {csv_path}. '
            f'Expected substrings: {percentile_substrings}. '
            f'Available columns: {list(dataframe.columns)}'
        )

    # 0a — drop NaN
    n_raw       = len(dataframe)
    df_0a       = dataframe.dropna(subset=['peptide', 'allele', percentile_column]).copy()
    n_after_nan = len(df_0a)

    # 0b — apply percentile threshold
    df_0b               = df_0a[df_0a[percentile_column] <= threshold].copy()
    n_after_threshold   = len(df_0b)

    # Reference counts for comparison display (strong and weak binder breakpoints)
    n_strong_binders = int((df_0a[percentile_column] <= 0.5).sum())
    n_weak_binders   = int((df_0a[percentile_column] <= 2.0).sum())

    return {
        'df_0a':       df_0a,
        'df_0b':       df_0b,
        'pct_col':     percentile_column,
        'n_raw':       n_raw,
        'n_0a':        n_after_nan,
        'dropped_nan': n_raw - n_after_nan,
        'n_0b':        n_after_threshold,
        'dropped_thr': n_after_nan - n_after_threshold,
        'n_strong':    n_strong_binders,
        'n_weak':      n_weak_binders,
    }


def _consolidate(
    df: pd.DataFrame,
    pct_col: str,
    tool_prefix: str,
    metric_name: str,
) -> pd.DataFrame:
    """Phase 1 — collapses N rows per allele into 1 row per peptide.

    Output columns (prefix-only naming, contract for downstream steps):
      peptide, {tool_prefix}_best_allele, {tool_prefix}_{metric_name}_percentile (min across alleles),
      {tool_prefix}_alleles (alphabetical, ';'-joined), {tool_prefix}_{metric_name}_percentiles_all (same order),
      {tool_prefix}_num_alleles. tool_prefix ∈ {netmhcpan, mhcflurry}; metric_name ∈ {el, presentation}.
    """
    if df.empty:
        return pd.DataFrame()

    aggregated_rows = []

    for peptide_sequence, allele_group in df.groupby('peptide'):
        best_percentile_per_allele = (
            allele_group.groupby('allele')[pct_col]
            .min()
            .sort_index()
        )
        alleles_joined     = ';'.join(best_percentile_per_allele.index.tolist())
        percentiles_joined = ';'.join(f"{value:.2f}" for value in best_percentile_per_allele.values)
        number_of_alleles  = len(best_percentile_per_allele)

        best_row_index  = allele_group[pct_col].idxmin()
        best_percentile = round(float(allele_group.loc[best_row_index, pct_col]), 2)
        best_allele     = str(allele_group.loc[best_row_index, 'allele'])

        aggregated_rows.append({
            'peptide':                                          peptide_sequence,
            f'{tool_prefix}_best_allele':                      best_allele,
            f'{tool_prefix}_{metric_name}_percentile':         best_percentile,
            f'{tool_prefix}_alleles':                          alleles_joined,
            f'{tool_prefix}_{metric_name}_percentiles_all':    percentiles_joined,
            f'{tool_prefix}_num_alleles':                      number_of_alleles,
        })

    return pd.DataFrame(aggregated_rows)


def _intersect(net_consolidated: pd.DataFrame, flurry_consolidated: pd.DataFrame) -> dict:
    """Phase 2 — keeps peptides present in BOTH tools; returns the common-peptide DataFrame + audit counts."""
    if net_consolidated.empty or flurry_consolidated.empty:
        return {
            'common':       pd.DataFrame(columns=['peptide']),
            'net_only':     len(net_consolidated),
            'flurry_only':  len(flurry_consolidated),
            'common_count': 0,
        }

    peptides_in_net    = set(net_consolidated['peptide'])
    peptides_in_flurry = set(flurry_consolidated['peptide'])
    peptides_in_both   = peptides_in_net & peptides_in_flurry

    return {
        'common':       pd.DataFrame({'peptide': sorted(peptides_in_both)}),
        'net_only':     len(peptides_in_net - peptides_in_flurry),
        'flurry_only':  len(peptides_in_flurry - peptides_in_net),
        'common_count': len(peptides_in_both),
    }


def _finalize(
    common_peptides: pd.DataFrame,
    net_consolidated: pd.DataFrame,
    flurry_consolidated: pd.DataFrame,
) -> pd.DataFrame:
    """Phase 3 — filters each tool to common peptides and merges (column names already
    follow the prefix-only convention from _consolidate)."""
    if common_peptides.empty:
        return pd.DataFrame(columns=['peptide'])

    net_filtered    = net_consolidated[net_consolidated['peptide'].isin(common_peptides['peptide'])].copy()
    flurry_filtered = flurry_consolidated[flurry_consolidated['peptide'].isin(common_peptides['peptide'])].copy()

    merged = pd.merge(common_peptides, net_filtered,    on='peptide', how='left')
    merged = pd.merge(merged,          flurry_filtered, on='peptide', how='left')
    return merged


# ── Stage 2 — Calis 2013 ─────────────────────────────────────────────────────

# Allele → anchor positions (1-based), reproduced from predict_immunogenicity.py
# to resolve the mask without calling validate() (which uses sys.exit and reads a file).
_CALIS_ALLELES = {
    "H-2-Db":"2,5,9","H-2-Dd":"2,3,5","H-2-Kb":"2,3,9","H-2-Kd":"2,5,9",
    "H-2-Kk":"2,8,9","H-2-Ld":"2,5,9","HLA-A0101":"2,3,9","HLA-A0201":"1,2,9",
    "HLA-A0202":"1,2,9","HLA-A0203":"1,2,9","HLA-A0206":"1,2,9","HLA-A0211":"1,2,9",
    "HLA-A0301":"1,2,9","HLA-A1101":"1,2,9","HLA-A2301":"2,7,9","HLA-A2402":"2,7,9",
    "HLA-A2601":"1,2,9","HLA-A2902":"2,7,9","HLA-A3001":"1,3,9","HLA-A3002":"2,7,9",
    "HLA-A3101":"1,2,9","HLA-A3201":"1,2,9","HLA-A3301":"1,2,9","HLA-A6801":"1,2,9",
    "HLA-A6802":"1,2,9","HLA-A6901":"1,2,9","HLA-B0702":"1,2,9","HLA-B0801":"2,5,9",
    "HLA-B1501":"1,2,9","HLA-B1502":"1,2,9","HLA-B1801":"1,2,9","HLA-B2705":"2,3,9",
    "HLA-B3501":"1,2,9","HLA-B3901":"1,2,9","HLA-B4001":"1,2,9","HLA-B4002":"1,2,9",
    "HLA-B4402":"2,3,9","HLA-B4403":"2,3,9","HLA-B4501":"1,2,9","HLA-B4601":"1,2,9",
    "HLA-B5101":"1,2,9","HLA-B5301":"1,2,9","HLA-B5401":"1,2,9","HLA-B5701":"1,2,9",
    "HLA-B5801":"1,2,9",
}


def _score_calis(peptide_list: list, allele_imgt: str | None) -> dict:
    """Scores peptides via Calis 2013 (predict_immunogenicity.Prediction); returns {peptide: score}.
    Captures stdout because Prediction.predict() prints instead of returning.
    allele_imgt: IMGT format (e.g. 'HLA-A*02:01') or None → default anchor mask."""
    from .predict_immunogenicity import Prediction

    allele_calis = None
    anchor_mask  = None
    if allele_imgt:
        allele_normalized = allele_imgt.replace('*', '').replace(':', '')
        if allele_normalized in _CALIS_ALLELES:
            allele_calis = allele_normalized
            anchor_mask  = _CALIS_ALLELES[allele_normalized]

    captured_output = io.StringIO()
    with contextlib.redirect_stdout(captured_output):
        Prediction().predict((list(peptide_list), anchor_mask, allele_calis))

    # Parse output lines with format "PEPTIDE,length,score"
    peptide_scores = {}
    for output_line in captured_output.getvalue().splitlines():
        output_line = output_line.strip()
        if not output_line or ',' not in output_line:
            continue
        if output_line.startswith(('peptide,', 'masking:', 'masked', 'allele:')):
            continue
        line_parts = output_line.split(',')
        if len(line_parts) == 3:
            try:
                peptide_scores[line_parts[0]] = float(line_parts[2])
            except ValueError:
                continue
    return peptide_scores


def _apply_calis(consensus_df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Applies Calis 2013 grouped by netmhcpan_best_allele; returns (df filtered to score > 0, audit)."""
    if consensus_df.empty:
        return consensus_df.assign(calis_score=pd.Series(dtype=float)), {
            'input': 0, 'scored': 0, 'survived': 0, 'unsupported_allele': 0,
        }

    peptide_score_map  = {}
    peptide_allele_map = {}
    unsupported_count  = 0

    calis_groups_by_allele = list(
        consensus_df.groupby('netmhcpan_best_allele', dropna=False)
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(
            style="yellow",
            complete_style="green",
            finished_style="green",
            pulse_style="bold yellow",
        ),
        MofNCompleteColumn(),
        console=console,
        transient=True,
    ) as calis_progress_bar:
        calis_task_id = calis_progress_bar.add_task(
            "  Scoring Calis 2013 immunogenicity (per allele)",
            total=len(calis_groups_by_allele),
        )
        for allele_name, allele_group_df in calis_groups_by_allele:
            peptides_in_allele_group = allele_group_df['peptide'].tolist()
            if pd.isna(allele_name) or not allele_name:
                calis_score_result      = _score_calis(peptides_in_allele_group, None)
                allele_label_for_audit  = 'default_mask'
                unsupported_count      += len(peptides_in_allele_group)
            else:
                allele_normalized = allele_name.replace('*', '').replace(':', '')
                if allele_normalized not in _CALIS_ALLELES:
                    calis_score_result      = _score_calis(peptides_in_allele_group, None)
                    allele_label_for_audit  = f'{allele_name} (default_mask)'
                    unsupported_count      += len(peptides_in_allele_group)
                else:
                    calis_score_result      = _score_calis(peptides_in_allele_group, allele_name)
                    allele_label_for_audit  = allele_name

            for peptide_sequence in peptides_in_allele_group:
                peptide_score_map[peptide_sequence]  = calis_score_result.get(peptide_sequence)
                peptide_allele_map[peptide_sequence] = allele_label_for_audit

            calis_progress_bar.advance(calis_task_id)

    scored_df = consensus_df.copy()
    scored_df['calis_score']       = scored_df['peptide'].map(peptide_score_map)
    scored_df['calis_allele_used'] = scored_df['peptide'].map(peptide_allele_map)

    immunogenic_df = scored_df[
        scored_df['calis_score'].notna() & (scored_df['calis_score'] > 0)
    ].copy()

    return immunogenic_df, {
        'input':               len(consensus_df),
        'scored':              int(scored_df['calis_score'].notna().sum()),
        'survived':            len(immunogenic_df),
        'unsupported_allele':  unsupported_count,
    }


