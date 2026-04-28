"""
Phase functions for the consensus filter step.

Adapted from existing_scripts/analise_consenso_predicao.py — same 4-phase
structure, same Portuguese-leaning column names (HLAs_agregados, Num_HLAs,
melhor_allele) so the column-naming contract carries through to step
"select_representatives" untouched. The PROC_NET / PROC_FLURRY paths from the
original were dropped (we use only NetMHCpan EL and MHCflurry presentation).
The Phase 0b threshold filter is new — the original consumed pre-filtered
tables, but we run on raw step "predict_binding" output.

Phases:
  0a  load + drop NaN in essential columns (peptide, allele, percentile)
  0b  apply percentile threshold (≤ user-chosen value)
  1   consolidate per peptide: aggregate alleles, count, pick best percentile row
  2   intersect peptides between the two tools
  3   finalize: append _net_pred / _flurry_pred suffixes and merge

Each phase returns audit info so the orchestrator can render a control panel
showing how many peptides survived each filter — a scientific necessity since
the user wants to verify nothing was silently lost.
"""

from pathlib import Path

import pandas as pd


# ── Helpers ───────────────────────────────────────────────────────────────────

def find_column_by_substring_match(
    dataframe: pd.DataFrame,
    required_substrings: list,
) -> str | None:
    """
    Returns the first column name in `dataframe` whose name contains EVERY
    substring in `required_substrings`. Equivalent to find_column_name() in
    analise_consenso_predicao.py — substring matching, not exact match, because
    the actual column names emitted by step "predict_binding" depend on the
    tool versions and the IEDB API response shape.
    """
    for column_name in dataframe.columns:
        if all(substring in column_name for substring in required_substrings):
            return column_name
    return None


def _normalize_column_name(raw_name: str) -> str:
    """
    Lowercase, replace spaces with underscores, '#' with 'num'.
    Equivalent to limpar_nome_coluna() in analise_consenso_predicao.py.
    """
    cleaned = raw_name.strip()
    if cleaned.lower() == 'seq #':
        return 'sequence_number'
    return cleaned.replace(' ', '_').replace('#', 'num').lower()


def load_prediction_csv_with_decimal_normalization(csv_path: Path) -> pd.DataFrame:
    """
    Loads a prediction CSV from step "predict_binding". Behavior preserved from
    carregar_e_limpar_csv() in analise_consenso_predicao.py:
      - autodetect ',' vs ';' separator
      - lowercase / normalize column names
      - strip 'peptide' values
      - convert decimal commas to dots in numeric-looking columns
    """
    raw_dataframe = pd.read_csv(
        csv_path, sep=',', dtype=str, keep_default_na=False,
        na_values=['', 'NA', 'N/A'], engine='python',
    )
    if raw_dataframe.shape[1] == 1:
        raw_dataframe = pd.read_csv(
            csv_path, sep=';', dtype=str, keep_default_na=False,
            na_values=['', 'NA', 'N/A'], engine='python',
        )

    raw_dataframe.columns = [_normalize_column_name(c) for c in raw_dataframe.columns]

    if 'peptide' in raw_dataframe.columns:
        raw_dataframe['peptide'] = raw_dataframe['peptide'].str.strip()

    numeric_like_columns = [
        column_name for column_name in raw_dataframe.columns
        if 'percentile' in column_name
        or 'score' in column_name
        or 'ic50' in column_name
        or 'rank' in column_name
        or 'affinity' in column_name
    ]
    for column_name in numeric_like_columns:
        raw_dataframe[column_name] = (
            raw_dataframe[column_name].str.replace(',', '.', regex=False)
        )
        raw_dataframe[column_name] = pd.to_numeric(
            raw_dataframe[column_name], errors='coerce',
        )

    return raw_dataframe


# ── PHASE 0 — Load + Clean + Threshold ────────────────────────────────────────

def phase0_clean_and_threshold(
    raw_dataframe: pd.DataFrame,
    source_label: str,
    percentile_column_substrings: list,
    percentile_threshold_max: float,
) -> dict:
    """
    Phase 0 in two sub-phases:
      0a) Drop rows with NaN in essential columns (peptide, allele, percentile).
      0b) Keep only rows where percentile ≤ percentile_threshold_max.

    Args:
        raw_dataframe:                  the prediction table from step "predict_binding"
        source_label:                   'NetMHCpan' or 'MHCflurry' (for error msg)
        percentile_column_substrings:   passed to find_column_by_substring_match
        percentile_threshold_max:       e.g. 2.0 (weak+strong) or 0.5 (strong only)

    Returns a dict with both intermediate dataframes and per-sub-phase counts.
    """
    percentile_column_name = find_column_by_substring_match(
        raw_dataframe, percentile_column_substrings,
    )
    if percentile_column_name is None:
        raise ValueError(
            f'[{source_label}] No column matched substrings '
            f'{percentile_column_substrings}. '
            f'Available columns: {list(raw_dataframe.columns)}'
        )

    essential_column_names = ['peptide', 'allele', percentile_column_name]
    missing_essential_columns = [
        column for column in essential_column_names
        if column not in raw_dataframe.columns
    ]
    if missing_essential_columns:
        raise ValueError(
            f'[{source_label}] Missing essential columns: {missing_essential_columns}. '
            f'Available: {list(raw_dataframe.columns)}'
        )

    # Sub-phase 0a — drop NaN rows
    phase0a_input_count = len(raw_dataframe)
    phase0a_dataframe = raw_dataframe.dropna(subset=essential_column_names).copy()
    phase0a_output_count = len(phase0a_dataframe)
    phase0a_dropped_count = phase0a_input_count - phase0a_output_count

    # Sub-phase 0b — apply percentile threshold
    phase0b_dataframe = phase0a_dataframe[
        phase0a_dataframe[percentile_column_name] <= percentile_threshold_max
    ].copy()
    phase0b_output_count = len(phase0b_dataframe)
    phase0b_dropped_count = phase0a_output_count - phase0b_output_count

    return {
        'phase0a_dataframe':         phase0a_dataframe,
        'phase0b_dataframe':         phase0b_dataframe,
        'percentile_column_name':    percentile_column_name,
        'phase0a_input_count':       phase0a_input_count,
        'phase0a_dropped_count':     phase0a_dropped_count,
        'phase0a_output_count':      phase0a_output_count,
        'phase0b_dropped_count':     phase0b_dropped_count,
        'phase0b_output_count':      phase0b_output_count,
    }


# ── PHASE 1 — Consolidation per peptide ──────────────────────────────────────

def phase1_consolidate_per_peptide(
    thresholded_dataframe: pd.DataFrame,
    percentile_column_name: str,
) -> pd.DataFrame:
    """
    Phase 1 — collapse multiple per-allele rows into one row per peptide.

    Adapted from consolidar_epitopos() in analise_consenso_predicao.py.
    Same column naming convention. Difference: the original dropped 'allele'
    entirely; we rename it to 'melhor_allele' so Stage 2 (Calis) knows which
    allele to use for the anchor-position mask.

    Output columns:
      peptide                        — group key
      melhor_allele                  — allele that produced the best (lowest) percentile
      <percentile_column_name>       — that best percentile value
      <other source columns>         — from the row of melhor_allele
      HLAs_agregados                 — semicolon-joined unique HLAs for this peptide
      Num_HLAs                       — count of unique HLAs
    """
    if thresholded_dataframe.empty:
        return pd.DataFrame()

    grouped_by_peptide = thresholded_dataframe.groupby('peptide')

    hla_aggregated_per_peptide = grouped_by_peptide['allele'].apply(
        lambda alleles_for_peptide: ';'.join(sorted(alleles_for_peptide.unique()))
    ).reset_index().rename(columns={'allele': 'HLAs_agregados'})
    hla_aggregated_per_peptide['Num_HLAs'] = (
        hla_aggregated_per_peptide['HLAs_agregados']
        .apply(lambda joined: len(joined.split(';')))
    )

    best_percentile_row_indices = grouped_by_peptide[percentile_column_name].idxmin()
    best_percentile_rows = thresholded_dataframe.loc[best_percentile_row_indices].copy()
    best_percentile_rows = best_percentile_rows.rename(
        columns={'allele': 'melhor_allele'}
    )

    consolidated_dataframe = pd.merge(
        best_percentile_rows,
        hla_aggregated_per_peptide,
        on='peptide',
        how='left',
    )
    return consolidated_dataframe


# ── PHASE 2 — Intersection ───────────────────────────────────────────────────

def phase2_intersect_common_peptides(
    netmhcpan_consolidated: pd.DataFrame,
    mhcflurry_consolidated: pd.DataFrame,
) -> dict:
    """
    Phase 2 — peptides present in BOTH consolidated tables.

    Returns a dict carrying the result dataframe and audit counts (how many
    peptides only one tool found — useful for the control panel).
    """
    if netmhcpan_consolidated.empty or mhcflurry_consolidated.empty:
        return {
            'common_peptides_dataframe': pd.DataFrame(columns=['peptide']),
            'netmhcpan_only_count':      len(netmhcpan_consolidated),
            'mhcflurry_only_count':      len(mhcflurry_consolidated),
            'common_count':              0,
        }

    netmhcpan_peptide_set = set(netmhcpan_consolidated['peptide'])
    mhcflurry_peptide_set = set(mhcflurry_consolidated['peptide'])
    common_peptide_set    = netmhcpan_peptide_set & mhcflurry_peptide_set

    return {
        'common_peptides_dataframe': pd.DataFrame({
            'peptide': sorted(common_peptide_set),
        }),
        'netmhcpan_only_count': len(netmhcpan_peptide_set - mhcflurry_peptide_set),
        'mhcflurry_only_count': len(mhcflurry_peptide_set - netmhcpan_peptide_set),
        'common_count':         len(common_peptide_set),
    }


# ── PHASE 3 — Final table with source-of-origin prefixes ─────────────────────

def phase3_finalize_with_source_prefixes(
    common_peptides_dataframe: pd.DataFrame,
    netmhcpan_consolidated: pd.DataFrame,
    mhcflurry_consolidated: pd.DataFrame,
) -> pd.DataFrame:
    """
    Phase 3 — assemble the final consensus table.

    For each common peptide, pull the per-source consolidated row, append the
    suffix _net_pred / _flurry_pred to every column except 'peptide' (so the
    same column names from the two tools never collide), and merge on peptide.

    This naming convention is the contract carried into step
    "select_representatives" — keep it stable.

    Output columns:
      peptide
      melhor_allele_net_pred,    HLAs_agregados_net_pred,    Num_HLAs_net_pred,    <other>_net_pred
      melhor_allele_flurry_pred, HLAs_agregados_flurry_pred, Num_HLAs_flurry_pred, <other>_flurry_pred
    """
    if common_peptides_dataframe.empty:
        return pd.DataFrame(columns=['peptide'])

    netmhcpan_in_common = netmhcpan_consolidated[
        netmhcpan_consolidated['peptide'].isin(common_peptides_dataframe['peptide'])
    ].copy()
    mhcflurry_in_common = mhcflurry_consolidated[
        mhcflurry_consolidated['peptide'].isin(common_peptides_dataframe['peptide'])
    ].copy()

    netmhcpan_in_common = netmhcpan_in_common.rename(columns={
        column_name: f'{column_name}_net_pred'
        for column_name in netmhcpan_in_common.columns if column_name != 'peptide'
    })
    mhcflurry_in_common = mhcflurry_in_common.rename(columns={
        column_name: f'{column_name}_flurry_pred'
        for column_name in mhcflurry_in_common.columns if column_name != 'peptide'
    })

    final_dataframe = pd.merge(
        common_peptides_dataframe, netmhcpan_in_common, on='peptide', how='left',
    )
    final_dataframe = pd.merge(
        final_dataframe, mhcflurry_in_common, on='peptide', how='left',
    )
    return final_dataframe
