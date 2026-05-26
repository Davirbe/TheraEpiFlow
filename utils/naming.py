"""
Standardized file and column naming conventions for TheraEpiFlow.

File naming pattern:  {STEP}_{TOOL}_{TRACK_ID}.{ext}
Column naming pattern: {tool}_{metric}  (prefix-only — no redundant tool suffix)

Examples:
  Files:   PRED_NET_HPV16_E1.csv
           CONSENSUS_IMMUNOGENIC_HPV16_E1.csv
           CLUSTER_HPV16_E1.csv
           TOXICITY_SAFE_HPV16_E1.csv

  Columns: netmhcpan_el_percentile          (not netmhcpan_el_percentile_net_pred)
           mhcflurry_presentation_percentile (not mhcflurry_presentation_percentile_flurry_pred)
           netmhcpan_alleles                 (not HLAs_agregados_net_pred)
           netmhcpan_el_percentiles_all      (all allele EL%, semicolon-separated)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


# ── NetMHCpan output columns ──────────────────────────────────────────────────

COLUMN_PEPTIDE                         = "peptide"
COLUMN_NETMHC_BEST_ALLELE              = "netmhcpan_best_allele"
COLUMN_NETMHC_EL_PERCENTILE            = "netmhcpan_el_percentile"
COLUMN_NETMHC_ALLELES                  = "netmhcpan_alleles"
COLUMN_NETMHC_NUM_ALLELES              = "netmhcpan_num_alleles"
COLUMN_NETMHC_EL_PERCENTILES_ALL       = "netmhcpan_el_percentiles_all"

# ── MHCFlurry output columns ──────────────────────────────────────────────────

COLUMN_FLURRY_BEST_ALLELE              = "mhcflurry_best_allele"
COLUMN_FLURRY_PERCENTILE               = "mhcflurry_presentation_percentile"
COLUMN_FLURRY_ALLELES                  = "mhcflurry_alleles"
COLUMN_FLURRY_NUM_ALLELES              = "mhcflurry_num_alleles"
COLUMN_FLURRY_PERCENTILES_ALL          = "mhcflurry_presentation_percentiles_all"

# ── Cross-tool aggregated columns ─────────────────────────────────────────────

COLUMN_ALLELES_UNITED                  = "alleles_united"
COLUMN_NUM_ALLELES_UNITED              = "num_alleles_united"
COLUMN_BEST_REPRESENTATIVE             = "BEST_REPRESENTATIVE"

# Marker written into the BEST_REPRESENTATIVE column by select_representatives
# and consumed as a filter by every downstream step (analyze_conservation,
# population_coverage, predict_murine, ...). Centralised here so the literal
# character lives in one place — a future refactor may replace it with a
# boolean column to eliminate any CSV encoding risk.
STAR_MARKER                            = "★"


# ── Track ID ──────────────────────────────────────────────────────────────────

def build_track_id(organism_label: str, protein_label: str) -> str:
    """
    Builds a standardized track ID from organism and protein labels.

    Examples:
        build_track_id("HPV16", "E1")       → "HPV16_E1"
        build_track_id("ZIKV", "ENVELOPE")  → "ZIKV_ENVELOPE"
        build_track_id("SARS-CoV-2", "S")   → "SARS-CoV-2_S"
    """
    return f"{organism_label.upper()}_{protein_label.upper()}"


def parse_track_id(track_id: str, project_config: dict) -> tuple[str, str]:
    """
    Returns (organism_label, protein_label) for a track by reading project_config.

    Reads the canonical labels stored at create-project time rather than splitting
    the string — naive split breaks on organism labels containing dashes
    (e.g. SARS-CoV-2_S would split as ['SARS', 'CoV', '2_S']).

    Falls back to the *_name fields if *_label is absent, and finally to the
    raw track_id pieces if neither is present.
    """
    track_meta = project_config.get('tracks', {}).get(track_id, {})
    organism = track_meta.get('organism_label') or track_meta.get('organism_name')
    protein  = track_meta.get('protein_label')  or track_meta.get('protein_name')
    if organism and protein:
        return organism, protein

    fallback_pieces = track_id.rsplit('_', 1)
    if len(fallback_pieces) == 2:
        return organism or fallback_pieces[0], protein or fallback_pieces[1]
    return organism or track_id, protein or ''


# ── File name builders ────────────────────────────────────────────────────────

def get_prediction_filename(tool: str, track_id: str) -> str:
    """
    Returns the standardized filename for prediction output files.

    Args:
        tool:     "NET_PRED", "FLURRY_PRED", "NET_PROC", or "FLURRY_PROC"
        track_id: e.g. "HPV16_E1"

    Examples:
        get_prediction_filename("NET_PRED", "HPV16_E1")    → "PRED_NET_HPV16_E1.csv"
        get_prediction_filename("FLURRY_PRED", "HPV16_E1") → "PRED_FLURRY_HPV16_E1.csv"
    """
    tool_map = {
        "NET_PRED":    "PRED_NET",
        "FLURRY_PRED": "PRED_FLURRY",
        "NET_PROC":    "PROC_NET",
        "FLURRY_PROC": "PROC_FLURRY",
    }
    prefix = tool_map.get(tool.upper(), tool.upper())
    return f"{prefix}_{track_id}.csv"


def get_step_filename(step: str, track_id: str, tool: str = "", ext: str = "csv") -> str:
    """
    Returns the standardized filename for any pipeline step output.

    Args:
        step:     Pipeline step label (e.g. "CONSENSUS_IMMUNOGENIC", "CLUSTER", "TOXICITY_SAFE")
        track_id: e.g. "HPV16_E1"
        tool:     Optional tool label inserted between step and track_id
        ext:      File extension (default: "csv")

    Examples:
        get_step_filename("CONSENSUS_IMMUNOGENIC", "HPV16_E1")    → "CONSENSUS_IMMUNOGENIC_HPV16_E1.csv"
        get_step_filename("CLUSTER", "HPV16_E1")                  → "CLUSTER_HPV16_E1.csv"
        get_step_filename("CLUSTER_REPR", "HPV16_E1", ext="xlsx") → "CLUSTER_REPR_HPV16_E1.xlsx"
        get_step_filename("TOXICITY", "HPV16_E1", tool="SAFE")    → "TOXICITY_SAFE_HPV16_E1.csv"
    """
    filename_parts = [step.upper()]
    if tool:
        filename_parts.append(tool.upper())
    filename_parts.append(track_id)
    return f"{'_'.join(filename_parts)}.{ext}"


# ── HLA validation + normalization ────────────────────────────────────────────
#
# The IEDB classic NetMHCpan API and the MHCFlurry presentation predictor both
# expect the standard IMGT format `HLA-A*02:01`. Users sometimes type the
# allele without the asterisk (`HLA-A02:01`) or in lowercase. These functions
# detect and auto-correct those input forms so a typo does not silently
# corrupt the run.

import re as _re

_HLA_PATTERN = _re.compile(
    r"^HLA-([ABCDEG])\*?(\d{2,3}):(\d{2,3})$",
    _re.IGNORECASE,
)


def parse_hla_allele(raw: str) -> tuple[str | None, str | None]:
    """
    Validates and normalizes a single HLA Class I allele to IMGT format.

    Returns:
        (normalized, correction_note) where:
        - normalized:      the IMGT-formatted allele (e.g. 'HLA-A*02:01'),
                           or None when the input does not look like an HLA allele.
        - correction_note: None when the input was already in canonical form,
                           or a short human-readable string describing what
                           was changed (e.g. 'added missing "*"').

    Murine alleles (e.g. 'H-2Db') return (None, None) so the caller knows they
    are not HLA — murine handling lives in predict_murine.

    Examples:
        parse_hla_allele('HLA-A*02:01') → ('HLA-A*02:01', None)
        parse_hla_allele('HLA-A02:01')  → ('HLA-A*02:01', 'added missing "*"')
        parse_hla_allele('hla-a*02:01') → ('HLA-A*02:01', 'uppercased')
        parse_hla_allele('foo')         → (None, None)
    """
    if not raw:
        return None, None
    cleaned = raw.strip()
    match = _HLA_PATTERN.match(cleaned)
    if match is None:
        return None, None

    locus_letter, family, member = match.group(1), match.group(2), match.group(3)
    normalized = f"HLA-{locus_letter.upper()}*{family}:{member}"

    if cleaned == normalized:
        return normalized, None

    # Build a tiny audit string so the user can see what was changed.
    corrections: list[str] = []
    if "*" not in cleaned:
        corrections.append('added missing "*"')
    if cleaned != cleaned.upper() and cleaned.upper() != normalized:
        # uppercase only matters when it actually flipped the case of letters
        pass
    if cleaned.upper() != cleaned and "*" in cleaned:
        corrections.append("uppercased")
    elif cleaned.upper() != cleaned:
        corrections.append("uppercased")
    note = " + ".join(corrections) if corrections else "normalized"
    return normalized, note


# ── Column name resolution ────────────────────────────────────────────────────

def find_column_name(dataframe: pd.DataFrame, possible_names: list[str]) -> str | None:
    """
    Finds the first matching column name from a list of candidates.

    Returns the first matching column name found, or None if none match.

    Example:
        col = find_column_name(df, ["EL_Rank", "el_rank", "%Rank_EL"])
    """
    for candidate_name in possible_names:
        if candidate_name in dataframe.columns:
            return candidate_name
    return None
