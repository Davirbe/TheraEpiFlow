"""
Standardized file and column naming conventions for TheraEPIflow.

File naming pattern:  {STEP}_{TOOL}_{TRACK_ID}.{ext}
Column naming pattern: {value}_{tool}_{type}

Examples:
  Files:   PRED_NET_HPV16_E1.csv
           CONSENSUS_HPV16_E1.csv
           CLUSTER_IEDB_HPV16_E1.csv
           TOXICITY_SAFE_HPV16_E1.csv

  Columns: netmhcpan_el_percentile_net_pred
           HLAs_agregados_flurry_pred
           Num_HLAs_unidos_net_flurry_pred
"""

from __future__ import annotations

from pathlib import Path


# ── Column suffixes ───────────────────────────────────────────────────────────
# These suffixes identify the origin of each column across all pipeline steps.

SUFFIX_NET_PRED    = "_net_pred"       # NetMHCpan — EL prediction
SUFFIX_FLURRY_PRED = "_flurry_pred"    # MHCFlurry — EL prediction
SUFFIX_NET_PROC    = "_net_proc"       # NetMHCpan — BA processing
SUFFIX_FLURRY_PROC = "_flurry_proc"    # MHCFlurry — BA processing


# ── Critical column names (fixed anchors used across steps) ───────────────────

COLUMN_PEPTIDE                    = "peptide"
COLUMN_NETMHC_EL_PERCENTILE       = f"netmhcpan_el_percentile{SUFFIX_NET_PRED}"
COLUMN_FLURRY_PERCENTILE          = f"mhcflurry_percentile{SUFFIX_FLURRY_PRED}"
COLUMN_HLAS_NET                   = f"HLAs_agregados{SUFFIX_NET_PRED}"
COLUMN_HLAS_FLURRY                = f"HLAs_agregados{SUFFIX_FLURRY_PRED}"
COLUMN_NUM_HLAS_NET               = f"Num_HLAs{SUFFIX_NET_PRED}"
COLUMN_NUM_HLAS_FLURRY            = f"Num_HLAs{SUFFIX_FLURRY_PRED}"
COLUMN_HLAS_UNITED                = f"HLAs_agregados_united_net_flurry_pred"
COLUMN_NUM_HLAS_UNITED            = f"Num_HLAs_united_net_flurry_pred"
COLUMN_BEST_REPRESENTATIVE        = "BEST_REPRESENTATIVE"   # '★' for selected epitopes


# ── Track ID ──────────────────────────────────────────────────────────────────

def build_track_id(genotype: str, protein: str) -> str:
    """
    Builds a standardized track ID from genotype and protein.

    Examples:
        build_track_id("HPV16", "E1")       → "HPV16_E1"
        build_track_id("ZIKV", "ENVELOPE")  → "ZIKV_ENVELOPE"
        build_track_id("SARS-CoV-2", "S")   → "SARS-CoV-2_S"
    """
    return f"{genotype.upper()}_{protein.upper()}"


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
        get_prediction_filename("NET_PROC", "HPV16_E1")    → "PROC_NET_HPV16_E1.csv"
        get_prediction_filename("FLURRY_PROC", "HPV16_E1") → "PROC_FLURRY_HPV16_E1.csv"
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
        step:     Pipeline step label (e.g. "CONSENSUS", "CLUSTER_IEDB", "TOXICITY_SAFE")
        track_id: e.g. "HPV16_E1"
        tool:     Optional tool label inserted between step and track_id
        ext:      File extension (default: "csv")

    Examples:
        get_step_filename("CONSENSUS", "HPV16_E1")              → "CONSENSUS_HPV16_E1.csv"
        get_step_filename("CLUSTER", "HPV16_E1", tool="IEDB")   → "CLUSTER_IEDB_HPV16_E1.csv"
        get_step_filename("CLUSTER_REPR", "HPV16_E1", ext="xlsx")→ "CLUSTER_REPR_HPV16_E1.xlsx"
        get_step_filename("TOXICITY", "HPV16_E1", tool="SAFE")  → "TOXICITY_SAFE_HPV16_E1.csv"
        get_step_filename("VARIANTS", "HPV16_E1", ext="fasta")  → "VARIANTS_HPV16_E1.fasta"
    """
    parts = [step.upper()]
    if tool:
        parts.append(tool.upper())
    parts.append(track_id)
    return f"{'_'.join(parts)}.{ext}"


# ── HLA format conversion ─────────────────────────────────────────────────────

def allele_to_netmhcpan_format(allele: str) -> str:
    """
    Converts standard IMGT allele format to NetMHCpan format.

    NetMHCpan does not accept '*' or ':' in allele names.

    Examples:
        allele_to_netmhcpan_format("HLA-A*02:01") → "HLA-A0201"
        allele_to_netmhcpan_format("HLA-B*07:02") → "HLA-B0702"
        allele_to_netmhcpan_format("H-2Db")       → "H-2Db"  (murine, unchanged)
    """
    return allele.replace("*", "").replace(":", "")


def alleles_to_netmhcpan_format(alleles: list[str]) -> list[str]:
    """Converts a list of alleles to NetMHCpan format."""
    return [allele_to_netmhcpan_format(a) for a in alleles]


# ── Column name resolution ────────────────────────────────────────────────────
# Used across steps when column names may vary between prediction tool outputs.

def find_column_name(dataframe: pd.DataFrame, possible_names: list[str]) -> str | None:
    """
    Finds the first matching column name from a list of candidates.

    Returns the first matching column name found, or None if none match.

    Example:
        col = find_column_name(df, ["EL_Rank", "el_rank", "%Rank_EL"])
    """
    for name in possible_names:
        if name in dataframe.columns:
            return name
    return None


def require_column(dataframe: pd.DataFrame, possible_names: list[str], context: str = "") -> str:
    """
    Like find_column_name, but raises ValueError if no match is found.

    Args:
        context: Optional label for the error message (e.g. step name).
    """
    result = find_column_name(dataframe, possible_names)
    if result is None:
        prefix = f"[{context}] " if context else ""
        raise ValueError(
            f"{prefix}None of the expected columns found: {possible_names}. "
            f"Available columns: {list(dataframe.columns)}"
        )
    return result
