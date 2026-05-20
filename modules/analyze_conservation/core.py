"""
Conservation domain logic for analyze_conservation: BLOSUM62 substitution
scoring, MHC-I anchor verdicts, sliding-window identity, per-position stats and
the IEDB-style summary builder. Pure compute — no Rich, no openpyxl, no prompts.
"""

from collections import Counter

import pandas as pd
from Bio.Align import substitution_matrices

import config
from utils.naming import (
    COLUMN_ALLELES_UNITED,
    COLUMN_NUM_ALLELES_UNITED,
    COLUMN_PEPTIDE,
)

# Accept variants whose length is within ±_LENGTH_TOLERANCE of the reference.
_LENGTH_TOLERANCE = 0.25

# Max mutations to qualify a (peptide, variant) pair for the MUTATIONS XLSX
# and to qualify an epitope for inclusion in the heatmap PNG.
_MAX_MUTATIONS_FOR_HEATMAP = 2

# BLOSUM62 substitution matrix from Biopython. Score >= 0 = conservative.
_BLOSUM62_MATRIX = substitution_matrices.load("BLOSUM62")
# ── Substitution scoring + anchor helpers ─────────────────────────────────────

def _blosum62_score(ref_aa: str, var_aa: str) -> int:
    """BLOSUM62 substitution score. >=0 = conservative substitution."""
    try:
        return int(_BLOSUM62_MATRIX[(ref_aa, var_aa)])
    except (KeyError, IndexError):
        return -4


def _anchor_positions(peptide_length: int) -> tuple[int, int]:
    """Returns (p2_idx, pomega_idx) 0-based for MHC-I anchor positions."""
    return 1, peptide_length - 1


def _classify_mhc_verdict(n_mutations: int, anchor_hit: bool, blosum62_min: int) -> str:
    """
    Heuristic verdict for variant tolerance based on anchor + BLOSUM62.

    likely_lost     → any mutation in P2 or PΩ
    excellent_match → 1 mutation, non-anchor, BLOSUM62 ≥ 0
    tolerated       → 1-2 non-anchor mutations otherwise
    """
    if anchor_hit:
        return "likely_lost"
    if n_mutations == 1 and blosum62_min >= 0:
        return "excellent_match"
    return "tolerated"


def _format_pct_fraction(numerator: int, denominator: int) -> str:
    """Returns "XX.XX% (n/N)". Empty denominator → "—"."""
    if denominator <= 0:
        return "—"
    pct = numerator / denominator * 100
    return f"{pct:.2f}% ({numerator}/{denominator})"


# ── Core computation ──────────────────────────────────────────────────────────

def compute_max_window_identity(peptide: str, sequence: str) -> tuple[float, str]:
    """
    Slides a window of len(peptide) over sequence, returning the maximum
    per-position identity (0.0–1.0) and the corresponding window string.
    Returns (0.0, "") when the sequence is shorter than the peptide.
    """
    peptide_length = len(peptide)
    best_identity  = 0.0
    best_window    = ""
    for start in range(len(sequence) - peptide_length + 1):
        window   = sequence[start : start + peptide_length]
        identity = sum(ref == var for ref, var in zip(peptide, window)) / peptide_length
        if identity > best_identity:
            best_identity = identity
            best_window   = window
    return best_identity, best_window


def classify_conservation_label(mean_max_identity: float, n_variants_total: int) -> str:
    if n_variants_total == 0:
        return "conservation_unknown"
    if mean_max_identity == 1.0:
        return "perfect"
    if mean_max_identity >= config.CONSERVATION_HIGH_THRESHOLD:
        return "high"
    if mean_max_identity >= config.CONSERVATION_MODERATE_THRESHOLD:
        return "moderate"
    return "low"


def _variant_label(record) -> str:
    """Returns a short 'accession(organism)' label from a SeqRecord."""
    desc_parts = record.description.strip().split(" ", 1)
    organism_short = (
        desc_parts[1].split("|")[0].strip()[:30]
        if len(desc_parts) > 1
        else record.id
    )
    return f"{record.id}({organism_short})"


def compute_epitope_conservation(
    peptide: str,
    records: list,
    analysis_threshold: float,
) -> tuple[dict, list[tuple[str, float, str]]]:
    """
    Computes all conservation metrics for a single peptide against every
    variant record. Returns (metrics_dict, alignment_tuples).

    alignment_tuples: [(variant_label, identity, best_window), ...] for every
    record, re-used downstream for mutation records and the position heatmap.

    Metrics carry RAW COUNTS (not pre-formatted percentages). The CSV/XLSX
    writers compose "XX.XX% (n/N)" strings from these.
    """
    n_total: int = len(records)
    alignment_tuples: list[tuple[str, float, str]] = []

    if n_total == 0:
        empty_metrics = {
            "n_variants_used":      0,
            "n_passed_threshold":   0,
            "n_exact_match":        0,
            "n_identity_90":        0,
            "n_identity_80":        0,
            "analysis_threshold":   analysis_threshold,
            "mean_max_identity":    0.0,
            "min_max_identity":     0.0,
            "max_max_identity":     0.0,
            "conservation_label":   "conservation_unknown",
        }
        return empty_metrics, alignment_tuples

    per_variant_identities: list[float] = []
    for record in records:
        sequence = str(record.seq).upper()
        identity, window = compute_max_window_identity(peptide, sequence)
        per_variant_identities.append(identity)
        alignment_tuples.append((_variant_label(record), identity, window))

    mean_id  = sum(per_variant_identities) / n_total
    min_id   = min(per_variant_identities)
    max_id   = max(per_variant_identities)
    n_exact  = sum(1 for i in per_variant_identities if i == 1.0)
    n_90pct  = sum(1 for i in per_variant_identities if i >= 0.90)
    n_80pct  = sum(1 for i in per_variant_identities if i >= 0.80)
    n_passed = sum(1 for i in per_variant_identities if i >= analysis_threshold)

    metrics = {
        "n_variants_used":      n_total,
        "n_passed_threshold":   n_passed,
        "n_exact_match":        n_exact,
        "n_identity_90":        n_90pct,
        "n_identity_80":        n_80pct,
        "analysis_threshold":   analysis_threshold,
        "mean_max_identity":    round(mean_id, 4),
        "min_max_identity":     round(min_id, 4),
        "max_max_identity":     round(max_id, 4),
        "conservation_label":   classify_conservation_label(mean_id, n_total),
    }
    return metrics, alignment_tuples


def compute_mutation_records(
    peptide: str,
    alignment_tuples: list[tuple[str, float, str]],
    alleles_united: str,
) -> list[dict]:
    """
    For each (variant, peptide) pair where the best window has 1 or 2
    substitutions, builds a record capturing position(s), BLOSUM62 score, and
    the MHC verdict. Pairs with 0 mutations or >2 mutations are skipped.
    """
    if not peptide or not alignment_tuples:
        return []

    peptide_length      = len(peptide)
    p2_idx, pomega_idx  = _anchor_positions(peptide_length)
    anchor_position_set = {p2_idx, pomega_idx}

    records = []
    for variant_label, _identity, window in alignment_tuples:
        if len(window) != peptide_length:
            continue

        mutation_triples = [
            (pos, ref_aa, var_aa)
            for pos, (ref_aa, var_aa) in enumerate(zip(peptide, window))
            if ref_aa != var_aa
        ]
        n_mutations = len(mutation_triples)
        if n_mutations < 1 or n_mutations > _MAX_MUTATIONS_FOR_HEATMAP:
            continue

        anchor_hit    = any(pos in anchor_position_set for pos, _, _ in mutation_triples)
        blosum_scores = [_blosum62_score(ref, var) for _, ref, var in mutation_triples]
        blosum_min    = min(blosum_scores)
        verdict       = _classify_mhc_verdict(n_mutations, anchor_hit, blosum_min)

        accession      = variant_label.split("(")[0]
        mutations_str  = "; ".join(f"P{pos + 1}:{ref}→{var}" for pos, ref, var in mutation_triples)
        positions_str  = ",".join(str(pos + 1) for pos, _, _ in mutation_triples)
        n_match        = peptide_length - n_mutations

        records.append({
            "peptide":            peptide,
            "length":             peptide_length,
            "variant_accession":  accession,
            "variant_window":     window,
            "identity":           _format_pct_fraction(n_match, peptide_length),
            "n_mutations":        n_mutations,
            "mutations":          mutations_str,
            "mutation_positions": positions_str,
            "anchor_hit":         "Y" if anchor_hit else "N",
            "blosum62_min":       blosum_min,
            "mhc_verdict":        verdict,
            "alleles_united":     alleles_united,
        })

    return records


def _has_variant_within_max_mut(peptide: str, alignment_tuples: list) -> bool:
    """True when at least one variant has ≤ _MAX_MUTATIONS_FOR_HEATMAP differences."""
    peptide_length = len(peptide)
    for _, _, window in alignment_tuples:
        if len(window) != peptide_length:
            continue
        n_mut = sum(a != b for a, b in zip(peptide, window))
        if n_mut <= _MAX_MUTATIONS_FOR_HEATMAP:
            return True
    return False


def _align_position_stats_to_anchor(
    pos_stats: list[dict],
    peptide_length: int,
    max_length: int,
) -> tuple[list, list]:
    """
    Maps a peptide's per-position stats into a max_length-wide grid by
    anchor-aligning both ends: N-terminal half left-justified, C-terminal
    half right-justified. Middle positions (the MHC-I "bulge") become None.

    Guarantees: P2 of every peptide lands at column 1, PΩ at column
    max_length - 1 — so anchor borders line up across mixed lengths.
    """
    aligned_pcts:  list = [None] * max_length
    aligned_annot: list = [""]   * max_length
    n_term_half = peptide_length // 2
    for pep_idx in range(peptide_length):
        col = (
            pep_idx
            if pep_idx < n_term_half
            else max_length - (peptide_length - pep_idx)
        )
        pct = pos_stats[pep_idx]["conservation_pct"]
        aligned_pcts[col] = pct
        if pct < 100.0:
            aligned_annot[col] = (
                f"{int(round(pct))}%\n{pos_stats[pep_idx]['top_mutation']}"
            )
    return aligned_pcts, aligned_annot


def _build_anchor_aligned_xtick_labels(max_length: int) -> list[str]:
    """
    X-axis labels matching the layout produced by
    _align_position_stats_to_anchor: N-terminal half uses absolute
    numbering (P1, P2, ...); C-terminal half uses offset-from-end
    (..., PΩ-1, PΩ).
    """
    n_term_half = max_length // 2
    n_term_labels = [f"P{i + 1}" for i in range(n_term_half)]
    c_term_count  = max_length - n_term_half
    c_term_labels = [
        ("PΩ" if (c_term_count - 1 - k) == 0 else f"PΩ-{c_term_count - 1 - k}")
        for k in range(c_term_count)
    ]
    return n_term_labels + c_term_labels


def _compute_position_stats(peptide: str, alignment_tuples: list[tuple]) -> list[dict]:
    """
    For each position in the peptide returns:
      conservation_pct  — % of variants with the same AA as reference (0–100)
      top_mutation      — most common substitution string e.g. "Q→R", or "" if 100%
    """
    n_variants = len(alignment_tuples)
    stats = []
    for pos_idx, ref_aa in enumerate(peptide):
        if n_variants == 0:
            stats.append({"conservation_pct": 0.0, "top_mutation": ""})
            continue
        n_match: int = 0
        mutation_counter: Counter = Counter()
        for _, _, window in alignment_tuples:
            if len(window) > pos_idx:
                var_aa = window[pos_idx]
                if var_aa == ref_aa:
                    n_match += 1
                else:
                    mutation_counter[f"{ref_aa}→{var_aa}"] += 1
        pct     = n_match / n_variants * 100
        top_mut = mutation_counter.most_common(1)[0][0] if mutation_counter else ""
        stats.append({"conservation_pct": pct, "top_mutation": top_mut})
    return stats


# ── IEDB-style summary table builder ──────────────────────────────────────────

_SUMMARY_COLUMNS = [
    "#",
    "peptide",
    "length",
    "pct_identity_100",
    "pct_identity_90",
    "pct_identity_80",
    "pct_identity_threshold",
    "min_identity",
    "max_identity",
    "avg_identity",
    "conservation_label",
    "n_excellent_match",
    "n_tolerated",
    "variants_exact_match",
    "variants_tolerable",
    "alleles_united",
    "num_alleles_united",
]


def _build_iedb_summary_rows(result_df: pd.DataFrame) -> pd.DataFrame:
    """
    Builds the IEDB-style presentation DataFrame. Tier columns are formatted
    as "XX.XX% (n/N)" strings; min/max/avg identity as "XX.XX%". Every other
    field passes through.
    """
    rows = []
    for idx, row in enumerate(result_df.itertuples(index=False), start=1):
        n_total      = int(getattr(row, "n_variants_used", 0))
        peptide      = getattr(row, COLUMN_PEPTIDE)
        num_alleles  = getattr(row, COLUMN_NUM_ALLELES_UNITED, 0)
        num_alleles  = int(num_alleles) if pd.notna(num_alleles) else 0

        rows.append({
            "#":                      idx,
            "peptide":                peptide,
            "length":                 len(peptide),
            "pct_identity_100":       _format_pct_fraction(int(getattr(row, "n_exact_match",      0)), n_total),
            "pct_identity_90":        _format_pct_fraction(int(getattr(row, "n_identity_90",      0)), n_total),
            "pct_identity_80":        _format_pct_fraction(int(getattr(row, "n_identity_80",      0)), n_total),
            "pct_identity_threshold": _format_pct_fraction(int(getattr(row, "n_passed_threshold", 0)), n_total),
            "min_identity":           f"{float(getattr(row, 'min_max_identity',  0.0)) * 100:.2f}%",
            "max_identity":           f"{float(getattr(row, 'max_max_identity',  0.0)) * 100:.2f}%",
            "avg_identity":           f"{float(getattr(row, 'mean_max_identity', 0.0)) * 100:.2f}%",
            "conservation_label":     getattr(row, "conservation_label", "conservation_unknown"),
            "n_excellent_match":      int(getattr(row, "n_excellent_match", 0)),
            "n_tolerated":            int(getattr(row, "n_tolerated", 0)),
            "variants_exact_match":   getattr(row, "variants_exact_match", "") or "",
            "variants_tolerable":     getattr(row, "variants_tolerable",   "") or "",
            "alleles_united":         getattr(row, COLUMN_ALLELES_UNITED, "") or "",
            "num_alleles_united":     num_alleles,
        })
    return pd.DataFrame(rows, columns=_SUMMARY_COLUMNS)


