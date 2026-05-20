"""
Domain logic for population_coverage: database loading, allele-frequency maps
and the IEDB diploid coverage math. Pure functions — no terminal, no file I/O.
"""

import pickle
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd

# ── Database location ─────────────────────────────────────────────────────────

_DB_PATH    = Path(__file__).parent / "data" / "population_genotype_map.p"
_MHC_CLASS  = "I"  # this pipeline is MHC-I only

# Coverage thresholds for colour bands (purely presentational).
_COVERAGE_HIGH_BAND     = 50.0
_COVERAGE_MODERATE_BAND = 10.0


# ── Database loading ──────────────────────────────────────────────────────────

def _load_population_database() -> dict:
    """
    Loads the vendored IEDB allele-frequency pickle, returning only the
    population_coverage dict:

        { mhc_class: { population: { locus: [(allele, freq), ...] } } }

    The two trailing tables (country_ethnicity, ethnicity) are sequenced in
    the same pickle stream and must be read to reach them, then discarded.
    """
    with open(_DB_PATH, "rb") as fh:
        population_db    = pickle.load(fh)
        _country_eth_map = pickle.load(fh)  # noqa: F841 — sequenced in pickle stream
        _ethnicity_map   = pickle.load(fh)  # noqa: F841
    return population_db


def _list_available_populations(population_db: dict) -> list[str]:
    """Sorted list of populations available for class I."""
    return sorted(population_db.get(_MHC_CLASS, {}).keys())


# ── Frequency map + coverage math ─────────────────────────────────────────────

def _get_locus(hla_allele: str) -> str:
    """Returns the locus prefix (e.g. 'HLA-A' from 'HLA-A*01:01')."""
    return hla_allele.split("*", 1)[0]


def _build_population_freq_map(population_db: dict, population: str) -> dict[str, float]:
    """
    Builds a {allele: per-locus-normalized frequency} map for the given
    population. If a locus's raw frequencies sum to > 1 in the source data,
    each frequency is divided by the locus sum (so the locus totals to 1).
    Otherwise raw frequencies are kept (already gametic).
    """
    population_block = population_db.get(_MHC_CLASS, {}).get(population, {})
    freq_map: dict[str, float] = {}
    for _locus, allele_freq_pairs in population_block.items():
        locus_total = sum(freq for _, freq in allele_freq_pairs)
        for allele, freq in allele_freq_pairs:
            freq_map[allele] = (freq / locus_total) if locus_total > 1 else freq
    return freq_map


def _compute_epitope_coverage(
    epitope_alleles: set[str],
    freq_map:        dict[str, float],
) -> tuple[float, set[str]]:
    """
    Returns (coverage_pct, matched_alleles).

    coverage_pct is in [0, 100].
    matched_alleles is the subset of epitope_alleles found in freq_map.
    """
    locus_frequency_sum: dict[str, float] = defaultdict(float)
    matched_alleles: set[str] = set()
    for allele in epitope_alleles:
        if allele in freq_map:
            locus_frequency_sum[_get_locus(allele)] += freq_map[allele]
            matched_alleles.add(allele)

    if not locus_frequency_sum:
        return 0.0, matched_alleles

    not_covered_product = 1.0
    for q in locus_frequency_sum.values():
        q = min(q, 1.0)
        p_locus              = 1.0 - (1.0 - q) ** 2
        not_covered_product *= (1.0 - p_locus)

    coverage_pct = round((1.0 - not_covered_product) * 100.0, 2)
    return coverage_pct, matched_alleles


def _parse_alleles_united(raw_value: object) -> set[str]:
    """Splits the semicolon-joined alleles_united field. Returns an empty set
    when the value is missing or NaN."""
    if pd.isna(raw_value):
        return set()
    return {token.strip() for token in str(raw_value).split(";") if token.strip()}


def _safe_filename_token(value: str) -> str:
    """Sanitises a population name for use in a filename."""
    cleaned: list[str] = []
    for character in value.strip():
        if character.isalnum() or character in {"_", "-"}:
            cleaned.append(character)
        elif character == " ":
            cleaned.append("_")
        # every other character is dropped
    safe_value = "".join(cleaned)
    return safe_value or "Population"


def _coverage_band(coverage_pct: float, cutoff: Optional[float]) -> str:
    """Returns one of 'high', 'moderate', 'low', 'unknown' for colouring."""
    if coverage_pct <= 0.0:
        return "unknown"
    if coverage_pct >= _COVERAGE_HIGH_BAND:
        return "high"
    effective_min = cutoff if cutoff is not None else _COVERAGE_MODERATE_BAND
    if coverage_pct >= effective_min:
        return "moderate"
    return "low"
