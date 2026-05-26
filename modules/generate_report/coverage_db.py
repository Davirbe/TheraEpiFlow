"""generate_report — IEDB pickle reduction for in-browser coverage math.

The HTML report needs to compute cumulative population coverage every time
the user selects or deselects an epitope. Doing that on every click via HTTP
would defeat the offline guarantee. The cheap alternative: bake a small JSON
inline into the HTML and let JS run the diploid formula directly.

The shape is `{population: {locus: {allele: frequency}}}` — alleles grouped
by their HLA locus (HLA-A, HLA-B, HLA-C, ...). Grouping is required because
the diploid coverage formula sums frequencies WITHIN each locus first
(`q_locus = Σ f_allele`), applies `p_locus = 1 - (1 - q_locus)²`, then
combines across loci (`1 - Π_locus (1 - p_locus)`). A flat
`{allele: freq}` map cannot reproduce this — see
`modules/population_coverage/core.py:_compute_epitope_coverage` for the
reference implementation.

The pickle has thousands of allele entries per population; the project
only needs the ~27 it actually configured. Drops the embedded JSON size
from ~3 MB to a few KB.
"""

from __future__ import annotations

from modules.population_coverage.core import (
    _build_population_freq_map,
    _get_locus,
    _load_population_database,
)


def build_coverage_db(
    populations: list[str],
    alleles_in_project: set[str],
) -> dict[str, dict[str, dict[str, float]]]:
    """Returns `{population: {locus: {allele: per-locus-normalized frequency}}}`.

    The per-locus normalization happens inside `_build_population_freq_map`
    (so each locus's frequencies sum to ≤ 1). Here we just re-bucket the
    flat output by locus before returning, and drop alleles not present in
    the project.
    """
    population_db = _load_population_database()
    coverage_db: dict[str, dict[str, dict[str, float]]] = {}
    for population in populations:
        full_freq_map = _build_population_freq_map(population_db, population)
        per_locus: dict[str, dict[str, float]] = {}
        for allele, frequency in full_freq_map.items():
            if allele not in alleles_in_project:
                continue
            locus = _get_locus(allele)
            per_locus.setdefault(locus, {})[allele] = float(frequency)
        coverage_db[population] = per_locus
    return coverage_db
