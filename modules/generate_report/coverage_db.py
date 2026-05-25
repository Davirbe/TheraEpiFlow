"""generate_report — IEDB pickle reduction for in-browser coverage math.

The vaxbuilder HTML needs to compute cumulative population coverage every
time the user selects/deselects an epitope. Doing that on every click via
HTTP would defeat the offline guarantee. The cheap alternative: bake a
small `{population: {allele: frequency}}` JSON inline into the HTML and let
JS run the `1 - Π(1 - f_i)` formula directly.

This module loads the vendored IEDB pickle (3.6 MB) and reduces it to that
small JSON, restricted to (a) the populations selected in `project_config`
and (b) the alleles actually present in the project's peptides.
"""

from __future__ import annotations

from modules.population_coverage.core import (
    _build_population_freq_map,
    _load_population_database,
)


def build_coverage_db(
    populations: list[str],
    alleles_in_project: set[str],
) -> dict[str, dict[str, float]]:
    """Returns {population: {allele: per-locus-normalized frequency}}.

    Only the alleles in `alleles_in_project` are kept — the IEDB pickle has
    thousands of allele entries per population; the project only needs the
    ~27 it actually configured. Drops the size from ~3 MB to ~5 KB inline.

    Reuses `modules.population_coverage.core` for the loader and the
    per-locus normalization, so any future correction to the IEDB math
    propagates automatically.
    """
    population_db = _load_population_database()
    coverage_db: dict[str, dict[str, float]] = {}
    for population in populations:
        full_freq_map = _build_population_freq_map(population_db, population)
        restricted = {
            allele: float(frequency)
            for allele, frequency in full_freq_map.items()
            if allele in alleles_in_project
        }
        coverage_db[population] = restricted
    return coverage_db
