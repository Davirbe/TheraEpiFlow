"""
population_coverage step.

Computes, for each ★ representative from select_representatives, the fraction
of one or more human populations carrying at least one of the epitope's HLA
alleles. Uses the vendored IEDB allele-frequency pickle — no network calls.

Diploid coverage model (IEDB Population Coverage tool, Bui 2006), per epitope:

    For each locus the epitope binds:
        q       = Σ frequency(allele)  over the epitope's HLAs in that locus
        p_locus = 1 - (1 - q)²              (probability of cover, diploid)
    Combine independent loci:
        coverage_pct = (1 - Π (1 - p_locus)) · 100

QUALITATIVE step — no epitopes are removed; coverage is recorded on every ★
representative for downstream filtering in the HTML report.

Module layout (split by responsibility):
    step.py     — PopulationCoverageStep orchestration (preflight/run/postflight)
    core.py     — database loading + coverage math (pure functions)
    prompts.py  — interactive population / cutoff selection
    io.py       — CSV / XLSX writers
    charts.py   — hit-chart + matrix PNGs
    render.py   — Rich console summary
"""

from .step import PopulationCoverageStep

__all__ = ["PopulationCoverageStep"]
