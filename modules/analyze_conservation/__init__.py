"""
analyze_conservation step.

For each ★ representative, measures how conserved the epitope is across the
track's variants (sliding-window identity) and classifies the MHC-I mutation
tolerance (P2 + PΩ anchors + BLOSUM62). Qualitative — never removes epitopes.

Module layout (split by responsibility):
    step.py     — AnalyzeConservationStep orchestration (preflight/run/postflight)
    core.py     — BLOSUM62 scoring, anchor verdicts, identity, position stats
    io.py       — FASTA loading + XLSX writers (conservation palette)
    charts.py   — dual-panel conservation PNG
    prompts.py  — threshold + local-FASTA override prompts
    render.py   — Rich conservation + FASTA-status tables
"""

from .step import AnalyzeConservationStep

__all__ = ["AnalyzeConservationStep"]
