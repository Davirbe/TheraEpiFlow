"""
fetch_sequences step.

Resolves a protein track to its sequence: searches UniProt (or loads a local
FASTA), lets the user pick among candidates, validates and writes the per-track
input FASTA plus a registry entry.

Module layout (split by responsibility):
    step.py     — FetchSequencesStep orchestration (run / describe_outputs)
    core.py     — organism normalization, UniProt search/scoring, validation
    io.py       — UniProt FASTA download + local-FASTA loading
    prompts.py  — interactive candidate selection
    render.py   — Rich candidates table

ORGANISM_ALIASES is re-exported here because utils.project_manager imports it.
"""

from .core import ORGANISM_ALIASES
from .step import FetchSequencesStep

__all__ = ["FetchSequencesStep", "ORGANISM_ALIASES"]
