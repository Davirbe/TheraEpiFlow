"""
search_variants step.

Retrieves protein-sequence variants of a track's reference from UniProt
(species-wide or genus-wide scope, optional host filter), validates them and
writes the variants FASTA + metadata used by analyze_conservation.

Module layout (split by responsibility):
    step.py     — SearchVariantsStep orchestration (run / describe_outputs)
    core.py     — HTTP helper, taxonomy, UniProt search, identity, record building
    io.py       — reference loading + empty-output writer
    prompts.py  — cache-redo, scope/host filter, candidate multi-selection
    render.py   — Rich candidates table
"""

from .step import SearchVariantsStep

__all__ = ["SearchVariantsStep"]
