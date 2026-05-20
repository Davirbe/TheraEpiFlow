"""
consensus_filter step.

Receives prediction CSVs from predict_binding and applies two filters: a
Stage-1 presentation filter (NetMHCpan EL %Rank AND MHCFlurry presentation
%ile ≤ threshold) and a Stage-2 Calis 2013 immunogenicity score (> 0).

Module layout (split by responsibility):
    step.py                  — ConsensusFilterStep orchestration (run / describe_outputs)
    core.py                  — column resolution, Stage-1 filter, Stage-2 Calis scoring
    prompts.py               — interactive threshold selection
    render.py                — Rich progressive per-stage tables
    predict_immunogenicity.py — vendored Calis 2013 implementation (DO NOT MODIFY)
"""

from .step import ConsensusFilterStep

__all__ = ["ConsensusFilterStep"]
