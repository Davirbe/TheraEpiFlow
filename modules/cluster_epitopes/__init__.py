"""
cluster_epitopes step.

Groups safe epitopes into sequence-similarity clusters with NetworkX, then
writes a CLUSTER CSV + XLSX with one row per epitope and its cluster id.

Module layout (split by responsibility):
    step.py     — ClusterEpitopesStep orchestration (run / describe_outputs)
    core.py     — input resolution, similarity matrix, clustering algorithms
    io.py       — CLUSTER XLSX writer
    prompts.py  — interactive clustering-parameter selection
"""

from .step import ClusterEpitopesStep

__all__ = ["ClusterEpitopesStep"]
