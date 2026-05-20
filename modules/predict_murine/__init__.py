"""
predict_murine step.

Re-runs NetMHCpan EL + MHCFlurry on the ★ representatives, this time against
murine H-2 alleles, to flag epitopes that are translatable to mouse models.

Module layout (split by responsibility):
    step.py     — PredictMurineStep orchestration (run / describe_outputs)
    core.py     — ★ loading, synthetic records, tier labels, aggregation
    prompts.py  — interactive murine-strain selection
"""

from .step import PredictMurineStep

__all__ = ["PredictMurineStep"]
