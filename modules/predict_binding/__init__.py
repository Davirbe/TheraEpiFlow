"""
predict_binding step.

Runs MHC-I binding predictions on every sequence in a track using two tools:
NetMHCpan 4.1 EL (IEDB HTTP API) and MHCFlurry 2.0 (local). Writes one CSV per
tool plus an audit JSON. Requires `mhcflurry-downloads fetch` once.

Module layout (split by responsibility):
    step.py     — PredictBindingStep orchestration (run / describe_outputs)
    core.py     — NetMHCpan + MHCFlurry runners, binder-count summaries
    io.py       — FASTA loading / peptide expansion
    prompts.py  — interactive allele + peptide-length selection

The two runners are re-exported here because predict_murine reuses them.
"""

from .core import _run_mhcflurry_with_progress, _run_netmhcpan_iedb_silent
from .step import PredictBindingStep

__all__ = [
    "PredictBindingStep",
    "_run_mhcflurry_with_progress",
    "_run_netmhcpan_iedb_silent",
]
