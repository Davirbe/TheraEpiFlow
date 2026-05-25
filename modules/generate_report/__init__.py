"""generate_report — global step facade.

Re-exports GenerateReportStep so STEP_REGISTRY can keep importing the module
without caring about the SRP-split internals (step / core / io / coverage_db).
"""

from .step import GenerateReportStep

__all__ = ['GenerateReportStep']
