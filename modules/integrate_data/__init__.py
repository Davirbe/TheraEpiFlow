"""integrate_data — global step facade.

Re-exports IntegrateDataStep so STEP_REGISTRY can keep importing the module
without caring about the SRP-split internals (step / core / io / prompts).
"""

from .step import IntegrateDataStep

__all__ = ['IntegrateDataStep']
