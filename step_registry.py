"""
Step registry — the single source of truth for the pipeline's steps and order.

`STEP_REGISTRY` is an ordered dict (insertion order = pipeline order). Each entry
maps step_name → (module_path, class_name, step_type), where step_type is 'track'
(runs once per sequence track) or 'global' (runs once after all per-track steps).

To add/remove/reorder steps, edit ONLY this dict. step_name is also the folder
name under modules/ and the key in pipeline.json. The dynamic loaders below make
the split into per-role packages transparent — they import the package and pull
the class off it (the package __init__ re-exports it).
"""

from typing import Optional


STEP_REGISTRY: dict = {
    # Per-track steps (run once per sequence track)
    'fetch_sequences':         ('modules.fetch_sequences',         'FetchSequencesStep',         'track'),
    'predict_binding':         ('modules.predict_binding',         'PredictBindingStep',         'track'),
    'consensus_filter':        ('modules.consensus_filter',        'ConsensusFilterStep',        'track'),
    'screen_toxicity':         ('modules.screen_toxicity',         'ScreenToxicityStep',         'track'),
    'cluster_epitopes':        ('modules.cluster_epitopes',        'ClusterEpitopesStep',        'track'),
    'select_representatives':  ('modules.select_representatives',  'SelectRepresentativesStep', 'track'),
    'search_variants':         ('modules.search_variants',         'SearchVariantsStep',         'track'),
    'analyze_conservation':    ('modules.analyze_conservation',    'AnalyzeConservationStep',    'track'),
    'population_coverage':     ('modules.population_coverage',     'PopulationCoverageStep',     'track'),
    'predict_murine':          ('modules.predict_murine',          'PredictMurineStep',          'track'),
    'curate_murine':           ('modules.curate_murine',           'CurateMurineStep',           'track'),
    # Global steps (run once after all per-track steps complete on every track)
    'integrate_data':          ('modules.integrate_data',          'IntegrateDataStep',          'global'),
    'generate_report':         ('modules.generate_report',         'GenerateReportStep',         'global'),
    # NOTE: export_bundle was removed from STEP_REGISTRY in 2026-05.
    # Its tar.gz packaging logic lives at `utils.download_ui.offer_download_menu`
    # and is reached from the project REPL via the [z] shortcut. It is not a
    # pipeline step (computes nothing new — only packages files already on disk).
}

TRACK_STEPS:  list = [name for name, entry in STEP_REGISTRY.items() if entry[2] == 'track']
GLOBAL_STEPS: list = [name for name, entry in STEP_REGISTRY.items() if entry[2] == 'global']

# ── Step loading ──────────────────────────────────────────────────────────────

def _import_step_class(step_name: str):
    """
    Dynamically imports and returns the step class for a given step name.
    Returns None if the module exists but the class is not yet implemented
    (empty __init__.py) or the step name is not in the registry.
    """
    if step_name not in STEP_REGISTRY:
        return None
    module_path, class_name, _ = STEP_REGISTRY[step_name]
    try:
        import importlib
        module = importlib.import_module(module_path)
        return getattr(module, class_name, None)
    except ImportError:
        return None


def _step_is_implemented(step_name: str) -> bool:
    return _import_step_class(step_name) is not None


def _resolve_step_name_from_user_input(user_text: str) -> Optional[str]:
    """
    Allows users to type a unique prefix instead of the full step name.
    Returns the matched step_name, or None if no match / ambiguous.
    """
    candidates = [s for s in STEP_REGISTRY if s.startswith(user_text)]
    if len(candidates) == 1:
        return candidates[0]
    if user_text in STEP_REGISTRY:
        return user_text
    return None

