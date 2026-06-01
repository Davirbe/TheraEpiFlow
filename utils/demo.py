"""Guided single-protein demo, launched from the main menu.

Seeds a throwaway HPV16_E7 project (98 aa, the fastest realistic target) and
runs the whole pipeline end to end, so a new user can watch every step produce
output without filling in the project wizard.

It doubles as an environment check: each step is run in order and, if one fails
(prediction models not downloaded, a service unreachable, the toxicity model
missing, the environment not activated), the run stops at that exact step and
points to TROUBLESHOOTING.md instead of failing silently further down.
"""

from __future__ import annotations

import json

import config
from step_registry import STEP_REGISTRY, GLOBAL_STEPS, _import_step_class
from utils.console import console, confirm
from utils.naming import build_track_id
from utils.project_manager import (
    PROJECTS_DIR,
    TRACK_INTERMEDIATE_FOLDERS,
    create_project,
    delete_project,
    load_project_config,
    save_project_config,
)

DEMO_PROJECT_NAME = "demo_hpv16_e7"

# One didactic track: the HPV16 E7 oncoprotein (UniProt P03129, 98 aa). Short
# enough to run the full pipeline in a couple of minutes, and well characterised
# in the literature (the run recovers the known epitope LLMGTLGIV).
_DEMO_TRACK = {
    "organism_name":        "HPV16",
    "organism_label":       "H16",
    "protein_name":         "E7",
    "protein_label":        "E7",
    "seed_accession":       "P03129",
    "seed_size":            98,
    "tax_id":               333760,
    "variants_scope":       "interspecific",
    "variants_host_filter": "Homo sapiens",
    "variants_family_taxid": 2169595,
    "input_source":         "uniprot",
    "local_file_path":      None,
}

# When a step fails, point the user at the relevant section of the troubleshooting
# guide rather than dumping a raw traceback. Keyed by step name.
_FAILURE_HINTS: dict[str, str] = {
    "fetch_sequences":  "Could not reach UniProt. Check your internet connection.",
    "search_variants":  "Could not reach UniProt. Check your internet connection.",
    "predict_binding":  "Binding prediction failed. Run 'mhcflurry-downloads fetch' once, "
                        "and confirm the IEDB NetMHCpan endpoint is reachable.",
    "predict_murine":   "Murine prediction failed. Run 'mhcflurry-downloads fetch' and "
                        "confirm the IEDB NetMHCpan endpoint is reachable.",
    "screen_toxicity":  "Toxicity screening failed. Confirm the ToxinPred3 model files are present.",
}
_DEFAULT_HINT = "See TROUBLESHOOTING.md for the common install and connectivity problems."


def seed_demo_project(overwrite: bool = True) -> str:
    """Create projects/demo_hpv16_e7 fully configured for a headless run."""
    if (PROJECTS_DIR / DEMO_PROJECT_NAME).exists() and overwrite:
        delete_project(DEMO_PROJECT_NAME)

    create_project(project_name=DEMO_PROJECT_NAME, description="Guided pipeline demo (HPV16 E7).")

    track_id = build_track_id(_DEMO_TRACK["organism_label"], _DEMO_TRACK["protein_label"])
    project_config = load_project_config(DEMO_PROJECT_NAME)
    project_config.update({
        "tracks":                 {track_id: dict(_DEMO_TRACK)},
        "hla_alleles":            list(config.DEFAULT_HLA_ALLELES),
        "peptide_lengths":        list(config.DEFAULT_PEPTIDE_LENGTHS),
        "consensus_threshold":    config.CONSENSUS_NETMHCPAN_EL_RANK_MAX_PERCENT,
        "toxicity_threshold":     config.TOXICITY_SCORE_THRESHOLD,
        "cluster_method":         "cluster_break",
        "cluster_threshold":      config.CLUSTER_IDENTITY_CUTOFF,
        "conservation_threshold": config.CONSERVATION_SIMILARITY_DEFAULT,
        "coverage_populations":   ["World"],
        "coverage_minimum_pct":   None,
        "murine_strain":          "default",
        "step_overrides":         {"integrate_data": {"view_columns": "default"}},
    })
    save_project_config(DEMO_PROJECT_NAME, project_config)

    # Reproduce the folder layout and pipeline state the wizard would create.
    project_dir = PROJECTS_DIR / DEMO_PROJECT_NAME
    (project_dir / "data" / "input" / track_id).mkdir(parents=True, exist_ok=True)
    for sub in TRACK_INTERMEDIATE_FOLDERS:
        (project_dir / "data" / "intermediate" / track_id / sub).mkdir(parents=True, exist_ok=True)

    pipeline_state = {
        "tracks":       {track_id: {"current_step": 0, "steps": {}}},
        "global_steps": {name: {"status": "pending"} for name in GLOBAL_STEPS},
    }
    (project_dir / "pipeline.json").write_text(
        json.dumps(pipeline_state, indent=2, ensure_ascii=False), encoding="utf-8",
    )
    return DEMO_PROJECT_NAME


def _run_one_step(step_name: str, scope: str, project_name: str, project_config: dict,
                  track_ids: list[str]) -> bool:
    """Run a single registry step (per track, or once for globals). Return True on success."""
    cls = _import_step_class(step_name)
    if cls is None:
        console.print(f"[red]Step '{step_name}' could not be imported.[/red]")
        return False

    runs = (
        [cls(project_name=project_name, project_config=project_config, track_id=t) for t in track_ids]
        if scope == "track"
        else [cls(project_name=project_name, project_config=project_config)]
    )
    for instance in runs:
        outcome = instance.execute(force_rerun=False, reconfigure=False)
        if outcome.get("status") == "error":
            hint = _FAILURE_HINTS.get(step_name, _DEFAULT_HINT)
            console.print(
                f"\n[bold red]The demo stopped at '{step_name}'.[/bold red]\n"
                f"[red]{outcome.get('error_message', 'Unknown error.')}[/red]\n"
                f"[yellow]{hint}[/yellow]"
            )
            return False
    return True


def run_demo() -> str | None:
    """Seed and run the demo end to end. Return the project name on success, else None."""
    console.print(
        "\n[bold cyan]Guided demo: HPV16 E7[/bold cyan]\n"
        "Runs the full pipeline on a single 98 aa protein so you can see every "
        "step produce output. It needs internet access (UniProt + IEDB) and the "
        "MHCFlurry models ('mhcflurry-downloads fetch'); it takes a couple of "
        "minutes. If a step fails it stops there and tells you what to fix.\n"
    )
    if not confirm("Run the demo now?", default=True):
        return None

    project_name = seed_demo_project(overwrite=True)
    project_config = load_project_config(project_name)
    track_ids = list(project_config["tracks"].keys())

    for step_index, (step_name, (_, _, scope)) in enumerate(STEP_REGISTRY.items(), start=1):
        console.rule(f"[bold]Step {step_index}/{len(STEP_REGISTRY)}: {step_name}[/bold]")
        if not _run_one_step(step_name, scope, project_name, project_config, track_ids):
            return None

    console.print(
        f"\n[bold green]Demo finished.[/bold green] The project '{project_name}' is ready to "
        "explore: the master table and the interactive HTML report are in its output folder."
    )
    return project_name
