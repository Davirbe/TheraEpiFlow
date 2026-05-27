"""Disposable-project seeder for the validation suite.

Creates a fully configured TheraEpiFlow project programmatically so that
every wizard is skipped at run time. Reuses the production helpers from
`utils.project_manager` — only the interactive prompts are bypassed.

Each replicate gets its own project folder under `projects/bench_*` so that
the cache from one replicate cannot bleed into another (and so wall-clock
timing reflects a cold pipeline run).

Usage as a library:
    from tests.validation.seed_project import seed_project, TrackSpec
    seed_project("bench_hpv16_e7_tr1_r01", [
        TrackSpec(organism_label="H16", protein_label="E7",
                  organism_name="HPV16", protein_name="E7",
                  seed_accession="P03129", seed_size=98, tax_id=333760,
                  variants_scope="interspecific",
                  variants_host_filter="Homo sapiens",
                  variants_family_taxid=2169595),
    ])

Usage on the command line:
    python -m tests.validation.seed_project --preset hpv16_e7_tr1 --label bench_hpv16_e7_tr1_r01
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional

# Allow `python -m tests.validation.seed_project` to import the top-level package
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from step_registry import GLOBAL_STEPS  # noqa: E402
from utils.project_manager import (  # noqa: E402
    PROJECTS_DIR,
    TRACK_INTERMEDIATE_FOLDERS,
    create_project,
    delete_project,
    load_project_config,
    save_project_config,
)
from utils.naming import build_track_id  # noqa: E402


# ── Defaults aligned with existing_scripts/validation_tools.md ────────────────

DEFAULT_HLA_PANEL: list[str] = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*27:05",
    "HLA-B*39:01", "HLA-B*40:01", "HLA-B*58:01", "HLA-B*15:01",
]
DEFAULT_PEPTIDE_LENGTHS: list[int] = [9]
DEFAULT_CONSENSUS_THRESHOLD: float = 2.0
DEFAULT_TOXICITY_THRESHOLD: float = 0.38
DEFAULT_CLUSTER_METHOD: str = "cluster_break"
DEFAULT_CLUSTER_THRESHOLD: float = 0.8
DEFAULT_CONSERVATION_THRESHOLD: float = 1.0
DEFAULT_COVERAGE_POPULATIONS: list[str] = ["World"]
DEFAULT_MURINE_STRAIN: str = "default"  # = C57BL/6 + BALB/c (see predict_murine prompts)


@dataclass
class TrackSpec:
    """Everything a project_config['tracks'][track_id] entry needs to skip every wizard."""

    organism_name:           str
    organism_label:          str
    protein_name:            str
    protein_label:           str
    seed_accession:          Optional[str] = None
    seed_size:               Optional[int] = None
    tax_id:                  Optional[int] = None
    variants_scope:          Optional[str] = None       # "intraspecific" | "interspecific"
    variants_host_filter:    Optional[str] = None       # e.g. "Homo sapiens" or None
    variants_family_taxid:   Optional[int] = None       # required when scope == "interspecific"
    input_source:            str = "uniprot"
    local_file_path:         Optional[str] = None

    def as_track_entry(self) -> dict:
        """Mirror the schema seen in projects/hpv16/project_config.json (tracks[track_id])."""
        return {k: v for k, v in asdict(self).items() if v is not None or k in (
            "local_file_path", "variants_host_filter",
        )}


# ── Track presets (one per benchmark project layout) ──────────────────────────

def _hpv16_e7_track(letter_suffix: str = "") -> TrackSpec:
    label = f"H16{letter_suffix}"
    return TrackSpec(
        organism_name="HPV16",
        organism_label=label,
        protein_name="E7",
        protein_label="E7",
        seed_accession="P03129",
        seed_size=98,
        tax_id=333760,
        variants_scope="interspecific",
        variants_host_filter="Homo sapiens",
        variants_family_taxid=2169595,
    )


def _sars_ncap_track(letter_suffix: str = "") -> TrackSpec:
    label = f"SARS2{letter_suffix}"
    return TrackSpec(
        organism_name="SARS-CoV-2",
        organism_label=label,
        protein_name="Nucleocapsid",
        protein_label="N",
        seed_accession="P0DTC9",
        seed_size=419,
        tax_id=2697049,
        variants_scope="intraspecific",
        variants_host_filter=None,
        variants_family_taxid=None,
    )


PRESETS: dict[str, list[TrackSpec]] = {
    # Experiment 1
    "hpv16_e7_tr1":       [_hpv16_e7_track()],
    "hpv16_e7_tr3":       [_hpv16_e7_track("A"), _hpv16_e7_track("B"), _hpv16_e7_track("C")],
    "sars_nucleo_tr1":    [_sars_ncap_track()],
    "sars_nucleo_tr3":    [_sars_ncap_track("A"), _sars_ncap_track("B"), _sars_ncap_track("C")],
    # Experiment 2 (clinical validation)
    "clinic":             [_hpv16_e7_track(), _sars_ncap_track()],
}


# ── Project creation ──────────────────────────────────────────────────────────

def _materialize_track_folders(project_name: str, track_ids: list[str]) -> None:
    """Reproduce the folder layout normally created by setup_project_tracks_interactive()."""
    project_dir = PROJECTS_DIR / project_name
    for track_id in track_ids:
        (project_dir / "data" / "input" / track_id).mkdir(parents=True, exist_ok=True)
        for sub in TRACK_INTERMEDIATE_FOLDERS:
            (project_dir / "data" / "intermediate" / track_id / sub).mkdir(parents=True, exist_ok=True)


def _seed_pipeline_state(project_name: str, track_ids: list[str]) -> None:
    """Initialise pipeline.json with one empty entry per track + every global step pending."""
    pipeline_path = PROJECTS_DIR / project_name / "pipeline.json"
    pipeline_state = {
        "tracks": {tid: {"current_step": 0, "steps": {}} for tid in track_ids},
        "global_steps": {name: {"status": "pending"} for name in GLOBAL_STEPS},
    }
    pipeline_path.write_text(json.dumps(pipeline_state, indent=2, ensure_ascii=False))


def seed_project(
    project_name:           str,
    tracks:                 list[TrackSpec],
    *,
    description:            str = "",
    hla_alleles:            Optional[list[str]] = None,
    peptide_lengths:        Optional[list[int]] = None,
    consensus_threshold:    Optional[float] = None,
    toxicity_threshold:     Optional[float] = None,
    cluster_method:         Optional[str] = None,
    cluster_threshold:      Optional[float] = None,
    conservation_threshold: Optional[float] = None,
    coverage_populations:   Optional[list[str]] = None,
    murine_strain:          Optional[str] = None,
    step_overrides:         Optional[dict] = None,
    overwrite:              bool = False,
) -> str:
    """Create projects/{project_name} fully configured for a headless run.

    Returns the project name. If `overwrite` is True and the project already
    exists, it is deleted first (use cautiously — irreversible).
    """
    if not tracks:
        raise ValueError("seed_project requires at least one TrackSpec.")

    if (PROJECTS_DIR / project_name).exists():
        if not overwrite:
            raise FileExistsError(
                f"Project '{project_name}' already exists. Pass overwrite=True to replace."
            )
        delete_project(project_name)

    create_project(project_name=project_name, description=description)

    # Build the tracks dict using the same track_id convention as the wizard.
    tracks_payload: dict[str, dict] = {}
    seen_ids: set[str] = set()
    for spec in tracks:
        track_id = build_track_id(spec.organism_label, spec.protein_label)
        if track_id in seen_ids:
            raise ValueError(
                f"Duplicate track_id '{track_id}' — use a different organism_label suffix "
                f"(e.g. H16A / H16B / H16C) for replicates of the same organism+protein."
            )
        seen_ids.add(track_id)
        tracks_payload[track_id] = spec.as_track_entry()

    project_config = load_project_config(project_name)
    project_config.update({
        "tracks":                 tracks_payload,
        "hla_alleles":            list(hla_alleles)            if hla_alleles            is not None else list(DEFAULT_HLA_PANEL),
        "peptide_lengths":        list(peptide_lengths)        if peptide_lengths        is not None else list(DEFAULT_PEPTIDE_LENGTHS),
        "consensus_threshold":    consensus_threshold          if consensus_threshold    is not None else DEFAULT_CONSENSUS_THRESHOLD,
        "toxicity_threshold":     toxicity_threshold           if toxicity_threshold     is not None else DEFAULT_TOXICITY_THRESHOLD,
        "cluster_method":         cluster_method               if cluster_method         is not None else DEFAULT_CLUSTER_METHOD,
        "cluster_threshold":      cluster_threshold            if cluster_threshold      is not None else DEFAULT_CLUSTER_THRESHOLD,
        "conservation_threshold": conservation_threshold       if conservation_threshold is not None else DEFAULT_CONSERVATION_THRESHOLD,
        "coverage_populations":   list(coverage_populations)   if coverage_populations   is not None else list(DEFAULT_COVERAGE_POPULATIONS),
        "coverage_minimum_pct":   None,
        "murine_strain":          murine_strain                if murine_strain          is not None else DEFAULT_MURINE_STRAIN,
        "step_overrides":         step_overrides               if step_overrides         is not None else {"integrate_data": {"view_columns": "default"}},
    })
    save_project_config(project_name, project_config)

    _materialize_track_folders(project_name, list(tracks_payload.keys()))
    _seed_pipeline_state(project_name, list(tracks_payload.keys()))
    return project_name


def destroy_project(project_name: str, missing_ok: bool = True) -> None:
    """Permanent removal — wrapper around utils.project_manager.delete_project()."""
    if not (PROJECTS_DIR / project_name).exists():
        if missing_ok:
            return
        raise FileNotFoundError(f"Project '{project_name}' not found.")
    delete_project(project_name)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preset", required=True, choices=sorted(PRESETS.keys()),
        help="Which track layout to use.",
    )
    parser.add_argument(
        "--label", required=True,
        help="Project name (folder under projects/). Convention: bench_<preset>_r<NN>.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Delete and recreate the project if it already exists.",
    )
    parser.add_argument(
        "--destroy", action="store_true",
        help="Instead of creating, permanently delete the named project.",
    )
    args = parser.parse_args()

    if args.destroy:
        destroy_project(args.label, missing_ok=True)
        print(f"Deleted projects/{args.label} (if it existed).")
        return

    name = seed_project(
        project_name=args.label,
        tracks=PRESETS[args.preset],
        overwrite=args.overwrite,
    )
    cfg = load_project_config(name)
    track_count = len(cfg["tracks"])
    print(f"Seeded projects/{name} — preset '{args.preset}' with {track_count} track(s).")


if __name__ == "__main__":
    main()
