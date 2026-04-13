"""
Project lifecycle management.

A project contains one or more sequence tracks (e.g. HPV16_E1, HPV18_E6).
All tracks share the same pipeline configuration (alleles, thresholds, etc.)
but each track progresses independently through steps 01-12.
Steps 13-14 are global and run once after all tracks complete.

File structure per project:
  projects/{project_name}/
    project_config.json   → shared settings for all tracks
    pipeline.json         → per-track progress + global step states
    data/
      input/              → one FASTA per track: {TRACK_ID}.fasta
      intermediate/
        {track_id}/       → all intermediate files for that track
          predictions/    → PRED_NET, PRED_FLURRY, PROC_NET, PROC_FLURRY CSVs
          consensus/      → CONSENSUS CSV
          clusters/       → CLUSTER_IEDB CSV + CLUSTER_REPR XLSX
          toxicity/       → TOXICITY_ALL + TOXICITY_SAFE CSVs
          variants/       → VARIANTS FASTA (for conservation)
          conservation/   → CONSERVATION CSV
          coverage/       → COVERAGE XLSX
          murine/         → MURINE_PRED CSV + MURINE_CURATED XLSX
      output/             → master_table.xlsx, report.html
"""

import json
import shutil
import datetime
from pathlib import Path

from utils.naming import build_track_id

PROJECTS_DIR = Path("projects")
REGISTRY_FILE = PROJECTS_DIR / "projects_registry.json"

TRACK_INTERMEDIATE_FOLDERS = [
    "predictions",
    "consensus",
    "clusters",
    "toxicity",
    "variants",
    "conservation",
    "coverage",
    "murine",
]

TOTAL_TRACK_STEPS  = 12
TOTAL_GLOBAL_STEPS = 2


# ── Registry helpers ──────────────────────────────────────────────────────────

def _ensure_projects_dir():
    """Creates projects/ and projects_registry.json if they don't exist."""
    PROJECTS_DIR.mkdir(exist_ok=True)
    if not REGISTRY_FILE.exists():
        with open(REGISTRY_FILE, "w") as f:
            json.dump({"projects": {}}, f, indent=2)


def _load_registry() -> dict:
    _ensure_projects_dir()
    with open(REGISTRY_FILE) as f:
        return json.load(f)


def _save_registry(registry: dict):
    with open(REGISTRY_FILE, "w") as f:
        json.dump(registry, f, indent=2, ensure_ascii=False)


# ── Project creation ──────────────────────────────────────────────────────────

def create_project(
    project_name: str,
    sequences: list[dict],
    alleles: list[str],
    peptide_lengths: list[int],
    description: str = "",
    entrez_email: str = "",
    include_murine: bool = False,
    murine_strain: str = "complete",
) -> Path:
    """
    Creates a new project with all required folders and config files.

    Each entry in `sequences` must contain:
      - genotype: str  (e.g. "HPV16", "ZIKV", "SARS-CoV-2")
      - protein:  str  (e.g. "E1", "ENVELOPE", "SPIKE")
      - accession: str (e.g. "NP_041323.1")  — optional if uploading FASTA

    track_id is auto-generated as "{GENOTYPE}_{PROTEIN}".

    Returns the project directory path.
    """
    _ensure_projects_dir()

    project_dir = PROJECTS_DIR / project_name
    if project_dir.exists():
        raise FileExistsError(f"Project '{project_name}' already exists.")

    # Build track list with auto-generated track_ids
    tracks = []
    for seq in sequences:
        track_id = build_track_id(seq["genotype"], seq["protein"])
        tracks.append({
            "track_id":  track_id,
            "genotype":  seq["genotype"].upper(),
            "protein":   seq["protein"].upper(),
            "accession": seq.get("accession", ""),
        })

    # Create shared folders
    for folder in ["data/input", "data/output"]:
        (project_dir / folder).mkdir(parents=True, exist_ok=True)

    # Create per-track intermediate folders
    for track in tracks:
        for folder in TRACK_INTERMEDIATE_FOLDERS:
            (project_dir / "data" / "intermediate" / track["track_id"] / folder).mkdir(
                parents=True, exist_ok=True
            )

    # Write project_config.json
    project_config = {
        "project_name":  project_name,
        "description":   description,
        "created_at":    datetime.datetime.now().isoformat(),
        "entrez_email":  entrez_email,
        "sequences":     tracks,
        "prediction": {
            "alleles":                 alleles,
            "peptide_lengths":         peptide_lengths,
            "netmhc_el_rank_cutoff":   0.5,
            "flurry_percentile_cutoff": 2.0,
        },
        "clustering": {
            "identity_cutoff": 0.9,
        },
        "toxicity": {
            "threshold": 0.6,
        },
        "conservation": {
            "similarity_cutoff": 1.0,
            "search_types": ["strain", "isolate", "genotype"],
        },
        "coverage": {
            "populations": ["World"],
            "mhc_class":   "I",
        },
        "murine": {
            "enabled": include_murine,
            "strain":  murine_strain,
        },
    }

    with open(project_dir / "project_config.json", "w") as f:
        json.dump(project_config, f, indent=2, ensure_ascii=False)

    # Write pipeline.json
    pipeline_state = {
        "tracks": {
            track["track_id"]: {
                "current_step": 0,
                "steps": {},
            }
            for track in tracks
        },
        "global_steps": {
            "step13_integrate_data":  {"status": "pending"},
            "step14_generate_report": {"status": "pending"},
        },
    }

    with open(project_dir / "pipeline.json", "w") as f:
        json.dump(pipeline_state, f, indent=2, ensure_ascii=False)

    # Register project
    registry = _load_registry()
    registry["projects"][project_name] = {
        "created_at":  project_config["created_at"],
        "last_used":   project_config["created_at"],
        "description": description,
        "track_count": len(tracks),
        "status":      "in_progress",
    }
    _save_registry(registry)

    return project_dir


# ── Project loading ───────────────────────────────────────────────────────────

def load_project_config(project_name: str) -> dict:
    """Returns the project_config for an existing project."""
    config_path = PROJECTS_DIR / project_name / "project_config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Project '{project_name}' not found.")
    with open(config_path) as f:
        return json.load(f)


def save_project_config(project_name: str, config: dict):
    """Saves updated project_config (e.g. after user changes a setting)."""
    config_path = PROJECTS_DIR / project_name / "project_config.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)


def update_last_used(project_name: str):
    registry = _load_registry()
    if project_name in registry["projects"]:
        registry["projects"][project_name]["last_used"] = (
            datetime.datetime.now().isoformat()
        )
        _save_registry(registry)


# ── Project listing ───────────────────────────────────────────────────────────

def list_projects() -> list[dict]:
    """Returns all projects with metadata for display in the TUI."""
    registry = _load_registry()
    projects = []

    for name, meta in registry["projects"].items():
        pipeline_path = PROJECTS_DIR / name / "pipeline.json"
        completed_tracks = 0
        total_tracks = meta.get("track_count", 0)

        if pipeline_path.exists():
            with open(pipeline_path) as f:
                state = json.load(f)
            completed_tracks = sum(
                1 for track in state.get("tracks", {}).values()
                if track.get("current_step", 0) >= TOTAL_TRACK_STEPS
            )

        projects.append({
            "name":             name,
            "description":      meta.get("description", ""),
            "status":           meta.get("status", "unknown"),
            "created_at":       meta.get("created_at", ""),
            "last_used":        meta.get("last_used", ""),
            "track_count":      total_tracks,
            "completed_tracks": completed_tracks,
        })

    return projects


# ── Track utilities ───────────────────────────────────────────────────────────

def get_track_dir(project_name: str, track_id: str) -> Path:
    """Returns the intermediate data directory for a specific track."""
    return PROJECTS_DIR / project_name / "data" / "intermediate" / track_id


def get_tracks_by_protein(project_name: str, protein: str) -> list[dict]:
    """
    Returns all tracks in a project that share the same protein type.
    Useful for conservation: compare HPV16_E1 against all other E1 tracks.
    """
    config = load_project_config(project_name)
    return [
        seq for seq in config["sequences"]
        if seq["protein"].upper() == protein.upper()
    ]


# ── Project deletion ──────────────────────────────────────────────────────────

def delete_project(project_name: str):
    """Permanently deletes a project folder and removes it from the registry."""
    project_dir = PROJECTS_DIR / project_name
    if not project_dir.exists():
        raise FileNotFoundError(f"Project '{project_name}' not found.")

    shutil.rmtree(project_dir)

    registry = _load_registry()
    registry["projects"].pop(project_name, None)
    _save_registry(registry)
