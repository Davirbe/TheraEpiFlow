"""tar.gz archive helpers for shipping projects (or step outputs) off-machine.

Two entry points, both using `tarfile` from the stdlib (no extra dep):

  - `archive_project(name, dst_dir, include_predictions=False)`
        bundles `projects/{name}/` into `{name}_full_{stamp}.tar.gz`.
        By default excludes `data/intermediate/*/predictions/` because raw
        prediction CSVs are large and reproducible — the user can opt in
        with `include_predictions=True` if they need a forensic copy.

  - `archive_step(name, step_name, dst_dir)`
        bundles a single step's outputs across every track of the project
        into `{name}_{step_name}_{stamp}.tar.gz`. Looks at the canonical
        sub-folder for that step inside each track plus, for global steps,
        the `data/output/` files that carry the step prefix.

Both return the `Path` of the archive they wrote.
"""

from __future__ import annotations

import datetime
import tarfile
from pathlib import Path


_PROJECTS_DIR = Path("projects")

# Per-track sub-folder owned by each track-level step. Used by archive_step
# to know what to grab. Keys match the canonical step_name values from
# step_registry.STEP_REGISTRY.
_STEP_TO_TRACK_SUBFOLDER: dict[str, str] = {
    "fetch_sequences":        "../input",        # special: pulls from data/input/{track}
    "predict_binding":        "predictions",
    "consensus_filter":       "consensus",
    "screen_toxicity":        "toxicity",
    "cluster_epitopes":       "clusters",
    "select_representatives": "clusters",
    "search_variants":        "variants",
    "analyze_conservation":   "conservation",
    "population_coverage":    "coverage",
    "predict_murine":         "murine",
    "curate_murine":          "murine",
}

# Global steps emit files into data/output/ matching one of these filename prefixes.
_GLOBAL_STEP_OUTPUT_PREFIXES: dict[str, list[str]] = {
    "integrate_data":  ["MASTER_TABLE_FULL_", "MASTER_TABLE_VIEW_", "MASTER_TABLE_AUDIT_"],
    "generate_report": ["REPORT_"],
}


def _timestamp() -> str:
    """Compact, sortable timestamp for filenames (no colons → safe on Windows)."""
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def _project_dir(project_name: str) -> Path:
    project_dir = _PROJECTS_DIR / project_name
    if not project_dir.exists():
        raise FileNotFoundError(f"Project not found: {project_dir}")
    return project_dir


def archive_project(
    project_name:        str,
    destination_dir:     Path,
    include_predictions: bool = False,
) -> Path:
    """Bundles `projects/{project_name}/` into a tar.gz under `destination_dir`.

    By default the heavy `data/intermediate/*/predictions/` folders are
    excluded (raw NetMHCpan/MHCflurry outputs are large and trivially
    re-derivable). Set `include_predictions=True` to keep them.
    """
    project_dir       = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path      = destination_dir / f"{project_name}_full_{_timestamp()}.tar.gz"

    def _filter(tarinfo: tarfile.TarInfo) -> tarfile.TarInfo | None:
        if not include_predictions and "/predictions/" in ("/" + tarinfo.name + "/"):
            return None
        return tarinfo

    with tarfile.open(archive_path, "w:gz") as archive_file:
        # arcname keeps the bundle rooted under the project name so the user
        # can simply `tar xzf … && cd {project_name}/` later.
        archive_file.add(project_dir, arcname=project_name, filter=_filter)

    return archive_path


def archive_step(
    project_name:    str,
    step_name:       str,
    destination_dir: Path,
) -> Path:
    """Bundles a single step's outputs across all tracks into a tar.gz.

    For per-track steps: grabs each track's canonical sub-folder for that
    step (e.g. `consensus_filter` → `data/intermediate/{track}/consensus/`).
    For global steps: grabs the matching files from `data/output/`.
    """
    project_dir       = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path      = destination_dir / f"{project_name}_{step_name}_{_timestamp()}.tar.gz"

    track_subfolder   = _STEP_TO_TRACK_SUBFOLDER.get(step_name)
    global_prefixes   = _GLOBAL_STEP_OUTPUT_PREFIXES.get(step_name, [])

    with tarfile.open(archive_path, "w:gz") as archive_file:
        if track_subfolder is not None:
            intermediate_root = project_dir / "data" / "intermediate"
            input_root        = project_dir / "data" / "input"
            for track_dir in sorted(intermediate_root.glob("*")):
                if not track_dir.is_dir():
                    continue
                if track_subfolder == "../input":
                    candidate = input_root / track_dir.name
                else:
                    candidate = track_dir / track_subfolder
                if candidate.exists():
                    archive_file.add(
                        candidate,
                        arcname=f"{project_name}/{step_name}/{track_dir.name}",
                    )

        for prefix in global_prefixes:
            output_root = project_dir / "data" / "output"
            for matching_file in sorted(output_root.glob(f"{prefix}*")):
                archive_file.add(
                    matching_file,
                    arcname=f"{project_name}/{step_name}/{matching_file.name}",
                )

    return archive_path
