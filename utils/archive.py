"""Archive helpers for shipping project outputs off-machine.

Two pairs of entry points — one pair for tar.gz (Linux/native), one for zip
(Windows/WSL). Both pairs share the same structure:

  archive_project(name, dst_dir, include_predictions=False)  → tar.gz
  archive_project_zip(name, dst_dir, include_predictions=False) → .zip

  archive_step(name, step_name, dst_dir)  → tar.gz
  archive_step_zip(name, step_name, dst_dir)  → .zip

The download_ui selects the format: zip when running under WSL (so the file
is natively openable in Windows Explorer), tar.gz on pure Linux/macOS.

Both formats exclude `data/intermediate/*/predictions/` by default because
raw prediction CSVs are large and trivially re-derivable.
"""

from __future__ import annotations

import datetime
import tarfile
import zipfile
from pathlib import Path


_PROJECTS_DIR = Path("projects")

# Per-track sub-folder owned by each track-level step.
_STEP_TO_TRACK_SUBFOLDER: dict[str, str] = {
    "fetch_sequences":        "../input",
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

_GLOBAL_STEP_OUTPUT_PREFIXES: dict[str, list[str]] = {
    "integrate_data":  ["MASTER_TABLE_FULL_", "MASTER_TABLE_VIEW_", "MASTER_TABLE_AUDIT_"],
    "generate_report": ["REPORT_"],
}


def _timestamp() -> str:
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def _project_dir(project_name: str) -> Path:
    project_dir = _PROJECTS_DIR / project_name
    if not project_dir.exists():
        raise FileNotFoundError(f"Project not found: {project_dir}")
    return project_dir


def _is_predictions_path(path: Path) -> bool:
    """True if any component of the path is named 'predictions'."""
    return "predictions" in path.parts


# ── tar.gz ────────────────────────────────────────────────────────────────────

def archive_project(
    project_name:        str,
    destination_dir:     Path,
    include_predictions: bool = False,
) -> Path:
    """Bundle `projects/{project_name}/` into a tar.gz."""
    project_dir  = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path = destination_dir / f"{project_name}_full_{_timestamp()}.tar.gz"

    def _filter(tarinfo: tarfile.TarInfo) -> tarfile.TarInfo | None:
        if not include_predictions and "/predictions/" in ("/" + tarinfo.name + "/"):
            return None
        return tarinfo

    with tarfile.open(archive_path, "w:gz") as tf:
        tf.add(project_dir, arcname=project_name, filter=_filter)
    return archive_path


def archive_step(
    project_name:    str,
    step_name:       str,
    destination_dir: Path,
) -> Path:
    """Bundle a single step's outputs into a tar.gz."""
    project_dir    = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path   = destination_dir / f"{project_name}_{step_name}_{_timestamp()}.tar.gz"

    track_subfolder = _STEP_TO_TRACK_SUBFOLDER.get(step_name)
    global_prefixes = _GLOBAL_STEP_OUTPUT_PREFIXES.get(step_name, [])

    with tarfile.open(archive_path, "w:gz") as tf:
        if track_subfolder is not None:
            intermediate_root = project_dir / "data" / "intermediate"
            input_root        = project_dir / "data" / "input"
            for track_dir in sorted(intermediate_root.glob("*")):
                if not track_dir.is_dir():
                    continue
                candidate = (input_root / track_dir.name
                             if track_subfolder == "../input"
                             else track_dir / track_subfolder)
                if candidate.exists():
                    tf.add(candidate,
                           arcname=f"{project_name}/{step_name}/{track_dir.name}")
        for prefix in global_prefixes:
            for f in sorted((project_dir / "data" / "output").glob(f"{prefix}*")):
                tf.add(f, arcname=f"{project_name}/{step_name}/{f.name}")
    return archive_path


# ── zip (Windows-friendly) ────────────────────────────────────────────────────

def _add_dir_to_zip(zf: zipfile.ZipFile, src_dir: Path, arcname_prefix: str,
                    include_predictions: bool) -> None:
    """Recursively add all files in src_dir into the zipfile."""
    for file_path in sorted(src_dir.rglob("*")):
        if not file_path.is_file():
            continue
        if not include_predictions and _is_predictions_path(
                file_path.relative_to(src_dir)):
            continue
        rel = file_path.relative_to(src_dir)
        zf.write(file_path, arcname=f"{arcname_prefix}/{rel}")


def archive_project_zip(
    project_name:        str,
    destination_dir:     Path,
    include_predictions: bool = False,
) -> Path:
    """Bundle `projects/{project_name}/` into a .zip (Windows-friendly)."""
    project_dir  = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path = destination_dir / f"{project_name}_full_{_timestamp()}.zip"

    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        _add_dir_to_zip(zf, project_dir, project_name, include_predictions)
    return archive_path


def archive_step_zip(
    project_name:    str,
    step_name:       str,
    destination_dir: Path,
) -> Path:
    """Bundle a single step's outputs into a .zip (Windows-friendly)."""
    project_dir    = _project_dir(project_name)
    destination_dir.mkdir(parents=True, exist_ok=True)
    archive_path   = destination_dir / f"{project_name}_{step_name}_{_timestamp()}.zip"

    track_subfolder = _STEP_TO_TRACK_SUBFOLDER.get(step_name)
    global_prefixes = _GLOBAL_STEP_OUTPUT_PREFIXES.get(step_name, [])

    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        if track_subfolder is not None:
            intermediate_root = project_dir / "data" / "intermediate"
            input_root        = project_dir / "data" / "input"
            for track_dir in sorted(intermediate_root.glob("*")):
                if not track_dir.is_dir():
                    continue
                candidate = (input_root / track_dir.name
                             if track_subfolder == "../input"
                             else track_dir / track_subfolder)
                if candidate.exists():
                    arc_prefix = f"{project_name}/{step_name}/{track_dir.name}"
                    for file_path in sorted(candidate.rglob("*")):
                        if file_path.is_file():
                            rel = file_path.relative_to(candidate)
                            zf.write(file_path, arcname=f"{arc_prefix}/{rel}")

        for prefix in global_prefixes:
            for f in sorted((project_dir / "data" / "output").glob(f"{prefix}*")):
                zf.write(f, arcname=f"{project_name}/{step_name}/{f.name}")

    return archive_path
