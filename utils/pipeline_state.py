"""
Helpers for reading and writing pipeline.json.

pipeline.json tracks progress at two levels:
  - Per track: each sequence track has its own dict of {step_name: status_record}
  - Global:    one dict for steps that run after all tracks have completed

Step keys are the bare step name (e.g. "fetch_sequences"). Older pipeline.json
files (written before 27/04/2026) used "stepNN_<name>" prefixes — those are
silently migrated on load by `_migrate_legacy_step_keys`.
"""

import json
import re
import datetime
from pathlib import Path

PROJECTS_DIR = Path("projects")

_LEGACY_STEP_KEY_PATTERN = re.compile(r'^step\d{2}_(.+)$')


# ── Load / Save ───────────────────────────────────────────────────────────────

def load_pipeline_state(project_name: str) -> dict:
    path = PROJECTS_DIR / project_name / "pipeline.json"
    with open(path) as f:
        state = json.load(f)
    return _migrate_legacy_step_keys(state)


def save_pipeline_state(project_name: str, state: dict):
    path = PROJECTS_DIR / project_name / "pipeline.json"
    with open(path, "w") as f:
        json.dump(state, f, indent=2, ensure_ascii=False)


def _migrate_legacy_step_keys(state: dict) -> dict:
    """
    Strips 'stepNN_' prefix from any legacy step keys in the loaded state.

    Pre-27/04/2026 pipeline.json files stored keys like 'step01_fetch_sequences'
    or 'step13_integrate_data'. Today step keys are just the bare step name.
    This shim is a read-time migration: when the migrated state is later saved
    via set_track_step_status / set_global_step_status, the new (bare-name)
    keys are persisted, so legacy keys disappear naturally.
    """
    for track_state in state.get("tracks", {}).values():
        legacy_steps_dict = track_state.get("steps", {})
        if any(_LEGACY_STEP_KEY_PATTERN.match(key) for key in legacy_steps_dict):
            track_state["steps"] = {
                _strip_legacy_prefix(key): value
                for key, value in legacy_steps_dict.items()
            }

    legacy_global_dict = state.get("global_steps", {})
    if legacy_global_dict and any(
        _LEGACY_STEP_KEY_PATTERN.match(key) for key in legacy_global_dict
    ):
        state["global_steps"] = {
            _strip_legacy_prefix(key): value
            for key, value in legacy_global_dict.items()
        }

    return state


def _strip_legacy_prefix(possibly_legacy_key: str) -> str:
    match = _LEGACY_STEP_KEY_PATTERN.match(possibly_legacy_key)
    return match.group(1) if match else possibly_legacy_key


# ── Track step helpers ────────────────────────────────────────────────────────

def get_track_step_status(project_name: str, track_id: str, step_key: str) -> str:
    """Returns status of a step for a specific track: 'done', 'error', or 'pending'."""
    state = load_pipeline_state(project_name)
    return (
        state.get("tracks", {})
             .get(track_id, {})
             .get("steps", {})
             .get(step_key, {})
             .get("status", "pending")
    )


def set_track_step_status(
    project_name: str,
    track_id: str,
    step_key: str,
    status: str,
    output: str = None,
    error: str = None,
):
    """Updates the status of a step for a specific track."""
    state = load_pipeline_state(project_name)
    step_entry = {"status": status}

    if status == "done":
        step_entry["completed_at"] = datetime.datetime.now().isoformat()
        if output:
            step_entry["output"] = output

    if status == "error":
        step_entry["failed_at"] = datetime.datetime.now().isoformat()
        if error:
            step_entry["error"] = error

    state["tracks"][track_id]["steps"][step_key] = step_entry

    # Track which step name was completed most recently — purely informational.
    if status == "done":
        state["tracks"][track_id]["last_completed_step"] = step_key

    save_pipeline_state(project_name, state)


def reset_track_step(project_name: str, track_id: str, step_key: str):
    """Resets a track step so it runs again on next execute()."""
    state = load_pipeline_state(project_name)
    state["tracks"][track_id]["steps"].pop(step_key, None)
    save_pipeline_state(project_name, state)


def count_done_track_steps(project_name: str, track_id: str) -> int:
    """Counts how many track steps are marked 'done' for a given track."""
    state = load_pipeline_state(project_name)
    track_steps = state.get("tracks", {}).get(track_id, {}).get("steps", {})
    return sum(1 for step_record in track_steps.values()
               if step_record.get("status") == "done")


# ── Global step helpers ───────────────────────────────────────────────────────

def get_global_step_status(project_name: str, step_key: str) -> str:
    """Returns status of a global step."""
    state = load_pipeline_state(project_name)
    return state.get("global_steps", {}).get(step_key, {}).get("status", "pending")


def set_global_step_status(
    project_name: str,
    step_key: str,
    status: str,
    output: str = None,
    error: str = None,
):
    """Updates the status of a global step."""
    state = load_pipeline_state(project_name)
    step_entry = {"status": status}

    if status == "done":
        step_entry["completed_at"] = datetime.datetime.now().isoformat()
        if output:
            step_entry["output"] = output

    if status == "error":
        step_entry["failed_at"] = datetime.datetime.now().isoformat()
        if error:
            step_entry["error"] = error

    state["global_steps"][step_key] = step_entry
    save_pipeline_state(project_name, state)


# ── Progress summary ──────────────────────────────────────────────────────────

def get_project_progress(project_name: str) -> dict:
    """
    Returns a summary of progress across all tracks. Useful for displays.
    """
    state = load_pipeline_state(project_name)
    summary = {}

    for track_id, track_data in state.get("tracks", {}).items():
        steps_dict = track_data.get("steps", {})
        summary[track_id] = {
            "last_completed_step": track_data.get("last_completed_step", None),
            "done_steps": [
                name for name, value in steps_dict.items()
                if value.get("status") == "done"
            ],
            "error_steps": [
                name for name, value in steps_dict.items()
                if value.get("status") == "error"
            ],
        }

    return summary


def all_tracks_completed(project_name: str, expected_track_step_names: list) -> bool:
    """
    Returns True when every track has 'done' status for every expected step name.

    The caller passes the canonical list of per-track step names (typically
    main.TRACK_STEPS) — pipeline_state.py stays free of any step registry knowledge.
    """
    state = load_pipeline_state(project_name)
    expected_set = set(expected_track_step_names)
    for track_data in state.get("tracks", {}).values():
        done_set = {
            name for name, value in track_data.get("steps", {}).items()
            if value.get("status") == "done"
        }
        if not expected_set.issubset(done_set):
            return False
    return True
