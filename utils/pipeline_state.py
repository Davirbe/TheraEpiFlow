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

    track_entry = state["tracks"].setdefault(track_id, {"steps": {}})
    track_entry.setdefault("steps", {})[step_key] = step_entry

    # Track which step name was completed most recently — purely informational.
    if status == "done":
        track_entry["last_completed_step"] = step_key

    save_pipeline_state(project_name, state)


def reset_track_step(project_name: str, track_id: str, step_key: str):
    """Resets a track step so it runs again on next execute()."""
    state = load_pipeline_state(project_name)
    state["tracks"].get(track_id, {}).get("steps", {}).pop(step_key, None)
    save_pipeline_state(project_name, state)


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

    state.setdefault("global_steps", {})[step_key] = step_entry
    save_pipeline_state(project_name, state)


