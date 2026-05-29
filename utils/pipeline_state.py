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
    with open(path, encoding="utf-8") as f:
        state = json.load(f)
    return _migrate_legacy_step_keys(state)


def save_pipeline_state(project_name: str, state: dict):
    path = PROJECTS_DIR / project_name / "pipeline.json"
    with open(path, "w", encoding="utf-8") as f:
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


def reset_track_steps_from(project_name: str, track_id: str, step_keys: list[str]):
    """Resets several track steps in one load/save (used by the cascade rewind).

    Pops every key in `step_keys` from the track's step map, so the next run
    re-executes them. Missing keys are ignored. Batching avoids one file write
    per step when rewinding a long downstream tail."""
    state = load_pipeline_state(project_name)
    track_steps = state.get("tracks", {}).get(track_id, {}).get("steps", {})
    for step_key in step_keys:
        track_steps.pop(step_key, None)
    save_pipeline_state(project_name, state)


# ── Intro-page memory (hybrid pre-step page) ──────────────────────────────────
#
# Each project tracks which steps the user has already seen the full pre-step
# page for. First time a step runs in a project → full intro renders.
# Subsequent runs → compact intro renders. Independent of step type (track vs
# global) because the user reads the intro once per (project, step) and the
# REPL flow is the same.

def has_seen_step_intro(project_name: str, step_key: str) -> bool:
    """True when the full pre-step page has been rendered at least once for
    this (project, step) pair. False otherwise (first-time user)."""
    state = load_pipeline_state(project_name)
    seen_steps = state.get("intro_seen", [])
    return step_key in seen_steps


def mark_step_intro_seen(project_name: str, step_key: str) -> None:
    """Persist that the user has now seen the full intro for this step.
    Idempotent. Adds the key to `intro_seen` in pipeline.json."""
    state = load_pipeline_state(project_name)
    seen_steps = state.setdefault("intro_seen", [])
    if step_key not in seen_steps:
        seen_steps.append(step_key)
        save_pipeline_state(project_name, state)


def reset_step_intro(project_name: str, step_key: str) -> None:
    """Drop the seen-intro flag so the full pre-step page renders again on the
    next run. Used by `--help` interaction tests and by users who want to
    re-read the introduction."""
    state = load_pipeline_state(project_name)
    seen_steps = state.get("intro_seen", [])
    if step_key in seen_steps:
        seen_steps.remove(step_key)
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


def reset_global_step(project_name: str, step_key: str):
    """Resets a global step so it runs again on next execute().

    Mirror of `reset_track_step` for the `global_steps` map. Used whenever a
    per-track rewind (or a track edit) invalidates the aggregated outputs that
    `integrate_data` / `generate_report` produce from every track."""
    state = load_pipeline_state(project_name)
    state.get("global_steps", {}).pop(step_key, None)
    save_pipeline_state(project_name, state)


