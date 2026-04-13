"""
Helpers for reading and writing pipeline.json.

pipeline.json tracks progress at two levels:
  - Per track (steps 01-12): each sequence has its own step states
  - Global (steps 13-14): run once after all tracks complete
"""

import json
import datetime
from pathlib import Path

PROJECTS_DIR = Path("projects")


# ── Load / Save ───────────────────────────────────────────────────────────────

def load_pipeline_state(project_name: str) -> dict:
    path = PROJECTS_DIR / project_name / "pipeline.json"
    with open(path) as f:
        return json.load(f)


def save_pipeline_state(project_name: str, state: dict):
    path = PROJECTS_DIR / project_name / "pipeline.json"
    with open(path, "w") as f:
        json.dump(state, f, indent=2, ensure_ascii=False)


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

    # Update current_step to highest completed step number
    done_steps = [
        int(k.split("_")[0].replace("step", ""))
        for k, v in state["tracks"][track_id]["steps"].items()
        if v.get("status") == "done"
    ]
    if done_steps:
        state["tracks"][track_id]["current_step"] = max(done_steps)

    save_pipeline_state(project_name, state)


def reset_track_step(project_name: str, track_id: str, step_key: str):
    """Resets a track step so it runs again on next execute()."""
    state = load_pipeline_state(project_name)
    state["tracks"][track_id]["steps"].pop(step_key, None)
    save_pipeline_state(project_name, state)


# ── Global step helpers ───────────────────────────────────────────────────────

def get_global_step_status(project_name: str, step_key: str) -> str:
    """Returns status of a global step (step13, step14)."""
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
    Returns a summary of progress across all tracks.
    Useful for displaying status in the TUI.
    """
    state = load_pipeline_state(project_name)
    summary = {}

    for track_id, track_data in state.get("tracks", {}).items():
        summary[track_id] = {
            "current_step": track_data.get("current_step", 0),
            "done_steps": [
                k for k, v in track_data.get("steps", {}).items()
                if v.get("status") == "done"
            ],
            "error_steps": [
                k for k, v in track_data.get("steps", {}).items()
                if v.get("status") == "error"
            ],
        }

    return summary


def all_tracks_completed(project_name: str, total_track_steps: int = 12) -> bool:
    """Returns True when all tracks have completed all per-track steps."""
    state = load_pipeline_state(project_name)
    return all(
        track.get("current_step", 0) >= total_track_steps
        for track in state.get("tracks", {}).values()
    )
