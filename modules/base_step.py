"""
Base classes for all pipeline steps.

Steps 01-12 are per-track (BaseTrackStep): they run once for each sequence track.
Steps 13-14 are global (BaseGlobalStep): they run once after all tracks complete.

Contract: execute() NEVER raises. It returns a status dict so the caller
(the interactive REPL in main.py) can decide what to do next — retry, skip,
go back, etc. Exceptions are caught, logged in pipeline.json as 'error',
and surfaced via the returned dict.

A step marked 'error' in pipeline.json is NOT treated as blocking — the next
run will try it again. Only 'done' is honored as cache.
"""

import traceback
from abc import ABC, abstractmethod
from pathlib import Path

from rich.console import Console

from utils.pipeline_state import (
    get_track_step_status,
    set_track_step_status,
    reset_track_step,
    get_global_step_status,
    set_global_step_status,
)

console = Console()


# ── Return-value builders ─────────────────────────────────────────────────────

def _build_done_outcome(step_output_data) -> dict:
    return {'status': 'done', 'result': step_output_data}


def _build_error_outcome(error_message: str, error_traceback_text: str = '') -> dict:
    return {
        'status':           'error',
        'error_message':    error_message,
        'error_traceback':  error_traceback_text,
    }


def _build_skipped_outcome(skip_reason: str) -> dict:
    return {'status': 'skipped', 'reason': skip_reason}


# ── Per-track steps (01-12) ───────────────────────────────────────────────────

class BaseTrackStep(ABC):
    """
    Base for per-track steps. Runs once per sequence track with the same project config.

    Step identity is the `step_name` only — there is no numeric prefix anywhere
    (folder name, pipeline.json key, registry key all use the bare name). This
    makes adding/removing/reordering steps a registry-only edit, not a renaming
    cascade across the codebase.
    """
    step_name: str = ""

    def __init__(self, project_name: str, project_config: dict, track_id: str):
        self.project_name    = project_name
        self.project_config  = project_config
        self.track_id        = track_id
        self.project_dir     = Path("projects") / project_name
        self.track_dir       = self.project_dir / "data" / "intermediate" / track_id
        self.input_dir       = self.project_dir / "data" / "input"
        self.output_dir      = self.project_dir / "data" / "output"

    @property
    def step_key(self) -> str:
        return self.step_name

    @abstractmethod
    def run(self, input_data=None):
        """Core logic for this step. Must be implemented by each subclass."""
        pass

    def execute(self, input_data=None, force_rerun: bool = False) -> dict:
        """
        Runs the step with cache checking and state tracking.

        Returns a dict with at minimum a 'status' key:
          {'status': 'done',    'result': <whatever run() returned>}
          {'status': 'skipped', 'reason': 'already_done'}
          {'status': 'error',   'error_message': str, 'error_traceback': str}

        Never raises — callers (main.py REPL) use the status to decide next action.

        Args:
          input_data  — forwarded to run()
          force_rerun — if True, ignores cached 'done' status and runs again
        """
        cached_step_status = get_track_step_status(
            self.project_name, self.track_id, self.step_key,
        )

        if cached_step_status == 'done' and not force_rerun:
            console.print(
                f"[dim]⏭  [{self.track_id}] [{self.step_key}] Already done — skipping.[/dim]"
            )
            return _build_skipped_outcome('already_done')

        # When forcing a rerun, clear any previous state first so run() starts clean
        if force_rerun and cached_step_status == 'done':
            reset_track_step(self.project_name, self.track_id, self.step_key)

        console.print(
            f"[bold cyan]▶  [{self.track_id}] [{self.step_key}] Running...[/bold cyan]"
        )

        try:
            run_return_value = self.run(input_data)
        except Exception as caught_run_exception:
            formatted_traceback_text = traceback.format_exc()
            set_track_step_status(
                self.project_name, self.track_id, self.step_key,
                status='error', error=str(caught_run_exception),
            )
            console.print(
                f"[bold red]✗  [{self.track_id}] [{self.step_key}] "
                f"Failed: {caught_run_exception}[/bold red]"
            )
            return _build_error_outcome(
                error_message=str(caught_run_exception),
                error_traceback_text=formatted_traceback_text,
            )

        set_track_step_status(
            self.project_name, self.track_id, self.step_key,
            status='done',
            output=str(run_return_value) if run_return_value is not None else None,
        )
        console.print(
            f"[bold green]✓  [{self.track_id}] [{self.step_key}] Done.[/bold green]"
        )
        return _build_done_outcome(run_return_value)


# ── Global steps (13-14) ──────────────────────────────────────────────────────

class BaseGlobalStep(ABC):
    """
    Base for global steps that run once after all tracks complete, combining all results.
    Same naming policy as BaseTrackStep — step_name only, no numbers.
    """
    step_name: str = ""

    def __init__(self, project_name: str, project_config: dict):
        self.project_name     = project_name
        self.project_config   = project_config
        self.project_dir      = Path("projects") / project_name
        self.intermediate_dir = self.project_dir / "data" / "intermediate"
        self.output_dir       = self.project_dir / "data" / "output"

    @property
    def step_key(self) -> str:
        return self.step_name

    @abstractmethod
    def run(self, input_data=None):
        pass

    def execute(self, input_data=None, force_rerun: bool = False) -> dict:
        """Same contract as BaseTrackStep.execute() — returns a status dict, never raises."""
        cached_step_status = get_global_step_status(self.project_name, self.step_key)

        if cached_step_status == 'done' and not force_rerun:
            console.print(f"[dim]⏭  [{self.step_key}] Already done — skipping.[/dim]")
            return _build_skipped_outcome('already_done')

        console.print(f"[bold cyan]▶  [{self.step_key}] Running...[/bold cyan]")

        try:
            run_return_value = self.run(input_data)
        except Exception as caught_run_exception:
            formatted_traceback_text = traceback.format_exc()
            set_global_step_status(
                self.project_name, self.step_key,
                status='error', error=str(caught_run_exception),
            )
            console.print(
                f"[bold red]✗  [{self.step_key}] Failed: {caught_run_exception}[/bold red]"
            )
            return _build_error_outcome(
                error_message=str(caught_run_exception),
                error_traceback_text=formatted_traceback_text,
            )

        set_global_step_status(
            self.project_name, self.step_key,
            status='done',
            output=str(run_return_value) if run_return_value is not None else None,
        )
        console.print(f"[bold green]✓  [{self.step_key}] Done.[/bold green]")
        return _build_done_outcome(run_return_value)
