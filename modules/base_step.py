"""
Base class for all pipeline steps.

Steps 01-12 are per-track: they run once for each sequence track.
Steps 13-14 are global: they run once after all tracks complete.

Subclasses implement run() and set step_number + step_name.
BaseStep handles cache checking, state saving, and logging.
"""

import datetime
from abc import ABC, abstractmethod
from pathlib import Path

from rich.console import Console

from utils.pipeline_state import (
    get_track_step_status,
    set_track_step_status,
    get_global_step_status,
    set_global_step_status,
)

console = Console()


class BaseTrackStep(ABC):
    """
    Base for per-track steps (01-12).
    Runs once per sequence track with the same project config.
    """
    step_number: int = 0
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
        return f"step{self.step_number:02d}_{self.step_name}"

    @abstractmethod
    def run(self, input_data=None):
        """Core logic for this step. Must be implemented by each subclass."""
        pass

    def execute(self, input_data=None):
        """
        Runs the step with cache checking and state tracking.
        Skips if already completed. Saves state on success or failure.
        """
        status = get_track_step_status(self.project_name, self.track_id, self.step_key)

        if status == "done":
            console.print(
                f"[dim]⏭  [{self.track_id}] [{self.step_key}] Already done — skipping.[/dim]"
            )
            return None

        console.print(
            f"[bold cyan]▶  [{self.track_id}] [{self.step_key}] Running...[/bold cyan]"
        )

        try:
            result = self.run(input_data)
            set_track_step_status(
                self.project_name, self.track_id, self.step_key,
                status="done", output=str(result) if result is not None else None,
            )
            console.print(
                f"[bold green]✅ [{self.track_id}] [{self.step_key}] Done.[/bold green]"
            )
            return result

        except Exception as error:
            set_track_step_status(
                self.project_name, self.track_id, self.step_key,
                status="error", error=str(error),
            )
            console.print(
                f"[bold red]✗  [{self.track_id}] [{self.step_key}] Failed: {error}[/bold red]"
            )
            raise


class BaseGlobalStep(ABC):
    """
    Base for global steps (13-14).
    Runs once after all tracks complete, combining all results.
    """
    step_number: int = 0
    step_name: str = ""

    def __init__(self, project_name: str, project_config: dict):
        self.project_name   = project_name
        self.project_config = project_config
        self.project_dir    = Path("projects") / project_name
        self.intermediate_dir = self.project_dir / "data" / "intermediate"
        self.output_dir     = self.project_dir / "data" / "output"

    @property
    def step_key(self) -> str:
        return f"step{self.step_number:02d}_{self.step_name}"

    @abstractmethod
    def run(self, input_data=None):
        pass

    def execute(self, input_data=None):
        status = get_global_step_status(self.project_name, self.step_key)

        if status == "done":
            console.print(f"[dim]⏭  [{self.step_key}] Already done — skipping.[/dim]")
            return None

        console.print(f"[bold cyan]▶  [{self.step_key}] Running...[/bold cyan]")

        try:
            result = self.run(input_data)
            set_global_step_status(
                self.project_name, self.step_key,
                status="done", output=str(result) if result is not None else None,
            )
            console.print(f"[bold green]✅ [{self.step_key}] Done.[/bold green]")
            return result

        except Exception as error:
            set_global_step_status(
                self.project_name, self.step_key,
                status="error", error=str(error),
            )
            console.print(f"[bold red]✗  [{self.step_key}] Failed: {error}[/bold red]")
            raise
