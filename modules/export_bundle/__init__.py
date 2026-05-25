"""export_bundle step.

Final pipeline step. Packages the project (or one earlier step's outputs)
into a tar.gz archive the user can hand off. Wraps `utils.archive` and
asks the user interactively for the scope (full project vs single step)
and the destination folder (in-project / ~/Downloads / Windows Downloads
when running under WSL).

Implemented as a BaseGlobalStep so it appears alongside the other steps
in `[s]` status, in `[Enter]` / `[a]` auto-run, and in the file browser
("Downloads/" section).
"""

from __future__ import annotations

import getpass
import time
from pathlib import Path
from typing import ClassVar

from rich import box
from rich.panel import Panel
from rich.text import Text

from modules.base_step import BaseGlobalStep
from utils.archive import archive_project, archive_step
from utils.console import ask, confirm, console, is_interactive_session
from utils.step_summary import print_step_summary
from step_registry import STEP_REGISTRY


def _is_running_under_wsl() -> bool:
    """Detects WSL by inspecting /proc/version for the 'microsoft' marker."""
    proc_version_path = Path('/proc/version')
    if not proc_version_path.exists():
        return False
    try:
        return 'microsoft' in proc_version_path.read_text(errors='ignore').lower()
    except OSError:
        return False


def _windows_downloads_under_wsl() -> Path | None:
    """Tries to locate the Windows-side Downloads folder when running in WSL.

    Returns None when the path cannot be resolved (no WSL, no /mnt/c, or no
    Downloads folder for the current Windows user).
    """
    if not _is_running_under_wsl():
        return None

    mnt_c_users = Path('/mnt/c/Users')
    if not mnt_c_users.is_dir():
        return None

    # First try the matching Linux username; many WSL setups mirror it.
    linux_username = getpass.getuser()
    candidate = mnt_c_users / linux_username / 'Downloads'
    if candidate.is_dir():
        return candidate

    # Otherwise pick the only human-looking user under /mnt/c/Users (skip
    # the standard Windows pseudo-accounts).
    pseudo_accounts = {'Public', 'Default', 'Default User', 'All Users', 'WsiAccount'}
    real_user_dirs  = [
        user_dir for user_dir in mnt_c_users.iterdir()
        if user_dir.is_dir() and user_dir.name not in pseudo_accounts
    ]
    if len(real_user_dirs) == 1:
        candidate = real_user_dirs[0] / 'Downloads'
        if candidate.is_dir():
            return candidate

    return None


def _candidate_destinations(project_name: str) -> list[tuple[str, Path]]:
    """Builds the ordered list of (label, path) destination options.

    The in-project destination is always present; the others appear only
    when their target folder actually exists, so the user does not see
    a Windows-Downloads option on a pure Linux box.
    """
    in_project_downloads = Path('projects') / project_name / 'downloads'
    candidates: list[tuple[str, Path]] = [
        ('In-project folder (always available)', in_project_downloads),
    ]

    user_downloads = Path.home() / 'Downloads'
    if user_downloads.is_dir():
        candidates.append((f'Your Downloads folder ({user_downloads})', user_downloads))

    windows_downloads = _windows_downloads_under_wsl()
    if windows_downloads is not None:
        candidates.append((
            f'Windows Downloads folder ({windows_downloads})', windows_downloads,
        ))

    return candidates


def _prompt_destination(project_name: str) -> Path | None:
    """Returns the chosen destination Path, or None if the user cancelled."""
    candidates = _candidate_destinations(project_name)

    if not is_interactive_session() or len(candidates) == 1:
        return candidates[0][1]

    console.print('\n[bold]Where should the archive go?[/bold]')
    for one_based_index, (label, path) in enumerate(candidates, start=1):
        console.print(f'  [cyan]{one_based_index}.[/cyan] {label}')

    raw_choice = ask(
        f'Pick a destination (1–{len(candidates)})',
        default='1',
    )
    raw_choice = raw_choice.strip()
    if not raw_choice.isdigit() or not (1 <= int(raw_choice) <= len(candidates)):
        console.print(f'[yellow]"{raw_choice}" is not a valid number — cancelled.[/yellow]')
        return None
    return candidates[int(raw_choice) - 1][1]


def _prompt_scope() -> str | None:
    """Returns 'full', 'step', or None to cancel."""
    if not is_interactive_session():
        return 'full'
    choice = ask(
        '\nArchive scope — type [f]ull project, [s]tep, or [c]ancel',
        default='f',
        choices=['f', 's', 'c'],
    )
    if choice == 'c':
        return None
    return 'full' if choice == 'f' else 'step'


def _prompt_step_name() -> str | None:
    """Lets the user pick which earlier step to archive. Returns None on cancel."""
    available_step_names = [name for name in STEP_REGISTRY if name != 'export_bundle']
    console.print('\n[bold]Available steps[/bold]')
    for one_based_index, step_name in enumerate(available_step_names, start=1):
        console.print(f'  [cyan]{one_based_index:>2}.[/cyan] {step_name}')

    raw_selection = ask('\nPick step (number)', default='').strip()
    if not raw_selection:
        return None
    if not (raw_selection.isdigit() and 1 <= int(raw_selection) <= len(available_step_names)):
        console.print(f'[yellow]"{raw_selection}" is not a valid number — cancelled.[/yellow]')
        return None
    return available_step_names[int(raw_selection) - 1]


class ExportBundleStep(BaseGlobalStep):
    step_name = 'export_bundle'

    description = (
        "Packages the project into a tar.gz the user can hand off — the "
        "pipeline's final step. Offers in-project, ~/Downloads, or the "
        "Windows Downloads folder (auto-detected under WSL) as destinations."
    )
    long_description: ClassVar[str] = (
        "Wraps `utils.archive` so the bundling logic stays in one place. "
        "On every run the step asks interactively for the scope (the whole "
        "project, or a single earlier step's outputs across all tracks) "
        "and the destination folder. Heavy `predictions/` folders are "
        "excluded by default and can be opted in.\n\n"
        "Re-running this step always generates a new timestamped archive — "
        "the previous archives stay on disk so the user can hand off the "
        "exact build they tested."
    )
    methodology: ClassVar[str] = (
        "1. Asks the user the scope (full / per-step / cancel).\n"
        "2. Asks where to save it. Always offers the in-project "
        "`downloads/` folder; appends `~/Downloads` when it exists, and "
        "the Windows-side Downloads folder when running under WSL "
        "(detected via /proc/version + /mnt/c/Users).\n"
        "3. Calls `archive_project` or `archive_step` from `utils.archive`, "
        "which writes a gzipped tar built via the stdlib `tarfile` module."
    )
    references: ClassVar[list] = []
    data_format: ClassVar[str] = (
        "No file inputs — reads the project tree directly from "
        "`projects/{project}/`. Two interactive prompts (scope + "
        "destination) and one optional confirmation for the heavy "
        "`predictions/` opt-in."
    )
    outputs_overview: ClassVar[str] = (
        "[bold]{project}_full_{stamp}.tar.gz[/bold] — full project bundle "
        "(default scope). Contains every output, intermediate, and config "
        "file in the project, excluding `predictions/` unless opted in.\n"
        "[bold]{project}_{step}_{stamp}.tar.gz[/bold] — single-step bundle "
        "(when the per-step scope is chosen). Contains only the canonical "
        "output folder for that step across every track, plus the global "
        "files for that step in `data/output/`."
    )
    tips: ClassVar[list] = [
        "Pick the in-project destination if you only need to keep a "
        "snapshot inside the repo (the standard place colleagues will find).",
        "Pick ~/Downloads or the Windows Downloads folder when you want the "
        "archive in the same place your browser puts files.",
        "Use the per-step scope to share one stage's outputs with a "
        "collaborator who doesn't need the whole pipeline.",
        "Heavy `predictions/` folders are excluded by default — opt in only "
        "when you need the raw NetMHCpan/MHCflurry CSVs for a forensic copy.",
    ]

    _last_archive_path: Path | None = None

    def describe_outputs(self) -> dict:
        if self._last_archive_path is not None and self._last_archive_path.exists():
            return {
                self._last_archive_path:
                    "tar.gz bundle just written. Use [b] → Downloads later to find earlier archives.",
            }
        return {}

    def run(self, input_data=None):
        run_start_time = time.time()

        chosen_scope = _prompt_scope()
        if chosen_scope is None:
            raise RuntimeError('Archive cancelled by user.')

        if chosen_scope == 'full':
            include_predictions = confirm(
                'Include the heavy `predictions/` folder '
                '(raw NetMHCpan/MHCflurry CSVs)?',
                default=False,
            )
            chosen_step_name = None
        else:
            include_predictions = False
            chosen_step_name = _prompt_step_name()
            if chosen_step_name is None:
                raise RuntimeError('Archive cancelled by user.')

        destination_dir = _prompt_destination(self.project_name)
        if destination_dir is None:
            raise RuntimeError('Archive cancelled by user.')

        console.print(Panel(
            Text.from_markup(
                f"[bold]Scope:[/bold] "
                f"{'full project' if chosen_scope == 'full' else 'step ' + chosen_step_name}\n"
                f"[bold]Destination:[/bold] {destination_dir}\n"
                f"[bold]Include predictions/:[/bold] "
                f"{'yes' if include_predictions else 'no'}"
            ),
            title='Export bundle',
            border_style='cyan',
            box=box.ROUNDED,
        ))

        console.print('[dim]Creating archive…[/dim]')
        if chosen_scope == 'full':
            archive_path = archive_project(
                project_name        = self.project_name,
                destination_dir     = destination_dir,
                include_predictions = include_predictions,
            )
        else:
            archive_path = archive_step(
                project_name    = self.project_name,
                step_name       = chosen_step_name,
                destination_dir = destination_dir,
            )

        self._last_archive_path = archive_path
        archive_size_bytes      = archive_path.stat().st_size

        elapsed_seconds = time.time() - run_start_time

        narrative_lines = [
            f"[bold]Scope:[/bold] "
            f"{'full project' if chosen_scope == 'full' else 'step ' + chosen_step_name}",
            f"[bold]Destination folder:[/bold] {destination_dir}",
            f"[bold]Archive size:[/bold] "
            f"{archive_size_bytes / 1024 / 1024:.1f} MB"
            f" [dim]({archive_size_bytes:,} bytes)[/dim]",
            f"[bold]Includes predictions/:[/bold] "
            f"{'yes' if include_predictions else 'no'}",
        ]
        if not destination_dir.samefile(Path('projects') / self.project_name / 'downloads'):
            narrative_lines.append(
                "[dim]Tip:[/dim] this archive lives outside the project "
                "folder — re-runs will not touch it."
            )

        print_step_summary(
            step_title      = f'Archive written for {self.project_name}',
            elapsed_seconds = elapsed_seconds,
            narrative_lines = narrative_lines,
            output_files    = [archive_path],
        )

        return {
            'archive_path':     str(archive_path),
            'scope':            chosen_scope,
            'step_name':        chosen_step_name,
            'destination_dir':  str(destination_dir),
            'size_bytes':       int(archive_size_bytes),
            'predictions_kept': bool(include_predictions),
        }

