"""Interactive download menu for a finished project.

This used to be the final pipeline step (`modules/export_bundle/`). It was
moved here in 2026-05 because the action is *not* a step that computes
anything — it just packages files that already exist on disk. Treating it as
a 14th step put it in the progress bar and made the pipeline feel unfinished
even after `generate_report` had produced every paper artifact.

The menu is reached from the project REPL via the `[z]` shortcut. It asks
the user for scope (full project / single earlier step), an optional opt-in
for the heavy `predictions/` folder, and a destination (in-project / Linux
~/Downloads / Windows-side Downloads under WSL — auto-detected via
/proc/version + /mnt/c/Users). Heavy lifting (the gzipped tar itself) is
delegated to `utils.archive`.
"""

from __future__ import annotations

import getpass
import time
from pathlib import Path
from typing import Optional

from rich import box
from rich.panel import Panel
from rich.text import Text

from step_registry import STEP_REGISTRY
from utils.archive import archive_project, archive_step
from utils.console import ask, confirm, console, is_interactive_session


# ── WSL detection ─────────────────────────────────────────────────────────────

def _is_running_under_wsl() -> bool:
    """True when /proc/version contains 'microsoft' (the WSL marker)."""
    proc_version_path = Path("/proc/version")
    if not proc_version_path.exists():
        return False
    try:
        return "microsoft" in proc_version_path.read_text(errors="ignore").lower()
    except OSError:
        return False


def _windows_downloads_under_wsl() -> Optional[Path]:
    """Locate the Windows-side Downloads folder when running in WSL.

    Returns None when the path cannot be resolved (no WSL, no /mnt/c, no
    Downloads folder for the current Windows user).
    """
    if not _is_running_under_wsl():
        return None

    mnt_c_users = Path("/mnt/c/Users")
    if not mnt_c_users.is_dir():
        return None

    # First try the matching Linux username — many WSL setups mirror it.
    linux_username = getpass.getuser()
    candidate = mnt_c_users / linux_username / "Downloads"
    if candidate.is_dir():
        return candidate

    # Otherwise pick the only human-looking user (skip Windows pseudo-accounts).
    pseudo_accounts = {"Public", "Default", "Default User", "All Users", "WsiAccount"}
    real_user_dirs = [
        user_dir for user_dir in mnt_c_users.iterdir()
        if user_dir.is_dir() and user_dir.name not in pseudo_accounts
    ]
    if len(real_user_dirs) == 1:
        candidate = real_user_dirs[0] / "Downloads"
        if candidate.is_dir():
            return candidate

    return None


# ── Destination + scope prompts ───────────────────────────────────────────────

def _candidate_destinations(project_name: str) -> list[tuple[str, Path]]:
    """Ordered list of (label, path) destinations.

    The in-project folder is always present. ~/Downloads and the WSL-side
    Windows Downloads only appear when they actually exist, so a Linux box
    doesn't see a Windows option.
    """
    in_project_downloads = Path("projects") / project_name / "downloads"
    candidates: list[tuple[str, Path]] = [
        ("In-project folder (always available)", in_project_downloads),
    ]

    user_downloads = Path.home() / "Downloads"
    if user_downloads.is_dir():
        candidates.append((f"Your Downloads folder ({user_downloads})", user_downloads))

    windows_downloads = _windows_downloads_under_wsl()
    if windows_downloads is not None:
        candidates.append((
            f"Windows Downloads folder ({windows_downloads})", windows_downloads,
        ))
    return candidates


def _prompt_destination(project_name: str) -> Optional[Path]:
    candidates = _candidate_destinations(project_name)
    if not is_interactive_session() or len(candidates) == 1:
        return candidates[0][1]

    console.print("\n[bold]Where should the archive go?[/bold]")
    for one_based_index, (label, _path) in enumerate(candidates, start=1):
        console.print(f"  [cyan]{one_based_index}.[/cyan] {label}")

    raw_choice = ask(
        f"Pick a destination (1–{len(candidates)})",
        default="1",
    ).strip()
    if not raw_choice.isdigit() or not (1 <= int(raw_choice) <= len(candidates)):
        console.print(f'[yellow]"{raw_choice}" is not a valid number — cancelled.[/yellow]')
        return None
    return candidates[int(raw_choice) - 1][1]


def _prompt_scope() -> Optional[str]:
    """Returns 'full', 'step', or None to cancel."""
    if not is_interactive_session():
        return "full"
    choice = ask(
        "\nArchive scope — type [f]ull project, [s]pecific step, or [c]ancel",
        default="f",
        choices=["f", "s", "c"],
    )
    if choice == "c":
        return None
    return "full" if choice == "f" else "step"


def _prompt_step_name() -> Optional[str]:
    """Lets the user pick which earlier step to archive. Returns None on cancel."""
    # Every step in the registry is now a candidate (export_bundle no longer
    # appears here since it was removed from STEP_REGISTRY).
    available_step_names = list(STEP_REGISTRY.keys())
    console.print("\n[bold]Available steps[/bold]")
    for one_based_index, step_name in enumerate(available_step_names, start=1):
        console.print(f"  [cyan]{one_based_index:>2}.[/cyan] {step_name}")

    raw_selection = ask("\nPick step (number)", default="").strip()
    if not raw_selection:
        return None
    if not (raw_selection.isdigit() and 1 <= int(raw_selection) <= len(available_step_names)):
        console.print(f'[yellow]"{raw_selection}" is not a valid number — cancelled.[/yellow]')
        return None
    return available_step_names[int(raw_selection) - 1]


# ── Public API ────────────────────────────────────────────────────────────────

def offer_download_menu(project_name: str) -> Optional[Path]:
    """Run the interactive download UI for the given project.

    Returns the Path of the archive that was written, or None when the user
    cancelled. Side effects are all visible in the Rich console.
    """
    run_start_time = time.time()

    console.print(Panel(
        Text.from_markup(
            f"[bold]Project:[/bold] {project_name}\n"
            f"Bundle the project (or a single step's outputs) into a tar.gz "
            f"archive you can hand off, archive, or send to a collaborator."
        ),
        title="[bold]Download project[/bold]",
        border_style="cyan",
        box=box.ROUNDED,
    ))

    chosen_scope = _prompt_scope()
    if chosen_scope is None:
        console.print("[dim]Download cancelled.[/dim]")
        return None

    if chosen_scope == "full":
        include_predictions = confirm(
            "Include the heavy `predictions/` folder "
            "(raw NetMHCpan/MHCflurry CSVs)?",
            default=False,
        )
        chosen_step_name = None
    else:
        include_predictions = False
        chosen_step_name = _prompt_step_name()
        if chosen_step_name is None:
            console.print("[dim]Download cancelled.[/dim]")
            return None

    destination_dir = _prompt_destination(project_name)
    if destination_dir is None:
        console.print("[dim]Download cancelled.[/dim]")
        return None

    console.print(Panel(
        Text.from_markup(
            f"[bold]Scope:[/bold] "
            f"{'full project' if chosen_scope == 'full' else 'step ' + chosen_step_name}\n"
            f"[bold]Destination:[/bold] {destination_dir}\n"
            f"[bold]Include predictions/:[/bold] "
            f"{'yes' if include_predictions else 'no'}"
        ),
        title="Confirm archive",
        border_style="cyan",
        box=box.ROUNDED,
    ))

    console.print("[dim]Creating archive…[/dim]")
    if chosen_scope == "full":
        archive_path = archive_project(
            project_name        = project_name,
            destination_dir     = destination_dir,
            include_predictions = include_predictions,
        )
    else:
        archive_path = archive_step(
            project_name    = project_name,
            step_name       = chosen_step_name,
            destination_dir = destination_dir,
        )

    elapsed_seconds = time.time() - run_start_time
    archive_size_bytes = archive_path.stat().st_size

    console.print(Panel(
        Text.from_markup(
            f"[bold green]✓ Archive written[/bold green]  "
            f"[dim]({elapsed_seconds:.1f}s)[/dim]\n\n"
            f"[bold]File:[/bold] {archive_path.name}\n"
            f"[bold]Folder:[/bold] {archive_path.parent}\n"
            f"[bold]Size:[/bold] {archive_size_bytes / 1024 / 1024:.1f} MB "
            f"[dim]({archive_size_bytes:,} bytes)[/dim]\n"
            f"[bold]Includes predictions/:[/bold] "
            f"{'yes' if include_predictions else 'no'}"
        ),
        title=f"Download ready — {project_name}",
        border_style="green",
        box=box.HEAVY,
    ))

    return archive_path
