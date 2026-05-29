"""Interactive download menu for a finished project.

This used to be the final pipeline step (`modules/export_bundle/`). It was
moved here in 2026-05 because the action is *not* a step that computes
anything — it just packages files that already exist on disk.

The menu is reached from the project REPL via the `[z]` shortcut. It asks
for scope (full project / single step), an optional opt-in for the heavy
`predictions/` folder, and a destination (in-project / Linux ~/Downloads /
Windows-side Downloads — auto-detected via /proc/version + WSLENV + /mnt/c).

Format selection is automatic:
  - WSL (running under Windows Subsystem for Linux) → .zip (natively openable
    in Windows Explorer with a double-click)
  - Pure Linux / macOS → .tar.gz
"""

from __future__ import annotations

import getpass
import os
import time
from pathlib import Path
from typing import Optional

from rich import box
from rich.panel import Panel
from rich.text import Text

from step_registry import STEP_REGISTRY
from utils.archive import (
    archive_project, archive_project_zip,
    archive_step,    archive_step_zip,
)
from utils.console import ask, confirm, console, is_interactive_session


# ── Environment detection ──────────────────────────────────────────────────────

def _is_running_under_wsl() -> bool:
    """True when the process is running inside WSL (WSL1 or WSL2).

    Three independent checks — any one suffices:
      1. WSLENV environment variable is set (WSL2 always sets it).
      2. /proc/version contains 'microsoft' (reliable on both WSL1 and WSL2).
      3. /proc/sys/fs/binfmt_misc/WSLInterop exists (WSL2 interop marker).
    """
    if os.environ.get("WSLENV") is not None:
        return True
    proc_version = Path("/proc/version")
    if proc_version.exists():
        try:
            if "microsoft" in proc_version.read_text(errors="ignore").lower():
                return True
        except OSError:
            pass
    if Path("/proc/sys/fs/binfmt_misc/WSLInterop").exists():
        return True
    return False


def _windows_downloads_under_wsl() -> Optional[Path]:
    """Locate the Windows-side Downloads folder when running in WSL.

    Returns None when the path cannot be resolved.
    """
    if not _is_running_under_wsl():
        return None

    mnt_c_users = Path("/mnt/c/Users")
    if not mnt_c_users.is_dir():
        return None

    # Try matching Linux username first (common WSL setup).
    candidate = mnt_c_users / getpass.getuser() / "Downloads"
    if candidate.is_dir():
        return candidate

    # Scan for a single real Windows user account.
    pseudo = {"Public", "Default", "Default User", "All Users", "WsiAccount"}
    real_users = [d for d in mnt_c_users.iterdir()
                  if d.is_dir() and d.name not in pseudo]
    if len(real_users) == 1:
        candidate = real_users[0] / "Downloads"
        if candidate.is_dir():
            return candidate

    # Multiple Windows accounts — list them so the user can pick manually.
    return None


# ── Destination + scope prompts ───────────────────────────────────────────────

def _candidate_destinations(project_name: str) -> list[tuple[str, Path]]:
    in_project = Path("projects") / project_name / "downloads"
    candidates: list[tuple[str, Path]] = [
        ("In-project folder (always available)", in_project),
    ]

    linux_dl = Path.home() / "Downloads"
    if linux_dl.is_dir():
        candidates.append((f"Linux Downloads ({linux_dl})", linux_dl))

    win_dl = _windows_downloads_under_wsl()
    if win_dl is not None:
        candidates.append((
            f"Windows Downloads ({win_dl})  [recommended on WSL]", win_dl,
        ))
    return candidates


def _prompt_destination(project_name: str) -> Optional[Path]:
    candidates = _candidate_destinations(project_name)

    # In non-interactive mode or only one option → pick automatically.
    if not is_interactive_session() or len(candidates) == 1:
        return candidates[0][1]

    # Default to Windows Downloads when available (last entry if WSL).
    default_idx = len(candidates)  # 1-based

    console.print("\n[bold]Where should the archive go?[/bold]")
    for i, (label, _path) in enumerate(candidates, start=1):
        marker = " [dim](default)[/dim]" if i == default_idx else ""
        console.print(f"  [cyan]{i}.[/cyan] {label}{marker}")

    raw = ask(
        f"Pick a destination (1–{len(candidates)})",
        default=str(default_idx),
    ).strip()
    if not raw.isdigit() or not (1 <= int(raw) <= len(candidates)):
        console.print(f'[yellow]"{raw}" is not a valid number; cancelled.[/yellow]')
        return None
    return candidates[int(raw) - 1][1]


def _prompt_scope() -> Optional[str]:
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
    available = list(STEP_REGISTRY.keys())
    console.print("\n[bold]Available steps[/bold]")
    for i, name in enumerate(available, start=1):
        console.print(f"  [cyan]{i:>2}.[/cyan] {name}")

    raw = ask("\nPick step (number)", default="").strip()
    if not raw:
        return None
    if not (raw.isdigit() and 1 <= int(raw) <= len(available)):
        console.print(f'[yellow]"{raw}" is not a valid number; cancelled.[/yellow]')
        return None
    return available[int(raw) - 1]


# ── Public API ────────────────────────────────────────────────────────────────

def offer_download_menu(project_name: str) -> Optional[Path]:
    """Run the interactive download UI.

    Returns the Path of the archive written, or None on cancel.
    Format: .zip when running under WSL (Windows-friendly), .tar.gz otherwise.
    """
    run_start = time.time()
    use_zip   = _is_running_under_wsl()
    fmt_label = ".zip (Windows)" if use_zip else ".tar.gz (Linux)"

    console.print(Panel(
        Text.from_markup(
            f"[bold]Project:[/bold] {project_name}\n"
            f"Bundle the project (or a single step's outputs) into an archive "
            f"you can hand off, archive, or send to a collaborator.\n"
            f"[dim]Format: {fmt_label}[/dim]"
        ),
        title="[bold]Download project[/bold]",
        border_style="cyan",
        box=box.ROUNDED,
    ))

    scope = _prompt_scope()
    if scope is None:
        console.print("[dim]Download cancelled.[/dim]")
        return None

    if scope == "full":
        include_predictions = confirm(
            "Include the heavy `predictions/` folder "
            "(raw NetMHCpan/MHCflurry CSVs)?",
            default=False,
        )
        step_name = None
    else:
        include_predictions = False
        step_name = _prompt_step_name()
        if step_name is None:
            console.print("[dim]Download cancelled.[/dim]")
            return None

    dest_dir = _prompt_destination(project_name)
    if dest_dir is None:
        console.print("[dim]Download cancelled.[/dim]")
        return None

    console.print(Panel(
        Text.from_markup(
            f"[bold]Scope:[/bold] "
            f"{'full project' if scope == 'full' else 'step ' + step_name}\n"
            f"[bold]Destination:[/bold] {dest_dir}\n"
            f"[bold]Format:[/bold] {fmt_label}\n"
            f"[bold]Include predictions/:[/bold] "
            f"{'yes' if include_predictions else 'no'}"
        ),
        title="Confirm archive",
        border_style="cyan",
        box=box.ROUNDED,
    ))

    console.print("[dim]Creating archive…[/dim]")

    if scope == "full":
        archive_fn = archive_project_zip if use_zip else archive_project
        archive_path = archive_fn(
            project_name=project_name,
            destination_dir=dest_dir,
            include_predictions=include_predictions,
        )
    else:
        archive_fn = archive_step_zip if use_zip else archive_step
        archive_path = archive_fn(
            project_name=project_name,
            step_name=step_name,
            destination_dir=dest_dir,
        )

    elapsed      = time.time() - run_start
    size_bytes   = archive_path.stat().st_size

    console.print(Panel(
        Text.from_markup(
            f"[bold green]✓ Archive written[/bold green]  "
            f"[dim]({elapsed:.1f}s)[/dim]\n\n"
            f"[bold]File:[/bold] {archive_path.name}\n"
            f"[bold]Folder:[/bold] {archive_path.parent}\n"
            f"[bold]Size:[/bold] {size_bytes / 1024 / 1024:.1f} MB "
            f"[dim]({size_bytes:,} bytes)[/dim]\n"
            f"[bold]Includes predictions/:[/bold] "
            f"{'yes' if include_predictions else 'no'}"
        ),
        title=f"Download ready: {project_name}",
        border_style="green",
        box=box.HEAVY,
    ))

    return archive_path
