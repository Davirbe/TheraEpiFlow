"""Standardized end-of-step summary rendering.

Every implemented step finishes with a short narrative panel that explains
what just happened in plain language, followed by a list of files written
to disk. The previous ad-hoc `Table(..., show_header=False)` key/value
layout was hard to scan; this module replaces it with a consistent voice.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from rich import box
from rich.console import Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from utils.console import console, format_file_size_human


def _build_files_table(file_paths: Iterable[Path]) -> Table | None:
    """Return a Rich table with file name + size, or None if no files exist on disk."""
    files_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim", pad_edge=False)
    files_table.add_column("File", style="cyan", overflow="fold")
    files_table.add_column("Size", justify="right", style="dim")

    has_any_existing_file: bool = False
    for file_path in file_paths:
        if not file_path.exists():
            continue
        has_any_existing_file = True
        files_table.add_row(file_path.name, format_file_size_human(file_path.stat().st_size))

    return files_table if has_any_existing_file else None


def print_step_summary(
    step_title: str,
    elapsed_seconds: float,
    narrative_lines: list[str],
    output_files: list[Path] | None = None,
) -> None:
    """Render a uniform end-of-step summary panel.

    Parameters
    ----------
    step_title:
        Short headline of what completed (e.g. "Binding prediction complete").
    elapsed_seconds:
        Wall time of the step. Formatted as `1m23s` / `45s`.
    narrative_lines:
        Bullet-style sentences describing what happened. One sentence per
        item. Plain text or Rich markup are both accepted.
    output_files:
        Optional list of files written to disk. Only existing files are
        listed (missing entries are silently skipped, since some steps emit
        optional artifacts).
    """
    elapsed_label: str
    if elapsed_seconds >= 60:
        elapsed_minutes = int(elapsed_seconds // 60)
        elapsed_remainder = int(elapsed_seconds % 60)
        elapsed_label = f"{elapsed_minutes}m{elapsed_remainder:02d}s"
    else:
        elapsed_label = f"{elapsed_seconds:.1f}s"

    headline_text = Text.from_markup(f"[bold green]✓[/bold green] {step_title}  [dim]({elapsed_label})[/dim]")

    narrative_block_lines: list[Text] = []
    for narrative_line in narrative_lines:
        narrative_block_lines.append(Text.from_markup(f"  • {narrative_line}"))

    panel_renderables: list = [headline_text]
    if narrative_block_lines:
        panel_renderables.append(Text(""))
        panel_renderables.extend(narrative_block_lines)

    if output_files:
        files_table = _build_files_table(output_files)
        if files_table is not None:
            panel_renderables.append(Text(""))
            panel_renderables.append(Text.from_markup("[bold]Files written:[/bold]"))
            panel_renderables.append(files_table)

    console.print(Panel(Group(*panel_renderables), box=box.ROUNDED, border_style="green", padding=(1, 2)))
