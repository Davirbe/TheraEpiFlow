"""Centralized Rich Console and TTY detection.

Every module should import the shared `console` from here instead of
constructing its own `Console()`. This guarantees consistent width and
makes it easy to inject a different console (e.g., for tests) later.
"""

from __future__ import annotations

import sys
from typing import Optional

from rich.console import Console
from rich.prompt import Confirm, Prompt

DEFAULT_TERMINAL_WIDTH = 120

console: Console = Console(width=DEFAULT_TERMINAL_WIDTH)

# Subtle zebra striping for data tables (alternating rows readable but not noisy).
# Apply by passing `row_styles=DEFAULT_TABLE_ROW_STYLES` to `rich.table.Table(...)`.
# Skip this for key/value summary tables (`show_header=False` with 2–6 rows) —
# zebra striping only helps when scanning long lists of similar rows.
DEFAULT_TABLE_ROW_STYLES: list[str] = ["", "on grey7"]


def format_file_size_human(num_bytes: int) -> str:
    """Return a short human-readable representation of a file size."""
    size_in_bytes: float = float(num_bytes)
    for unit_label in ("B", "KB", "MB", "GB"):
        if size_in_bytes < 1024.0 or unit_label == "GB":
            return f"{size_in_bytes:,.0f} {unit_label}" if unit_label == "B" else f"{size_in_bytes:.1f} {unit_label}"
        size_in_bytes /= 1024.0
    return f"{size_in_bytes:.1f} GB"


def is_interactive_session() -> bool:
    """Return True when stdin is attached to a terminal.

    Used to gate interactive prompts and post-step menus so they never fire
    in CI, redirected stdin, or one-shot `--step` invocations that the user
    explicitly did not want to interact with.
    """
    return sys.stdin.isatty()


def ask(
    text: str,
    default: Optional[str] = None,
    choices: Optional[list[str]] = None,
    required: bool = False,
) -> str:
    """Rich-styled prompt with non-interactive fallback.

    Falls back to `default` (or "" when default is None) when stdin is not a
    TTY, so any caller can use the same prompt API in CI / piped runs without
    crashing on EOFError.

    When `required=True` and there is no `default`, an empty answer is rejected
    and the prompt is re-shown until a non-empty value is provided.
    """
    if not is_interactive_session():
        return default if default is not None else ""

    while True:
        response = Prompt.ask(
            text,
            default=default,
            choices=choices,
            show_default=default is not None,
            console=console,
        )
        if required and default is None and not str(response).strip():
            console.print('[red]Empty input not allowed. Please type a value.[/red]')
            continue
        return response


def confirm(text: str, default: bool = False) -> bool:
    """Rich-styled yes/no prompt with non-interactive fallback."""
    if not is_interactive_session():
        return default
    return Confirm.ask(text, default=default, console=console)


def press_enter_to_continue(text: str = 'Press Enter to continue') -> None:
    """Pauses execution until the user presses Enter. No-op in non-interactive mode."""
    if not is_interactive_session():
        return
    try:
        input(f'\n  {text}…')
    except EOFError:
        return


def flush_stdin() -> None:
    """Discards any pending keystrokes from the terminal input buffer.

    Call after long-running operations where the user may have impatiently
    pressed Enter thinking the program froze. Without this, those buffered
    keystrokes would be consumed by the next `input()` and silently skip
    prompts. POSIX-only (Linux/macOS/WSL); no-op on platforms without termios.
    """
    if not is_interactive_session():
        return
    try:
        import termios
    except ImportError:
        return
    try:
        termios.tcflush(sys.stdin.fileno(), termios.TCIFLUSH)
    except (termios.error, OSError):
        return


def confirm_value(label: str, value: str, indent: str = '') -> bool:
    """Echoes a value the user just typed and asks for explicit y/n confirmation.

    Used on critical wizard fields (organism, protein, FASTA path, HLA list,
    threshold, cluster method) so that a typo can be corrected with a single
    keystroke instead of having to re-run the whole phase. Returns True when
    the user accepts the value, False to signal a re-prompt is needed.

    Empty / whitespace-only `value` always returns False (a non-empty value
    is a prerequisite of asking for confirmation).
    """
    if not is_interactive_session():
        return True
    if not str(value).strip():
        return False
    console.print(
        f"{indent}[dim]→[/dim] [bold]{label}:[/bold] [cyan]{value}[/cyan]"
    )
    while True:
        try:
            answer = input(f"{indent}  [y] confirm  [n] re-enter: ").strip().lower()
        except EOFError:
            return True
        if answer in ('', 'y', 'yes'):
            return True
        if answer in ('n', 'no'):
            return False
        console.print(f"{indent}  [red]Type y or n.[/red]")


def show_recap_and_confirm(
    title: str,
    fields: list[tuple[str, str]],
    proceed_label: str = 'Proceed with these values',
) -> bool:
    """Renders a recap Panel of all collected values and asks one final y/n.

    `fields` is an ordered list of (label, value) pairs. Returns True when the
    user accepts the recap, False to signal the caller should restart the
    surrounding phase. Non-interactive sessions always return True.

    Designed for the moment right before a wizard saves anything to disk —
    the caller is responsible for treating False as "do not persist".
    """
    from rich import box
    from rich.panel import Panel
    from rich.table import Table

    if not is_interactive_session():
        return True

    recap_table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1), expand=False)
    recap_table.add_column(style='bold', no_wrap=True, justify='right')
    recap_table.add_column(style='cyan')
    for label, value in fields:
        recap_table.add_row(f"{label}:", str(value) if value else '[dim]—[/dim]')

    console.print()
    console.print(Panel(
        recap_table,
        title=f'[bold cyan]{title}[/bold cyan]',
        title_align='left',
        box=box.HEAVY_EDGE,
        border_style='cyan',
        padding=(1, 2),
    ))
    while True:
        try:
            answer = input(f"  {proceed_label}? [y] yes  [n] restart: ").strip().lower()
        except EOFError:
            return True
        if answer in ('', 'y', 'yes'):
            return True
        if answer in ('n', 'no'):
            return False
        console.print("  [red]Type y or n.[/red]")
