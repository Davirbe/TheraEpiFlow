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
