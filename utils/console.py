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
) -> str:
    """Rich-styled prompt with non-interactive fallback.

    Falls back to `default` (or "" when default is None) when stdin is not a
    TTY, so any caller can use the same prompt API in CI / piped runs without
    crashing on EOFError.
    """
    if not is_interactive_session():
        return default if default is not None else ""
    return Prompt.ask(
        text,
        default=default,
        choices=choices,
        show_default=default is not None,
        console=console,
    )


def confirm(text: str, default: bool = False) -> bool:
    """Rich-styled yes/no prompt with non-interactive fallback."""
    if not is_interactive_session():
        return default
    return Confirm.ask(text, default=default, console=console)
