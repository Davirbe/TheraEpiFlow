"""Centralized Rich Console and TTY detection.

Every module should import the shared `console` from here instead of
constructing its own `Console()`. This guarantees consistent width and
makes it easy to inject a different console (e.g., for tests) later.
"""

from __future__ import annotations

import sys

from rich.console import Console

DEFAULT_TERMINAL_WIDTH = 120

console: Console = Console(width=DEFAULT_TERMINAL_WIDTH)


def is_interactive_session() -> bool:
    """Return True when stdin is attached to a terminal.

    Used to gate interactive prompts and post-step menus so they never fire
    in CI, redirected stdin, or one-shot `--step` invocations that the user
    explicitly did not want to interact with.
    """
    return sys.stdin.isatty()
