"""Small text-formatting helpers for the CLI render layer.

These exist so every step that prints tables / panels stays visually consistent
and never lets a single value blow up the layout. Used heavily by render-side
modules (`*/render.py`); pure functions, no Rich or pandas dependency.
"""

from __future__ import annotations


def compact_num(value: object) -> str:
    """Format a number with thousands separator (English locale).

    Non-numeric input passes through as ``str(value)`` so it's safe to use as
    a one-liner in table builders that mix numeric and text cells.

    Examples:
        compact_num(24300)   →  "24,300"
        compact_num(0.85)    →  "0.85"
        compact_num("n/a")   →  "n/a"
    """
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,d}"
    if isinstance(value, float):
        if value.is_integer():
            return f"{int(value):,d}"
        return f"{value:,.4g}"
    return str(value)


def truncate_with_ellipsis(text: object, max_len: int = 60) -> str:
    """Cap a string at `max_len` characters, replacing the tail with `…`.

    Returns the input unchanged if it already fits. Always renders something
    even when the input is None / non-string.
    """
    if text is None:
        return ""
    string_value = str(text)
    if len(string_value) <= max_len:
        return string_value
    if max_len <= 1:
        return "…"
    return string_value[: max_len - 1] + "…"
