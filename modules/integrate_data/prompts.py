"""integrate_data — VIEW customization UI.

Rich-based numbered-list toggle prompt. No external deps beyond Rich.

Output: a list of column internal names (in catalog order) representing the
user's chosen VIEW projection. Persisted to project_config under
`step_overrides.integrate_data.view_columns`.
"""

from __future__ import annotations

from rich import box
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from utils.console import ask, confirm, console, is_interactive_session

from .core import default_column_selection


_GROUP_HEADERS = {
    'identity':     'Identity (required)',
    'binding':      'Binding & breadth',
    'conservation': 'Conservation',
    'coverage':     'Population coverage',
    'murine':       'Murine model',
}
_GROUP_ORDER = ['identity', 'binding', 'conservation', 'coverage', 'murine']


def _render_catalog(catalog: list[dict], selected: set[str]) -> None:
    """Pretty-prints the catalog as a numbered table grouped by category."""
    table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim", pad_edge=False)
    table.add_column("#",         style="dim",  justify="right")
    table.add_column("",          style="bold")   # checkbox
    table.add_column("Column",    style="cyan")
    table.add_column("Header",    style="white")
    table.add_column("Group",     style="dim")

    for one_based_index, entry in enumerate(catalog, start=1):
        is_selected = entry['name'] in selected
        is_required = entry['required']
        if is_required:
            check_glyph = "[bold green]■[/bold green]"
        elif is_selected:
            check_glyph = "[green]■[/green]"
        else:
            check_glyph = "[dim]□[/dim]"
        group_label = _GROUP_HEADERS.get(entry['group'], entry['group'])
        table.add_row(
            str(one_based_index),
            check_glyph,
            entry['name'],
            entry['header'],
            group_label,
        )
    console.print(table)


def prompt_view_customization(catalog: list[dict]) -> list[str]:
    """Interactive flow: ask Y/n for default; on n, render checkboxes + take toggles.

    In non-interactive sessions (no TTY), returns the default selection silently —
    this lets `--step integrate_data` pipelines work without prompting.
    """
    if not is_interactive_session():
        return default_column_selection(catalog)

    use_default = confirm("Use default VIEW columns?", default=True)
    if use_default:
        return default_column_selection(catalog)

    selected: set[str] = set(default_column_selection(catalog))

    console.print()
    console.print(Panel.fit(
        Text.from_markup(
            "[bold]Customize VIEW columns[/bold]\n"
            "[dim]Required identity columns are always on. Toggle others by typing "
            "their numbers (space-separated). Press Enter on an empty line to confirm.[/dim]"
        ),
        box=box.ROUNDED,
        border_style="cyan",
    ))

    while True:
        _render_catalog(catalog, selected)
        raw_input_text = ask(
            "Toggle (e.g. '4 9 12'), or Enter to confirm",
            default="",
        )
        if not raw_input_text.strip():
            break

        try:
            requested_indices = [int(token) for token in raw_input_text.split()]
        except ValueError:
            console.print("[yellow]Could not parse: type space-separated numbers like '4 9 12'.[/yellow]")
            continue

        for one_based_index in requested_indices:
            if not (1 <= one_based_index <= len(catalog)):
                console.print(f"[yellow]Skipped {one_based_index}: out of range.[/yellow]")
                continue
            entry = catalog[one_based_index - 1]
            if entry['required']:
                console.print(f"[yellow]Skipped {one_based_index}: '{entry['name']}' is required.[/yellow]")
                continue
            if entry['name'] in selected:
                selected.discard(entry['name'])
            else:
                selected.add(entry['name'])

    # Return in catalog order so downstream rendering is stable.
    return [entry['name'] for entry in catalog if entry['name'] in selected]
