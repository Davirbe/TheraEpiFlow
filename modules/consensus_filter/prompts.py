"""Interactive presentation-percentile threshold prompt for consensus_filter."""

from rich import box
from rich.panel import Panel

from utils.console import console, is_interactive_session
from utils.project_manager import save_project_config

# ── Threshold ─────────────────────────────────────────────────────────────────

def _ask_threshold(project_name: str, project_config: dict, is_rerun: bool = False) -> float:
    """Returns the consensus threshold from project_config or prompts once and caches it.
    On a rerun (is_rerun) the saved threshold is offered for editing instead of reused silently."""
    existing_threshold = project_config.get('consensus_threshold')
    if existing_threshold is not None:
        if not (is_rerun and is_interactive_session()):
            return float(existing_threshold)
        console.print(
            f"[dim]  Consensus threshold (saved): ≤ {existing_threshold}[/dim]\n"
            "  [cyan][1][/cyan] Keep saved   [cyan][2][/cyan] Change threshold"
        )
        try:
            edit_choice = input("> ").strip()
        except EOFError:
            edit_choice = "1"
        if edit_choice != "2":
            return float(existing_threshold)

    console.print(Panel.fit(
        '[bold cyan]Consensus threshold[/bold cyan]\n'
        '[dim]Keeps peptides with percentile rank ≤ threshold in BOTH tools[/dim]\n\n'
        '  [cyan][1][/cyan] Strong + weak binders  (≤ 2.0)  [dim]default[/dim]\n'
        '  [cyan][2][/cyan] Strong binders only    (≤ 0.5)\n'
        '  [cyan][3][/cyan] Custom value',
        box=box.ROUNDED, border_style='cyan',
    ))

    while True:
        try:
            user_choice = input('> ').strip()
        except EOFError:
            user_choice = '1'
        if user_choice in ('', '1'):
            chosen_threshold = 2.0
            break
        if user_choice == '2':
            chosen_threshold = 0.5
            break
        if user_choice == '3':
            try:
                raw_threshold = input('Enter threshold (e.g. 1.5): ').strip().replace(',', '.')
            except EOFError:
                raw_threshold = '2.0'
            try:
                chosen_threshold = float(raw_threshold)
                if chosen_threshold <= 0:
                    raise ValueError
                break
            except ValueError:
                console.print('[red]Invalid value — enter a number greater than 0.[/red]')
                continue
        console.print('[red]Choose 1, 2 or 3.[/red]')

    project_config['consensus_threshold'] = chosen_threshold
    save_project_config(project_name, project_config)
    console.print(f'[dim]Threshold set to ≤ {chosen_threshold} — saved to project_config.[/dim]')
    return chosen_threshold


