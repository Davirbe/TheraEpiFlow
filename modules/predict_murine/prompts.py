"""Interactive murine-strain selection prompt for predict_murine."""

from rich.prompt import Prompt

import config
from utils.console import console, is_interactive_session
from utils.project_manager import save_project_config

_STRAIN_MENU_ORDER = ['default', 'C57BL/6', 'BALB/c', 'CBA/C3H/AKR', 'all']
# ── Strain prompt ─────────────────────────────────────────────────────────────

def _ask_murine_strain(
    project_name: str, project_config: dict, is_rerun: bool = False,
) -> tuple[str, list[str]]:
    """Returns (strain_group_name, list_of_h2_alleles).
    On a rerun (is_rerun) the saved strain is offered for editing instead of reused silently."""
    saved_strain_name = project_config.get('murine_strain')
    if saved_strain_name and saved_strain_name in config.MURINE_ALLELES:
        saved_h2_alleles = list(config.MURINE_ALLELES[saved_strain_name])
        console.print(
            f"[dim]  Strain (saved): {saved_strain_name} "
            f"({len(saved_h2_alleles)} alleles)[/dim]"
        )
        if not (is_rerun and is_interactive_session()):
            return saved_strain_name, saved_h2_alleles
        console.print(
            "  [cyan][1][/cyan] Keep saved strain   [cyan][2][/cyan] Pick a different strain"
        )
        try:
            edit_choice = input("> ").strip()
        except EOFError:
            edit_choice = "1"
        if edit_choice != "2":
            return saved_strain_name, saved_h2_alleles
        console.print("[dim]  Re-selecting murine strain…[/dim]")

    console.print("\n[bold]Murine MHC-I strain group[/bold]")
    menu_rows = [
        ("default",     "C57BL/6 + BALB/c",        5, "[ENTER]"),
        ("C57BL/6",     "",                         2, ""),
        ("BALB/c",      "",                         3, ""),
        ("CBA/C3H/AKR", "",                         2, ""),
        ("all",         "every consolidated H-2",   7, ""),
    ]
    for menu_index, (group_label, group_description, num_alleles, tail) in enumerate(menu_rows, start=1):
        console.print(
            f"  [cyan][{menu_index}][/cyan] "
            f"{group_label:<14} "
            f"{group_description:<26} "
            f"({num_alleles} alleles)  {tail}"
        )

    try:
        chosen_index_text = Prompt.ask("Choice", default="1")
    except EOFError:
        chosen_index_text = "1"

    try:
        chosen_index = int(chosen_index_text.strip())
    except ValueError:
        chosen_index = 1
    if chosen_index < 1 or chosen_index > len(_STRAIN_MENU_ORDER):
        console.print("[yellow]  Invalid choice — falling back to default.[/yellow]")
        chosen_index = 1

    chosen_strain_name = _STRAIN_MENU_ORDER[chosen_index - 1]
    chosen_h2_alleles  = list(config.MURINE_ALLELES[chosen_strain_name])

    project_config['murine_strain'] = chosen_strain_name
    save_project_config(project_name, project_config)

    console.print(
        f"[dim]  → {chosen_strain_name}: {len(chosen_h2_alleles)} alleles "
        f"({', '.join(chosen_h2_alleles)})[/dim]"
    )
    return chosen_strain_name, chosen_h2_alleles


