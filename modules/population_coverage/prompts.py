"""
Interactive preflight prompts for population_coverage: population selection and
the optional informational coverage cutoff. The only layer that calls input();
selections are persisted into project_config.
"""

from typing import Optional

from rich import box
from rich.panel import Panel

from utils.console import console, is_interactive_session
from utils.project_manager import save_project_config

from .core import _list_available_populations

# Curated short list of populations for the interactive picker. The user can
# always type a free-form name that is validated against the pickle.
_SUGGESTED_POPULATIONS = [
    "World",
    "Europe",
    "East Asia",
    "South Asia",
    "Africa",
    "North America (Native American)",
    "South America",
    "Brazil",
]


def _prompt_populations(project_name: str, project_config: dict, population_db: dict) -> list[str]:
    """Returns the list of populations to compute. Persists into project_config."""
    saved_populations = project_config.get("coverage_populations")

    if not is_interactive_session():
        return list(saved_populations) if saved_populations else ["World"]

    available_populations = _list_available_populations(population_db)

    if saved_populations:
        console.print(Panel(
            f"[bold]Population coverage — saved selection[/bold]\n\n"
            f"Current populations: [cyan]{', '.join(saved_populations)}[/cyan]\n\n"
            "  [cyan][1][/cyan] Keep current selection\n"
            "  [cyan][2][/cyan] Change selection",
            box=box.ROUNDED, title="Setup: population_coverage", title_align="left",
        ))
        while True:
            try:
                choice = input("> ").strip()
            except EOFError:
                choice = "1"
            if choice in ("1", ""):
                return list(saved_populations)
            if choice == "2":
                break
            console.print("[dim]Type 1 or 2.[/dim]")

    console.print(Panel(
        "[bold]Population coverage — pick one or more populations[/bold]\n\n"
        "[dim]The IEDB allele-frequency database covers 239 populations.\n"
        "Type a comma-separated list of numbers (e.g. '1,3,8') or 'other'\n"
        "to type a free-form name validated against the full list.[/dim]\n\n"
        + "\n".join(
            f"  [cyan][{i + 1}][/cyan] {pop_name}"
            for i, pop_name in enumerate(_SUGGESTED_POPULATIONS)
        )
        + "\n  [cyan][other][/cyan] enter a custom population name",
        box=box.ROUNDED, title="Setup: population_coverage", title_align="left",
    ))

    selected_populations: list[str] = []
    while not selected_populations:
        try:
            raw_input_text = input("> ").strip().lower()
        except EOFError:
            raw_input_text = "1"

        if raw_input_text == "other":
            try:
                custom_name = input("Population name (exact): ").strip()
            except EOFError:
                custom_name = ""
            if custom_name in available_populations:
                selected_populations = [custom_name]
            else:
                console.print(
                    f"[red]'{custom_name}' is not in the database.[/red]\n"
                    f"[dim]Hint: search is case-sensitive. Examples:\n"
                    f"{', '.join(available_populations[:8])}, …[/dim]"
                )
            continue

        try:
            indices = [int(token.strip()) for token in raw_input_text.split(",") if token.strip()]
        except ValueError:
            console.print("[red]Invalid input — type numbers separated by commas, or 'other'.[/red]")
            continue

        candidate_populations: list[str] = []
        for idx in indices:
            if 1 <= idx <= len(_SUGGESTED_POPULATIONS):
                candidate = _SUGGESTED_POPULATIONS[idx - 1]
                if candidate in available_populations:
                    candidate_populations.append(candidate)
                else:
                    console.print(
                        f"[yellow]'{candidate}' not in this version of the IEDB pickle, skipping.[/yellow]"
                    )
            else:
                console.print(f"[yellow]Index {idx} out of range, skipping.[/yellow]")

        if candidate_populations:
            seen: set[str] = set()
            selected_populations = [
                pop for pop in candidate_populations if not (pop in seen or seen.add(pop))
            ]

    project_config["coverage_populations"] = selected_populations
    save_project_config(project_name, project_config)
    console.print(
        f"[dim]Populations saved to project_config: "
        f"{', '.join(selected_populations)}[/dim]"
    )
    return selected_populations


def _prompt_coverage_cutoff(project_name: str, project_config: dict) -> Optional[float]:
    """Returns an optional informational cutoff percentage. Persisted."""
    saved_cutoff = project_config.get("coverage_minimum_pct")

    if not is_interactive_session():
        return float(saved_cutoff) if saved_cutoff is not None else None

    if saved_cutoff is not None:
        console.print(
            f"[dim]Coverage cutoff saved earlier: "
            f"[cyan]{saved_cutoff}%[/cyan] (informational only).[/dim]"
        )
        return float(saved_cutoff)

    console.print(Panel(
        "[bold]Coverage cutoff (optional)[/bold]\n\n"
        "[dim]A reference threshold used only to flag low-coverage epitopes\n"
        "in the final HTML report. It does NOT remove anything here.[/dim]\n\n"
        "  [cyan][1][/cyan] Skip (no cutoff)\n"
        "  [cyan][2][/cyan] Set a value",
        box=box.ROUNDED, title="Setup: population_coverage", title_align="left",
    ))
    while True:
        try:
            choice = input("> ").strip()
        except EOFError:
            choice = "1"
        if choice in ("1", ""):
            project_config["coverage_minimum_pct"] = None
            save_project_config(project_name, project_config)
            return None
        if choice == "2":
            try:
                raw_value = input("Cutoff % (e.g. 10): ").strip().replace(",", ".")
            except EOFError:
                raw_value = ""
            try:
                cutoff_value = float(raw_value)
                if 0.0 <= cutoff_value <= 100.0:
                    project_config["coverage_minimum_pct"] = cutoff_value
                    save_project_config(project_name, project_config)
                    return cutoff_value
            except ValueError:
                pass
            console.print("[red]Invalid value. Use a number between 0 and 100.[/red]")
        else:
            console.print("[dim]Type 1 or 2.[/dim]")
