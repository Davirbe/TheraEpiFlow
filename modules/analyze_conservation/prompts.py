"""Interactive prompts for analyze_conservation: the identity threshold and the
preflight local-FASTA override prompts."""

from pathlib import Path

from rich import box
from rich.panel import Panel

from utils.console import console, is_interactive_session
from utils.project_manager import save_project_config

# ── Threshold prompt ──────────────────────────────────────────────────────────

def prompt_analysis_threshold(project_name: str, project_config: dict, is_rerun: bool = False) -> float:
    """
    Returns the analysis threshold for this run.
    Interactive rerun: shows the saved value and lets the user change it.
    First run / non-interactive / downstream cascade re-run: uses the saved value
    (or default 1.0) without prompting, so only an explicit redo re-asks.
    Persists to project_config['conservation_threshold'] after any change.
    """
    saved = project_config.get("conservation_threshold")

    if saved is not None and not (is_rerun and is_interactive_session()):
        return float(saved)

    if not is_interactive_session():
        return float(saved) if saved is not None else 1.0

    if saved is not None:
        console.print(Panel(
            f"[bold]Conservation analysis threshold[/bold]\n\n"
            f"Current threshold: [cyan]{int(round(float(saved) * 100))}%[/cyan]\n\n"
            "[dim]Determines which variants are shown as passed/failed.\n"
            "Summary counts (100%, 90%, 80%) and row colours are fixed\n"
            "regardless of this value.[/dim]\n\n"
            "  [cyan][1][/cyan] Keep current value\n"
            "  [cyan][2][/cyan] Change threshold",
            box=box.ROUNDED, title="Setup: analyze_conservation", title_align="left",
        ))
        while True:
            try:
                choice = input("> ").strip()
            except EOFError:
                choice = "1"
            if choice in ("1", ""):
                return float(saved)
            if choice == "2":
                break
            console.print("[dim]Type 1 or 2.[/dim]")

    console.print(Panel(
        "[bold]Conservation analysis threshold[/bold]\n\n"
        "[dim]Variants with identity >= threshold appear as 'passed'.\n"
        "Failed variants show their actual best-matching window and mutations.[/dim]\n\n"
        "  [cyan][1][/cyan] 1.00 — exact match only (default)\n"
        "  [cyan][2][/cyan] 0.90 — 90% identity\n"
        "  [cyan][3][/cyan] 0.80 — 80% identity\n"
        "  [cyan][4][/cyan] Custom value (e.g. 0.75, 0.85, 0.65)",
        box=box.ROUNDED, title="Setup: analyze_conservation", title_align="left",
    ))

    preset_map = {"1": 1.0, "2": 0.90, "3": 0.80}
    while True:
        try:
            choice = input("> ").strip()
        except EOFError:
            choice = "1"
        if choice in preset_map:
            threshold = preset_map[choice]
            break
        if choice == "4":
            while True:
                try:
                    raw = input("Threshold (e.g. 0.75): ").strip().replace(",", ".")
                except EOFError:
                    raw = "1.0"
                try:
                    threshold = float(raw)
                    if 0.0 < threshold <= 1.0:
                        break
                    console.print("[red]Value must be between 0 (exclusive) and 1 (inclusive).[/red]")
                except ValueError:
                    console.print("[red]Invalid value. Use a dot as decimal separator (e.g. 0.75).[/red]")
            break
        console.print("[dim]Invalid option. Type 1, 2, 3 or 4.[/dim]")

    project_config["conservation_threshold"] = threshold
    save_project_config(project_name, project_config)
    console.print(f"[dim]Threshold set to {threshold:.2f}, saved to project_config.[/dim]")
    return threshold


def prompt_length_filter(is_rerun: bool = False) -> bool:
    """Whether to apply the ±length filter to variants. Default ON (returns True).

    On a rerun (interactive) it offers to disable the filter — useful when a track
    errored or returned too few variants because partial/fragment sequences were dropped.
    First runs and non-interactive sessions keep the default (filter ON)."""
    if not (is_rerun and is_interactive_session()):
        return True

    console.print(
        "\n[bold]Variant length filter[/bold] "
        "[dim](±length tolerance vs. reference)[/dim]"
    )
    console.print("  [cyan]1[/cyan] — Keep filter ON (default — drops partial/fragment variants)")
    console.print("  [cyan]2[/cyan] — Turn filter OFF (include every variant regardless of length)")
    try:
        choice = input("> ").strip()
    except EOFError:
        choice = "1"
    if choice == "2":
        console.print("[dim]→ Length filter OFF — keeping all variants.[/dim]")
        return False
    console.print("[dim]→ Length filter ON.[/dim]")
    return True


def _ask_for_local_fasta_overrides(track_status_rows: list[dict]) -> dict:
    """
    Loop interactively letting the user attach a local FASTA path to one or
    more track ids. Returns {track_id: path_string}. Empty dict if the user
    declines or supplies nothing usable.
    """
    fasta_overrides: dict[str, str] = {}
    track_id_set = {row["track_id"] for row in track_status_rows}

    while True:
        try:
            track_input = input("  Track id to attach FASTA to (or blank to finish): ").strip()
        except EOFError:
            track_input = ""

        if not track_input:
            break
        if track_input not in track_id_set:
            console.print(f"  [yellow]Unknown track id '{track_input}'. Pick one from the table above.[/yellow]")
            continue

        try:
            path_input = input(f"  Local FASTA path for {track_input}: ").strip()
        except EOFError:
            path_input = ""

        candidate_path = Path(path_input).expanduser() if path_input else None
        if candidate_path is None or not candidate_path.exists() or candidate_path.stat().st_size == 0:
            console.print("  [red]File not found or empty. Skipping.[/red]")
            continue

        fasta_overrides[track_input] = str(candidate_path)
        console.print(f"  [green]✓ Will use {candidate_path} for {track_input}.[/green]")

    return fasta_overrides


