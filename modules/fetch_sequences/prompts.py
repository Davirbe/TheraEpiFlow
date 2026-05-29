"""Interactive selection prompt for fetch_sequences."""

from pathlib import Path

from utils.console import console, is_interactive_session
from utils.input_validation import validate_local_path

# ── Local-FASTA-first prompt ────────────────────────────────────────────────────

def _ask_local_fasta_first(track_id: str) -> str | None:
    """Asks up front whether the user has a local FASTA for this track.

    Returns a validated path string to use as the reference, or None to fall back to
    the UniProt search. Non-interactive sessions skip this (None)."""
    if not is_interactive_session():
        return None

    console.print(
        f"[bold]Do you have a local FASTA for {track_id}?[/bold] "
        "[dim](else it will be fetched from UniProt)[/dim]"
    )
    try:
        answer = input("  Use a local FASTA? [y/N]: ").strip().lower()
    except EOFError:
        answer = ""
    if answer not in {"y", "yes"}:
        return None

    try:
        path_input = input("  Path to local FASTA: ").strip()
    except EOFError:
        path_input = ""
    path_result = validate_local_path(path_input) if path_input else None
    if path_result is None or not path_result.ok:
        if path_result is not None:
            console.print(f"  [yellow]{path_result.error} Falling back to UniProt search.[/yellow]")
        return None
    candidate = Path(path_result.value).expanduser()
    if not candidate.exists() or candidate.stat().st_size == 0:
        console.print("  [yellow]File not found or empty; falling back to UniProt search.[/yellow]")
        return None

    console.print(f"  [green]✓ Using local FASTA: {candidate}[/green]")
    return str(candidate)


# ── Selection prompt ───────────────────────────────────────────────────────────

def _prompt_selection(hits: list[dict], non_interactive: bool = False) -> dict:
    """Prompts user to select a hit by row number.
    Non-interactive: Swiss-Prot (any) > non-flagged TrEMBL > first hit.
    Interactive Enter: first non-flagged hit, or first hit if all flagged."""
    if non_interactive:
        # Swiss-Prot is human-curated — always prefer it, even if flagged
        best = (
            next((h for h in hits if h['reviewed']), None)
            or next((h for h in hits if not h.get('flag')), hits[0])
        )
        console.print(
            f'[dim]→ Auto-selected: {best["accession"]} '
            f'{best["protein_name"]}[/dim]'
        )
        return best

    best_default = next((h for h in hits if not h.get('flag')), hits[0])

    console.print(
        f'\n[bold]Select sequence[/bold] '
        f'[dim](1–{len(hits)}, or Enter to accept best candidate)[/dim]'
    )
    try:
        raw = input('> ').strip()
    except EOFError:
        raw = ''

    if not raw:
        console.print(
            f'[dim]→ Selected: {best_default["accession"]} '
            f'{best_default["protein_name"]}[/dim]'
        )
        return best_default

    try:
        idx = int(raw) - 1
        if 0 <= idx < len(hits):
            return hits[idx]
    except ValueError:
        pass

    console.print('[yellow]Invalid selection, using best candidate.[/yellow]')
    return best_default


