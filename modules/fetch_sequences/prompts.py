"""Interactive selection prompt for fetch_sequences."""

from utils.console import console

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
            f'— {best["protein_name"]}[/dim]'
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
            f'— {best_default["protein_name"]}[/dim]'
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


