"""Rich console rendering for search_variants: the variant candidates table."""

from rich import box
from rich.table import Table

from utils.console import console

# ── Rich table display ────────────────────────────────────────────────────────

def _display_variants_table(candidates: list[dict]):
    """Renders a Rich table of variant candidates colour-coded by identity."""
    table = Table(
        box=box.ROUNDED, show_header=True, header_style="bold white",
        title=f"UniProt Variants — {len(candidates)} candidates",
        title_style="bold cyan",
    )
    table.add_column("#",           no_wrap=True, justify="right", min_width=3)
    table.add_column("Accession",   no_wrap=True, style="cyan",    min_width=12)
    table.add_column("Organism",    no_wrap=False,                  min_width=30, max_width=45)
    table.add_column("% Identity",  no_wrap=True, justify="right", min_width=10)
    table.add_column("Length",      no_wrap=True, justify="right", min_width=8)
    table.add_column("Status",      no_wrap=True,                  min_width=12)

    for i, c in enumerate(candidates, start=1):
        identity = c["identity"]

        if identity is None or identity < 50:
            row_style = "dim red"
        elif identity < 80:
            row_style = "dim"
        else:
            row_style = ""

        identity_str = f"{identity:.1f}%" if identity is not None else "—"
        status_str   = "[bold yellow]★ Swiss-Prot[/bold yellow]" if c["reviewed"] else "[dim]TrEMBL[/dim]"

        table.add_row(
            str(i),
            c["accession"],
            c["organism"][:45],
            identity_str,
            f"{c['length']} aa",
            status_str,
            style=row_style,
        )

    console.print(table)


