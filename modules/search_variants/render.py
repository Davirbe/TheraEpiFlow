"""Rich console rendering for search_variants: the variant candidates table."""

from rich import box
from rich.table import Table

from utils.console import console
from .core import _genotype_label, _group_candidates_by_genotype, _pick_best_per_genotype

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


def _display_genotype_grouped_table(candidates: list[dict]) -> list[dict]:
    """Renders one row per genotype with its best representative (★) and a count.

    Returns the ordered list of best-per-genotype representatives so the caller can
    map selection indices back to candidates. Used in interspecific mode to surface
    every genotype instead of burying low-identity types under the reference's isolates."""
    best_per_genotype = _pick_best_per_genotype(candidates)
    counts_by_genotype = {
        tax_id: len(bucket)
        for tax_id, bucket in _group_candidates_by_genotype(candidates).items()
    }

    table = Table(
        box=box.ROUNDED, show_header=True, header_style="bold white",
        title=f"Variants by genotype — {len(best_per_genotype)} genotypes "
              f"({len(candidates)} candidates)",
        title_style="bold cyan",
    )
    table.add_column("#",          no_wrap=True, justify="right", min_width=3)
    table.add_column("Genotype",   no_wrap=False,                 min_width=30, max_width=45)
    table.add_column("Best",       no_wrap=True, style="cyan",    min_width=12)
    table.add_column("% Identity", no_wrap=True, justify="right", min_width=10)
    table.add_column("Status",     no_wrap=True,                  min_width=12)
    table.add_column("n",          no_wrap=True, justify="right", min_width=4)

    for i, rep in enumerate(best_per_genotype, start=1):
        identity = rep["identity"]
        if identity is None or identity < 50:
            row_style = "dim red"
        elif identity < 80:
            row_style = "dim"
        else:
            row_style = ""

        identity_str = f"{identity:.1f}%" if identity is not None else "—"
        status_str   = "[bold yellow]★ Swiss-Prot[/bold yellow]" if rep["reviewed"] else "[dim]TrEMBL[/dim]"

        table.add_row(
            str(i),
            _genotype_label(rep)[:45],
            rep["accession"],
            identity_str,
            status_str,
            str(counts_by_genotype.get(rep.get("tax_id", 0), 1)),
            style=row_style,
        )

    console.print(table)
    return best_per_genotype


