"""Rich console rendering for search_variants: the variant candidates table."""

from rich import box
from rich.table import Table

from utils.console import console
from .core import _genotype_label, _group_candidates_by_genotype, _pick_best_per_genotype

# Human-readable, colour-coded labels for the advisory identity flags.
_FLAG_LABELS = {
    "possibly_unrelated": "[red]<30% unrelated?[/red]",
    "near_identical":     "[yellow]≥99% near-ident.[/yellow]",
    "uncut_polyprotein":  "[magenta]uncut polyprotein?[/magenta]",
}


def _format_flags(candidate: dict) -> str:
    """Renders a candidate's advisory flags for the table (empty string = none)."""
    return " ".join(_FLAG_LABELS.get(flag, flag) for flag in candidate.get("flags", []))


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
    table.add_column("Flags",       no_wrap=True,                  min_width=14)

    for row_number, candidate in enumerate(candidates, start=1):
        identity_percent = candidate["identity"]

        if identity_percent is None or identity_percent < 50:
            row_style = "dim red"
        elif identity_percent < 80:
            row_style = "dim"
        else:
            row_style = ""

        identity_text = f"{identity_percent:.1f}%" if identity_percent is not None else "—"
        status_text   = "[bold yellow]★ Swiss-Prot[/bold yellow]" if candidate["reviewed"] else "[dim]TrEMBL[/dim]"

        table.add_row(
            str(row_number),
            candidate["accession"],
            candidate["organism"][:45],
            identity_text,
            f"{candidate['length']} aa",
            status_text,
            _format_flags(candidate),
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
    table.add_column("Flags",      no_wrap=True,                  min_width=14)
    table.add_column("n",          no_wrap=True, justify="right", min_width=4)

    for row_number, genotype_representative in enumerate(best_per_genotype, start=1):
        identity_percent = genotype_representative["identity"]
        if identity_percent is None or identity_percent < 50:
            row_style = "dim red"
        elif identity_percent < 80:
            row_style = "dim"
        else:
            row_style = ""

        identity_text = f"{identity_percent:.1f}%" if identity_percent is not None else "—"
        status_text   = (
            "[bold yellow]★ Swiss-Prot[/bold yellow]"
            if genotype_representative["reviewed"] else "[dim]TrEMBL[/dim]"
        )

        table.add_row(
            str(row_number),
            _genotype_label(genotype_representative)[:45],
            genotype_representative["accession"],
            identity_text,
            status_text,
            _format_flags(genotype_representative),
            str(counts_by_genotype.get(genotype_representative.get("tax_id", 0), 1)),
            style=row_style,
        )

    console.print(table)
    return best_per_genotype


