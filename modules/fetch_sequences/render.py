"""Rich console rendering for fetch_sequences: the UniProt candidates table."""

from typing import Optional

from rich import box
from rich.table import Table

from utils.console import console

# ── Display table ──────────────────────────────────────────────────────────────

def _display_uniprot_table(hits: list[dict], organism: str, protein_name: Optional[str]):
    """Renders Rich selection table with Swiss-Prot / TrEMBL status and size flags."""
    reviewed_count   = sum(1 for h in hits if h['reviewed'])
    unreviewed_count = len(hits) - reviewed_count

    table = Table(
        box=box.ROUNDED, show_header=True, header_style='bold white',
        title=(
            f'UniProt — {organism} / {protein_name or "all proteins"}  '
            f'([bold yellow]{reviewed_count}[/bold yellow] Swiss-Prot  '
            f'[dim]{unreviewed_count} TrEMBL[/dim])'
        ),
        title_style='bold cyan',
    )
    table.add_column('#',         no_wrap=True, justify='right', min_width=3)
    table.add_column('Status',    no_wrap=True, min_width=30)
    table.add_column('Accession', no_wrap=True, style='cyan', min_width=12)
    table.add_column('Protein',   no_wrap=False, min_width=30, max_width=45)
    table.add_column('Length',    no_wrap=True, justify='right', min_width=12)

    for i, hit in enumerate(hits, start=1):
        if hit['reviewed']:
            status_str = '[bold yellow]⭐[/bold yellow] Reviewed (Swiss-Prot)'
            row_style  = 'yellow'
        else:
            status_str = '[dim]⚠️  Unreviewed (TrEMBL)[/dim]'
            row_style  = ''

        size_str = f'{hit["length"]} aa'
        if hit.get('flag'):
            size_str += f'  [dim red]{hit["flag"]}[/dim red]'
            row_style = 'dim'

        table.add_row(
            str(i),
            status_str,
            hit['accession'],
            hit['protein_name'][:45],
            size_str,
            style=row_style,
        )

    console.print(table)


