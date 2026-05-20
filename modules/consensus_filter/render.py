"""Rich progressive display for consensus_filter: the four per-stage tables."""

import pandas as pd
from rich import box
from rich.table import Table

from utils.console import console

# ── Progressive display ───────────────────────────────────────────────────────

def _print_stage1_filtering(net_data: dict, flu_data: dict, threshold: float):
    """Stage 1 — loading, NaN removal, and threshold with a cutoff comparison."""
    display_table = Table(
        box=box.SIMPLE, show_header=True, header_style='bold',
        title='[bold cyan]Stage 1 — Presentation threshold filter[/bold cyan]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('NetMHCpan', justify='right')
    display_table.add_column('MHCFlurry', justify='right')

    display_table.add_row('Prediction rows (input)',  str(net_data['n_raw']),        str(flu_data['n_raw']))
    display_table.add_row('  removed (NaN)',          f'-{net_data["dropped_nan"]}', f'-{flu_data["dropped_nan"]}')
    display_table.add_row('  remaining',              str(net_data['n_0a']),         str(flu_data['n_0a']))
    display_table.add_section()

    def _mark_active(cutoff):
        return ' [bold yellow]★[/bold yellow]' if cutoff == threshold else ''

    display_table.add_row(
        f'≤ 0.5  (strong binders only){_mark_active(0.5)}',
        str(net_data['n_strong']), str(flu_data['n_strong'])
    )
    display_table.add_row(
        f'≤ 2.0  (strong + weak binders){_mark_active(2.0)}',
        str(net_data['n_weak']), str(flu_data['n_weak'])
    )
    if threshold not in (0.5, 2.0):
        display_table.add_row(
            f'≤ {threshold}  (custom) [bold yellow]★[/bold yellow]',
            str(net_data['n_0b']), str(flu_data['n_0b'])
        )

    console.print(display_table)


def _print_stage2_consolidation(
    net_data: dict,
    flu_data: dict,
    net_consolidated: pd.DataFrame,
    flurry_consolidated: pd.DataFrame,
):
    """Stage 2 — consolidation of allele rows into unique peptides."""
    n_net_peptides    = len(net_consolidated)
    n_flurry_peptides = len(flurry_consolidated)

    display_table = Table(
        box=box.SIMPLE, show_header=True, header_style='bold',
        title='[bold green]Stage 2 — Consolidation (allele×peptide rows → unique peptides)[/bold green]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('NetMHCpan', justify='right')
    display_table.add_column('MHCFlurry', justify='right')

    display_table.add_row('Rows after threshold',     str(net_data['n_0b']),                     str(flu_data['n_0b']))
    display_table.add_row('Unique peptides',          str(n_net_peptides),                       str(n_flurry_peptides))
    display_table.add_row('Collapsed rows',           str(net_data['n_0b'] - n_net_peptides),    str(flu_data['n_0b'] - n_flurry_peptides))

    console.print(display_table)


def _print_stage3_intersection(intersection_data: dict):
    """Stage 3 — intersection between NetMHCpan and MHCFlurry."""
    display_table = Table(
        box=box.SIMPLE, show_header=False,
        title='[bold yellow]Stage 3 — Tool intersection[/bold yellow]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('', justify='right')

    display_table.add_row('NetMHCpan only',  str(intersection_data['net_only']))
    display_table.add_row('MHCFlurry only',  str(intersection_data['flurry_only']))
    display_table.add_row(
        'In consensus (both tools) [bold yellow]★[/bold yellow]',
        f'[bold green]{intersection_data["common_count"]}[/bold green]'
    )

    console.print(display_table)


def _print_stage4_immunogenicity(n_input: int, n_survivors: int):
    """Stage 4 — result of the Calis 2013 immunogenicity filter."""
    n_discarded  = n_input - n_survivors
    pct_survived = (n_survivors / n_input * 100) if n_input > 0 else 0.0
    pct_discarded = 100.0 - pct_survived

    display_table = Table(
        box=box.SIMPLE, show_header=False,
        title='[bold magenta]Stage 4 — Calis 2013 immunogenicity (score > 0)[/bold magenta]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('', justify='right')

    display_table.add_row('Input (consensus)',  str(n_input))
    display_table.add_row(
        'Survivors',
        f'[bold green]{n_survivors}[/bold green]  ({pct_survived:.0f}%)'
    )
    display_table.add_row(
        'Discarded',
        f'[dim]{n_discarded}[/dim]  ({pct_discarded:.0f}%)'
    )

    console.print(display_table)
    console.print(f'\n[bold green]FINAL RESULT: {n_survivors} immunogenic epitopes[/bold green]\n')


