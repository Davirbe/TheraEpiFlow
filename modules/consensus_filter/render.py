"""Rich progressive display for consensus_filter: the four per-stage panels.

Each stage renders as a Panel containing a short narrative + reference + a
small data table. The narrative is intentionally explicit so a researcher
who has never read the code can follow what the step actually did.
"""

import pandas as pd
from rich import box
from rich.console import Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from utils.console import DEFAULT_TABLE_ROW_STYLES, console
from utils.text_format import compact_num


_NARRATIVE_BLANK = "[dim]_[/dim]"   # placeholder so empty cells stay aligned


def _data_table(*, headers: list[str], rows: list[list[object]]) -> Table:
    """Internal helper: a zebra-striped Table with the consensus header style."""
    table = Table(
        box=box.SIMPLE,
        show_header=True,
        header_style="bold",
        row_styles=DEFAULT_TABLE_ROW_STYLES,
        pad_edge=False,
    )
    table.add_column(headers[0], style="bold cyan", no_wrap=True)
    for column_label in headers[1:]:
        table.add_column(column_label, justify="right")
    for row_cells in rows:
        table.add_row(*(str(cell) for cell in row_cells))
    return table


def _stage_panel(
    *,
    title:        str,
    what_happened: str,
    reference:    str | None,
    data_table:   Table,
    example:      str | None = None,
    border_style: str = "cyan",
) -> Panel:
    """Composes a stage panel: title + narrative + (optional example) + (optional ref) + table."""
    body_blocks: list = []
    body_blocks.append(Text.from_markup(f"[bold]What happened:[/bold] {what_happened}"))
    if example:
        body_blocks.append(Text(""))
        body_blocks.append(Text.from_markup(f"[bold]Example:[/bold] {example}"))
    if reference:
        body_blocks.append(Text(""))
        body_blocks.append(Text.from_markup(f"[bold dim]Reference:[/bold dim] [dim]{reference}[/dim]"))
    body_blocks.append(Text(""))
    body_blocks.append(data_table)

    return Panel(
        Group(*body_blocks),
        title=f"[bold]{title}[/bold]",
        title_align="left",
        border_style=border_style,
        box=box.ROUNDED,
        padding=(1, 2),
    )


# ── Stage renderers ───────────────────────────────────────────────────────────

def _print_stage1_filtering(net_data: dict, flu_data: dict, threshold: float):
    """Stage 1 — affinity threshold filter on the raw predictions."""
    base_rows: list[list[object]] = [
        ["Input rows",         compact_num(net_data["n_raw"]),         compact_num(flu_data["n_raw"])],
        ["  removed (NaN)",    f'-{compact_num(net_data["dropped_nan"])}', f'-{compact_num(flu_data["dropped_nan"])}'],
        ["  remaining",        compact_num(net_data["n_0a"]),          compact_num(flu_data["n_0a"])],
    ]

    def _mark(active_cutoff: float) -> str:
        return "  [bold yellow]★[/bold yellow]" if active_cutoff == threshold else ""

    band_rows: list[list[object]] = [
        [f'≤ 0.5  (strong binders only){_mark(0.5)}', compact_num(net_data["n_strong"]), compact_num(flu_data["n_strong"])],
        [f'≤ 2.0  (strong + weak){_mark(2.0)}',       compact_num(net_data["n_weak"]),   compact_num(flu_data["n_weak"])],
    ]
    if threshold not in (0.5, 2.0):
        band_rows.append([
            f'≤ {threshold}  (custom)  [bold yellow]★[/bold yellow]',
            compact_num(net_data["n_0b"]), compact_num(flu_data["n_0b"]),
        ])

    data_table = _data_table(
        headers = ["", "NetMHCpan", "MHCflurry"],
        rows    = base_rows + band_rows,
    )

    panel = _stage_panel(
        title          = "Stage 1 — Affinity threshold filter (%ile ≤ threshold)",
        what_happened  = (
            "Each tool produces one row per [italic]peptide × allele[/italic] prediction with a "
            "percentile score (lower = better binder). We kept only rows whose percentile is "
            f"≤ {threshold} (the threshold marked with ★)."
        ),
        reference      = "Reynisson 2020 (NetMHCpan-4.1, doi:10.1093/nar/gkaa379) · O'Donnell 2020 (MHCflurry 2.0, doi:10.1016/j.cels.2020.06.010)",
        data_table     = data_table,
        border_style   = "cyan",
    )
    console.print(panel)


def _print_stage2_consolidation(
    net_data:           dict,
    flu_data:           dict,
    net_consolidated:   pd.DataFrame,
    flurry_consolidated: pd.DataFrame,
):
    """Stage 2 — consolidation: per-allele rows collapse into one row per peptide."""
    n_net_peptides    = len(net_consolidated)
    n_flurry_peptides = len(flurry_consolidated)

    data_table = _data_table(
        headers = ["", "NetMHCpan", "MHCflurry"],
        rows = [
            ["Rows after threshold",  compact_num(net_data["n_0b"]),                       compact_num(flu_data["n_0b"])],
            ["Unique peptides",       compact_num(n_net_peptides),                         compact_num(n_flurry_peptides)],
            ["Rows collapsed",        compact_num(net_data["n_0b"] - n_net_peptides),      compact_num(flu_data["n_0b"] - n_flurry_peptides)],
        ],
    )

    panel = _stage_panel(
        title          = "Stage 2 — Consolidation (one row per peptide)",
        what_happened  = (
            "The same peptide can be predicted against multiple HLAs, so it appears in several "
            "rows. Here we collapse to one row per peptide, with the bound alleles aggregated "
            "into the `alleles_united` column (semicolon-separated)."
        ),
        example        = (
            "[mono]AETETAHAL[/mono] × {B*18:01, B*39:01, B*40:01, B*44:02} → 4 rows in, 1 row out "
            "(now carries `alleles_united = HLA-B*18:01;HLA-B*39:01;HLA-B*40:01;HLA-B*44:02`)."
        ),
        reference      = None,
        data_table     = data_table,
        border_style   = "green",
    )
    console.print(panel)


def _print_stage3_intersection(intersection_data: dict):
    """Stage 3 — consensus: keep only peptides flagged by BOTH tools."""
    data_table = _data_table(
        headers = ["", "Peptides"],
        rows = [
            ["NetMHCpan only",                   compact_num(intersection_data["net_only"])],
            ["MHCflurry only",                   compact_num(intersection_data["flurry_only"])],
            ["In consensus (both tools)  ★",     f'[bold green]{compact_num(intersection_data["common_count"])}[/bold green]'],
        ],
    )

    panel = _stage_panel(
        title          = "Stage 3 — Tool consensus (NetMHCpan ∩ MHCflurry)",
        what_happened  = (
            "A peptide is kept only if BOTH tools called it a binder. Anything flagged by only "
            "one tool is dropped — this is what makes the result a [italic]consensus[/italic] "
            "and the main reason this pipeline tends to produce a small, conservative shortlist."
        ),
        reference      = None,
        data_table     = data_table,
        border_style   = "yellow",
    )
    console.print(panel)


def _print_stage4_immunogenicity(n_input: int, n_survivors: int):
    """Stage 4 — Calis 2013 immunogenicity score filter."""
    n_discarded     = n_input - n_survivors
    pct_survived    = (n_survivors / n_input * 100) if n_input > 0 else 0.0
    pct_discarded   = 100.0 - pct_survived

    data_table = _data_table(
        headers = ["", "Peptides"],
        rows = [
            ["Input (consensus)",  compact_num(n_input)],
            ["Survivors  ★",       f'[bold green]{compact_num(n_survivors)}[/bold green]  ({pct_survived:.0f}%)'],
            ["Discarded",          f'[dim]{compact_num(n_discarded)}[/dim]  ({pct_discarded:.0f}%)'],
        ],
    )

    panel = _stage_panel(
        title          = "Stage 4 — Immunogenicity filter (Calis 2013, score > 0)",
        what_happened  = (
            "Applied the Calis 2013 immunogenicity score (based on amino-acid properties at the "
            "anchor positions of MHC-I epitopes). Kept only peptides predicted to be T-cell "
            "immunogenic (score > 0)."
        ),
        reference      = "Calis JJ et al., PLoS Comput Biol 2013, doi:10.1371/journal.pcbi.1003266",
        data_table     = data_table,
        border_style   = "magenta",
    )
    console.print(panel)

    console.print(
        Panel(
            Text.from_markup(f"[bold green]✓ Final result: {compact_num(n_survivors)} immunogenic consensus peptides[/bold green]"),
            box=box.HEAVY,
            border_style="green",
            padding=(0, 2),
        )
    )


# Suppress flake8/pyflakes false positive: _NARRATIVE_BLANK kept for future use
_ = _NARRATIVE_BLANK
