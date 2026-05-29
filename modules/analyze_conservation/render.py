"""Rich console rendering for analyze_conservation: the per-track conservation
table and the preflight FASTA-status table."""

import pandas as pd
from rich import box
from rich.table import Table

from utils.console import console
from utils.naming import COLUMN_PEPTIDE

# ── Rich console output ───────────────────────────────────────────────────────

def print_conservation_rich_table(
    df: pd.DataFrame, track_id: str, analysis_threshold: float
):
    thr_pct = int(round(analysis_threshold * 100))
    t = Table(
        box=box.SIMPLE,
        title=f"Conservation: {track_id}  (threshold={thr_pct}%)",
        show_lines=False,
    )
    t.add_column("Peptide",         style="bold", no_wrap=True)
    t.add_column("Variants",        justify="right")
    t.add_column(f"≥{thr_pct}%",    justify="right")
    t.add_column("100%",            justify="right")
    t.add_column("≥90%",            justify="right")
    t.add_column("≥80%",            justify="right")
    t.add_column("Avg ID",          justify="right")
    t.add_column("Label",           no_wrap=True)

    label_rich = {
        "perfect":              "[bold green]★ perfect[/bold green]",
        "high":                 "[green]  high[/green]",
        "moderate":             "[yellow]  moderate[/yellow]",
        "low":                  "[red]  low[/red]",
        "conservation_unknown": "[dim]  unknown[/dim]",
    }

    for _, row in df.iterrows():
        label   = row.get("conservation_label", "conservation_unknown")
        n_total = int(row.get("n_variants_used", 0))

        def fmt(num: int) -> str:
            return f"{num}/{n_total}" if n_total else "-"

        t.add_row(
            str(row.get(COLUMN_PEPTIDE, "")),
            str(n_total),
            fmt(int(row.get("n_passed_threshold", 0))),
            fmt(int(row.get("n_exact_match",      0))),
            fmt(int(row.get("n_identity_90",      0))),
            fmt(int(row.get("n_identity_80",      0))),
            f"{float(row.get('mean_max_identity', 0.0)) * 100:.1f}%",
            label_rich.get(label, label),
        )

    console.print(t)


def _render_track_status_table(track_status_rows: list[dict]) -> Table:
    """Build the consolidated Rich Table shown by `preflight()`."""
    status_table = Table(box=box.SIMPLE, header_style="bold cyan")
    status_table.add_column("Track", style="cyan")
    status_table.add_column("Ref length", justify="right")
    status_table.add_column("Variants file", overflow="fold")
    status_table.add_column("n total", justify="right")
    status_table.add_column("n after ±25% filter", justify="right")
    status_table.add_column("Status")

    status_color_map = {
        "ok":              "[green]ok[/green]",
        "missing_fasta":   "[yellow]missing variants FASTA[/yellow]",
        "no_reference":    "[red]no reference FASTA[/red]",
        "length_mismatch": "[yellow]all variants filtered out (length)[/yellow]",
    }

    for status_row in track_status_rows:
        status_table.add_row(
            status_row["track_id"],
            str(status_row["ref_length"]) if status_row["ref_length"] is not None else "-",
            status_row["fasta_path"].name if status_row["fasta_path"] else "-",
            str(status_row["n_records_raw"])  if status_row["fasta_exists"] else "-",
            str(status_row["n_records_kept"]) if status_row["fasta_exists"] else "-",
            status_color_map.get(status_row["status"], status_row["status"]),
        )

    return status_table


