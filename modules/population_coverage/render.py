"""Rich console output for population_coverage: the per-population summary table."""

from typing import Optional

import pandas as pd
from rich import box
from rich.table import Table

from utils.console import console


def _print_console_summary(
    track_id:    str,
    populations: list[str],
    summary_rows: list[dict],
    n_epitopes:   int,
    cutoff:       Optional[float],
):
    summary_df = pd.DataFrame(summary_rows)

    table = Table(
        title=f"Coverage: {track_id}  ({n_epitopes} ★ epitopes)",
        box=box.SIMPLE,
    )
    table.add_column("Population", style="bold")
    table.add_column("Mean coverage", justify="right")
    table.add_column("Max",           justify="right")
    table.add_column("Min",           justify="right")
    table.add_column("≥ cutoff",      justify="right")
    for population in populations:
        rows_for_pop  = summary_df[summary_df["population"] == population]
        coverage_vals = rows_for_pop["coverage_pct"].astype(float)
        if cutoff is not None:
            above_cutoff_text = f"{int((coverage_vals >= cutoff).sum())}/{len(coverage_vals)}"
        else:
            above_cutoff_text = "-"
        table.add_row(
            population,
            f"{coverage_vals.mean():.2f}%" if len(coverage_vals) else "-",
            f"{coverage_vals.max():.2f}%"  if len(coverage_vals) else "-",
            f"{coverage_vals.min():.2f}%"  if len(coverage_vals) else "-",
            above_cutoff_text,
        )
    console.print(table)
