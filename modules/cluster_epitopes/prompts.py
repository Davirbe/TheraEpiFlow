"""Interactive clustering-parameter prompt for cluster_epitopes."""

from rich import box
from rich.panel import Panel

import config
from utils.console import console
from utils.project_manager import save_project_config

# ── Parameter prompt ──────────────────────────────────────────────────────────

def _ask_clustering_params(project_name: str, project_config: dict) -> tuple:
    """Returns (identity_threshold, clustering_method), prompting once and caching on project_config."""
    saved_threshold = project_config.get("cluster_threshold")
    saved_method    = project_config.get("cluster_method")
    if saved_threshold is not None and saved_method is not None:
        console.print(
            f"[dim]  Clustering params (saved): "
            f"method={saved_method}, threshold={saved_threshold}[/dim]"
        )
        return float(saved_threshold), str(saved_method)

    default_threshold = config.CLUSTER_IDENTITY_CUTOFF

    console.print(Panel(
        "[bold]Step 1 of 2 — Identity threshold[/bold]\n\n"
        "[dim]Epitopes with pairwise sequence identity ≥ threshold are grouped "
        "into the same cluster. One representative is selected per cluster.[/dim]",
        box=box.ROUNDED, title="Setup: cluster_epitopes", title_align="left",
    ))

    while True:
        try:
            raw_threshold = input(
                f"  Identity threshold (0.0–1.0) [{default_threshold}]: "
            ).strip()
            identity_threshold = float(raw_threshold) if raw_threshold else default_threshold
            if 0.0 < identity_threshold <= 1.0:
                break
            console.print("[yellow]  Must be > 0 and ≤ 1.[/yellow]")
        except (ValueError, EOFError):
            identity_threshold = default_threshold
            break

    # Method panel rendered AFTER threshold is confirmed — bundling them confused users.
    console.print(Panel(
        "[bold]Step 2 of 2 — Clustering method[/bold]\n\n"
        "  [cyan][1][/cyan] cluster_break   — cohesive clusters, no bridge groupings [bold](recommended)[/bold]\n"
        "  [cyan][2][/cyan] single_linkage  — most permissive (IEDB-compatible)\n"
        "  [cyan][3][/cyan] clique          — all members pairwise similar (strictest)",
        box=box.ROUNDED, title="Setup: cluster_epitopes", title_align="left",
    ))

    method_choices = {"1": "cluster_break", "2": "single_linkage", "3": "clique"}
    try:
        raw_method = input("  Method [1]: ").strip()
        clustering_method = method_choices.get(raw_method, "cluster_break")
    except EOFError:
        clustering_method = "cluster_break"

    project_config["cluster_threshold"] = identity_threshold
    project_config["cluster_method"]    = clustering_method
    save_project_config(project_name, project_config)

    return identity_threshold, clustering_method


