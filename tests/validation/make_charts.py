"""Chart generators for the validation suite.

Reads the per-replicate JSONs produced by `run_experiment_1.py` and
`run_experiment_2.py` and renders publication-ready PNGs into the figures
directory (tracked in git so the paper has stable artifacts).

Charts produced:
  Experiment 1 (one PNG per preset):
    exp1_loss_by_stage_{preset}.png     — line plot, 10 reps overlaid + median bold
    exp1_time_by_stage_{preset}.png     — stacked bar, x=rep, y=seconds, color=step
    exp1_consistency_matrix_{preset}.png — heatmap, rows=★ peptides, cols=rep
    exp1_venn_{preset}.png              — 3-circle Venn of consensus stages (averaged)
  Experiment 1 (cross-preset):
    exp1_time_tr1_vs_tr3.png            — paired bars, mean total time
    exp1_time_protein_size.png          — paired bars across organism sizes
  Experiment 2 (one PNG per track):
    exp2_iedb_counts_{track}.png        — horizontal bars, IEDB hits per ★ peptide
    exp2_neg_control_comparison.png     — violin/box, ★ vs negative control

Usage:
    python -m tests.validation.make_charts --exp1-root tests/validation/results/exp1
    python -m tests.validation.make_charts --exp2-root tests/validation/results/exp2
    python -m tests.validation.make_charts \\
        --exp1-root tests/validation/results/exp1 \\
        --exp2-root tests/validation/results/exp2 \\
        --figures-root tests/validation/figures
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    from matplotlib_venn import venn2, venn2_circles  # type: ignore
    _VENN_AVAILABLE = True
except ImportError:
    _VENN_AVAILABLE = False

sns.set_style("whitegrid")
plt.rcParams.update({
    "figure.dpi":      120,
    "savefig.dpi":     200,
    "font.size":       10,
    "axes.titlesize":  12,
    "axes.labelsize":  10,
    "figure.figsize":  (9, 5),
})


_STAGE_DISPLAY_ORDER: list[tuple[str, str]] = [
    ("raw_netmhcpan",          "raw"),
    ("netmhcpan_survivors",    "netMHCpan"),
    ("mhcflurry_survivors",    "MHCflurry"),
    ("consensus_intersection", "consensus"),
    ("immunogenic_calis",      "calis"),
    ("toxicity_safe",          "toxicity"),
    ("cluster_reps_star",      "★ reps"),
]

# Loss-by-stage chart uses display labels suited for the paper and drops the
# `raw_netmhcpan` stage (the user wants the funnel to start at the first filter).
_LOSS_STAGE_DISPLAY_ORDER: list[tuple[str, str]] = [
    ("netmhcpan_survivors",    "Prediction netMHCpan"),
    ("mhcflurry_survivors",    "Prediction MHCflurry"),
    ("consensus_intersection", "Consensus"),
    ("immunogenic_calis",      "Immunogenecity"),
    ("toxicity_safe",          "Toxicity"),
    ("cluster_reps_star",      "Representatives of the clusters"),
]

# Per-track step canonical order (used by time_by_stage to drive per-track
# subplots in multi-track presets). Mirrors the per-track block in
# step_registry.STEP_REGISTRY. Global steps (integrate_data, generate_report,
# export_bundle) are EXCLUDED here on purpose: when we split the chart into
# one subplot per track, the globals (which run once per replicate, not per
# track) get a small footnote instead of a misleading per-track bar.
_TRACK_STEP_ORDER: list[str] = [
    "fetch_sequences",
    "predict_binding",
    "consensus_filter",
    "screen_toxicity",
    "cluster_epitopes",
    "select_representatives",
    "search_variants",
    "analyze_conservation",
    "population_coverage",
    "predict_murine",
    "curate_murine",
]
_GLOBAL_STEPS: set[str] = {"integrate_data", "generate_report", "export_bundle"}

# Fixed colour per pipeline step so the same step is always the same colour
# across every chart that renders timings (tr1 single-panel, tr3 multi-panel,
# any future protein-size or summary chart). Without this dict, tr1 picks
# colours from the larger step set (incl. globals) and tr3 picks from the
# smaller per-track set, so `predict_binding` would be a different shade
# in each chart and the eye loses the comparison.
_STEP_COLOR_MAP: dict[str, str] = {
    "fetch_sequences":        "#1f77b4",
    "predict_binding":        "#ff7f0e",
    "consensus_filter":       "#2ca02c",
    "screen_toxicity":        "#d62728",
    "cluster_epitopes":       "#9467bd",
    "select_representatives": "#8c564b",
    "search_variants":        "#e377c2",
    "analyze_conservation":   "#7f7f7f",
    "population_coverage":    "#bcbd22",
    "predict_murine":         "#17becf",
    "curate_murine":          "#aec7e8",
    "integrate_data":         "#ffbb78",
    "generate_report":        "#98df8a",
    "export_bundle":          "#ff9896",
}


# ── Generic IO ────────────────────────────────────────────────────────────────

def _load_replicate_jsons(exp_root: Path, preset_subdir: Optional[str]) -> list[dict]:
    """Load every `metrics.json` + `timings.json` pair under exp_root[/preset_subdir]/repNN/."""
    base = exp_root if preset_subdir is None else exp_root / preset_subdir
    if not base.exists():
        return []
    records: list[dict] = []
    for rep_dir in sorted(base.glob("rep*")):
        if not rep_dir.is_dir():
            continue
        metrics_path = rep_dir / "metrics.json"
        timings_path = rep_dir / "timings.json"
        metrics_payload = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.exists() else {}
        timings_payload = json.loads(timings_path.read_text(encoding="utf-8")) if timings_path.exists() else {}
        records.append({
            "rep_dir":   rep_dir,
            "rep_label": rep_dir.name,
            "metrics":   metrics_payload,
            "timings":   timings_payload,
        })
    return records


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


# ── Experiment 1 per-preset charts ────────────────────────────────────────────

def plot_loss_by_stage(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """Single funnel line per panel — peptide count at each filter stage.

    The pipeline is fully deterministic per track when the network succeeds:
    every non-failed replicate yields the *same* stage_counts vector. Showing
    10 overlapping semi-transparent lines plus an aggregate ("mean") just adds
    visual noise. Instead, plot one solid orange line with markers at the six
    funnel stages, using the canonical (= first non-failed) replicate's vector.

    For multi-track presets (tr3), the 3 subplots share both axes so the
    reader can compare the three predictions on the same y-scale.
    """
    series: dict[str, list[list[int]]] = {}  # track_id -> [stage_counts_per_rep]
    for record in replicates:
        for track_id, track_payload in record["metrics"].get("tracks", {}).items():
            counts = track_payload.get("stage_counts", {})
            row = [counts.get(stage_key, 0) for stage_key, _ in _LOSS_STAGE_DISPLAY_ORDER]
            series.setdefault(track_id, []).append(row)

    if not series:
        return

    fig, axes = plt.subplots(
        1, len(series),
        figsize=(5.5 * len(series), 4.5),
        sharex=True, sharey=True,
    )
    if len(series) == 1:
        axes = [axes]

    stage_labels = [label for _, label in _LOSS_STAGE_DISPLAY_ORDER]
    track_order = list(series.keys())
    for subplot_idx, (ax, track_id) in enumerate(zip(axes, track_order)):
        rep_rows = series[track_id]
        rep_array = np.array(rep_rows, dtype=float)
        x = np.arange(len(stage_labels))

        # Canonical funnel: the first non-failed replicate's vector. With a
        # deterministic pipeline, this equals every successful replicate.
        canonical_row = next((r for r in rep_array if r.sum() > 0), None)
        if canonical_row is not None:
            ax.plot(
                x, canonical_row,
                color="darkorange", linewidth=2.5, marker="o", markersize=8,
            )
            for x_value, y_value in zip(x, canonical_row):
                ax.annotate(
                    f"{int(y_value)}",
                    xy=(x_value, y_value),
                    xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, color="#444",
                )

        ax.set_xticks(x)
        ax.set_xticklabels(stage_labels, rotation=30, ha="right")
        ax.set_title(f"Prediction {subplot_idx + 1}", fontweight="bold")
        ax.set_xlabel("Steps", fontweight="bold")
        if subplot_idx == 0:
            ax.set_ylabel("peptides", fontweight="bold")
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_time_by_stage(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """Stacked bar of per-step wall time per replicate.

    - Single-track presets (tr1): one panel, all steps including globals stacked.
    - Multi-track presets (tr3): N panels (one per track), each labeled
      `Prediction 1/2/3`. Global steps are excluded from these panels (they run
      once per replicate, not per track) and noted in a footnote below the
      figure to avoid misleading the reader.
    """
    # Collect per-(rep, track, step) seconds AND track-list ordering.
    long_rows: list[dict] = []
    track_ids_seen: list[str] = []
    track_seen_set: set[str] = set()
    for record in replicates:
        for step_record in record["timings"].get("steps", []):
            step_name = step_record["step_name"]
            for run in step_record.get("runs", []):
                if run.get("elapsed_seconds") is None:
                    continue
                track_id = run.get("track_id")  # None for global steps
                if track_id is not None and track_id not in track_seen_set:
                    track_seen_set.add(track_id)
                    track_ids_seen.append(track_id)
                long_rows.append({
                    "rep":      record["rep_label"],
                    "step":     step_name,
                    "track_id": track_id,
                    "seconds":  float(run["elapsed_seconds"]),
                })
    if not long_rows:
        return

    df = pd.DataFrame(long_rows)
    has_multiple_tracks = len(track_ids_seen) > 1

    # Stable step ordering — full registry order, then any unknown steps trailing.
    def _ordered_columns(columns: list[str]) -> list[str]:
        registry_then_globals = _TRACK_STEP_ORDER + [
            "integrate_data", "generate_report", "export_bundle",
        ]
        ordered = [s for s in registry_then_globals if s in columns]
        leftover = [s for s in columns if s not in ordered]
        return ordered + leftover

    def _colors_for(columns: list[str]) -> list[str]:
        # _STEP_COLOR_MAP guarantees the same step → same colour across tr1 and tr3.
        return [_STEP_COLOR_MAP.get(step_name, "#cccccc") for step_name in columns]

    if not has_multiple_tracks:
        # Single panel — all steps stacked (track-level + globals).
        pivot = df.pivot_table(index="rep", columns="step", values="seconds", aggfunc="sum", fill_value=0)
        pivot = pivot.reindex(sorted(pivot.index))
        pivot = pivot.reindex(columns=_ordered_columns(list(pivot.columns)))

        fig, ax = plt.subplots(figsize=(11, 5))
        pivot.plot.bar(
            stacked=True, ax=ax,
            color=_colors_for(list(pivot.columns)),
            width=0.85,
        )
        ax.set_xlabel("replicate", fontweight="bold")
        ax.set_ylabel("seconds",   fontweight="bold")
        for tick in ax.get_xticklabels():
            tick.set_rotation(0)
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, title="step")
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        return

    # Multi-track preset: one panel per track, globals excluded with footnote.
    per_track_df = df[df["track_id"].notna()].copy()

    fig, axes = plt.subplots(1, len(track_ids_seen), figsize=(5.5 * len(track_ids_seen), 5), sharey=True)
    if len(track_ids_seen) == 1:
        axes = [axes]

    for subplot_idx, (ax, track_id) in enumerate(zip(axes, track_ids_seen)):
        track_df = per_track_df[per_track_df["track_id"] == track_id]
        if track_df.empty:
            continue
        pivot = track_df.pivot_table(
            index="rep", columns="step", values="seconds", aggfunc="sum", fill_value=0,
        )
        pivot = pivot.reindex(sorted(pivot.index))
        pivot = pivot.reindex(columns=_ordered_columns(list(pivot.columns)))

        pivot.plot.bar(
            stacked=True, ax=ax,
            color=_colors_for(list(pivot.columns)),
            width=0.85, legend=(subplot_idx == 0),
        )
        ax.set_title(f"Prediction {subplot_idx + 1}", fontweight="bold")
        ax.set_xlabel("replicate", fontweight="bold")
        ax.set_ylabel("seconds" if subplot_idx == 0 else "", fontweight="bold" if subplot_idx == 0 else "normal")
        for tick in ax.get_xticklabels():
            tick.set_rotation(0)
        if subplot_idx == 0:
            ax.legend(loc="center left", bbox_to_anchor=(-0.6, 0.5), fontsize=8, title="step")

    fig.text(
        0.5, 0.005,
        "Global steps (integrate_data, generate_report, export_bundle) ran once per replicate (<1s total) and are excluded.",
        ha="center", fontsize=8, style="italic", color="#555",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_consistency_matrix(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """Heatmap: rows = union of ★ peptides across reps, cols = rep index.

    Cells: green = peptide present in that replicate, light grey = peptide absent,
    red = the replicate's pipeline run produced ZERO ★ peptides for this track
    (per-track failure flag — covers both `raw==0` and the anomalous case where
    NetMHCpan returned data but the intersection collapsed to zero).
    """
    for_each_track: dict[str, dict[str, list[bool]]] = {}
    for rep_index, record in enumerate(replicates, start=1):
        for track_id, track_payload in record["metrics"].get("tracks", {}).items():
            star_peptides = set(track_payload.get("stages", {}).get("cluster_reps_star", []))
            container = for_each_track.setdefault(track_id, {})
            for peptide in star_peptides:
                container.setdefault(peptide, [False] * len(replicates))
                container[peptide][rep_index - 1] = True
            for peptide_row in container.values():
                while len(peptide_row) < len(replicates):
                    peptide_row.append(False)

    if not for_each_track:
        return

    # Per-track failure flag: this replicate produced zero ★ peptides for this track.
    # Catches both the "raw NetMHCpan empty" case (r04/r05 of hpv16_e7_tr1) AND the
    # anomalous "raw=721 but intersection=0" case (r09 H16A_E7 of hpv16_e7_tr3).
    per_track_failed: dict[str, set[int]] = {tid: set() for tid in for_each_track}
    for rep_index, record in enumerate(replicates, start=1):
        for track_id, track_payload in record["metrics"].get("tracks", {}).items():
            counts = track_payload.get("stage_counts", {})
            if counts.get("cluster_reps_star", 0) == 0:
                per_track_failed.setdefault(track_id, set()).add(rep_index)

    fig, axes = plt.subplots(
        1, len(for_each_track),
        figsize=(
            0.5 * len(replicates) * len(for_each_track) + 2,
            max(4, 0.25 * max((len(t) for t in for_each_track.values()), default=10)),
        ),
        squeeze=False,
    )

    track_order = list(for_each_track.keys())
    for ax_idx, track_id in enumerate(track_order):
        peptide_presence = for_each_track[track_id]
        sorted_peptides = sorted(peptide_presence.keys())
        failed_for_this_track = per_track_failed.get(track_id, set())

        # Cells: 0 = absent (light grey), 1 = present (green), 2 = failed (red).
        matrix = np.zeros((len(sorted_peptides), len(replicates)), dtype=int)
        for row_idx, peptide in enumerate(sorted_peptides):
            for col_idx in range(len(replicates)):
                if (col_idx + 1) in failed_for_this_track:
                    matrix[row_idx, col_idx] = 2
                elif peptide_presence[peptide][col_idx]:
                    matrix[row_idx, col_idx] = 1
        ax = axes[0][ax_idx]
        sns.heatmap(
            matrix, ax=ax, cmap=["#f5f5f5", "#2e7d32", "#c62828"], cbar=False, vmin=0, vmax=2,
            xticklabels=[f"r{r+1:02d}" for r in range(len(replicates))],
            yticklabels=sorted_peptides, linewidths=0.3, linecolor="white",
        )
        ax.set_title(f"Prediction {ax_idx + 1}", fontweight="bold")
        ax.set_xlabel("replicate", fontweight="bold")
        ax.set_ylabel("peptide",   fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_venn_consensus(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """2-circle Venn: NetMHCpan survivors ∩ MHCflurry survivors.

    By definition, the consensus filter output equals NetMHCpan ∩ MHCflurry, so
    the central overlap region IS the Consensus set. We label it explicitly as
    "Consensus = N" instead of forcing a third circle that would be entirely
    contained in the overlap and add no information.
    """
    if not _VENN_AVAILABLE:
        return

    net_sizes:        list[int] = []
    flurry_sizes:     list[int] = []
    intersection_sizes: list[int] = []

    for record in replicates:
        for track_payload in record["metrics"].get("tracks", {}).values():
            stages = track_payload.get("stages", {})
            net    = set(stages.get("netmhcpan_survivors", []))
            flurry = set(stages.get("mhcflurry_survivors", []))
            net_sizes.append(len(net))
            flurry_sizes.append(len(flurry))
            intersection_sizes.append(len(net & flurry))

    if not net_sizes:
        return

    mean_net      = float(np.mean(net_sizes))
    mean_flurry   = float(np.mean(flurry_sizes))
    mean_consensus = float(np.mean(intersection_sizes))

    # Venn-2 sub-regions: (10 = Net only, 01 = Flurry only, 11 = both)
    net_only    = max(mean_net    - mean_consensus, 0)
    flurry_only = max(mean_flurry - mean_consensus, 0)
    overlap     = max(mean_consensus, 0)

    # Semantic colours — the legend patches MUST equal the actual circle
    # colours, otherwise the reader can't trace which circle is which set.
    color_netmhcpan = "#1f78b4"  # blue
    color_mhcflurry = "#33a02c"  # green
    color_consensus = "#999999"  # gray (overlap region = Consensus)

    fig, ax = plt.subplots(figsize=(7, 6))
    diagram = venn2(
        subsets=(net_only, flurry_only, overlap),
        set_labels=("NetMHCpan ≤ 2%", "MHCflurry ≤ 2%"),
        ax=ax,
    )

    # Override matplotlib_venn's default colours so each region matches its
    # legend swatch. set_color() also resets alpha; reapply afterwards.
    region_colors = {"10": color_netmhcpan, "01": color_mhcflurry, "11": color_consensus}
    for region_id, color in region_colors.items():
        patch = diagram.get_patch_by_id(region_id) if diagram else None
        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(0.55)
            patch.set_edgecolor("none")

    # Format value labels: one decimal + relabel the centre as "Consensus = N".
    for region_id in ("10", "01", "11"):
        label = diagram.get_label_by_id(region_id) if diagram else None
        if label is None:
            continue
        try:
            value_text = f"{float(label.get_text()):.1f}"
        except (TypeError, ValueError):
            value_text = label.get_text()
        if region_id == "11":
            label.set_text(f"Consensus = {value_text}")
        else:
            label.set_text(value_text)

    venn2_circles(
        subsets=(net_only, flurry_only, overlap), linewidth=0.8, ax=ax,
    )

    # Legend uses the same colours we forced on the patches.
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=color_netmhcpan, edgecolor="none", alpha=0.55, label=f"NetMHCpan: {mean_net:.1f}"),
        Patch(facecolor=color_mhcflurry, edgecolor="none", alpha=0.55, label=f"MHCflurry: {mean_flurry:.1f}"),
        Patch(facecolor=color_consensus, edgecolor="none", alpha=0.55, label=f"Consensus: {mean_consensus:.1f}"),
    ]
    ax.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=3, fontsize=9)
    ax.set_title(f"Consensus funnel — {preset_key}")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ── Experiment 1 cross-preset charts ──────────────────────────────────────────

def plot_time_tr1_vs_tr3(exp1_root: Path, figures_root: Path) -> None:
    """Paired bars: mean total seconds, tr1 vs tr3 per organism."""
    pairs = [
        ("HPV16_E7",  "hpv16_e7_tr1",  "hpv16_e7_tr3"),
        ("SARS_NCAP", "sars_nucleo_tr1", "sars_nucleo_tr3"),
    ]
    rows: list[dict] = []
    for organism_label, tr1_key, tr3_key in pairs:
        for layout_label, preset_key in (("1 track", tr1_key), ("3 tracks", tr3_key)):
            replicates = _load_replicate_jsons(exp1_root, preset_key)
            totals = [
                rec["timings"].get("total_elapsed_seconds")
                for rec in replicates
                if rec["timings"].get("total_elapsed_seconds") is not None
            ]
            if not totals:
                continue
            rows.append({
                "organism": organism_label,
                "layout":   layout_label,
                "mean_s":   float(np.mean(totals)),
                "stddev_s": float(np.std(totals, ddof=1)) if len(totals) > 1 else 0.0,
            })
    if not rows:
        return

    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.barplot(
        data=df, x="organism", y="mean_s", hue="layout",
        palette={"1 track": "#4c78a8", "3 tracks": "#f58518"}, ax=ax,
    )
    for idx, row in df.iterrows():
        ax.errorbar(
            x=idx // 2 + (-0.2 if row["layout"] == "1 track" else 0.2),
            y=row["mean_s"], yerr=row["stddev_s"], color="black", capsize=4, fmt="none",
        )
    ax.set_ylabel("mean total seconds")
    ax.set_xlabel("")
    ax.set_title("Pipeline wall time — tr1 vs tr3")
    fig.tight_layout()
    fig.savefig(figures_root / "exp1_time_tr1_vs_tr3.png", bbox_inches="tight")
    plt.close(fig)


def plot_time_protein_size(exp1_root: Path, figures_root: Path) -> None:
    """Paired bars: HPV16_E7 (98 aa) vs SARS_NCAP (419 aa), for each track count."""
    pairs = [
        ("1 track",  "hpv16_e7_tr1",   "sars_nucleo_tr1"),
        ("3 tracks", "hpv16_e7_tr3",   "sars_nucleo_tr3"),
    ]
    rows: list[dict] = []
    for layout_label, hpv_key, sars_key in pairs:
        for protein_label, preset_key in (("HPV16_E7 (98 aa)", hpv_key), ("SARS_NCAP (419 aa)", sars_key)):
            replicates = _load_replicate_jsons(exp1_root, preset_key)
            totals = [
                rec["timings"].get("total_elapsed_seconds")
                for rec in replicates
                if rec["timings"].get("total_elapsed_seconds") is not None
            ]
            if not totals:
                continue
            rows.append({
                "layout":   layout_label,
                "protein":  protein_label,
                "mean_s":   float(np.mean(totals)),
                "stddev_s": float(np.std(totals, ddof=1)) if len(totals) > 1 else 0.0,
            })
    if not rows:
        return

    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.barplot(data=df, x="layout", y="mean_s", hue="protein", ax=ax)
    ax.set_ylabel("mean total seconds")
    ax.set_xlabel("")
    ax.set_title("Pipeline wall time — protein-size effect")
    fig.tight_layout()
    fig.savefig(figures_root / "exp1_time_protein_size.png", bbox_inches="tight")
    plt.close(fig)


# ── Experiment 2 charts ───────────────────────────────────────────────────────

def plot_iedb_counts(exp2_root: Path, figures_root: Path) -> None:
    """Horizontal bar plots — one per track."""
    iedb_path = exp2_root / "iedb_hits.json"
    if not iedb_path.exists():
        return
    payload = json.loads(iedb_path.read_text(encoding="utf-8"))

    for track_id, track_payload in payload.get("tracks", {}).items():
        rows = []
        for record in track_payload.get("star", []):
            rows.append({
                "peptide":           record["peptide"],
                "tcell_assay_count": record["tcell_assay_count"],
                "mhc_assay_count":   record["mhc_assay_count"],
                "total":             record["tcell_assay_count"] + record["mhc_assay_count"],
                "presence":          record.get("presence_count", 0),
            })
        if not rows:
            continue
        df = pd.DataFrame(rows).sort_values("total", ascending=True)

        fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(df))))
        ax.barh(df["peptide"], df["tcell_assay_count"], color="#1f77b4", label="T-cell assays")
        ax.barh(df["peptide"], df["mhc_assay_count"], left=df["tcell_assay_count"], color="#ff7f0e", label="MHC-I assays")
        ax.set_xlabel("Assays in IEDB", fontweight="bold")
        ax.set_ylabel("peptide", fontweight="bold")
        ax.set_title(f"Literature support per ★ peptide — {track_id}")
        ax.legend(loc="lower right", fontsize=9)
        ax.grid(False)
        # Denser, integer X ticks so the reader can read each peptide's count
        # without squinting (default matplotlib picks ~5-6 ticks; ~15 is more useful).
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=15))
        sns.despine(ax=ax)
        fig.tight_layout()
        fig.savefig(figures_root / f"exp2_iedb_counts_{track_id}.png", bbox_inches="tight")
        plt.close(fig)


def plot_neg_control_comparison(exp2_root: Path, figures_root: Path) -> None:
    """Violin/box plot of literature hits, ★ vs negative-control."""
    iedb_path = exp2_root / "iedb_hits.json"
    if not iedb_path.exists():
        return
    payload = json.loads(iedb_path.read_text(encoding="utf-8"))

    rows: list[dict] = []
    for track_id, track_payload in payload.get("tracks", {}).items():
        for record in track_payload.get("star", []):
            rows.append({
                "track":    track_id,
                "category": "★ final",
                "hits":     record["tcell_assay_count"] + record["mhc_assay_count"],
            })
        for record in track_payload.get("negative_control", []):
            rows.append({
                "track":    track_id,
                "category": "negative control",
                "hits":     record["tcell_assay_count"] + record["mhc_assay_count"],
            })
    if not rows:
        return

    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 5))

    category_palette = {"★ final": "#1f77b4", "negative control": "#d62728"}

    # Box plot summary + strip plot overlay so each individual peptide is visible.
    sns.boxplot(
        data=df, x="track", y="hits", hue="category",
        palette=category_palette, ax=ax,
        fliersize=0,            # outliers are shown by the strip plot
        linewidth=1.2,
    )
    sns.stripplot(
        data=df, x="track", y="hits", hue="category",
        dodge=True, alpha=0.7, size=4, color="black",
        ax=ax, legend=False,
    )

    ax.set_ylabel("IEDB assays (T-cell + MHC-I)", fontweight="bold")
    ax.set_xlabel("track", fontweight="bold")
    ax.set_title("Clinical literature support — ★ vs negative control")
    ax.set_ylim(bottom=0)
    ax.grid(False)
    sns.despine(ax=ax)

    # Drop the duplicated stripplot legend handles so the legend stays clean.
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        seen = set()
        deduped = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
        ax.legend([h for h, _ in deduped], [l for _, l in deduped], loc="upper right")

    fig.tight_layout()
    fig.savefig(figures_root / "exp2_neg_control_comparison.png", bbox_inches="tight")
    plt.close(fig)


# ── Orchestrator ──────────────────────────────────────────────────────────────

def render_exp1(exp1_root: Path, figures_root: Path) -> None:
    if not exp1_root.exists():
        return
    _ensure_dir(figures_root)
    preset_dirs = [d for d in exp1_root.iterdir() if d.is_dir() and d.name != "_adhoc"]
    for preset_dir in preset_dirs:
        preset_key = preset_dir.name
        replicates = _load_replicate_jsons(exp1_root, preset_key)
        if not replicates:
            continue
        plot_loss_by_stage(replicates,         preset_key, figures_root / f"exp1_loss_by_stage_{preset_key}.png")
        plot_time_by_stage(replicates,         preset_key, figures_root / f"exp1_time_by_stage_{preset_key}.png")
        plot_consistency_matrix(replicates,    preset_key, figures_root / f"exp1_consistency_matrix_{preset_key}.png")
        plot_venn_consensus(replicates,        preset_key, figures_root / f"exp1_venn_{preset_key}.png")
    plot_time_tr1_vs_tr3(exp1_root, figures_root)
    plot_time_protein_size(exp1_root, figures_root)


def render_exp2(exp2_root: Path, figures_root: Path) -> None:
    if not exp2_root.exists():
        return
    _ensure_dir(figures_root)
    plot_iedb_counts(exp2_root, figures_root)
    plot_neg_control_comparison(exp2_root, figures_root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp1-root", type=Path, default=None)
    parser.add_argument("--exp2-root", type=Path, default=None)
    parser.add_argument(
        "--figures-root", type=Path, default=Path("tests/validation/figures"),
        help="Destination for PNGs (kept under git so the paper has stable artifacts).",
    )
    args = parser.parse_args()

    if not (args.exp1_root or args.exp2_root):
        parser.error("Pass at least --exp1-root or --exp2-root.")

    if args.exp1_root:
        render_exp1(args.exp1_root, args.figures_root)
    if args.exp2_root:
        render_exp2(args.exp2_root, args.figures_root)

    print(f"Done. Figures written under {args.figures_root}")


if __name__ == "__main__":
    main()
