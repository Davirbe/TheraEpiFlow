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

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    from matplotlib_venn import venn3, venn3_circles  # type: ignore
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
    """Line plot: one (semi-transparent) line per replicate × track; bold median across reps."""
    series: dict[str, list[list[int]]] = {}  # track_id -> [stage_counts_per_rep]
    for record in replicates:
        for track_id, track_payload in record["metrics"].get("tracks", {}).items():
            counts = track_payload.get("stage_counts", {})
            row = [counts.get(stage_key, 0) for stage_key, _ in _STAGE_DISPLAY_ORDER]
            series.setdefault(track_id, []).append(row)

    if not series:
        return

    fig, axes = plt.subplots(1, len(series), figsize=(5.5 * len(series), 4.5), sharey=False)
    if len(series) == 1:
        axes = [axes]

    stage_labels = [label for _, label in _STAGE_DISPLAY_ORDER]
    for ax, (track_id, rep_rows) in zip(axes, series.items()):
        rep_array = np.array(rep_rows, dtype=float)
        x = np.arange(len(stage_labels))
        for r in rep_array:
            ax.plot(x, r, color="steelblue", alpha=0.35, linewidth=1.0)
        median_row = np.median(rep_array, axis=0)
        ax.plot(x, median_row, color="darkorange", linewidth=2.5, label="median")
        ax.set_yscale("symlog", linthresh=10)
        ax.set_xticks(x)
        ax.set_xticklabels(stage_labels, rotation=30, ha="right")
        ax.set_title(f"{preset_key} · {track_id}")
        ax.set_ylabel("peptides")
        ax.legend(loc="upper right", fontsize=9)

    fig.suptitle(f"Loss by stage — {preset_key}", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_time_by_stage(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """Stacked bar: x=replicate, y=seconds, color=step."""
    long_rows: list[dict] = []
    for record in replicates:
        for step_record in record["timings"].get("steps", []):
            step_name = step_record["step_name"]
            for run in step_record.get("runs", []):
                if run.get("elapsed_seconds") is None:
                    continue
                long_rows.append({
                    "rep":     record["rep_label"],
                    "step":    step_name,
                    "seconds": float(run["elapsed_seconds"]),
                })
    if not long_rows:
        return

    df = pd.DataFrame(long_rows)
    pivot = df.pivot_table(index="rep", columns="step", values="seconds", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(sorted(pivot.index))

    fig, ax = plt.subplots(figsize=(11, 5))
    pivot.plot.bar(stacked=True, ax=ax, colormap="tab20", width=0.85)
    ax.set_xlabel("replicate")
    ax.set_ylabel("seconds")
    ax.set_title(f"Per-step wall time — {preset_key}")
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, title="step")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_consistency_matrix(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """Heatmap: rows = union of ★ peptides across reps, cols = rep index, cells = present?"""
    for_each_track: dict[str, dict[str, list[bool]]] = {}
    for rep_index, record in enumerate(replicates, start=1):
        for track_id, track_payload in record["metrics"].get("tracks", {}).items():
            star_peptides = set(track_payload.get("stages", {}).get("cluster_reps_star", []))
            container = for_each_track.setdefault(track_id, {})
            for peptide in star_peptides:
                container.setdefault(peptide, [False] * len(replicates))
                container[peptide][rep_index - 1] = True
            # back-fill rows for peptides discovered later
            for peptide_row in container.values():
                while len(peptide_row) < len(replicates):
                    peptide_row.append(False)

    if not for_each_track:
        return

    fig, axes = plt.subplots(1, len(for_each_track), figsize=(0.5 * len(replicates) * len(for_each_track) + 2, max(4, 0.25 * max((len(t) for t in for_each_track.values()), default=10))), squeeze=False)

    for ax_idx, (track_id, peptide_presence) in enumerate(for_each_track.items()):
        sorted_peptides = sorted(peptide_presence.keys())
        matrix = np.array(
            [[1 if peptide_presence[p][r] else 0 for r in range(len(replicates))]
             for p in sorted_peptides],
            dtype=int,
        )
        ax = axes[0][ax_idx]
        sns.heatmap(
            matrix, ax=ax, cmap=["#f5f5f5", "#2e7d32"], cbar=False,
            xticklabels=[f"r{r+1:02d}" for r in range(len(replicates))],
            yticklabels=sorted_peptides, linewidths=0.3, linecolor="white",
        )
        ax.set_title(f"{preset_key} · {track_id}\n{matrix.shape[0]} unique ★ peptides")
        ax.set_xlabel("replicate")
        ax.set_ylabel("peptide")

    fig.suptitle(f"Replicate consistency of ★ peptides — {preset_key}", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_venn_consensus(replicates: list[dict], preset_key: str, out_path: Path) -> None:
    """3-circle Venn: NetMHCpan survivors, MHCflurry survivors, calis-passed.
    Uses the means across replicates for each set size."""
    if not _VENN_AVAILABLE:
        return

    net_sizes:    list[int] = []
    flurry_sizes: list[int] = []
    calis_sizes:  list[int] = []
    inter_net_flurry:  list[int] = []
    inter_net_calis:   list[int] = []
    inter_flurry_calis: list[int] = []
    inter_all:         list[int] = []

    for record in replicates:
        for track_payload in record["metrics"].get("tracks", {}).values():
            stages = track_payload.get("stages", {})
            net    = set(stages.get("netmhcpan_survivors", []))
            flurry = set(stages.get("mhcflurry_survivors", []))
            calis  = set(stages.get("immunogenic_calis", []))
            net_sizes.append(len(net))
            flurry_sizes.append(len(flurry))
            calis_sizes.append(len(calis))
            inter_net_flurry.append(len(net & flurry))
            inter_net_calis.append(len(net & calis))
            inter_flurry_calis.append(len(flurry & calis))
            inter_all.append(len(net & flurry & calis))

    if not net_sizes:
        return

    # Convert to the seven Venn-3 sub-regions:
    # (100=net_only, 010=flurry_only, 001=calis_only, 110, 101, 011, 111)
    mean_net    = np.mean(net_sizes)
    mean_flurry = np.mean(flurry_sizes)
    mean_calis  = np.mean(calis_sizes)
    mean_nf     = np.mean(inter_net_flurry)
    mean_nc     = np.mean(inter_net_calis)
    mean_fc     = np.mean(inter_flurry_calis)
    mean_all    = np.mean(inter_all)

    abc = max(mean_all, 0)
    ab  = max(mean_nf - abc, 0)
    ac  = max(mean_nc - abc, 0)
    bc  = max(mean_fc - abc, 0)
    a   = max(mean_net    - ab - ac - abc, 0)
    b   = max(mean_flurry - ab - bc - abc, 0)
    c   = max(mean_calis  - ac - bc - abc, 0)

    fig, ax = plt.subplots(figsize=(7, 6))
    venn3(
        subsets=(a, b, ab, c, ac, bc, abc),
        set_labels=("NetMHCpan ≤ 2%", "MHCflurry ≤ 2%", "Calis > 0"),
        ax=ax,
    )
    venn3_circles(
        subsets=(a, b, ab, c, ac, bc, abc), linewidth=0.8, ax=ax,
    )
    ax.set_title(f"Consensus funnel — {preset_key} (mean of {len(net_sizes)} track-replicate runs)")
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
        for index, peptide_row in enumerate(df.itertuples()):
            ax.text(
                peptide_row.total + 0.5, index, f"{peptide_row.presence}/10",
                va="center", fontsize=8, color="#444",
            )
        ax.set_xlabel("assays in IEDB (Human, MHC-I)")
        ax.set_title(f"Literature support per ★ peptide — {track_id}")
        ax.legend(loc="lower right", fontsize=9)
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
    sns.violinplot(
        data=df, x="track", y="hits", hue="category",
        split=True, inner="quartile",
        palette={"★ final": "#1f77b4", "negative control": "#d62728"},
        ax=ax, cut=0,
    )
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("IEDB assays (T-cell + MHC-I)")
    ax.set_xlabel("")
    ax.set_title("Clinical literature support — ★ vs negative control")
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
