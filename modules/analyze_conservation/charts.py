"""PNG generation for analyze_conservation: the dual-panel position-conservation
+ identity-tier figure. Matplotlib is imported lazily inside the writer."""

from pathlib import Path

from .core import (
    _align_position_stats_to_anchor,
    _anchor_positions,
    _build_anchor_aligned_xtick_labels,
    _compute_position_stats,
    _has_variant_within_max_mut,
)

_LABEL_HEX = {
    "perfect":              "#00B050",
    "high":                 "#92D050",
    "moderate":             "#FFFF99",
    "low":                  "#FF9999",
    "conservation_unknown": "#D9D9D9",
}
_ANCHOR_BORDER_HEX = "#2E75B6"
# ── Dual-panel PNG (position conservation + identity tiers) ──────────────────

def write_conservation_dual_panel_png(
    alignment_data: list[dict],
    track_id: str,
    n_variants_total: int,
    analysis_threshold: float,
    output_path: Path,
):
    """
    Single figure with three axes sharing the Y axis (peptides):
      [label sidebar] [position conservation P1..PΩ] [identity tiers]

    Filter: epitopes with ≥1 variant having ≤_MAX_MUTATIONS_FOR_HEATMAP
    substitutions. This includes 100%-conserved epitopes (0 muts ≤ 2).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    import seaborn as sns

    filtered = [
        item for item in alignment_data
        if _has_variant_within_max_mut(item["peptide"], item["alignment_tuples"])
    ]
    if not filtered or n_variants_total == 0:
        return

    enriched: list[dict] = []
    for item in filtered:
        peptide              = item["peptide"]
        pos_stats            = _compute_position_stats(peptide, item["alignment_tuples"])
        n_positions          = len(peptide)
        p2_idx, pomega_idx   = _anchor_positions(n_positions)
        p2_pct               = pos_stats[p2_idx]["conservation_pct"]     if len(pos_stats) > p2_idx     else 0.0
        pomega_pct           = pos_stats[pomega_idx]["conservation_pct"] if len(pos_stats) > pomega_idx else 0.0
        enriched.append({
            "peptide":            peptide,
            "pos_stats":          pos_stats,
            "anchor_score":       min(p2_pct, pomega_pct),
            "n_positions":        n_positions,
            "conservation_label": item["conservation_label"],
            "n_exact_match":      item["n_exact_match"],
            "n_identity_90":      item["n_identity_90"],
            "n_identity_80":      item["n_identity_80"],
            "n_passed_threshold": item["n_passed_threshold"],
        })

    enriched.sort(key=lambda x: x["anchor_score"], reverse=True)

    n_rows         = len(enriched)
    n_positions    = max(item["n_positions"] for item in enriched)
    p2_idx, pomega_idx = _anchor_positions(n_positions)
    threshold_pct  = int(round(analysis_threshold * 100))

    # Anchor-aligned matrix: every peptide's P2 sits at col 1 and PΩ at the
    # last col, regardless of its length. The middle bulge becomes NaN.
    pos_matrix   = np.full((n_rows, n_positions), np.nan, dtype=float)
    pos_annot    = np.empty((n_rows, n_positions), dtype=object)
    pos_annot[:] = ""
    for row_idx, item in enumerate(enriched):
        aligned_pcts, aligned_annot = _align_position_stats_to_anchor(
            item["pos_stats"], item["n_positions"], n_positions
        )
        for col_idx, pct in enumerate(aligned_pcts):
            if pct is not None:
                pos_matrix[row_idx, col_idx] = pct
            pos_annot[row_idx, col_idx] = aligned_annot[col_idx]

    tier_matrix = np.array(
        [
            [
                item["n_exact_match"]      / n_variants_total * 100,
                item["n_identity_90"]      / n_variants_total * 100,
                item["n_identity_80"]      / n_variants_total * 100,
                item["n_passed_threshold"] / n_variants_total * 100,
            ]
            for item in enriched
        ]
    )
    tier_annot = np.array(
        [
            [
                f"{item['n_exact_match']}/{n_variants_total}",
                f"{item['n_identity_90']}/{n_variants_total}",
                f"{item['n_identity_80']}/{n_variants_total}",
                f"{item['n_passed_threshold']}/{n_variants_total}",
            ]
            for item in enriched
        ]
    )

    peptides = [item["peptide"]            for item in enriched]
    labels   = [item["conservation_label"] for item in enriched]

    fig_height = max(4.0, min(22.0, 0.45 * n_rows + 2.0))
    fig_width  = 14.0

    fig, axes = plt.subplots(
        1, 3,
        figsize=(fig_width, fig_height),
        gridspec_kw={"width_ratios": [0.04, n_positions * 0.55, 4 * 0.7]},
    )

    # ─── Left: label sidebar ───
    ax_bar = axes[0]
    for i, lbl in enumerate(labels):
        ax_bar.add_patch(
            mpatches.Rectangle((0, i), 1, 1, color=_LABEL_HEX.get(lbl, "#D9D9D9"), linewidth=0)
        )
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n_rows)
    ax_bar.set_xticks([])
    ax_bar.set_yticks(np.arange(n_rows) + 0.5)
    ax_bar.set_yticklabels(peptides, fontsize=8, fontfamily="monospace")
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Label", fontsize=8)
    ax_bar.tick_params(axis="y", length=0)

    # ─── Middle: position conservation ───
    ax_pos = axes[1]
    show_annot = n_rows <= 60
    sns.heatmap(
        pos_matrix,
        ax=ax_pos,
        cmap="RdYlGn",
        vmin=0, vmax=100,
        annot=pos_annot if show_annot else False,
        fmt="",
        annot_kws={"size": 7},
        linewidths=0.4, linecolor="#cccccc",
        cbar=False,
        xticklabels=_build_anchor_aligned_xtick_labels(n_positions),
        yticklabels=False,
        mask=np.isnan(pos_matrix),
    )
    ax_pos.set_facecolor("#f0f0f0")
    ax_pos.set_xlabel("")
    ax_pos.set_title("Position conservation", fontsize=10, fontweight="bold", pad=8)
    ax_pos.tick_params(axis="x", labelsize=9)

    for anchor_idx in (p2_idx, pomega_idx):
        ax_pos.add_patch(plt.Rectangle(
            (anchor_idx, 0), 1, n_rows,
            fill=False, edgecolor=_ANCHOR_BORDER_HEX, linewidth=2.5, clip_on=False,
        ))

    # ─── Right: identity tiers ───
    ax_tier = axes[2]
    tier_headers = ["≥100%", "≥90%", "≥80%", f"≥{threshold_pct}%"]
    sns.heatmap(
        tier_matrix,
        ax=ax_tier,
        cmap="RdYlGn",
        vmin=0, vmax=100,
        annot=tier_annot if show_annot else False,
        fmt="",
        annot_kws={"size": 7},
        linewidths=0.4, linecolor="#cccccc",
        cbar_kws={"shrink": 0.7, "label": "% of variants"},
        xticklabels=tier_headers,
        yticklabels=False,
    )
    ax_tier.set_xlabel("")
    ax_tier.set_title(f"Identity tiers (n={n_variants_total} variants)", fontsize=10, fontweight="bold", pad=8)
    ax_tier.tick_params(axis="x", labelsize=9)

    fig.suptitle(
        f"Conservation: {track_id} ({n_rows} epitopes shown, {n_variants_total} variants analysed)",
        fontsize=12, fontweight="bold", y=0.995,
    )

    legend_patches = [
        mpatches.Patch(color=_LABEL_HEX["perfect"],  label="★ Perfect"),
        mpatches.Patch(color=_LABEL_HEX["high"],     label="High"),
        mpatches.Patch(color=_LABEL_HEX["moderate"], label="Moderate"),
        mpatches.Patch(color=_LABEL_HEX["low"],      label="Low"),
        mpatches.Patch(facecolor="white", edgecolor=_ANCHOR_BORDER_HEX, linewidth=2, label="Anchor (P2, PΩ)"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=5,
        fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.005),
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


