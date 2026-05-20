"""
PNG generation for population_coverage: the IEDB-style per-population hit chart
and the comparative coverage heatmap. Matplotlib/seaborn are imported lazily
inside the writers so the rest of the step does not pay for them.
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from .core import _COVERAGE_HIGH_BAND, _COVERAGE_MODERATE_BAND, _MHC_CLASS

# IEDB-style hit chart palette (matches the user's prototype).
_HIT_CHART_HEADER_DARK  = "#1a3a5c"
_HIT_CHART_HEADER_MED   = "#2e6096"
_HIT_CHART_HEADER_LIGHT = "#3a7abf"
_HIT_CHART_PLUS_COLOR   = "#cc0000"
_HIT_CHART_MINUS_COLOR  = "#1a3a5c"
_HIT_CHART_ROW_EVEN     = "#f0f5fb"
_HIT_CHART_ROW_ODD      = "#ffffff"
_HIT_CHART_COV_HIGH     = "#d4edda"
_HIT_CHART_COV_MED      = "#fff3cd"
_HIT_CHART_COV_LOW      = "#f8d7da"
_HIT_CHART_HIT_BG       = "#e8f5e9"
_HIT_CHART_BORDER       = "#c0c0c0"
_HIT_CHART_TEXT_HEADER  = "white"


def write_hit_chart_png(
    population:        str,
    track_id:          str,
    epitope_results:   list[dict],
    hlas_with_freq:    list[tuple[str, float]],
    cutoff:            Optional[float],
    output_path:       Path,
):
    """English port of the user's prototype, matching IEDB layout."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    n_epitopes = len(epitope_results)
    n_hlas     = len(hlas_with_freq)
    if n_epitopes == 0 or n_hlas == 0:
        return

    effective_cutoff = cutoff if cutoff is not None else _COVERAGE_MODERATE_BAND

    col_epitope_width  = 1.6
    col_coverage_width = 0.9
    col_hla_width      = 0.72
    row_height         = 0.38
    header_height      = 0.9

    total_width  = col_epitope_width + col_coverage_width + n_hlas * col_hla_width + 0.5
    total_height = header_height + (n_epitopes + 1) * row_height + 1.0

    fig, ax = plt.subplots(figsize=(total_width, total_height))
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, total_height)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    def _draw_cell(
        x_pos: float, y_pos: float, w: float, h: float,
        text: str, background_color: str, text_color: str = "black",
        font_size: float = 7.5, bold: bool = False,
        h_align: str = "center", v_align: str = "center",
    ):
        ax.add_patch(mpatches.Rectangle(
            (x_pos, y_pos), w, h,
            facecolor=background_color, edgecolor=_HIT_CHART_BORDER,
            linewidth=0.5, zorder=1,
        ))
        ax.text(
            x_pos + w / 2, y_pos + h / 2, text,
            ha=h_align, va=v_align,
            fontsize=font_size, color=text_color,
            fontweight=("bold" if bold else "normal"),
            zorder=2, wrap=False, clip_on=True,
        )

    # Title bar
    title_height = 0.5
    y_title_base = total_height - title_height
    _draw_cell(
        0.1, y_title_base, total_width - 0.2, title_height - 0.05,
        f"{track_id}  —  Epitope Coverage  |  Population: {population}  |  Class: MHC {_MHC_CLASS}",
        _HIT_CHART_HEADER_DARK, _HIT_CHART_TEXT_HEADER, font_size=9, bold=True,
    )

    # HLA group sub-header
    y_hla_caption = y_title_base - 0.32
    x_hla_start   = 0.1 + col_epitope_width + col_coverage_width
    _draw_cell(
        x_hla_start, y_hla_caption,
        n_hlas * col_hla_width, 0.28,
        "HLA allele (genotypic frequency %)",
        _HIT_CHART_HEADER_MED, _HIT_CHART_TEXT_HEADER, font_size=8, bold=True,
    )

    # Column headers
    y_header = y_hla_caption - header_height
    x_cursor = 0.1
    _draw_cell(
        x_cursor, y_header, col_epitope_width, header_height,
        "Epitope", _HIT_CHART_HEADER_DARK, _HIT_CHART_TEXT_HEADER,
        font_size=8.5, bold=True,
    )
    x_cursor += col_epitope_width
    _draw_cell(
        x_cursor, y_header, col_coverage_width, header_height,
        f"Coverage\nClass {_MHC_CLASS}",
        _HIT_CHART_HEADER_DARK, _HIT_CHART_TEXT_HEADER, font_size=8, bold=True,
    )
    x_cursor += col_coverage_width
    for hla_index, (allele, frequency_pct) in enumerate(hlas_with_freq):
        short_label = allele.replace("HLA-", "")
        header_bg   = _HIT_CHART_HEADER_MED if hla_index % 2 == 0 else _HIT_CHART_HEADER_LIGHT
        _draw_cell(
            x_cursor, y_header, col_hla_width, header_height,
            f"{short_label}\n({frequency_pct}%)",
            header_bg, _HIT_CHART_TEXT_HEADER, font_size=6.5, bold=True,
        )
        x_cursor += col_hla_width

    # Epitope rows
    for row_index, result in enumerate(epitope_results):
        y_row    = y_header - (row_index + 1) * row_height
        x_cursor = 0.1
        row_bg   = _HIT_CHART_ROW_EVEN if row_index % 2 == 0 else _HIT_CHART_ROW_ODD

        _draw_cell(
            x_cursor, y_row, col_epitope_width, row_height,
            f"#{row_index + 1}: {result['peptide']}",
            row_bg, "black", font_size=7.5, h_align="center",
        )
        x_cursor += col_epitope_width

        coverage_pct = result["coverage_pct"]
        if coverage_pct >= _COVERAGE_HIGH_BAND:
            coverage_bg = _HIT_CHART_COV_HIGH
        elif coverage_pct >= effective_cutoff:
            coverage_bg = _HIT_CHART_COV_MED
        else:
            coverage_bg = _HIT_CHART_COV_LOW

        _draw_cell(
            x_cursor, y_row, col_coverage_width, row_height,
            f"{coverage_pct}%",
            coverage_bg, "black", font_size=8, bold=True,
        )
        x_cursor += col_coverage_width

        for allele, _ in hlas_with_freq:
            is_restricted = allele in result["matched_alleles"]
            symbol        = "+" if is_restricted else "-"
            text_color    = _HIT_CHART_PLUS_COLOR if is_restricted else _HIT_CHART_MINUS_COLOR
            cell_bg       = _HIT_CHART_HIT_BG if is_restricted else row_bg
            _draw_cell(
                x_cursor, y_row, col_hla_width, row_height,
                symbol, cell_bg, text_color, font_size=9, bold=is_restricted,
            )
            x_cursor += col_hla_width

    # Totals row — combined coverage = 1 - Π(1 - p_i) per epitope.
    y_totals = y_header - (n_epitopes + 1) * row_height
    x_cursor = 0.1
    _draw_cell(
        x_cursor, y_totals, col_epitope_width, row_height,
        "Epitope set", _HIT_CHART_HEADER_DARK, _HIT_CHART_TEXT_HEADER,
        font_size=7.5, bold=True,
    )
    x_cursor += col_epitope_width

    not_covered_product = 1.0
    for result in epitope_results:
        not_covered_product *= (1.0 - result["coverage_pct"] / 100.0)
    epitope_set_coverage = round((1.0 - not_covered_product) * 100.0, 2)

    _draw_cell(
        x_cursor, y_totals, col_coverage_width, row_height,
        f"{epitope_set_coverage}%",
        _HIT_CHART_HEADER_DARK, _HIT_CHART_TEXT_HEADER, font_size=7.5, bold=True,
    )
    x_cursor += col_coverage_width

    for allele, _ in hlas_with_freq:
        hit_count = sum(
            1 for result in epitope_results if allele in result["matched_alleles"]
        )
        _draw_cell(
            x_cursor, y_totals, col_hla_width, row_height,
            str(hit_count),
            _HIT_CHART_HEADER_MED, _HIT_CHART_TEXT_HEADER, font_size=7.5, bold=True,
        )
        x_cursor += col_hla_width

    # Legend
    y_legend_line = y_totals - 0.35
    ax.text(
        0.1, y_legend_line,
        "  +  restricted     -  not restricted     "
        f"■ ≥{_COVERAGE_HIGH_BAND:g}%     ■ ≥{effective_cutoff:g}%     ■ <{effective_cutoff:g}%",
        fontsize=7, color="#444444", va="top",
    )
    legend_entries = [
        (2.2, _HIT_CHART_COV_HIGH),
        (3.4, _HIT_CHART_COV_MED),
        (4.6, _HIT_CHART_COV_LOW),
    ]
    for legend_x, legend_color in legend_entries:
        ax.add_patch(mpatches.Rectangle(
            (legend_x, y_legend_line - 0.18), 0.18, 0.18,
            facecolor=legend_color, edgecolor="#888", linewidth=0.5,
        ))

    plt.tight_layout(pad=0.3)
    fig.savefig(str(output_path), dpi=180, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)


def write_coverage_matrix_png(
    track_id:           str,
    summary_rows:       list[dict],
    populations_order:  list[str],
    output_path:        Path,
):
    """
    Heatmap: rows = epitopes (sorted by mean coverage desc), columns = populations.
    Only generated when at least two populations are present.
    """
    if len(populations_order) < 2:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    summary_df = pd.DataFrame(summary_rows)
    pivot      = summary_df.pivot_table(
        index="peptide", columns="population", values="coverage_pct", aggfunc="first"
    ).reindex(columns=populations_order)

    pivot["__mean__"] = pivot[populations_order].mean(axis=1)
    pivot = pivot.sort_values("__mean__", ascending=False).drop(columns="__mean__")

    n_rows = len(pivot)
    if n_rows == 0:
        return

    fig_height = max(4.0, min(22.0, 0.32 * n_rows + 2.0))
    fig_width  = max(4.5, 1.4 * len(populations_order) + 2.0)
    fig, ax    = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        pivot.values,
        ax=ax,
        cmap="RdYlGn",
        vmin=0, vmax=100,
        annot=np.array([[f"{v:.1f}%" for v in row] for row in pivot.values]),
        fmt="",
        annot_kws={"size": 7},
        linewidths=0.4, linecolor="#cccccc",
        cbar_kws={"shrink": 0.7, "label": "Coverage (%)"},
        xticklabels=list(pivot.columns),
        yticklabels=list(pivot.index),
    )
    ax.set_title(f"Population coverage — {track_id}", fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=9)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
