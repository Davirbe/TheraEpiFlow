"""
analyze_conservation step.

Takes the ★ representatives from select_representatives and measures how
faithfully each epitope appears in the variant sequences produced by
search_variants. Conservation is quantified via sliding window identity:
for each variant sequence the best-matching window of len(peptide) residues
is found, and its per-position match ratio is recorded.

The analysis threshold (default 1.0 = exact match) is configurable per
project. It determines which variants are labelled "passed" vs "failed"
and which mutations are surfaced in position_mutation_profile. The
conservation_label and all row colours are based on mean_max_identity and
never change regardless of the chosen threshold.

The step is qualitative — no epitopes are removed.

Inputs (clusters/):
    CLUSTER_REPR_{track_id}.csv          select_representatives output

Inputs (variants/):
    VARIANTS_{track_id}.fasta            search_variants output or user-supplied

Outputs (conservation/):
    CONSERVATION_{track_id}.csv
    CONSERVATION_{track_id}.xlsx         row-coloured by conservation_label
    CONSERVATION_VISUAL_{track_id}.xlsx  position alignment matrix, one sheet per epitope
    CONSERVATION_HEATMAP_{track_id}.png  matplotlib/seaborn heatmap image
    CONSERVATION_AUDIT_{track_id}.json
"""

import datetime
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from rich import box
from rich.panel import Panel
from rich.table import Table

import config
from modules.base_step import BaseTrackStep
from utils.console import console, is_interactive_session
from utils.naming import (
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_PEPTIDE,
    get_step_filename,
)
from utils.project_manager import save_project_config

# ── Colour palette ────────────────────────────────────────────────────────────

_FILL_ORANGE    = PatternFill("solid", fgColor="FFD966")
_FILL_GRAY      = PatternFill("solid", fgColor="D9D9D9")
_FILL_PERFECT   = PatternFill("solid", fgColor="00B050")
_FILL_HIGH      = PatternFill("solid", fgColor="92D050")
_FILL_MODERATE  = PatternFill("solid", fgColor="FFFF99")
_FILL_LOW       = PatternFill("solid", fgColor="FF9999")
_FILL_UNKNOWN   = PatternFill("solid", fgColor="D9D9D9")

# Heatmap position colours (gradient by conservation %)
_FILL_POS_100    = PatternFill("solid", fgColor="00B050")   # 100%  — dark green
_FILL_POS_90     = PatternFill("solid", fgColor="92D050")   # ≥90%  — light green
_FILL_POS_70     = PatternFill("solid", fgColor="FFFF99")   # ≥70%  — yellow
_FILL_POS_50     = PatternFill("solid", fgColor="FFC000")   # ≥50%  — orange
_FILL_POS_LOW    = PatternFill("solid", fgColor="FF9999")   # <50%  — red
_FILL_ANCHOR_HDR = PatternFill("solid", fgColor="2E75B6")   # anchor column header — blue
_FILL_FOOTER     = PatternFill("solid", fgColor="BDD7EE")   # footer row — light blue

_LABEL_FILL = {
    "perfect":              _FILL_PERFECT,
    "high":                 _FILL_HIGH,
    "moderate":             _FILL_MODERATE,
    "low":                  _FILL_LOW,
    "conservation_unknown": _FILL_UNKNOWN,
}

_NUMERIC_COLS = frozenset({
    "n_variants_used",
    "pct_exact_match",
    "pct_identity_90",
    "pct_identity_80",
    "mean_max_identity",
})

# Accept variants whose length is within ±LENGTH_TOLERANCE of the reference.
# Sequences outside this window are scored via min(len) denominator inflation
# and produce artificially high or low identity scores.
_LENGTH_TOLERANCE = 0.25


# ── Core computation ──────────────────────────────────────────────────────────

def compute_max_window_identity(peptide: str, sequence: str) -> tuple[float, str]:
    """
    Slides a window of len(peptide) over sequence, returning the maximum
    per-position identity (0.0–1.0) and the corresponding window string.
    Returns (0.0, "") when the sequence is shorter than the peptide.
    """
    peptide_length = len(peptide)
    best_identity  = 0.0
    best_window    = ""
    for start in range(len(sequence) - peptide_length + 1):
        window   = sequence[start : start + peptide_length]
        identity = sum(ref == var for ref, var in zip(peptide, window)) / peptide_length
        if identity > best_identity:
            best_identity = identity
            best_window   = window
    return best_identity, best_window


def build_position_mutation_profile(
    peptide: str, all_windows: list[str], failed_windows: list[str]
) -> str:
    """
    Per-position conservation profile based on ALL variant windows.

    Format: "P1:47%(A→N) | P2:41%(I→E) | P3:100% | ..."
    - Conservation % = fraction of ALL variants that match the reference at that position.
    - Top mutation = most common substitution among variants that FAILED at this position.
    Returns "" when all_windows is empty.
    """
    if not all_windows:
        return ""

    n = len(all_windows)
    parts = []
    for pos_idx, ref_aa in enumerate(peptide):
        match_count = 0
        mutation_counter: Counter = Counter()
        for window in all_windows:
            if len(window) > pos_idx:
                var_aa = window[pos_idx]
                if var_aa == ref_aa:
                    match_count += 1
                else:
                    mutation_counter[f"{ref_aa}→{var_aa}"] += 1
        pct = round(match_count / n * 100)
        if mutation_counter:
            top_mut = mutation_counter.most_common(1)[0][0]
            parts.append(f"P{pos_idx + 1}:{pct}%({top_mut})")
        else:
            parts.append(f"P{pos_idx + 1}:{pct}%")
    return " | ".join(parts)


def classify_conservation_label(mean_max_identity: float, n_variants_total: int) -> str:
    if n_variants_total == 0:
        return "conservation_unknown"
    if mean_max_identity == 1.0:
        return "perfect"
    if mean_max_identity >= config.CONSERVATION_HIGH_THRESHOLD:
        return "high"
    if mean_max_identity >= config.CONSERVATION_MODERATE_THRESHOLD:
        return "moderate"
    return "low"


def build_conservation_summary(
    n_total: int,
    n_passed: int,
    analysis_threshold: float,
    n_identity_90pct: int,
    n_identity_80pct: int,
) -> str:
    """
    Builds a human-readable tier summary sorted from most to least stringent.
    The user's analysis_threshold is always included; 90% and 80% are added
    as fixed reference tiers when they differ from the chosen threshold.

    Examples:
        threshold=1.00  ->  "100%: 5/10 | 90%: 7/10 | 80%: 9/10"
        threshold=0.85  ->  "90%: 7/10 | 85%: 8/10 | 80%: 9/10"
        threshold=0.75  ->  "90%: 7/10 | 80%: 9/10 | 75%: 10/10"
        threshold=0.90  ->  "90%: 7/10 | 80%: 9/10"
    """
    tier_map: dict[float, int] = {
        0.90: n_identity_90pct,
        0.80: n_identity_80pct,
        analysis_threshold: n_passed,
    }
    sorted_tiers = sorted(tier_map.items(), reverse=True)
    parts = [
        f"{int(round(thr * 100))}%: {cnt}/{n_total}"
        for thr, cnt in sorted_tiers
    ]
    return " | ".join(parts)


def _variant_label(record) -> str:
    """Returns a short 'accession(organism)' label from a SeqRecord."""
    desc_parts = record.description.strip().split(" ", 1)
    organism_short = (
        desc_parts[1].split("|")[0].strip()[:30]
        if len(desc_parts) > 1
        else record.id
    )
    return f"{record.id}({organism_short})"


def compute_epitope_conservation(
    peptide: str,
    records: list,
    analysis_threshold: float,
) -> tuple[dict, list[tuple[str, float, str]]]:
    """
    Computes all conservation metrics for a single peptide against every
    variant record. Returns (metrics_dict, alignment_tuples).

    alignment_tuples: [(label, identity, best_window), ...] for every record,
    re-used by write_position_alignment_xlsx without re-computing identities.
    """
    n_total: int = len(records)
    alignment_tuples: list[tuple[str, float, str]] = []

    if n_total == 0:
        empty_metrics = {
            "n_variants_used":             0,
            "pct_passed_threshold":        0.0,
            "pct_exact_match":             0.0,
            "pct_identity_90":             0.0,
            "pct_identity_80":             0.0,
            "analysis_threshold":          analysis_threshold,
            "mean_max_identity":           0.0,
            "conservation_label":          "conservation_unknown",
            "position_mutation_profile":   "",
        }
        return empty_metrics, alignment_tuples

    per_variant_identities: list[float] = []
    passed_entries: list[str] = []
    failed_entries: list[str] = []
    failed_windows: list[str] = []

    for record in records:
        sequence = str(record.seq).upper()
        identity, window = compute_max_window_identity(peptide, sequence)
        label = _variant_label(record)

        per_variant_identities.append(identity)
        alignment_tuples.append((label, identity, window))

        if identity >= analysis_threshold:
            passed_entries.append(label)
        else:
            failed_entries.append(f"{label}:{window}")
            failed_windows.append(window)

    mean_max_identity = sum(per_variant_identities) / n_total
    n_exact  = sum(1 for i in per_variant_identities if i == 1.0)
    n_90pct  = sum(1 for i in per_variant_identities if i >= 0.90)
    n_80pct  = sum(1 for i in per_variant_identities if i >= 0.80)
    n_passed = sum(1 for i in per_variant_identities if i >= analysis_threshold)

    all_windows = [window for _, _, window in alignment_tuples]

    metrics = {
        "n_variants_used":             n_total,
        "pct_passed_threshold":        round(n_passed / n_total * 100, 1),
        "pct_exact_match":             round(n_exact  / n_total * 100, 1),
        "pct_identity_90":             round(n_90pct  / n_total * 100, 1),
        "pct_identity_80":             round(n_80pct  / n_total * 100, 1),
        "analysis_threshold":          analysis_threshold,
        "mean_max_identity":           round(mean_max_identity, 4),
        "conservation_label":          classify_conservation_label(mean_max_identity, n_total),
        "position_mutation_profile":   build_position_mutation_profile(
                                           peptide, all_windows, failed_windows
                                       ),
    }
    return metrics, alignment_tuples


# ── FASTA loading ─────────────────────────────────────────────────────────────

def load_fasta_sequences(fasta_path: Path, ref_length: int = 0) -> tuple[list, int]:
    """Returns (records, n_excluded). When ref_length > 0, filters to ±LENGTH_TOLERANCE."""
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if ref_length <= 0:
        return records, 0
    lo = max(1, int(ref_length * (1 - _LENGTH_TOLERANCE)))
    hi = int(ref_length * (1 + _LENGTH_TOLERANCE))
    kept = [r for r in records if lo <= len(r.seq) <= hi]
    return kept, len(records) - len(kept)


def _is_non_interactive() -> bool:
    return not sys.stdin.isatty()


# ── Threshold prompt ──────────────────────────────────────────────────────────

def prompt_analysis_threshold(project_name: str, project_config: dict) -> float:
    """
    Returns the analysis threshold for this run.
    Interactive: shows current value and lets the user change it.
    Non-interactive: uses saved value or default 1.0.
    Persists to project_config['conservation_threshold'] after any change.
    """
    saved = project_config.get("conservation_threshold")

    if _is_non_interactive():
        return float(saved) if saved is not None else 1.0

    if saved is not None:
        console.print(Panel(
            f"[bold]Conservation analysis threshold[/bold]\n\n"
            f"Current threshold: [cyan]{int(round(float(saved) * 100))}%[/cyan]\n\n"
            "[dim]Determines which variants are shown as passed/failed and which\n"
            "mutations are surfaced. Summary counts (100%, 90%, 80%) and row\n"
            "colours are fixed regardless of this value.[/dim]\n\n"
            "  [cyan][1][/cyan] Keep current value\n"
            "  [cyan][2][/cyan] Change threshold",
            box=box.ROUNDED, title="Setup: analyze_conservation", title_align="left",
        ))
        while True:
            try:
                choice = input("> ").strip()
            except EOFError:
                choice = "1"
            if choice in ("1", ""):
                return float(saved)
            if choice == "2":
                break
            console.print("[dim]Type 1 or 2.[/dim]")

    console.print(Panel(
        "[bold]Conservation analysis threshold[/bold]\n\n"
        "[dim]Variants with identity >= threshold appear as 'passed'.\n"
        "Failed variants show their actual best-matching window and mutations.[/dim]\n\n"
        "  [cyan][1][/cyan] 1.00 — exact match only (default)\n"
        "  [cyan][2][/cyan] 0.90 — 90% identity\n"
        "  [cyan][3][/cyan] 0.80 — 80% identity\n"
        "  [cyan][4][/cyan] Custom value (e.g. 0.75, 0.85, 0.65)",
        box=box.ROUNDED, title="Setup: analyze_conservation", title_align="left",
    ))

    preset_map = {"1": 1.0, "2": 0.90, "3": 0.80}
    while True:
        try:
            choice = input("> ").strip()
        except EOFError:
            choice = "1"
        if choice in preset_map:
            threshold = preset_map[choice]
            break
        if choice == "4":
            while True:
                try:
                    raw = input("Threshold (e.g. 0.75): ").strip().replace(",", ".")
                except EOFError:
                    raw = "1.0"
                try:
                    threshold = float(raw)
                    if 0.0 < threshold <= 1.0:
                        break
                    console.print("[red]Value must be between 0 (exclusive) and 1 (inclusive).[/red]")
                except ValueError:
                    console.print("[red]Invalid value. Use a dot as decimal separator (e.g. 0.75).[/red]")
            break
        console.print("[dim]Invalid option. Type 1, 2, 3 or 4.[/dim]")

    project_config["conservation_threshold"] = threshold
    save_project_config(project_name, project_config)
    console.print(f"[dim]Threshold set to {threshold:.2f}, saved to project_config.[/dim]")
    return threshold


# ── XLSX writers ──────────────────────────────────────────────────────────────

_XLSX_PRESENTATION_COLUMNS = [
    "peptide",
    "conservation_label",
    "mean_max_identity",
    "pct_passed_threshold",
    "pct_exact_match",
    "pct_identity_90",
    "pct_identity_80",
    "n_variants_used",
    "alleles_united",
    "num_alleles_united",
    "position_mutation_profile",
]

_XLSX_WIDE_COLUMNS = frozenset({"position_mutation_profile"})


def write_conservation_summary_xlsx(df: pd.DataFrame, output_path: Path):
    """
    Writes a user-facing XLSX with a curated subset of columns in logical order.
    The full-data CSV (all 38 columns) is kept separately for downstream use.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Conservation Summary"

    # Use only columns that exist in the dataframe
    columns = [col for col in _XLSX_PRESENTATION_COLUMNS if col in df.columns]

    for col_idx, col_name in enumerate(columns, start=1):
        cell           = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font      = Font(bold=True)
        cell.alignment = Alignment(horizontal="center")
        cell.fill      = _FILL_ORANGE if col_name in _NUMERIC_COLS else _FILL_GRAY

    for row_idx, (_, row) in enumerate(df.iterrows(), start=2):
        label    = row.get("conservation_label", "conservation_unknown")
        row_fill = _LABEL_FILL.get(label, _FILL_UNKNOWN)
        for col_idx, col_name in enumerate(columns, start=1):
            cell           = ws.cell(row=row_idx, column=col_idx, value=row[col_name])
            cell.fill      = row_fill
            cell.alignment = Alignment(wrap_text=(col_name in _XLSX_WIDE_COLUMNS))

    for col_idx, col_name in enumerate(columns, start=1):
        col_letter = ws.cell(row=1, column=col_idx).column_letter
        if col_name in _XLSX_WIDE_COLUMNS:
            ws.column_dimensions[col_letter].width = 50
        else:
            max_width = max(
                len(str(col_name)),
                *(len(str(v)) for v in df[col_name].tolist()),
            )
            ws.column_dimensions[col_letter].width = min(max_width + 2, 40)

    wb.save(str(output_path))


def _position_conservation_fill(pct: float) -> PatternFill:
    if pct >= 100.0:
        return _FILL_POS_100
    if pct >= 90.0:
        return _FILL_POS_90
    if pct >= 70.0:
        return _FILL_POS_70
    if pct >= 50.0:
        return _FILL_POS_50
    return _FILL_POS_LOW


def _compute_position_stats(
    peptide: str, alignment_tuples: list[tuple]
) -> list[dict]:
    """
    For each position in the peptide returns:
      conservation_pct  — % of variants with the same AA as reference (0–100)
      top_mutation      — most common substitution string e.g. "Q→R", or "" if 100%
    """
    n_variants = len(alignment_tuples)
    stats = []
    for pos_idx, ref_aa in enumerate(peptide):
        if n_variants == 0:
            stats.append({"conservation_pct": 0.0, "top_mutation": ""})
            continue
        n_match: int = 0
        mutation_counter: Counter = Counter()
        for _, _, window in alignment_tuples:
            if len(window) > pos_idx:
                var_aa = window[pos_idx]
                if var_aa == ref_aa:
                    n_match += 1
                else:
                    mutation_counter[f"{ref_aa}→{var_aa}"] += 1
        pct     = n_match / n_variants * 100
        top_mut = mutation_counter.most_common(1)[0][0] if mutation_counter else ""
        stats.append({"conservation_pct": pct, "top_mutation": top_mut})
    return stats


def write_position_heatmap_xlsx(
    alignment_data: list[dict],
    output_path: Path,
):
    """
    Writes the visual conservation heatmap XLSX — a single sheet, one row per
    epitope, one column per peptide position.

    Layout:
      Row 1        : headers — "Peptide" | P1 .. PN (anchors P2 and Pc in blue)
                               | "Anchor Score" | "Mean ID" | "Label"
      Rows 2..N+1  : one row per epitope, sorted by anchor_score descending
                     Each position cell shows:
                       Line 1 — conservation % (e.g. "67%")
                       Line 2 — top mutation if not 100% (e.g. "Q→R")
                     Cell background: gradient from green (100%) to red (<50%)
      Row N+2      : footer — average conservation per position across all epitopes

    Anchor positions for MHC-I: P2 (index 1) and Pc (last index).
    anchor_score = min(P2 conservation, Pc conservation).

    alignment_data entries:
        {
          "peptide":           str,
          "alignment_tuples":  [(label, identity, window), ...],
          "mean_max_identity": float,
          "conservation_label": str,
        }
    """
    if not alignment_data:
        return

    n_positions  = len(alignment_data[0]["peptide"])
    anchor_p2    = 1                   # 0-based index of P2
    anchor_pc    = n_positions - 1     # 0-based index of last position

    # ── Compute per-epitope stats ─────────────────────────────────────────────
    rows: list[dict] = []
    for item in alignment_data:
        peptide    = item["peptide"]
        pos_stats  = _compute_position_stats(peptide, item["alignment_tuples"])
        p2_pct     = pos_stats[anchor_p2]["conservation_pct"] if len(pos_stats) > anchor_p2 else 0.0
        pc_pct     = pos_stats[anchor_pc]["conservation_pct"] if len(pos_stats) > anchor_pc else 0.0
        anchor_score = min(p2_pct, pc_pct)
        rows.append({
            "peptide":           peptide,
            "pos_stats":         pos_stats,
            "anchor_score":      anchor_score,
            "mean_max_identity": item["mean_max_identity"],
            "conservation_label": item["conservation_label"],
        })

    # Sort by anchor_score descending (best vaccine candidates on top)
    rows.sort(key=lambda r: r["anchor_score"], reverse=True)

    # ── Build workbook ────────────────────────────────────────────────────────
    wb = Workbook()
    ws = wb.active
    ws.title = "Conservation Heatmap"

    position_headers = [f"P{i + 1}" for i in range(n_positions)]
    all_headers      = ["Peptide"] + position_headers + ["Anchor Score", "Mean ID", "Label"]
    n_cols           = len(all_headers)

    # Header row
    for col_idx, header in enumerate(all_headers, start=1):
        cell           = ws.cell(row=1, column=col_idx, value=header)
        cell.font      = Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(horizontal="center", vertical="center")

        # Anchor columns (P2 and Pc) get blue header; others get gray
        pos_0based = col_idx - 2  # col 2 = P1 = pos_0based 0
        if col_idx >= 2 and col_idx <= n_positions + 1:
            if pos_0based == anchor_p2 or pos_0based == anchor_pc:
                cell.fill = _FILL_ANCHOR_HDR
            else:
                cell.fill = _FILL_GRAY
        else:
            cell.fill = _FILL_GRAY

    ws.row_dimensions[1].height = 20

    # Data rows
    for row_idx, row_data in enumerate(rows, start=2):
        peptide        = row_data["peptide"]
        pos_stats      = row_data["pos_stats"]
        anchor_score   = row_data["anchor_score"]
        mean_id        = row_data["mean_max_identity"]
        label          = row_data["conservation_label"]

        # Peptide name cell
        name_cell           = ws.cell(row=row_idx, column=1, value=peptide)
        name_cell.font      = Font(bold=True)
        name_cell.alignment = Alignment(horizontal="left", vertical="center")
        name_cell.fill      = _LABEL_FILL.get(label, _FILL_UNKNOWN)

        # Position cells
        for pos_idx, stat in enumerate(pos_stats):
            col_idx   = pos_idx + 2
            pct       = stat["conservation_pct"]
            top_mut   = stat["top_mutation"]

            cell_text = f"{int(round(pct))}%"
            if top_mut:
                cell_text += f"\n{top_mut}"

            cell           = ws.cell(row=row_idx, column=col_idx, value=cell_text)
            cell.fill      = _position_conservation_fill(pct)
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            cell.font      = Font(bold=(pct == 100.0))

        # Anchor score, mean ID, label
        anchor_col = n_positions + 2
        mean_col   = n_positions + 3
        label_col  = n_positions + 4

        anchor_cell           = ws.cell(row=row_idx, column=anchor_col, value=f"{int(round(anchor_score))}%")
        anchor_cell.fill      = _position_conservation_fill(anchor_score)
        anchor_cell.alignment = Alignment(horizontal="center", vertical="center")
        anchor_cell.font      = Font(bold=True)

        mean_cell           = ws.cell(row=row_idx, column=mean_col, value=f"{mean_id:.3f}")
        mean_cell.alignment = Alignment(horizontal="center", vertical="center")
        mean_cell.fill      = _LABEL_FILL.get(label, _FILL_UNKNOWN)

        label_cell           = ws.cell(row=row_idx, column=label_col, value=label)
        label_cell.alignment = Alignment(horizontal="center", vertical="center")
        label_cell.fill      = _LABEL_FILL.get(label, _FILL_UNKNOWN)

        ws.row_dimensions[row_idx].height = 30  # space for 2-line cells

    # Footer: average conservation per position across all epitopes
    footer_row = len(rows) + 2
    footer_label_cell           = ws.cell(row=footer_row, column=1, value="Avg conservation")
    footer_label_cell.fill      = _FILL_FOOTER
    footer_label_cell.font      = Font(bold=True)
    footer_label_cell.alignment = Alignment(horizontal="left", vertical="center")

    for pos_idx in range(n_positions):
        col_idx  = pos_idx + 2
        avg_pct  = (
            sum(r["pos_stats"][pos_idx]["conservation_pct"] for r in rows) / len(rows)
            if rows else 0.0
        )
        cell           = ws.cell(row=footer_row, column=col_idx, value=f"{int(round(avg_pct))}%")
        cell.fill      = _FILL_FOOTER
        cell.font      = Font(bold=True)
        cell.alignment = Alignment(horizontal="center", vertical="center")

    for col_idx in range(n_positions + 2, n_cols + 1):
        cell      = ws.cell(row=footer_row, column=col_idx, value="")
        cell.fill = _FILL_FOOTER

    ws.row_dimensions[footer_row].height = 20

    # Column widths
    ws.column_dimensions["A"].width = 14
    for pos_idx in range(n_positions):
        col_letter = ws.cell(row=1, column=pos_idx + 2).column_letter
        ws.column_dimensions[col_letter].width = 9
    ws.column_dimensions[ws.cell(row=1, column=n_positions + 2).column_letter].width = 14
    ws.column_dimensions[ws.cell(row=1, column=n_positions + 3).column_letter].width = 10
    ws.column_dimensions[ws.cell(row=1, column=n_positions + 4).column_letter].width = 20

    wb.save(str(output_path))


# ── Position conservation heatmap PNG ────────────────────────────────────────

def write_position_conservation_png(
    alignment_data: list[dict],
    track_id: str,
    output_path: Path,
):
    """
    PNG heatmap: rows = epitopes, columns = peptide positions.
    Each cell = % of variants that conserve the reference AA at that position.
    Color: green (100%) → red (0%). Cells annotated with conservation %.

    Separate from write_conservation_heatmap_png (which shows epitope-level
    summary metrics). This image focuses on WHERE mutations occur per position.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    import numpy as np

    if not alignment_data:
        return

    n_positions = len(alignment_data[0]["peptide"])
    peptides    = [d["peptide"] for d in alignment_data]
    labels      = [d["conservation_label"] for d in alignment_data]

    # Build matrix: rows=epitopes, cols=positions, value=conservation% (0–100)
    matrix = []
    for item in alignment_data:
        pos_stats = _compute_position_stats(item["peptide"], item["alignment_tuples"])
        matrix.append([s["conservation_pct"] for s in pos_stats])

    matrix_df = pd.DataFrame(
        matrix,
        index=peptides,
        columns=[f"P{i+1}" for i in range(n_positions)],
    )

    # Sort by mean conservation descending
    matrix_df = matrix_df.loc[
        matrix_df.mean(axis=1).sort_values(ascending=False).index
    ]
    sorted_labels = [labels[peptides.index(p)] for p in matrix_df.index]

    n = len(matrix_df)
    row_height  = max(0.4, min(0.7, 12 / n))
    fig_height  = max(4, n * row_height + 2.5)
    fig, axes   = plt.subplots(
        1, 2, figsize=(max(8, n_positions * 1.2 + 2.5), fig_height),
        gridspec_kw={"width_ratios": [0.06, 1]},
    )

    # Left colour bar — label
    ax_bar = axes[0]
    for i, lbl in enumerate(sorted_labels):
        ax_bar.add_patch(mpatches.Rectangle(
            (0, i), 1, 1, color=_LABEL_HEX.get(lbl, "#D9D9D9"), linewidth=0
        ))
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n)
    ax_bar.set_xticks([])
    ax_bar.set_yticks(np.arange(n) + 0.5)
    ax_bar.set_yticklabels(matrix_df.index, fontsize=8, fontfamily="monospace")
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Label", fontsize=8)
    ax_bar.tick_params(axis="y", length=0)

    # Main heatmap — position conservation %
    ax_hm = axes[1]
    show_annot = n <= 40
    annot_data = matrix_df.map(lambda v: f"{int(round(v))}%") if show_annot else False
    sns.heatmap(
        matrix_df,
        ax=ax_hm,
        cmap="RdYlGn",
        vmin=0, vmax=100,
        annot=annot_data, fmt="",
        annot_kws={"size": 8},
        linewidths=0.4, linecolor="#cccccc",
        cbar_kws={"shrink": 0.6, "label": "Conservation %"},
        yticklabels=False,
    )
    ax_hm.set_xlabel("")
    ax_hm.set_ylabel("")
    ax_hm.set_title(
        f"Position conservation — {track_id}",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax_hm.tick_params(axis="x", labelsize=10)

    legend_patches = [
        mpatches.Patch(color=_LABEL_HEX["perfect"],  label="★ Perfect"),
        mpatches.Patch(color=_LABEL_HEX["high"],     label="High"),
        mpatches.Patch(color=_LABEL_HEX["moderate"], label="Moderate"),
        mpatches.Patch(color=_LABEL_HEX["low"],      label="Low"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=4,
        fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.005),
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Rich console output ───────────────────────────────────────────────────────

def print_conservation_rich_table(
    df: pd.DataFrame, track_id: str, analysis_threshold: float
):
    thr_pct = int(round(analysis_threshold * 100))
    t = Table(
        box=box.SIMPLE,
        title=f"Conservation analysis — {track_id}  (threshold={thr_pct}%)",
        show_lines=False,
    )
    t.add_column("Peptide",             style="bold", no_wrap=True)
    t.add_column("Variants",            justify="right")
    t.add_column(f"≥{thr_pct}%",        justify="right")
    t.add_column("100%",                justify="right")
    t.add_column("≥90%",                justify="right")
    t.add_column("≥80%",                justify="right")
    t.add_column("Mean ID",             justify="right")
    t.add_column("Label",               no_wrap=True)

    label_rich = {
        "perfect":              "[bold green]★ perfect[/bold green]",
        "high":                 "[green]  high[/green]",
        "moderate":             "[yellow]  moderate[/yellow]",
        "low":                  "[red]  low[/red]",
        "conservation_unknown": "[dim]  unknown[/dim]",
    }

    for _, row in df.iterrows():
        label = row.get("conservation_label", "conservation_unknown")
        n     = int(row.get("n_variants_used", 0))
        t.add_row(
            str(row.get(COLUMN_PEPTIDE, "")),
            str(n),
            f"{row.get('pct_passed_threshold', 0.0):.1f}%",
            f"{row.get('pct_exact_match', 0.0):.1f}%",
            f"{row.get('pct_identity_90', 0.0):.1f}%",
            f"{row.get('pct_identity_80', 0.0):.1f}%",
            f"{float(row.get('mean_max_identity', 0.0)):.3f}",
            label_rich.get(label, label),
        )

    console.print(t)


# ── Heatmap PNG ───────────────────────────────────────────────────────────────

_LABEL_HEX = {
    "perfect":              "#00B050",
    "high":                 "#92D050",
    "moderate":             "#FFFF99",
    "low":                  "#FF9999",
    "conservation_unknown": "#D9D9D9",
}


def write_conservation_heatmap_png(df: pd.DataFrame, track_id: str, output_path: Path):
    """Saves a seaborn heatmap PNG: peptides × 4 identity metrics, with a label colour bar."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    import numpy as np

    if df.empty:
        return

    df_plot = df.sort_values("mean_max_identity", ascending=False).reset_index(drop=True)
    n = len(df_plot)

    heatmap_data = pd.DataFrame({
        "Mean identity":  df_plot["mean_max_identity"].values,
        "Exact (100%)":   (df_plot["pct_exact_match"] / 100).values,
        "≥ 90% identity": (df_plot["pct_identity_90"]  / 100).values,
        "≥ 80% identity": (df_plot["pct_identity_80"]  / 100).values,
    }, index=df_plot["peptide"].values)

    label_colors = [
        _LABEL_HEX.get(lbl, "#D9D9D9")
        for lbl in df_plot["conservation_label"]
    ]

    row_height  = max(0.35, min(0.55, 14 / n))
    fig_height  = max(4, n * row_height + 2.5)
    fig, axes   = plt.subplots(
        1, 2,
        figsize=(10, fig_height),
        gridspec_kw={"width_ratios": [0.06, 1]},
    )

    # Left colour bar — conservation label
    ax_bar = axes[0]
    label_array = np.arange(n).reshape(-1, 1)
    for i, color in enumerate(label_colors):
        ax_bar.add_patch(mpatches.Rectangle((0, i), 1, 1, color=color, linewidth=0))
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n)
    ax_bar.set_xticks([])
    ax_bar.set_yticks(np.arange(n) + 0.5)
    ax_bar.set_yticklabels(df_plot["peptide"], fontsize=8, fontfamily="monospace")
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Label", fontsize=8)
    ax_bar.tick_params(axis="y", length=0)

    # Main heatmap
    ax_hm = axes[1]
    vmin, vmax = 0.0, 1.0
    show_annot = n <= 50
    annot_data = heatmap_data.map(lambda v: f"{v:.2f}") if show_annot else False
    sns.heatmap(
        heatmap_data,
        ax=ax_hm,
        cmap="RdYlGn",
        vmin=vmin, vmax=vmax,
        annot=annot_data, fmt="",
        annot_kws={"size": 7},
        linewidths=0.3, linecolor="#cccccc",
        cbar_kws={"shrink": 0.6, "label": "Fraction"},
        yticklabels=False,
    )
    ax_hm.set_xlabel("")
    ax_hm.set_ylabel("")
    ax_hm.set_title(f"Conservation — {track_id}", fontsize=11, fontweight="bold", pad=10)
    ax_hm.tick_params(axis="x", labelsize=9)

    # Legend
    legend_patches = [
        mpatches.Patch(color=_LABEL_HEX["perfect"],  label="★ Perfect (100%)"),
        mpatches.Patch(color=_LABEL_HEX["high"],     label="High (≥90%)"),
        mpatches.Patch(color=_LABEL_HEX["moderate"], label="Moderate (≥70%)"),
        mpatches.Patch(color=_LABEL_HEX["low"],      label="Low (<70%)"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=4,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Step ──────────────────────────────────────────────────────────────────────

def _inspect_track_fasta_status(
    project_name: str,
    track_id: str,
) -> dict:
    """
    Return a status dict describing the variant FASTA situation for one track.

    Used by `preflight()` to render the consolidated status table before the
    user is asked whether they want to provide local FASTA files. The dict
    has the shape:

        {
            "track_id":      str,
            "ref_length":    int | None,
            "fasta_path":    Path | None,
            "fasta_exists":  bool,
            "n_records_raw": int,
            "n_records_kept": int,
            "status":        str,   # "ok" | "missing_fasta" | "no_reference" | "length_mismatch"
        }
    """
    project_root = Path("projects") / project_name
    track_dir = project_root / "data" / "intermediate" / track_id
    input_track_dir = project_root / "data" / "input" / track_id

    reference_fasta_path = input_track_dir / get_step_filename("SEQUENCES", track_id, ext="fasta")
    reference_length: Optional[int] = None
    if reference_fasta_path.exists():
        reference_records = list(SeqIO.parse(str(reference_fasta_path), "fasta"))
        if reference_records:
            reference_length = len(reference_records[0].seq)

    variants_fasta_path = track_dir / "variants" / get_step_filename("VARIANTS", track_id, ext="fasta")
    fasta_exists = variants_fasta_path.exists() and variants_fasta_path.stat().st_size > 0

    n_records_raw = 0
    n_records_kept = 0
    if fasta_exists:
        kept_records, excluded_count = load_fasta_sequences(variants_fasta_path, reference_length or 0)
        n_records_kept = len(kept_records)
        n_records_raw = n_records_kept + excluded_count

    if not reference_length:
        status_label = "no_reference"
    elif not fasta_exists:
        status_label = "missing_fasta"
    elif n_records_kept == 0 and n_records_raw > 0:
        status_label = "length_mismatch"
    else:
        status_label = "ok"

    return {
        "track_id":       track_id,
        "ref_length":     reference_length,
        "fasta_path":     variants_fasta_path if fasta_exists else None,
        "fasta_exists":   fasta_exists,
        "n_records_raw":  n_records_raw,
        "n_records_kept": n_records_kept,
        "status":         status_label,
    }


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
            str(status_row["ref_length"]) if status_row["ref_length"] is not None else "—",
            status_row["fasta_path"].name if status_row["fasta_path"] else "—",
            str(status_row["n_records_raw"]) if status_row["fasta_exists"] else "—",
            str(status_row["n_records_kept"]) if status_row["fasta_exists"] else "—",
            status_color_map.get(status_row["status"], status_row["status"]),
        )

    return status_table


def _ask_for_local_fasta_overrides(track_status_rows: list[dict]) -> dict:
    """
    Loop interactively letting the user attach a local FASTA path to one or
    more track ids. Returns {track_id: path_string}. Empty dict if the user
    declines or supplies nothing usable.
    """
    fasta_overrides: dict[str, str] = {}
    track_id_set = {row["track_id"] for row in track_status_rows}

    while True:
        try:
            track_input = input("  Track id to attach FASTA to (or blank to finish): ").strip()
        except EOFError:
            track_input = ""

        if not track_input:
            break
        if track_input not in track_id_set:
            console.print(f"  [yellow]Unknown track id '{track_input}'. Pick one from the table above.[/yellow]")
            continue

        try:
            path_input = input(f"  Local FASTA path for {track_input}: ").strip()
        except EOFError:
            path_input = ""

        candidate_path = Path(path_input).expanduser() if path_input else None
        if candidate_path is None or not candidate_path.exists() or candidate_path.stat().st_size == 0:
            console.print("  [red]File not found or empty. Skipping.[/red]")
            continue

        fasta_overrides[track_input] = str(candidate_path)
        console.print(f"  [green]✓ Will use {candidate_path} for {track_input}.[/green]")

    return fasta_overrides


def _identify_problematic_tracks(
    project_name: str,
    track_outcomes: dict,
    minimum_acceptable_variants: int = 5,
) -> list[dict]:
    """
    Return tracks whose conservation analysis did not produce a meaningful
    result. A track is problematic when (a) it errored, (b) it ran but had
    fewer than `minimum_acceptable_variants` usable sequences, or (c) it ran
    with no variants at all (conservation_unknown).
    """
    problematic_track_rows: list[dict] = []

    for track_id, outcome_dict in track_outcomes.items():
        outcome_status = outcome_dict.get("status")

        if outcome_status == "error":
            problematic_track_rows.append({
                "track_id": track_id,
                "problem":  outcome_dict.get("error_message", "unknown error"),
                "n_used":   None,
            })
            continue

        # Read the audit JSON if it exists to learn how many variants were used.
        audit_path = (
            Path("projects") / project_name / "data" / "intermediate" / track_id /
            "conservation" / get_step_filename("CONSERVATION_AUDIT", track_id, ext="json")
        )
        if not audit_path.exists():
            continue

        try:
            audit_payload = json.loads(audit_path.read_text(encoding="utf-8"))
        except Exception:
            continue

        n_variants_used = int(audit_payload.get("n_variants_used", 0))
        if n_variants_used < minimum_acceptable_variants:
            problematic_track_rows.append({
                "track_id": track_id,
                "problem":  f"only {n_variants_used} variants used"
                            f" (label_counts={audit_payload.get('label_counts', {})})",
                "n_used":   n_variants_used,
            })

    return problematic_track_rows


class AnalyzeConservationStep(BaseTrackStep):
    step_name = "analyze_conservation"

    @classmethod
    def preflight(cls, project_name: str, project_config: dict, track_ids: list[str]) -> Optional[dict]:
        """
        Show a single consolidated table of variant FASTA status across all
        tracks, then ask the user once whether they want to attach local
        FASTA files to any of them. Returns {track_id: path_string} for
        every override the user provided, or None if they declined.
        """
        if not track_ids:
            return None

        track_status_rows = [
            _inspect_track_fasta_status(project_name, track_id)
            for track_id in track_ids
        ]

        console.print(Panel(
            _render_track_status_table(track_status_rows),
            title="[bold]Variant FASTA status across all tracks[/bold]",
            border_style="cyan", box=box.ROUNDED,
        ))

        if not is_interactive_session():
            return None

        try:
            response = input(
                "  Do you have a local FASTA you want to attach to any track? [y/N]: "
            ).strip().lower()
        except EOFError:
            response = ""

        if response not in {"y", "yes"}:
            return None

        console.print(
            "[dim]  Type the track id (must match the table above), then the FASTA path. "
            "Leave the track id blank to stop.[/dim]"
        )
        fasta_overrides = _ask_for_local_fasta_overrides(track_status_rows)
        return fasta_overrides if fasta_overrides else None

    @classmethod
    def postflight(cls, project_name: str, project_config: dict, track_outcomes: dict) -> None:
        """
        After all tracks have been analysed, list the ones that ended up
        without a meaningful result and offer to re-run them with a
        user-provided FASTA.
        """
        if not track_outcomes or not is_interactive_session():
            return

        problematic_track_rows = _identify_problematic_tracks(project_name, track_outcomes)
        if not problematic_track_rows:
            return

        recovery_table = Table(box=box.SIMPLE, header_style="bold yellow")
        recovery_table.add_column("Track", style="cyan")
        recovery_table.add_column("Variants used", justify="right")
        recovery_table.add_column("Problem")
        for problem_row in problematic_track_rows:
            recovery_table.add_row(
                problem_row["track_id"],
                "—" if problem_row["n_used"] is None else str(problem_row["n_used"]),
                problem_row["problem"],
            )

        console.print(Panel(
            recovery_table,
            title="[bold]Tracks with weak conservation results[/bold]",
            border_style="yellow", box=box.ROUNDED,
        ))

        try:
            response = input(
                "  Re-run any of these with a local FASTA? [y/N]: "
            ).strip().lower()
        except EOFError:
            response = ""
        if response not in {"y", "yes"}:
            return

        from utils.pipeline_state import reset_track_step

        problematic_track_id_set = {row["track_id"] for row in problematic_track_rows}
        while True:
            try:
                target_track = input("  Track id to re-run (or blank to finish): ").strip()
            except EOFError:
                target_track = ""
            if not target_track:
                break
            if target_track not in problematic_track_id_set:
                console.print(f"  [yellow]Not in the list above: '{target_track}'.[/yellow]")
                continue

            try:
                fasta_path_input = input(f"  Local FASTA path for {target_track}: ").strip()
            except EOFError:
                fasta_path_input = ""
            override_path = Path(fasta_path_input).expanduser() if fasta_path_input else None
            if override_path is None or not override_path.exists() or override_path.stat().st_size == 0:
                console.print("  [red]File not found or empty. Skipping.[/red]")
                continue

            reset_track_step(project_name, target_track, cls.step_name)
            rerun_step_instance = cls(
                project_name=project_name,
                project_config=project_config,
                track_id=target_track,
            )
            rerun_step_instance.preflight_config = {target_track: str(override_path)}
            rerun_step_instance.execute(force_rerun=False)

    def describe_outputs(self) -> dict[Path, str]:
        conservation_dir = self.track_dir / "conservation"
        return {
            conservation_dir / get_step_filename("CONSERVATION", self.track_id):
                "Per-epitope conservation metrics (mean identity, % at thresholds, position profile).",
            conservation_dir / get_step_filename("CONSERVATION", self.track_id, ext="xlsx"):
                "Same metrics with row colouring by conservation label.",
            conservation_dir / get_step_filename("CONSERVATION_VISUAL", self.track_id, ext="xlsx"):
                "Position-by-position alignment matrix, one sheet per epitope.",
            conservation_dir / get_step_filename("CONSERVATION_HEATMAP", self.track_id, ext="png"):
                "Heatmap of conservation metrics across all representatives.",
            conservation_dir / get_step_filename("CONSERVATION_POSITIONS", self.track_id, ext="png"):
                "Heatmap showing per-position conservation across all representatives.",
            conservation_dir / get_step_filename("CONSERVATION_AUDIT", self.track_id, ext="json"):
                "Run audit — threshold, FASTA source, variants used, label counts.",
        }

    def run(self, input_data=None):
        analysis_threshold = prompt_analysis_threshold(
            self.project_name, self.project_config
        )

        # ── Load ★ representatives ────────────────────────────────────────────
        clusters_dir        = self.track_dir / "clusters"
        representatives_csv = clusters_dir / get_step_filename("CLUSTER_REPR", self.track_id)
        if not representatives_csv.exists():
            raise FileNotFoundError(
                f"select_representatives output not found: {representatives_csv}\n"
                "Run 'select_representatives' before 'analyze_conservation'."
            )

        df_repr  = pd.read_csv(representatives_csv)
        df_stars = df_repr[df_repr[COLUMN_BEST_REPRESENTATIVE] == "★"].copy()
        if df_stars.empty:
            raise ValueError(
                f"No ★ representatives found in {representatives_csv.name}. "
                "Run 'select_representatives' first."
            )

        peptides = df_stars[COLUMN_PEPTIDE].tolist()

        # ── Reference length (for variant length filtering) ───────────────────
        input_dir  = self.track_dir.parent.parent / "input" / self.track_id
        ref_fasta  = input_dir / get_step_filename("SEQUENCES", self.track_id, ext="fasta")
        ref_length = 0
        if ref_fasta.exists():
            ref_recs = list(SeqIO.parse(str(ref_fasta), "fasta"))
            if ref_recs:
                ref_length = len(ref_recs[0].seq)
                console.print(
                    f"[dim]→ Reference length: {ref_length} aa "
                    f"(±{int(_LENGTH_TOLERANCE*100)}% filter applied)[/dim]"
                )

        # ── Resolve FASTA source (silent — no per-track prompts) ──────────────
        # Order of resolution:
        #   1. Override path passed in via preflight (user provided one upfront)
        #   2. Variants cache produced by search_variants
        #   3. None — analysis proceeds and every epitope is marked
        #      'conservation_unknown'. The user can supply a local FASTA in
        #      postflight() and re-run only the affected track(s).
        variants_dir = self.track_dir / "variants"
        fasta_path   = variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta")
        fasta_source = "none"
        records: list = []

        override_fasta_path: Optional[Path] = None
        if isinstance(self.preflight_config, dict):
            override_value = self.preflight_config.get(self.track_id)
            if override_value:
                override_fasta_path = Path(override_value).expanduser()

        if override_fasta_path is not None and override_fasta_path.exists() and override_fasta_path.stat().st_size > 0:
            records, n_excl = load_fasta_sequences(override_fasta_path, ref_length)
            fasta_source = "user_provided"
            msg = f"[dim]→ Using user-provided FASTA {override_fasta_path.name} ({len(records)} sequences"
            if n_excl:
                msg += f", {n_excl} excluded by length filter"
            console.print(msg + ")[/dim]")
        elif fasta_path.exists() and fasta_path.stat().st_size > 0:
            records, n_excl = load_fasta_sequences(fasta_path, ref_length)
            fasta_source    = "variants_cache"
            msg = f"[dim]→ Using {fasta_path.name} ({len(records)} sequences"
            if n_excl:
                msg += f", {n_excl} excluded by length filter"
            console.print(msg + ")[/dim]")
        else:
            console.print(
                f"[yellow]⚠ No variant FASTA available for {self.track_id} — "
                "epitopes will be labelled 'conservation_unknown'. Provide a "
                "local FASTA after the loop completes if needed.[/yellow]"
            )

        # ── Analyse each peptide ──────────────────────────────────────────────
        metrics_rows:   list[dict] = []
        alignment_data: list[dict] = []

        for peptide in peptides:
            metrics, alignment_tuples = compute_epitope_conservation(
                peptide, records, analysis_threshold
            )
            metrics_rows.append(metrics)
            alignment_data.append({
                "peptide":            peptide,
                "alignment_tuples":   alignment_tuples,
                "mean_max_identity":  metrics["mean_max_identity"],
                "conservation_label": metrics["conservation_label"],
            })

        # ── Build result DataFrame ────────────────────────────────────────────
        result_df = pd.concat(
            [df_stars.reset_index(drop=True), pd.DataFrame(metrics_rows)],
            axis=1,
        )

        # ── Save outputs ──────────────────────────────────────────────────────
        conservation_dir = self.track_dir / "conservation"
        conservation_dir.mkdir(parents=True, exist_ok=True)

        output_csv          = conservation_dir / get_step_filename("CONSERVATION", self.track_id)
        output_xlsx         = conservation_dir / get_step_filename("CONSERVATION", self.track_id, ext="xlsx")
        output_visual_xlsx  = conservation_dir / get_step_filename("CONSERVATION_VISUAL", self.track_id, ext="xlsx")
        output_heatmap_png  = conservation_dir / get_step_filename("CONSERVATION_HEATMAP", self.track_id, ext="png")
        output_pos_png      = conservation_dir / get_step_filename("CONSERVATION_POSITIONS", self.track_id, ext="png")
        audit_path          = conservation_dir / get_step_filename("CONSERVATION_AUDIT", self.track_id, ext="json")

        import csv as _csv
        _csv_cols = [c for c in _XLSX_PRESENTATION_COLUMNS if c in result_df.columns]
        result_df[_csv_cols].to_csv(output_csv, index=False, quoting=_csv.QUOTE_NONNUMERIC)
        write_conservation_summary_xlsx(result_df, output_xlsx)
        write_position_heatmap_xlsx(alignment_data, output_visual_xlsx)
        write_conservation_heatmap_png(result_df, self.track_id, output_heatmap_png)
        write_position_conservation_png(alignment_data, self.track_id, output_pos_png)

        # ── Console summary ───────────────────────────────────────────────────
        print_conservation_rich_table(result_df, self.track_id, analysis_threshold)

        n_analyzed   = len(result_df)
        label_counts = result_df["conservation_label"].value_counts().to_dict()
        mean_all     = round(float(result_df["mean_max_identity"].mean()), 4) if n_analyzed else 0.0

        console.print(
            f"\n[bold green]RESULT: {n_analyzed} representatives analysed "
            f"({fasta_source}, threshold={int(round(analysis_threshold * 100))}%).[/bold green]\n"
        )

        # ── Audit JSON ────────────────────────────────────────────────────────
        audit = {
            "timestamp":                         datetime.datetime.now().isoformat(),
            "track_id":                          self.track_id,
            "analysis_threshold":                analysis_threshold,
            "fasta_source":                      fasta_source,
            "n_variants_used":                   len(records),
            "ref_length_aa":                     ref_length,
            "length_filter_tolerance":           _LENGTH_TOLERANCE,
            "n_representatives_analyzed":        n_analyzed,
            "label_counts":                      label_counts,
            "mean_identity_across_all_epitopes": mean_all,
            "reference_thresholds": {
                "high":     config.CONSERVATION_HIGH_THRESHOLD,
                "moderate": config.CONSERVATION_MODERATE_THRESHOLD,
            },
            "outputs": {
                "csv":          str(output_csv),
                "summary_xlsx": str(output_xlsx),
                "visual_xlsx":  str(output_visual_xlsx),
                "heatmap_png":      str(output_heatmap_png),
                "positions_png":    str(output_pos_png),
                "audit":            str(audit_path),
            },
        }
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False))

        return {
            "output_csv":           str(output_csv),
            "output_xlsx":          str(output_xlsx),
            "output_visual_xlsx":   str(output_visual_xlsx),
            "output_heatmap_png":   str(output_heatmap_png),
            "output_positions_png": str(output_pos_png),
            "n_analyzed":           n_analyzed,
            "fasta_source":         fasta_source,
            "label_counts":         label_counts,
        }
