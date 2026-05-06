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
    CONSERVATION_AUDIT_{track_id}.json
"""

import datetime
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

import config
from modules.base_step import BaseTrackStep
from utils.naming import (
    COLUMN_BEST_REPRESENTATIVE,
    COLUMN_PEPTIDE,
    get_step_filename,
)
from utils.project_manager import save_project_config

console = Console(width=120)

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
    "n_variants_total",
    "n_variants_exact_match",
    "n_variants_identity_90pct",
    "n_variants_identity_80pct",
    "n_variants_passed_threshold",
    "analysis_threshold",
    "mean_max_identity",
})


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


def build_position_mutation_profile(peptide: str, failed_windows: list[str]) -> str:
    """
    Identifies which positions differ from the reference peptide across all
    failed windows, and returns a compact human-readable profile.

    Example: "pos3:L>I(2x),L>V(1x); pos6:Q>R(3x)"
    Returns "" when failed_windows is empty.
    """
    if not failed_windows:
        return ""

    position_counters: dict[int, Counter] = {}
    for window in failed_windows:
        for pos, (ref_aa, var_aa) in enumerate(zip(peptide, window), start=1):
            if ref_aa != var_aa:
                position_counters.setdefault(pos, Counter())[f"{ref_aa}>{var_aa}"] += 1

    if not position_counters:
        return ""

    parts = []
    for pos in sorted(position_counters):
        changes = position_counters[pos]
        change_strs = [
            f"{mut}({cnt}x)" if cnt > 1 else mut
            for mut, cnt in changes.most_common()
        ]
        parts.append(f"pos{pos}:{','.join(change_strs)}")
    return "; ".join(parts)


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
            "n_variants_total":            0,
            "n_variants_exact_match":      0,
            "n_variants_identity_90pct":   0,
            "n_variants_identity_80pct":   0,
            "analysis_threshold":          analysis_threshold,
            "n_variants_passed_threshold": 0,
            "conservation_summary":        "0/0",
            "mean_max_identity":           0.0,
            "conservation_label":          "conservation_unknown",
            "variants_passed_threshold":   "",
            "variants_failed_threshold":   "",
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

    metrics = {
        "n_variants_total":            n_total,
        "n_variants_exact_match":      n_exact,
        "n_variants_identity_90pct":   n_90pct,
        "n_variants_identity_80pct":   n_80pct,
        "analysis_threshold":          analysis_threshold,
        "n_variants_passed_threshold": n_passed,
        "conservation_summary":        build_conservation_summary(
                                           n_total, n_passed, analysis_threshold,
                                           n_90pct, n_80pct,
                                       ),
        "mean_max_identity":           round(mean_max_identity, 4),
        "conservation_label":          classify_conservation_label(mean_max_identity, n_total),
        "variants_passed_threshold":   ";".join(passed_entries),
        "variants_failed_threshold":   ";".join(failed_entries),
        "position_mutation_profile":   build_position_mutation_profile(peptide, failed_windows),
    }
    return metrics, alignment_tuples


# ── FASTA loading ─────────────────────────────────────────────────────────────

def load_fasta_sequences(fasta_path: Path) -> list:
    return list(SeqIO.parse(str(fasta_path), "fasta"))


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

def write_conservation_summary_xlsx(df: pd.DataFrame, output_path: Path):
    """Writes the summary table with row background colour by conservation_label."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Conservation Summary"
    columns = list(df.columns)

    for col_idx, col_name in enumerate(columns, start=1):
        cell           = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font      = Font(bold=True)
        cell.alignment = Alignment(horizontal="center")
        cell.fill      = _FILL_ORANGE if col_name in _NUMERIC_COLS else _FILL_GRAY

    for row_idx, (_, row) in enumerate(df.iterrows(), start=2):
        label    = row.get("conservation_label", "conservation_unknown")
        row_fill = _LABEL_FILL.get(label, _FILL_UNKNOWN)
        for col_idx, col_name in enumerate(columns, start=1):
            cell      = ws.cell(row=row_idx, column=col_idx, value=row[col_name])
            cell.fill = row_fill

    for col_idx, col_name in enumerate(columns, start=1):
        col_letter = ws.cell(row=1, column=col_idx).column_letter
        max_width  = max(
            len(str(col_name)),
            *(len(str(v)) for v in df[col_name].tolist()),
        )
        ws.column_dimensions[col_letter].width = min(max_width + 2, 60)

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
    t.add_column("100%",                justify="right")
    t.add_column("≥90%",                justify="right")
    t.add_column("≥80%",                justify="right")
    t.add_column(f"Passed {thr_pct}%",  justify="right")
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
        label   = row.get("conservation_label", "conservation_unknown")
        n_total = int(row.get("n_variants_total", 0))
        t.add_row(
            str(row.get(COLUMN_PEPTIDE, "")),
            str(n_total),
            f"{int(row.get('n_variants_exact_match', 0))}/{n_total}",
            f"{int(row.get('n_variants_identity_90pct', 0))}/{n_total}",
            f"{int(row.get('n_variants_identity_80pct', 0))}/{n_total}",
            f"{int(row.get('n_variants_passed_threshold', 0))}/{n_total}",
            f"{float(row.get('mean_max_identity', 0.0)):.3f}",
            label_rich.get(label, label),
        )

    console.print(t)


# ── Step ──────────────────────────────────────────────────────────────────────

class AnalyzeConservationStep(BaseTrackStep):
    step_name = "analyze_conservation"

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

        # ── Resolve FASTA source ──────────────────────────────────────────────
        variants_dir = self.track_dir / "variants"
        fasta_path   = variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta")
        fasta_source = "none"
        records: list = []

        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            records      = load_fasta_sequences(fasta_path)
            fasta_source = "variants_cache"
            console.print(f"[dim]→ Using {fasta_path.name} ({len(records)} sequences)[/dim]")
        elif _is_non_interactive():
            console.print(
                f"[yellow]⚠ No variant FASTA found for {self.track_id}. "
                "All epitopes will be labelled 'conservation_unknown'.[/yellow]"
            )
        else:
            console.print(Panel(
                f"[yellow]No variant FASTA found for {self.track_id}.[/yellow]\n\n"
                "  [cyan][1][/cyan] Continue without FASTA (label: conservation_unknown)\n"
                "  [cyan][2][/cyan] Provide a local FASTA file path",
                box=box.ROUNDED,
                title="Missing variants FASTA", title_align="left",
            ))
            while True:
                try:
                    choice = input("> ").strip()
                except EOFError:
                    choice = "1"
                if choice in ("1", ""):
                    break
                if choice == "2":
                    while True:
                        try:
                            raw_path = input("FASTA path: ").strip()
                        except EOFError:
                            raw_path = ""
                        user_path = Path(raw_path).expanduser()
                        if user_path.exists() and user_path.stat().st_size > 0:
                            records      = load_fasta_sequences(user_path)
                            fasta_source = "user_provided"
                            console.print(
                                f"[dim]→ Loaded {len(records)} sequences from {user_path}[/dim]"
                            )
                            break
                        console.print("[red]File not found or empty. Try again.[/red]")
                    break
                console.print("[dim]Type 1 or 2.[/dim]")

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

        output_csv         = conservation_dir / get_step_filename("CONSERVATION", self.track_id)
        output_xlsx        = conservation_dir / get_step_filename("CONSERVATION", self.track_id, ext="xlsx")
        output_visual_xlsx = conservation_dir / get_step_filename("CONSERVATION_VISUAL", self.track_id, ext="xlsx")
        audit_path         = conservation_dir / get_step_filename("CONSERVATION_AUDIT", self.track_id, ext="json")

        result_df.to_csv(output_csv, index=False)
        write_conservation_summary_xlsx(result_df, output_xlsx)
        write_position_heatmap_xlsx(alignment_data, output_visual_xlsx)

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
            "n_variants_total":                  len(records),
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
                "audit":        str(audit_path),
            },
        }
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False))

        return {
            "output_csv":         str(output_csv),
            "output_xlsx":        str(output_xlsx),
            "output_visual_xlsx": str(output_visual_xlsx),
            "n_analyzed":         n_analyzed,
            "fasta_source":       fasta_source,
            "label_counts":       label_counts,
        }
