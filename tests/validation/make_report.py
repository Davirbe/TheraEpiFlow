"""Consolidate the validation suite into a single markdown report.

Reads the per-experiment summaries (`results/exp1/summary.json`,
`results/exp2/summary.json`) and the hardware probe (`results/exp1/hardware.json`),
and renders `figures/report.md` with:

  - Hardware/environment header (Lenovo IdeaPad + Ryzen 5 + WSL2 + 8 GB)
  - Measurement caveat (serial mode, hardware-dependent timings)
  - Experiment 1 summary tables (median/quartile total times per preset)
  - Inline references to every PNG in `figures/`
  - Experiment 2 summary table (★ peptides per track, IEDB hit ranges)
  - Pointer to the supplementary `iedb_pmid_titles.csv`

The markdown is intentionally minimal — Davi expands into prose for the
paper itself. The file's value is keeping the numbers, figures and
provenance in one place.

Usage:
    python -m tests.validation.make_report \\
        --exp1-root tests/validation/results/exp1 \\
        --exp2-root tests/validation/results/exp2 \\
        --figures-root tests/validation/figures \\
        --out tests/validation/figures/report.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np


def _read_json_optional(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _summarize_replicate_totals(summary_payload: dict) -> dict:
    totals = [
        rec.get("total_seconds")
        for rec in summary_payload.get("replicates", [])
        if rec.get("total_seconds") is not None
    ]
    if not totals:
        return {"count": 0, "median_s": None, "p25_s": None, "p75_s": None, "min_s": None, "max_s": None}
    arr = np.array(totals)
    return {
        "count":     int(arr.size),
        "median_s":  float(np.median(arr)),
        "p25_s":     float(np.percentile(arr, 25)),
        "p75_s":     float(np.percentile(arr, 75)),
        "min_s":     float(arr.min()),
        "max_s":     float(arr.max()),
    }


def _format_seconds(value: Optional[float]) -> str:
    if value is None:
        return "—"
    if value < 60:
        return f"{value:.1f}s"
    return f"{value/60:.1f} min"


def _render_hardware_block(hardware_payload: Optional[dict]) -> str:
    if not hardware_payload:
        return "_Hardware snapshot missing — re-run `hardware_probe.py`._\n"

    host = hardware_payload.get("host", {})
    packages = hardware_payload.get("packages", {})
    tools = hardware_payload.get("tools", {})

    lines = [
        "### Execution environment\n",
        "| Field | Value |",
        "|---|---|",
        f"| Captured at | `{hardware_payload.get('captured_at', '?')}` |",
        f"| Platform | `{host.get('platform', '?')}` |",
        f"| WSL | `{host.get('is_wsl', '?')}` |",
        f"| Kernel | `{host.get('kernel', '?')}` |",
        f"| Python | `{host.get('python_version', '?')}` |",
        f"| netMHCpan | `{tools.get('netmhcpan_in_path', '?')}` |",
        "",
        "### Package versions\n",
        "| Package | Version |",
        "|---|---|",
    ]
    for name in sorted(packages.keys()):
        lines.append(f"| `{name}` | `{packages[name]}` |")
    lines.append("")
    return "\n".join(lines)


def _render_exp1_block(exp1_root: Path, figures_root: Path) -> str:
    summary = _read_json_optional(exp1_root / "summary.json")
    if not summary:
        return "## Experiment 1 — Time + Reproducibility\n\n_summary.json missing._\n"

    presets = sorted({r.get("preset_key") for r in summary.get("replicates", []) if r.get("preset_key")})
    grouped = {p: {"replicates": []} for p in presets}
    for replicate in summary.get("replicates", []):
        grouped.setdefault(replicate["preset_key"], {"replicates": []})["replicates"].append(replicate)

    lines = [
        "## Experiment 1 — Execution time + reproducibility\n",
        f"- Started: `{summary.get('started_at')}`",
        f"- Ended:   `{summary.get('ended_at')}`",
        f"- Total wall: **{_format_seconds(summary.get('wall_seconds'))}**\n",
        "### Replicate timing summary\n",
        "| Preset | n | median | IQR | min | max |",
        "|---|---|---|---|---|---|",
    ]
    for preset_key in presets:
        s = _summarize_replicate_totals(grouped[preset_key])
        iqr = (
            f"{_format_seconds(s['p25_s'])} – {_format_seconds(s['p75_s'])}"
            if s["count"] else "—"
        )
        lines.append(
            f"| `{preset_key}` | {s['count']} | "
            f"{_format_seconds(s['median_s'])} | {iqr} | "
            f"{_format_seconds(s['min_s'])} | {_format_seconds(s['max_s'])} |"
        )
    lines.append("")

    lines.append("### Figures\n")
    for preset_key in presets:
        for chart_kind in ("loss_by_stage", "time_by_stage", "consistency_matrix", "venn"):
            png_path = figures_root / f"exp1_{chart_kind}_{preset_key}.png"
            if png_path.exists():
                lines.append(f"![exp1 {chart_kind} {preset_key}]({png_path.relative_to(figures_root)})")
                lines.append("")
    for cross_chart in ("exp1_time_tr1_vs_tr3.png", "exp1_time_protein_size.png"):
        png_path = figures_root / cross_chart
        if png_path.exists():
            lines.append(f"![{cross_chart}]({png_path.relative_to(figures_root)})")
            lines.append("")

    return "\n".join(lines)


def _render_exp2_block(exp2_root: Path, figures_root: Path) -> str:
    summary  = _read_json_optional(exp2_root / "summary.json")
    iedb     = _read_json_optional(exp2_root / "iedb_hits.json")
    if not summary:
        return "## Experiment 2 — Clinical validation\n\n_summary.json missing._\n"

    timing = _summarize_replicate_totals(summary)
    lines = [
        "## Experiment 2 — Clinical validation (IEDB IQ-API)\n",
        f"- Started: `{summary.get('started_at')}`",
        f"- Ended:   `{summary.get('ended_at')}`",
        f"- Wall: **{_format_seconds(summary.get('wall_seconds'))}**",
        f"- Pipeline replicates: {timing['count']}",
        f"- Median pipeline wall per replicate: {_format_seconds(timing['median_s'])}",
        "",
    ]

    if iedb:
        lines.append("### Literature support summary\n")
        lines.append("| Track | ★ unique | ★ with hits | neg.control unique | neg.control with hits |")
        lines.append("|---|---|---|---|---|")
        for track_id, payload in iedb.get("tracks", {}).items():
            stars = payload.get("star", [])
            negatives = payload.get("negative_control", [])
            stars_with_hits = sum(
                1 for r in stars
                if (r.get("tcell_assay_count", 0) + r.get("mhc_assay_count", 0)) > 0
            )
            negatives_with_hits = sum(
                1 for r in negatives
                if (r.get("tcell_assay_count", 0) + r.get("mhc_assay_count", 0)) > 0
            )
            lines.append(
                f"| `{track_id}` | {len(stars)} | {stars_with_hits} | "
                f"{len(negatives)} | {negatives_with_hits} |"
            )
        lines.append("")

    lines.append("### Figures\n")
    for track_id in (iedb or {}).get("tracks", {}):
        png_path = figures_root / f"exp2_iedb_counts_{track_id}.png"
        if png_path.exists():
            lines.append(f"![exp2 IEDB counts {track_id}]({png_path.relative_to(figures_root)})")
            lines.append("")
    neg_png = figures_root / "exp2_neg_control_comparison.png"
    if neg_png.exists():
        lines.append(f"![exp2 negative-control comparison]({neg_png.relative_to(figures_root)})")
        lines.append("")

    pmid_csv = exp2_root / "iedb_pmid_titles.csv"
    if pmid_csv.exists():
        lines.append(
            f"_Supplementary publications table (PMID + title): "
            f"`{pmid_csv.relative_to(exp2_root.parent.parent.parent)}`_\n"
        )

    return "\n".join(lines)


def render_report(
    exp1_root: Optional[Path],
    exp2_root: Optional[Path],
    figures_root: Path,
    out_path: Path,
) -> None:
    sections: list[str] = []
    sections.append("# TheraEpiFlow — Validation Suite Report\n")
    sections.append(
        "> Timings in this report were measured on a Lenovo IdeaPad with AMD Ryzen 5 5500U, "
        "AMD Radeon iGPU, 8 GB RAM, running WSL2 (Ubuntu) on Windows. "
        "All replicates ran one-at-a-time on a quiet machine (no Chrome / IDE / Teams). "
        "Absolute seconds do not generalize to other hardware; compare _ratios_ (tr3 : tr1, "
        "SARS_NCAP : HPV16_E7) and the _shape_ of the loss curves.\n"
    )

    hardware_payload = (
        _read_json_optional(exp1_root / "hardware.json") if exp1_root else None
    ) or (
        _read_json_optional(exp2_root / "hardware.json") if exp2_root else None
    )
    sections.append(_render_hardware_block(hardware_payload))

    if exp1_root:
        sections.append(_render_exp1_block(exp1_root, figures_root))
    if exp2_root:
        sections.append(_render_exp2_block(exp2_root, figures_root))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(sections), encoding="utf-8")
    print(f"Report written to {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp1-root", type=Path, default=None)
    parser.add_argument("--exp2-root", type=Path, default=None)
    parser.add_argument("--figures-root", type=Path, default=Path("tests/validation/figures"))
    parser.add_argument("--out", type=Path, default=Path("tests/validation/figures/report.md"))
    args = parser.parse_args()

    if not (args.exp1_root or args.exp2_root):
        parser.error("Pass at least --exp1-root or --exp2-root.")

    render_report(args.exp1_root, args.exp2_root, args.figures_root, args.out)


if __name__ == "__main__":
    main()
