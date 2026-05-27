"""Harvest per-stage epitope sets + audit counts from a finished replicate.

After `run_replicate.py` writes a project's intermediate artifacts, this
module reads them back and emits a single `metrics.json` per replicate
containing the funnel: which peptides survived each stage of the pipeline.

The downstream `make_charts.py` consumes this output to build the loss-by-
stage curve, the per-replicate consistency matrix, and the consensus Venn
diagram.

Funnel stages collected per track:
  raw_netmhcpan             — every peptide predicted (any allele) — from `0a_PRED_NET_no_nan.csv`
  raw_mhcflurry             — same for MHCFlurry — from `0a_PRED_FLURRY_no_nan.csv`
  netmhcpan_survivors       — passed NetMHCpan %Rank ≤ threshold — from `1_PRED_NET_consolidated.csv`
  mhcflurry_survivors       — passed MHCFlurry %ile ≤ threshold — from `1_PRED_FLURRY_consolidated.csv`
  consensus_intersection    — present in BOTH per-tool survivors — from `2_intersection.csv`
  immunogenic_calis         — survived Calis 2013 (> 0) — from `CONSENSUS_IMMUNOGENIC_*.csv`
  toxicity_safe             — non-toxic per ToxinPred3 — from `TOXICITY_SAFE_*.csv`
  cluster_reps_star         — ★ representative per cluster — from `CLUSTER_REPR_*.csv` filtered

Usage:
    python -m tests.validation.collect_metrics --project bench_smoke_test --out metrics.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from utils.project_manager import PROJECTS_DIR, load_project_config  # noqa: E402


# Each stage records (csv-name-template, optional ★-filter)
_TRACK_STAGE_PLAN: list[tuple[str, str, bool]] = [
    # stage_name,                         file_template,                                  star_filter
    ("raw_netmhcpan",             "consensus/0a_PRED_NET_no_nan.csv",             False),
    ("raw_mhcflurry",             "consensus/0a_PRED_FLURRY_no_nan.csv",          False),
    ("netmhcpan_survivors",       "consensus/1_PRED_NET_consolidated.csv",        False),
    ("mhcflurry_survivors",       "consensus/1_PRED_FLURRY_consolidated.csv",     False),
    ("consensus_intersection",    "consensus/2_intersection.csv",                 False),
    ("immunogenic_calis",         "consensus/CONSENSUS_IMMUNOGENIC_{track}.csv",  False),
    ("toxicity_safe",             "toxicity/TOXICITY_SAFE_{track}.csv",           False),
    ("cluster_reps_star",         "clusters/CLUSTER_REPR_{track}.csv",            True),
]


def _detect_separator(csv_path: Path) -> str:
    """Sniff the field separator from the header line.

    Pipeline CSVs are mostly comma-separated (`write_user_facing_csv` in `utils/csv_write.py`)
    but the consensus_filter consolidated audits use semicolon because their cells contain
    `;`-joined allele lists.
    """
    first_line = csv_path.read_text(encoding="utf-8-sig").splitlines()[0]
    return ";" if first_line.count(";") > first_line.count(",") else ","


def _read_peptides(csv_path: Path, star_filter: bool) -> list[str]:
    """Return the peptide column from `csv_path`. Empty list if file missing."""
    if not csv_path.exists():
        return []
    separator = _detect_separator(csv_path)
    try:
        df = pd.read_csv(csv_path, sep=separator, dtype=str, encoding="utf-8-sig")
    except pd.errors.EmptyDataError:
        return []
    if "peptide" not in df.columns:
        return []
    if star_filter:
        if "BEST_REPRESENTATIVE" not in df.columns:
            return []
        df = df[df["BEST_REPRESENTATIVE"].astype(str).str.contains("★", na=False)]
    return df["peptide"].astype(str).tolist()


def _read_consensus_audit(consensus_dir: Path) -> Optional[dict]:
    audit_path = consensus_dir / "consensus_audit_summary.json"
    if not audit_path.exists():
        return None
    try:
        return json.loads(audit_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def collect_track_funnel(
    project_name: str, track_id: str,
) -> dict:
    """Return per-stage peptide sets + counts for one track of one replicate."""
    intermediate_dir = PROJECTS_DIR / project_name / "data" / "intermediate" / track_id

    stages: dict[str, list[str]] = {}
    for stage_name, file_template, star_filter in _TRACK_STAGE_PLAN:
        relative_path_text = file_template.format(track=track_id)
        absolute_csv_path = intermediate_dir / relative_path_text
        peptide_list = _read_peptides(absolute_csv_path, star_filter)
        # Deduplicate while preserving order (set drops order, which we want stable for JSON)
        seen: dict[str, None] = {}
        for peptide in peptide_list:
            seen.setdefault(peptide, None)
            seen[peptide] = None
        stages[stage_name] = list(seen.keys())

    consensus_audit = _read_consensus_audit(intermediate_dir / "consensus")
    stage_counts = {name: len(stages[name]) for name, _, _ in _TRACK_STAGE_PLAN}

    return {
        "track_id":         track_id,
        "stage_counts":     stage_counts,
        "stages":           stages,
        "consensus_audit":  consensus_audit,
    }


def collect_project_metrics(project_name: str) -> dict:
    """Walk every track of `project_name` and return the full metrics dict."""
    project_config = load_project_config(project_name)
    track_ids = list(project_config.get("tracks", {}).keys())

    tracks_payload: dict[str, dict] = {}
    for track_id in track_ids:
        tracks_payload[track_id] = collect_track_funnel(project_name, track_id)

    # Optional: project-level master table (after global step `integrate_data`).
    # The VIEW CSV uses display-cased column headers (e.g. "Peptide"), so match case-insensitively.
    master_view_path = PROJECTS_DIR / project_name / "data" / "output" / f"MASTER_TABLE_VIEW_{project_name}.csv"
    master_view_peptides: list[str] = []
    if master_view_path.exists():
        try:
            master_df = pd.read_csv(master_view_path, dtype=str, encoding="utf-8-sig")
            peptide_col = next(
                (col for col in master_df.columns if col.strip().lower() == "peptide"),
                None,
            )
            if peptide_col is not None:
                master_view_peptides = master_df[peptide_col].astype(str).tolist()
        except Exception:
            master_view_peptides = []

    return {
        "project_name":            project_name,
        "track_count":             len(track_ids),
        "tracks":                  tracks_payload,
        "master_view_peptides":    master_view_peptides,
        "master_view_count":       len(master_view_peptides),
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, help="Project folder under projects/.")
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Destination JSON. Default: tests/validation/results/_adhoc/{project}_metrics.json",
    )
    args = parser.parse_args()

    metrics = collect_project_metrics(args.project)
    out_path = args.out or Path("tests/validation/results/_adhoc") / f"{args.project}_metrics.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(metrics, indent=2, ensure_ascii=False))
    print(f"Wrote metrics for {len(metrics['tracks'])} track(s) to {out_path}")


if __name__ == "__main__":
    main()
