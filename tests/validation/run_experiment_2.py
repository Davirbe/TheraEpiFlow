"""Experiment 2 orchestrator — clinical literature validation via IEDB.

Runs one multi-organism project (HPV16_E7 + SARS_NCAP, 2 tracks, `clinic`
preset) 10 times serially. After the runs, harvests:

  - The ★ peptides at the end of the pipeline (per track), deduplicated
    across all 10 replicates.
  - A negative-control set per track: the 10 peptides with the worst
    `max(net_percentile, flurry_percentile)` among those rejected by the
    consensus filter's Stage-1 intersection. These are by definition the
    "least binder-like" peptides that the pipeline filtered out — they
    should accumulate few or no IEDB literature hits.

For each peptide in both sets, calls the IEDB IQ-API (`tests.validation.iedb_query`)
to count T-cell + MHC-Class-I assays and collect supporting PMIDs.

Usage:
    # Full Exp 2 (~6 hours pipeline + queries)
    python -m tests.validation.run_experiment_2

    # Mini smoke (1 replicate)
    python -m tests.validation.run_experiment_2 --replicates 1
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import gc
import json
import os
import sys
import time
import traceback
from pathlib import Path

if sys.stdin.isatty():
    sys.stdin = open(os.devnull, "r")

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from typing import Optional

import pandas as pd  # noqa: E402

from tests.validation.collect_metrics import collect_project_metrics  # noqa: E402
from tests.validation.hardware_probe import collect_hardware_snapshot  # noqa: E402
from tests.validation.iedb_query import query_peptides  # noqa: E402
from tests.validation.run_replicate import run_replicate  # noqa: E402
from tests.validation.seed_project import (  # noqa: E402
    PRESETS,
    destroy_project,
    seed_project,
)
from utils.archive import archive_project  # noqa: E402
from utils.project_manager import PROJECTS_DIR  # noqa: E402


CLINIC_PRESET_KEY = "clinic"
NEGATIVE_CONTROL_SIZE = 10


def _now_iso() -> str:
    return _datetime.datetime.now().isoformat(timespec="seconds")


def _settle_cpu() -> None:
    gc.collect()
    time.sleep(2)


def _replicate_label(replicate_index: int) -> str:
    return f"bench_{CLINIC_PRESET_KEY}_r{replicate_index:02d}"


def _read_csv_smart(csv_path: Path) -> pd.DataFrame:
    """Read CSV auto-detecting comma vs semicolon separator (consensus audit files use `;`)."""
    if not csv_path.exists():
        return pd.DataFrame()
    first_line = csv_path.read_text(encoding="utf-8-sig").splitlines()[0]
    separator = ";" if first_line.count(";") > first_line.count(",") else ","
    try:
        return pd.read_csv(csv_path, sep=separator, encoding="utf-8-sig")
    except Exception:
        return pd.DataFrame()


def collect_worst_eliminated_peptides(
    project_name: str, track_id: str, top_n: int,
) -> list[dict]:
    """Pull the `top_n` worst eliminated peptides for one track.

    A peptide is "eliminated" if it appears in the raw NetMHCpan or MHCFlurry
    audit (`0a_*_no_nan.csv`) but not in the consensus intersection
    (`2_intersection.csv`). Among the eliminated, we sort by
    `max(net_percentile, flurry_percentile)` descending — the highest
    percentile (= the weakest binder) wins the negative-control slot.

    Returns a list of dicts `[{peptide, net_pct, flurry_pct, max_pct}, …]`.
    Missing values (e.g. peptide rejected at the NaN stage) are filled with infinity
    so they sort to the front; that matches "worst binder" semantics.
    """
    consensus_dir = PROJECTS_DIR / project_name / "data" / "intermediate" / track_id / "consensus"

    intersection_df = _read_csv_smart(consensus_dir / "2_intersection.csv")
    if "peptide" not in intersection_df.columns:
        kept_set: set[str] = set()
    else:
        kept_set = set(intersection_df["peptide"].astype(str).tolist())

    net_df    = _read_csv_smart(consensus_dir / "0a_PRED_NET_no_nan.csv")
    flurry_df = _read_csv_smart(consensus_dir / "0a_PRED_FLURRY_no_nan.csv")

    if "peptide" not in net_df.columns and "peptide" not in flurry_df.columns:
        return []

    net_pct_col = next(
        (c for c in net_df.columns if "netmhcpan" in c.lower() and "percentile" in c.lower()),
        None,
    )
    flurry_pct_col = next(
        (c for c in flurry_df.columns if "mhcflurry" in c.lower() and "percentile" in c.lower()),
        None,
    )

    # Min percentile per peptide across all alleles (best binder allele wins).
    net_best    = (
        net_df.groupby("peptide")[net_pct_col].min().rename("net_pct")
        if net_pct_col is not None else pd.Series(dtype=float, name="net_pct")
    )
    flurry_best = (
        flurry_df.groupby("peptide")[flurry_pct_col].min().rename("flurry_pct")
        if flurry_pct_col is not None else pd.Series(dtype=float, name="flurry_pct")
    )

    universe = pd.concat([net_best, flurry_best], axis=1).reset_index()
    universe = universe[~universe["peptide"].isin(kept_set)]
    if universe.empty:
        return []

    universe["max_pct"] = universe[["net_pct", "flurry_pct"]].max(axis=1, skipna=False)
    universe["max_pct"] = universe["max_pct"].fillna(float("inf"))
    universe = universe.sort_values("max_pct", ascending=False).head(top_n)

    return [
        {
            "peptide":     str(row["peptide"]),
            "net_pct":     None if pd.isna(row.get("net_pct"))    else float(row["net_pct"]),
            "flurry_pct":  None if pd.isna(row.get("flurry_pct")) else float(row["flurry_pct"]),
            "max_pct":     None if not pd.notna(row.get("max_pct")) or row["max_pct"] == float("inf") else float(row["max_pct"]),
        }
        for _, row in universe.iterrows()
    ]


def run_single_replicate(
    replicate_index:      int,
    out_root:             Path,
    keep_project_folder:  bool,
    make_archive:         bool,
    include_predictions:  bool,
) -> dict:
    project_name  = _replicate_label(replicate_index)
    replicate_dir = out_root / f"rep{replicate_index:02d}"
    replicate_dir.mkdir(parents=True, exist_ok=True)
    error_payload: dict | None = None
    archive_path: Optional[Path] = None

    seed_project(
        project_name=project_name,
        tracks=PRESETS[CLINIC_PRESET_KEY],
        overwrite=True,
    )
    _settle_cpu()

    try:
        timing_report = run_replicate(project_name)
        (replicate_dir / "timings.json").write_text(
            json.dumps(timing_report, indent=2, ensure_ascii=False),
        )
    except Exception as run_exception:
        error_payload = {
            "phase":     "run_replicate",
            "message":   str(run_exception),
            "traceback": traceback.format_exc(),
        }
        timing_report = {"error": error_payload}

    try:
        metrics_report = collect_project_metrics(project_name)
        (replicate_dir / "metrics.json").write_text(
            json.dumps(metrics_report, indent=2, ensure_ascii=False),
        )
    except Exception as collect_exception:
        metrics_report = {"error": {
            "phase":     "collect_metrics",
            "message":   str(collect_exception),
            "traceback": traceback.format_exc(),
        }}
        (replicate_dir / "metrics.json").write_text(
            json.dumps(metrics_report, indent=2, ensure_ascii=False),
        )

    # Capture negative control BEFORE the project is destroyed.
    eliminated_per_track: dict[str, list[dict]] = {}
    for track_id in [
        f"H16_E7",   # HPV16_E7 — built by _hpv16_e7_track()
        f"SARS2_N",  # SARS-CoV-2 Nucleocapsid — built by _sars_ncap_track()
    ]:
        try:
            eliminated_per_track[track_id] = collect_worst_eliminated_peptides(
                project_name, track_id, NEGATIVE_CONTROL_SIZE,
            )
        except Exception as eliminate_exception:
            eliminated_per_track[track_id] = []
            error_payload = {
                "phase":     "negative_control",
                "track_id":  track_id,
                "message":   str(eliminate_exception),
                "traceback": traceback.format_exc(),
            }
    (replicate_dir / "negative_control.json").write_text(
        json.dumps(eliminated_per_track, indent=2, ensure_ascii=False),
    )

    if make_archive:
        try:
            archive_path = archive_project(
                project_name=project_name,
                destination_dir=replicate_dir / "bundle",
                include_predictions=include_predictions,
            )
        except Exception as archive_exception:
            error_payload = error_payload or {
                "phase":     "archive_project",
                "message":   str(archive_exception),
                "traceback": traceback.format_exc(),
            }

    if not keep_project_folder:
        destroy_project(project_name, missing_ok=True)

    return {
        "replicate_index": replicate_index,
        "project_name":    project_name,
        "ran_at":          _now_iso(),
        "total_seconds":   timing_report.get("total_elapsed_seconds") if isinstance(timing_report, dict) else None,
        "kept_on_disk":    keep_project_folder,
        "archive_path":    str(archive_path) if archive_path else None,
        "error":           error_payload,
    }


def aggregate_iedb_lookups(
    out_root: Path, replicate_count: int,
) -> dict:
    """After all replicates finish: union the ★ peptides + negative controls, query IEDB,
    write the aggregated payload + per-track CSVs."""

    aggregated_stars:    dict[str, set[str]] = {"H16_E7": set(), "SARS2_N": set()}
    aggregated_negatives: dict[str, set[str]] = {"H16_E7": set(), "SARS2_N": set()}
    presence_count:      dict[str, dict[str, int]] = {"H16_E7": {}, "SARS2_N": {}}

    for replicate_index in range(1, replicate_count + 1):
        replicate_dir = out_root / f"rep{replicate_index:02d}"
        metrics_path  = replicate_dir / "metrics.json"
        negctl_path   = replicate_dir / "negative_control.json"

        if metrics_path.exists():
            metrics_payload = json.loads(metrics_path.read_text(encoding="utf-8"))
            for track_id, track_payload in metrics_payload.get("tracks", {}).items():
                if track_id not in aggregated_stars:
                    continue
                star_peptides = track_payload.get("stages", {}).get("cluster_reps_star", [])
                for peptide in star_peptides:
                    aggregated_stars[track_id].add(peptide)
                    presence_count[track_id][peptide] = presence_count[track_id].get(peptide, 0) + 1

        if negctl_path.exists():
            negctl_payload = json.loads(negctl_path.read_text(encoding="utf-8"))
            for track_id, records in negctl_payload.items():
                if track_id not in aggregated_negatives:
                    continue
                for record in records:
                    aggregated_negatives[track_id].add(record["peptide"])

    iedb_payload: dict = {"queried_at": _now_iso(), "tracks": {}}

    for track_id in aggregated_stars:
        star_list      = sorted(aggregated_stars[track_id])
        negative_list  = sorted(aggregated_negatives[track_id])
        print(f"[exp2] querying IEDB for {track_id}: {len(star_list)} ★ + {len(negative_list)} negative-controls")

        star_records      = [r.to_dict() for r in query_peptides(star_list)]
        negative_records  = [r.to_dict() for r in query_peptides(negative_list)]

        for record in star_records:
            record["presence_count"] = presence_count[track_id].get(record["peptide"], 0)
            record["category"]       = "star"
        for record in negative_records:
            record["presence_count"] = 0
            record["category"]       = "negative_control"

        iedb_payload["tracks"][track_id] = {
            "star":             star_records,
            "negative_control": negative_records,
        }

    (out_root / "iedb_hits.json").write_text(
        json.dumps(iedb_payload, indent=2, ensure_ascii=False),
    )

    # Tidy per-track CSV: one row per peptide with the aggregated counts.
    for track_id, payload in iedb_payload["tracks"].items():
        rows = []
        for category_label, category_records in (("star", payload["star"]), ("negative_control", payload["negative_control"])):
            for record in category_records:
                rows.append({
                    "category":          category_label,
                    "peptide":           record["peptide"],
                    "presence_count":    record.get("presence_count", 0),
                    "tcell_assay_count": record["tcell_assay_count"],
                    "mhc_assay_count":   record["mhc_assay_count"],
                    "publication_count": record["publication_count"],
                    "errors":            "; ".join(record.get("errors", [])),
                })
        if rows:
            csv_path = out_root / f"iedb_hits_{track_id}.csv"
            pd.DataFrame(rows).to_csv(csv_path, index=False, encoding="utf-8-sig")

    # Long-format publication table (one row per PMID) for the paper's supplementary material.
    pmid_rows = []
    for track_id, payload in iedb_payload["tracks"].items():
        for category_label, category_records in (("star", payload["star"]), ("negative_control", payload["negative_control"])):
            for record in category_records:
                for publication in record.get("publications", []):
                    pmid_rows.append({
                        "track":   track_id,
                        "category": category_label,
                        "peptide": record["peptide"],
                        "pmid":    publication["pmid"],
                        "title":   publication["title"],
                    })
    if pmid_rows:
        pd.DataFrame(pmid_rows).to_csv(out_root / "iedb_pmid_titles.csv", index=False, encoding="utf-8-sig")

    return iedb_payload


def run_experiment(
    replicate_count:      int,
    out_root:             Path,
    keep_replicates:      int,
    make_archive:         bool,
    include_predictions:  bool,
    skip_iedb:            bool,
) -> dict:
    out_root.mkdir(parents=True, exist_ok=True)
    hardware_snapshot = collect_hardware_snapshot(_REPO_ROOT)
    (out_root / "hardware.json").write_text(
        json.dumps(hardware_snapshot, indent=2, ensure_ascii=False),
    )

    started_at = _now_iso()
    print(f"[exp2] started at {started_at}")
    print(f"[exp2] preset: {CLINIC_PRESET_KEY} (HPV16_E7 + SARS_NCAP)")
    print(f"[exp2] replicates: {replicate_count}")
    print(f"[exp2] keep first N replicate folders on disk: {keep_replicates}")
    print(f"[exp2] tar.gz bundle per replicate: {make_archive} (predictions={'included' if include_predictions else 'excluded'})")
    print(f"[exp2] out_root: {out_root}\n")

    replicate_records: list[dict] = []
    overall_started = time.perf_counter()

    for replicate_index in range(1, replicate_count + 1):
        print(f"  rep{replicate_index:02d} starting at {_now_iso()}…")
        replicate_started = time.perf_counter()
        keep_this_folder = replicate_index <= keep_replicates
        try:
            record = run_single_replicate(
                replicate_index=replicate_index,
                out_root=out_root,
                keep_project_folder=keep_this_folder,
                make_archive=make_archive,
                include_predictions=include_predictions,
            )
        except Exception as outer_exception:
            record = {
                "replicate_index": replicate_index,
                "project_name":    _replicate_label(replicate_index),
                "ran_at":          _now_iso(),
                "total_seconds":   None,
                "error":           {
                    "phase":     "orchestrator",
                    "message":   str(outer_exception),
                    "traceback": traceback.format_exc(),
                },
            }
        replicate_records.append(record)
        elapsed_replicate = time.perf_counter() - replicate_started
        status = "ERR" if record.get("error") else "ok "
        total_seconds = record.get("total_seconds")
        elapsed_inner = f"{total_seconds:.1f}s" if total_seconds is not None else "—"
        print(
            f"  rep{replicate_index:02d} [{status}] pipeline={elapsed_inner} "
            f"orchestration={elapsed_replicate:.1f}s"
        )

    iedb_payload = {} if skip_iedb else aggregate_iedb_lookups(out_root, replicate_count)

    overall_elapsed = time.perf_counter() - overall_started
    summary = {
        "started_at":           started_at,
        "ended_at":             _now_iso(),
        "wall_seconds":         round(overall_elapsed, 2),
        "replicate_count":      replicate_count,
        "keep_replicates":      keep_replicates,
        "make_archive":         make_archive,
        "include_predictions":  include_predictions,
        "skipped_iedb":         skip_iedb,
        "replicates":           replicate_records,
        "iedb_summary":    {
            track_id: {
                "star_count":              len(payload.get("star", [])),
                "negative_control_count":  len(payload.get("negative_control", [])),
            }
            for track_id, payload in iedb_payload.get("tracks", {}).items()
        } if iedb_payload else {},
    }
    (out_root / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
    )
    print(f"\n[exp2] done — wall {overall_elapsed/60:.1f} min — summary at {out_root / 'summary.json'}")
    return summary


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--replicates", type=int, default=10,
        help="Number of replicates (default: 10).",
    )
    parser.add_argument(
        "--out-root", type=Path, default=Path("tests/validation/results/exp2"),
        help="Where per-replicate JSONs, IEDB hits and the summary go.",
    )
    parser.add_argument(
        "--keep-replicates", type=int, default=1,
        help="How many full project folders to keep on disk (default: 1). The first N "
             "replicates stay under projects/bench_clinic_r*; the rest are deleted after "
             "their metrics + negative_control + tar.gz are written.",
    )
    parser.add_argument(
        "--no-bundle", action="store_true",
        help="Skip the tar.gz archive of each replicate (default: bundle every replicate).",
    )
    parser.add_argument(
        "--include-predictions", action="store_true",
        help="Include raw NetMHCpan/MHCFlurry CSVs in the tar.gz bundle "
             "(default: excluded — they are large and re-derivable).",
    )
    parser.add_argument(
        "--skip-iedb", action="store_true",
        help="Run the pipeline replicates but skip the final IEDB query step (faster for smoke).",
    )
    args = parser.parse_args()

    if args.replicates < 1:
        parser.error("--replicates must be ≥ 1")
    if args.keep_replicates < 0:
        parser.error("--keep-replicates must be ≥ 0")
    if args.keep_replicates > args.replicates:
        parser.error("--keep-replicates cannot exceed --replicates")

    run_experiment(
        replicate_count=args.replicates,
        out_root=args.out_root,
        keep_replicates=args.keep_replicates,
        make_archive=not args.no_bundle,
        include_predictions=args.include_predictions,
        skip_iedb=args.skip_iedb,
    )


if __name__ == "__main__":
    main()
