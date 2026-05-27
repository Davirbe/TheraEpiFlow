"""Experiment 1 orchestrator — execution time + reproducibility.

Walks four disposable projects (HPV16_E7 × {1 track, 3 tracks}, SARS_NCAP ×
{1 track, 3 tracks}), 10 replicates each, strictly serial. Captures wall time
per step (via `run_replicate`) and per-stage epitope sets (via `collect_metrics`).

By design, never adds parallelism: one replicate finishes before the next
starts, and one project finishes before the next starts. The pipeline's
internal NetMHCpan + MHCFlurry threading inside `predict_binding` is
production-equivalent and left alone.

Usage:
    # Full Exp 1 (~12 hours serial on a Ryzen 5 mobile)
    python -m tests.validation.run_experiment_1

    # Mini smoke (only 1 project, 1 replicate)
    python -m tests.validation.run_experiment_1 \\
        --projects hpv16_e7_tr1 --replicates 1

    # Custom output root
    python -m tests.validation.run_experiment_1 --out-root tests/validation/results/exp1
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import gc
import json
import sys
import time
import traceback
from pathlib import Path
from typing import Optional

# Reroute stdin BEFORE importing pipeline modules — keeps every wizard
# auto-defaulting via `utils.console.is_interactive_session()`.
import os
if sys.stdin.isatty():
    sys.stdin = open(os.devnull, "r")

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation.collect_metrics import collect_project_metrics  # noqa: E402
from tests.validation.hardware_probe import collect_hardware_snapshot  # noqa: E402
from tests.validation.run_replicate import run_replicate  # noqa: E402
from tests.validation.seed_project import (  # noqa: E402
    PRESETS,
    destroy_project,
    seed_project,
)
from utils.archive import archive_project  # noqa: E402


EXPERIMENT_PROJECTS: list[str] = [
    "hpv16_e7_tr1",
    "hpv16_e7_tr3",
    "sars_nucleo_tr1",
    "sars_nucleo_tr3",
]


def _now_iso() -> str:
    return _datetime.datetime.now().isoformat(timespec="seconds")


def _replicate_label(preset_key: str, replicate_index: int) -> str:
    return f"bench_{preset_key}_r{replicate_index:02d}"


def _settle_cpu() -> None:
    gc.collect()
    time.sleep(2)


def run_single_replicate(
    preset_key:           str,
    replicate_index:      int,
    out_root:             Path,
    keep_project_folder:  bool,
    make_archive:         bool,
    include_predictions:  bool,
) -> dict:
    project_name = _replicate_label(preset_key, replicate_index)
    replicate_dir = out_root / preset_key / f"rep{replicate_index:02d}"
    replicate_dir.mkdir(parents=True, exist_ok=True)

    error_payload: dict | None = None
    archive_path: Optional[Path] = None

    seed_project(
        project_name=project_name,
        tracks=PRESETS[preset_key],
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
        (replicate_dir / "timings.json").write_text(
            json.dumps({"error": error_payload}, indent=2, ensure_ascii=False),
        )

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
        "preset_key":      preset_key,
        "replicate_index": replicate_index,
        "project_name":    project_name,
        "ran_at":          _now_iso(),
        "total_seconds":   timing_report.get("total_elapsed_seconds"),
        "kept_on_disk":    keep_project_folder,
        "archive_path":    str(archive_path) if archive_path else None,
        "error":           error_payload,
    }


def run_experiment(
    presets:              list[str],
    replicate_count:      int,
    out_root:             Path,
    keep_replicates:      int,
    make_archive:         bool,
    include_predictions:  bool,
) -> dict:
    out_root.mkdir(parents=True, exist_ok=True)
    hardware_snapshot = collect_hardware_snapshot(_REPO_ROOT)
    (out_root / "hardware.json").write_text(
        json.dumps(hardware_snapshot, indent=2, ensure_ascii=False),
    )

    started_at = _now_iso()
    print(f"[exp1] started at {started_at}")
    print(f"[exp1] presets: {', '.join(presets)}")
    print(f"[exp1] replicates per preset: {replicate_count}")
    print(f"[exp1] keep first N replicate folders on disk: {keep_replicates}")
    print(f"[exp1] tar.gz bundle per replicate: {make_archive} (predictions={'included' if include_predictions else 'excluded'})")
    print(f"[exp1] out_root: {out_root}\n")

    replicate_records: list[dict] = []
    overall_started = time.perf_counter()

    for preset_key in presets:
        print(f"[exp1] ── preset '{preset_key}' ──")
        for replicate_index in range(1, replicate_count + 1):
            print(f"  rep{replicate_index:02d} starting at {_now_iso()}…")
            replicate_started = time.perf_counter()
            keep_this_folder = replicate_index <= keep_replicates
            try:
                record = run_single_replicate(
                    preset_key=preset_key,
                    replicate_index=replicate_index,
                    out_root=out_root,
                    keep_project_folder=keep_this_folder,
                    make_archive=make_archive,
                    include_predictions=include_predictions,
                )
            except Exception as outer_exception:
                record = {
                    "preset_key":      preset_key,
                    "replicate_index": replicate_index,
                    "project_name":    _replicate_label(preset_key, replicate_index),
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

    overall_elapsed = time.perf_counter() - overall_started
    summary = {
        "started_at":           started_at,
        "ended_at":             _now_iso(),
        "wall_seconds":         round(overall_elapsed, 2),
        "presets":              presets,
        "replicate_count":      replicate_count,
        "keep_replicates":      keep_replicates,
        "make_archive":         make_archive,
        "include_predictions":  include_predictions,
        "replicates":           replicate_records,
    }
    (out_root / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
    )
    print(f"\n[exp1] done — wall {overall_elapsed/60:.1f} min — summary at {out_root / 'summary.json'}")
    return summary


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--projects", nargs="+", default=EXPERIMENT_PROJECTS,
        choices=EXPERIMENT_PROJECTS,
        help="Which preset projects to run (default: all four).",
    )
    parser.add_argument(
        "--replicates", type=int, default=10,
        help="Number of replicates per preset (default: 10).",
    )
    parser.add_argument(
        "--out-root", type=Path, default=Path("tests/validation/results/exp1"),
        help="Where per-replicate JSONs and the summary go.",
    )
    parser.add_argument(
        "--keep-replicates", type=int, default=1,
        help="How many full project folders to keep on disk per preset (default: 1). "
             "The first N replicates per preset stay under projects/bench_*; the rest "
             "are deleted after their metrics + tar.gz are written.",
    )
    parser.add_argument(
        "--no-bundle", action="store_true",
        help="Skip the tar.gz archive of each replicate (default: bundle every replicate). "
             "Use only if disk space is tight AND you do not need archival copies.",
    )
    parser.add_argument(
        "--include-predictions", action="store_true",
        help="Include raw NetMHCpan/MHCFlurry CSVs in the tar.gz bundle "
             "(default: excluded — they are large and re-derivable).",
    )
    args = parser.parse_args()

    if args.replicates < 1:
        parser.error("--replicates must be ≥ 1")
    if args.keep_replicates < 0:
        parser.error("--keep-replicates must be ≥ 0")
    if args.keep_replicates > args.replicates:
        parser.error("--keep-replicates cannot exceed --replicates")

    run_experiment(
        presets=args.projects,
        replicate_count=args.replicates,
        out_root=args.out_root,
        keep_replicates=args.keep_replicates,
        make_archive=not args.no_bundle,
        include_predictions=args.include_predictions,
    )


if __name__ == "__main__":
    main()
