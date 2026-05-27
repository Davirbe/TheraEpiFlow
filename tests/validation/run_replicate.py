"""Run one replicate of the pipeline end-to-end and capture per-step wall time.

This is the unit benchmark inside the validation suite. It walks STEP_REGISTRY
in order, instantiates each step class for every track (or once for global
steps), wraps `execute()` with `time.perf_counter()`, and writes a single
`timings.json` describing the run.

Strictly serial:
  - No multiprocessing, no thread pool, no asyncio added here.
  - The pipeline's internal NetMHCpan + MHCFlurry threading inside
    `predict_binding` is preserved (it is production-equivalent).
  - stdin is rerouted to /dev/null on startup so every Rich prompt falls back
    to its default — wizards never block, never crash.

Usage:
    python -m tests.validation.run_replicate --project bench_smoke_test
    python -m tests.validation.run_replicate --project bench_smoke_test --seed-only
    python -m tests.validation.run_replicate --project bench_smoke_test --out timings.json
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import gc
import json
import os
import sys
import time
from pathlib import Path
from typing import Optional

# Reroute stdin BEFORE importing any pipeline module so every Rich helper that
# checks `sys.stdin.isatty()` resolves to "non-interactive" mode for this run.
if sys.stdin.isatty():
    sys.stdin = open(os.devnull, "r")

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from step_registry import STEP_REGISTRY, _import_step_class  # noqa: E402
from utils.project_manager import PROJECTS_DIR, load_project_config  # noqa: E402


def _now_iso() -> str:
    return _datetime.datetime.now().isoformat(timespec="seconds")


def _track_step_timing(
    cls, project_name: str, project_config: dict, track_id: str,
) -> dict:
    instance = cls(project_name=project_name, project_config=project_config, track_id=track_id)
    started_at = time.perf_counter()
    outcome = instance.execute(force_rerun=False, reconfigure=False)
    elapsed_seconds = time.perf_counter() - started_at
    return {
        "track_id":         track_id,
        "elapsed_seconds":  round(elapsed_seconds, 4),
        "status":           outcome.get("status"),
        "error":            outcome.get("error_message"),
    }


def _global_step_timing(cls, project_name: str, project_config: dict) -> dict:
    instance = cls(project_name=project_name, project_config=project_config)
    started_at = time.perf_counter()
    outcome = instance.execute(force_rerun=False, reconfigure=False)
    elapsed_seconds = time.perf_counter() - started_at
    return {
        "track_id":         None,
        "elapsed_seconds":  round(elapsed_seconds, 4),
        "status":           outcome.get("status"),
        "error":            outcome.get("error_message"),
    }


def run_replicate(project_name: str) -> dict:
    """Run every step in STEP_REGISTRY order; return a timing report dict."""
    project_config = load_project_config(project_name)
    track_ids = list(project_config.get("tracks", {}).keys())
    if not track_ids:
        raise ValueError(f"Project '{project_name}' has no tracks defined.")

    report: dict = {
        "project_name":        project_name,
        "track_count":         len(track_ids),
        "track_ids":           track_ids,
        "started_at":          _now_iso(),
        "ended_at":            None,
        "total_elapsed_seconds": None,
        "steps":               [],
    }

    overall_started = time.perf_counter()

    for step_name, (_, _, scope) in STEP_REGISTRY.items():
        cls = _import_step_class(step_name)
        if cls is None:
            report["steps"].append({
                "step_name":      step_name,
                "scope":          scope,
                "runs":           [],
                "import_error":   True,
            })
            continue

        runs: list[dict] = []
        if scope == "track":
            for track_id in track_ids:
                runs.append(_track_step_timing(cls, project_name, project_config, track_id))
        else:  # global
            runs.append(_global_step_timing(cls, project_name, project_config))

        report["steps"].append({
            "step_name":     step_name,
            "scope":         scope,
            "runs":          runs,
            "import_error":  False,
        })

    overall_elapsed = time.perf_counter() - overall_started
    report["ended_at"] = _now_iso()
    report["total_elapsed_seconds"] = round(overall_elapsed, 4)
    return report


def _validate_seed_only(project_name: str) -> None:
    """Smoke check: confirm the project exists, the config loads, every step class imports."""
    project_dir = PROJECTS_DIR / project_name
    if not project_dir.exists():
        raise SystemExit(f"Project folder does not exist: {project_dir}")

    project_config = load_project_config(project_name)
    track_ids = list(project_config.get("tracks", {}).keys())
    if not track_ids:
        raise SystemExit(f"Project '{project_name}' has no tracks defined.")

    print(f"Project:      {project_name}")
    print(f"Tracks:       {', '.join(track_ids)} ({len(track_ids)} total)")
    print(f"Steps in order:")
    missing_classes: list[str] = []
    for step_name, (mod_path, cls_name, scope) in STEP_REGISTRY.items():
        cls = _import_step_class(step_name)
        ok_flag = "OK " if cls is not None else "ERR"
        if cls is None:
            missing_classes.append(step_name)
        print(f"  [{ok_flag}] {scope:<6} {step_name:<24} ({mod_path}.{cls_name})")

    if missing_classes:
        raise SystemExit(f"Missing step classes: {', '.join(missing_classes)}")
    print("\nAll step classes importable. Ready for run_replicate().")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project", required=True,
        help="Project folder under projects/ (typically bench_<preset>_r<NN>).",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Destination JSON path for the timing report. "
             "Default: tests/validation/results/_adhoc/{project}_timings.json",
    )
    parser.add_argument(
        "--seed-only", action="store_true",
        help="Validate setup (project exists, classes import) without running any step.",
    )
    args = parser.parse_args()

    if args.seed_only:
        _validate_seed_only(args.project)
        return

    gc.collect()
    time.sleep(2)  # let the CPU settle before time.perf_counter() starts

    report = run_replicate(args.project)

    out_path = args.out or Path("tests/validation/results/_adhoc") / f"{args.project}_timings.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False))

    print(f"\nTotal: {report['total_elapsed_seconds']:.1f}s over {len(report['steps'])} steps.")
    print(f"Report: {out_path}")


if __name__ == "__main__":
    main()
