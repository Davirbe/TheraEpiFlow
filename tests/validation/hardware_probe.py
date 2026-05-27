"""Capture the machine and software environment of a benchmark run.

Writes a single JSON file with everything the reader of the paper needs to
contextualize the timing numbers: CPU model, RAM, kernel, WSL detection,
Python interpreter, package versions and external binary versions.

Usage:
    python -m tests.validation.hardware_probe --out tests/validation/results/exp1/hardware.json
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import json
import platform
import subprocess
import sys
from importlib import metadata as _metadata
from pathlib import Path


_PACKAGES_TO_PROBE = [
    "pandas", "numpy", "biopython", "mhcflurry", "tensorflow",
    "scikit-learn", "networkx", "matplotlib", "seaborn", "matplotlib-venn",
    "openpyxl", "xlsxwriter", "rich", "jinja2", "requests", "joblib",
    "toxinpred3",
]


def _run_capture(command: list[str], timeout: int = 5) -> str:
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return ""
    return (completed.stdout or "") + (completed.stderr or "")


def _detect_wsl() -> bool:
    proc_version = Path("/proc/version")
    if not proc_version.exists():
        return False
    try:
        text = proc_version.read_text(errors="ignore").lower()
    except OSError:
        return False
    return "microsoft" in text or "wsl" in text


def _package_version(distribution_name: str) -> str:
    try:
        return _metadata.version(distribution_name)
    except _metadata.PackageNotFoundError:
        return "not_installed"


def collect_hardware_snapshot(repo_root: Path) -> dict:
    cpu_info_raw    = _run_capture(["lscpu"])
    mem_info_raw    = _run_capture(["free", "-h"])
    disk_info_raw   = _run_capture(["df", "-h", str(repo_root)])
    nvidia_smi_raw  = _run_capture(["nvidia-smi", "-L"])

    return {
        "captured_at":     _datetime.datetime.now().isoformat(timespec="seconds"),
        "host":            {
            "platform":        platform.platform(),
            "machine":         platform.machine(),
            "processor":       platform.processor(),
            "python_version":  sys.version.split()[0],
            "python_compiler": platform.python_compiler(),
            "is_wsl":          _detect_wsl(),
            "kernel":          _run_capture(["uname", "-a"]).strip(),
        },
        "cpu_raw":         cpu_info_raw,
        "memory_raw":      mem_info_raw,
        "disk_raw":        disk_info_raw,
        "packages":        {
            name: _package_version(name) for name in _PACKAGES_TO_PROBE
        },
        "tools": {
            "netmhcpan_predictor":    "IEDB classic API (tools-cluster-interface.iedb.org)",
            "mhcflurry_predictor":    "local (models in ~/.local/share/mhcflurry/)",
            "nvidia_smi":             nvidia_smi_raw.strip() or "no_gpu_detected",
        },
        "execution_policy": {
            "harness_parallelism": "serial",
            "pipeline_internal_parallelism": "preserved (ThreadPoolExecutor in predict_binding)",
            "notes": "All replicates and projects run one at a time. The pipeline keeps its production-equivalent internal threading.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", required=True, type=Path,
        help="Destination JSON path (parent directories are created).",
    )
    parser.add_argument(
        "--repo-root", type=Path, default=Path.cwd(),
        help="Path used as the working tree for disk-usage capture.",
    )
    args = parser.parse_args()

    snapshot = collect_hardware_snapshot(args.repo_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(snapshot, indent=2, ensure_ascii=False))
    print(f"Hardware snapshot written to {args.out}")


if __name__ == "__main__":
    main()
