"""
screen_toxicity step.

Takes the immunogenic epitopes produced by consensus_filter and drops the
toxic ones using ToxinPred3 (Model 1: AAC + DPC features fed into Extra
Trees, default cutoff 0.38).

Inputs  (track_dir/consensus/):
    CONSENSUS_IMMUNOGENIC_{track_id}.csv

Outputs (track_dir/toxicity/):
    TOXICITY_ALL_{track_id}.csv      every peptide with score, ppv, label
    TOXICITY_SAFE_{track_id}.csv     only the Non-Toxin rows
    TOXICITY_AUDIT_{track_id}.json   counts and parameters used in the run

Prediction is done in memory: the model .pkl is loaded once with joblib and
the AAC+DPC vectors are built directly, so we avoid temporary files and
subprocess calls.
"""

import importlib.util
import json
import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rich import box
from utils.console import console
from rich.panel import Panel
from rich.table import Table

import config
from modules.base_step import BaseTrackStep
from utils.naming import get_step_filename, COLUMN_PEPTIDE
from utils.project_manager import save_project_config

_STD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

LABEL_TOXIC = "Toxin"
LABEL_SAFE = "Non-Toxin"


# ToxinPred3 core (in-memory, no temporary files)

def _get_model_path() -> Path:
    spec = importlib.util.find_spec("toxinpred3")
    if spec is None:
        raise ImportError(
            "toxinpred3 package not found. Install it with: pip install toxinpred3>=1.4"
        )
    pkg_dir = Path(spec.submodule_search_locations[0])
    return pkg_dir / "model" / "toxinpred3.0_model.pkl"


def _aac_vector(seq: str) -> list:
    n = len(seq)
    return [(seq.count(aa) / n) * 100.0 for aa in _STD_AAS]


def _dpc_vector(seq: str) -> list:
    n = len(seq)
    denom = n - 1
    result = []
    for aa1 in _STD_AAS:
        for aa2 in _STD_AAS:
            count = sum(1 for i in range(n - 1) if seq[i] == aa1 and seq[i + 1] == aa2)
            result.append(count / denom * 100.0 if denom > 0 else 0.0)
    return result


def _predict(peptides: list, model_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Returns (scores, ppv): two 1-D arrays, one value per peptide."""
    clf = joblib.load(str(model_path))
    X = np.array([_aac_vector(p) + _dpc_vector(p) for p in peptides], dtype=np.float32)
    scores = clf.predict_proba(X)[:, -1]
    # PPV calibration reproduced from upstream ToxinPred3
    # (toxinpred3/python_scripts/toxinpred3.py:346, Model 1 = AAC+DPC, no MERCI).
    # The authors fitted the linear coefficients on their validation set so
    # that running the standalone tool gives the same PPV column we produce.
    ppv = np.clip(scores * 1.2341 - 0.1182, 0.0, 1.0)
    return scores, ppv


# Threshold configuration

def _ask_threshold(project_name: str, project_config: dict) -> float:
    """
    Reads the threshold from project_config or asks the user once.
    The chosen value is saved back to project_config so future reruns
    do not ask again.
    """
    existing = project_config.get("toxicity_threshold")
    if existing is not None:
        return float(existing)

    console.print(Panel(
        "[bold]ToxinPred3 toxicity threshold[/bold]\n\n"
        "[dim]Peptides with score >= threshold are flagged as toxic and removed.[/dim]\n\n"
        f"  [cyan][1][/cyan] 0.38  (ToxinPred3 default, recommended)\n"
        "  [cyan][2][/cyan] 0.50  (more restrictive)\n"
        "  [cyan][3][/cyan] Type a custom value",
        box=box.ROUNDED, title="Setup: screen_toxicity", title_align="left",
    ))

    while True:
        try:
            choice = input("> ").strip()
        except EOFError:
            choice = "1"

        if choice == "1" or choice == "":
            threshold = config.TOXICITY_SCORE_THRESHOLD
            break
        if choice == "2":
            threshold = 0.50
            break
        if choice == "3":
            while True:
                try:
                    raw = input("Threshold (e.g. 0.45): ").strip().replace(",", ".")
                except EOFError:
                    raw = str(config.TOXICITY_SCORE_THRESHOLD)
                try:
                    threshold = float(raw)
                    if 0.0 < threshold <= 1.0:
                        break
                    console.print("[red]Value must be between 0 and 1.[/red]")
                except ValueError:
                    console.print("[red]Invalid value. Use a dot as decimal separator.[/red]")
            break
        console.print("[dim]Invalid option. Type 1, 2 or 3.[/dim]")

    project_config["toxicity_threshold"] = threshold
    save_project_config(project_name, project_config)
    console.print(f"[dim]Threshold set to >= {threshold}, saved to project_config.[/dim]")
    return threshold


# Display

def _print_result(n_input: int, n_toxic: int, n_safe: int, threshold: float):
    t = Table(box=box.SIMPLE, show_header=False,
              title=f"Toxicity screen with ToxinPred3  (threshold={threshold:.2f})")
    t.add_column("", style="bold cyan", no_wrap=True)
    t.add_column("", justify="right")

    t.add_row("Input (immunogenic epitopes)", str(n_input))
    t.add_row("Toxic peptides removed",
              f"[bold red]-{n_toxic}[/bold red]" if n_toxic else "[dim]0[/dim]")
    t.add_row("Safe [bold green]✓[/bold green]",
              f"[bold green]{n_safe}[/bold green]")

    console.print(t)
    console.print(
        f"\n[bold green]RESULT: {n_safe} non-toxic epitopes carried forward to the next step.[/bold green]\n"
    )


# Step

class ScreenToxicityStep(BaseTrackStep):
    step_name = "screen_toxicity"

    def describe_outputs(self) -> dict:
        toxicity_dir = self.track_dir / "toxicity"
        return {
            toxicity_dir / get_step_filename("TOXICITY_SAFE", self.track_id):
                "Non-toxic peptides only — score below threshold. Feeds cluster_epitopes.",
            toxicity_dir / get_step_filename("TOXICITY_ALL", self.track_id):
                "All peptides screened, with toxicity score and is_toxic flag.",
            toxicity_dir / get_step_filename("TOXICITY_AUDIT", self.track_id, ext="json"):
                "Run audit — threshold used, counts of toxic/safe.",
        }

    def run(self, input_data=None):
        threshold = _ask_threshold(self.project_name, self.project_config)

        consensus_dir = self.track_dir / "consensus"
        input_csv = consensus_dir / get_step_filename("CONSENSUS_IMMUNOGENIC", self.track_id)
        if not input_csv.exists():
            raise FileNotFoundError(
                f"consensus_filter output not found: {input_csv}\n"
                "Run the 'consensus_filter' step before 'screen_toxicity'."
            )

        df = pd.read_csv(input_csv)
        if COLUMN_PEPTIDE not in df.columns:
            raise ValueError(
                f"Column '{COLUMN_PEPTIDE}' not found in {input_csv.name}. "
                f"Available columns: {list(df.columns)}"
            )

        # Defensive cleanup. The AAC and DPC vectors only know the 20
        # canonical amino acids, so any peptide that is empty, NaN, or has
        # ambiguous residues (X, B, Z) would either crash or be silently
        # mis-encoded. We drop those rows and report the count instead.
        canonical = set(_STD_AAS)
        df = df.dropna(subset=[COLUMN_PEPTIDE]).copy()
        valid_mask = df[COLUMN_PEPTIDE].apply(
            lambda p: isinstance(p, str) and len(p) > 0 and set(p) <= canonical
        )
        n_dropped = int((~valid_mask).sum())
        if n_dropped:
            console.print(
                f"[yellow]Skipped {n_dropped} peptide(s) with empty or non-canonical residues.[/yellow]"
            )
        df = df[valid_mask].reset_index(drop=True)

        peptides = df[COLUMN_PEPTIDE].tolist()
        n_input = len(peptides)

        console.print("\n[bold]Loading ToxinPred3 model...[/bold]")
        model_path = _get_model_path()
        scores, ppv = _predict(peptides, model_path)

        df["toxinpred3_score"] = np.round(scores, 3)
        df["toxinpred3_ppv"] = np.round(ppv, 3)
        df["toxinpred3_label"] = [LABEL_TOXIC if s >= threshold else LABEL_SAFE for s in scores]

        out = self.track_dir / "toxicity"
        out.mkdir(parents=True, exist_ok=True)

        all_csv = out / get_step_filename("TOXICITY_ALL", self.track_id)
        safe_csv = out / get_step_filename("TOXICITY_SAFE", self.track_id)

        df.to_csv(all_csv, index=False)

        df_safe = df[df["toxinpred3_label"] == LABEL_SAFE].copy()
        df_safe.to_csv(safe_csv, index=False)

        n_toxic = int((df["toxinpred3_label"] == LABEL_TOXIC).sum())
        n_safe = len(df_safe)

        _print_result(n_input, n_toxic, n_safe, threshold)

        audit = {
            "timestamp":          datetime.datetime.now().isoformat(),
            "track_id":           self.track_id,
            "tool":               "ToxinPred3",
            "model":              "AAC+DPC Extra Trees (Model 1)",
            "threshold":          threshold,
            "n_input":            n_input,
            "n_dropped_invalid":  n_dropped,
            "n_toxic":            n_toxic,
            "n_safe":             n_safe,
            "pct_removed":        round(n_toxic / n_input * 100, 1) if n_input else 0,
            "output_all":         str(all_csv),
            "output_safe":        str(safe_csv),
        }
        audit_path = out / get_step_filename("TOXICITY_AUDIT", self.track_id, ext="json")
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False))

        return {
            "output_all":  str(all_csv),
            "output_safe": str(safe_csv),
            "n_input":     n_input,
            "n_toxic":     n_toxic,
            "n_safe":      n_safe,
        }
