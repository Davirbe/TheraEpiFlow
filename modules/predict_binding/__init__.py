"""
Step — Predict Binding

Runs MHC-I binding predictions for all sequences in a track using two tools:
  1. NetMHCpan 4.1 EL — via IEDB classic API (http://tools-cluster-interface.iedb.org/tools_api/mhci/)
  2. MHCFlurry 2.0    — via Python API (mhcflurry.Class1PresentationPredictor)

External dependencies:
  - IEDB MHC-I prediction endpoint (no local install required).
  - `mhcflurry >= 2.0` (PyPI). Requires a one-time
    `mhcflurry-downloads fetch` to install the trained models.

Citations:
  Reynisson B et al. NetMHCpan-4.1 and NetMHCIIpan-4.0. NAR. 2020;48:W449–W454.
  O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0. Cell Systems. 2020;11(1):42–48.e7.

Inputs:
  data/input/{track_id}/SEQUENCES_{track_id}.fasta  (from fetch_sequences)
  Fallback: user provides a local FASTA path, or pastes sequences in the terminal.

Parameters (asked once, saved to project_config):
  hla_alleles:     list of alleles in IMGT format, e.g. ["HLA-A*02:01", "HLA-B*07:02"]
  peptide_lengths: list of ints, e.g. [9] or [8, 9, 10]

Outputs (saved to data/intermediate/{track_id}/predictions/):
  PRED_NET_{track_id}.csv      — NetMHCpan 4.1 EL predictions
  PRED_FLURRY_{track_id}.csv   — MHCFlurry 2.0 predictions
  PREDICT_AUDIT_{track_id}.json

Column contract with consensus_filter:
  NET:   column containing ['netmhcpan_el', 'percentile'] → netmhcpan_el_percentile
  FLURRY: column containing ['mhcflurry', 'presentation', 'percentile'] → mhcflurry_presentation_percentile
  Both:  'peptide' and 'allele' columns required

Terminology note:
  "Strong binder"        = NetMHCpan EL %rank ≤ 0.5
  "Intermediate binder"  = 0.5 < NetMHCpan EL %rank ≤ 2.0  (a.k.a. "weak" in the NetMHC literature)
  This module shows "intermediate" to users to avoid the misleading "weak" label
  for peptides that frequently elicit real T-cell responses.
"""

# ── Silence TF / Keras / MHCFlurry noise BEFORE the libraries are imported ───
# These environment variables only take effect if set before tensorflow's first
# import. We set them at module load even though `tensorflow` and `mhcflurry`
# are imported lazily inside `_run_mhcflurry` — by the time a user invokes the
# step the env vars are already in place.
import logging as _logging_module
import os as _os_module
import warnings as _warnings_module

_os_module.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
_os_module.environ.setdefault("TF_ENABLE_ONEDNN_OPTS", "0")
_os_module.environ.setdefault("CUDA_VISIBLE_DEVICES", _os_module.environ.get("CUDA_VISIBLE_DEVICES", ""))
_warnings_module.filterwarnings("ignore", category=DeprecationWarning)
_warnings_module.filterwarnings("ignore", category=FutureWarning)
# UserWarning is broadened to catch noisy warnings from MHCFlurry's import
# chain (e.g. `pkg_resources is deprecated` from setuptools). Narrowing it to
# `module="tensorflow"` previously left the pkg_resources warning bleeding
# through, which collided with the Progress spinner and looked like a freeze.
_warnings_module.filterwarnings("ignore", category=UserWarning)
_logging_module.getLogger("tensorflow").setLevel(_logging_module.ERROR)
_logging_module.getLogger("mhcflurry").setLevel(_logging_module.ERROR)
_logging_module.getLogger("absl").setLevel(_logging_module.ERROR)
_logging_module.getLogger("h5py").setLevel(_logging_module.ERROR)


import datetime
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from rich import box
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
)
from rich.prompt import Prompt
from rich.text import Text

from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.fasta_utils import generate_peptides, is_valid_sequence
from utils.naming import get_prediction_filename, get_step_filename
from utils.output_capture import capture_fd_output
from utils.retry_helpers import retry_network_call
from utils.step_summary import print_step_summary

IEDB_MHCI_URL = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'

# Cutoffs used only for the narrative summary (the actual cutoffs that drive
# the consensus filter live in config.py and are not duplicated here).
SUMMARY_STRONG_BINDER_RANK_MAX = 0.5
SUMMARY_INTERMEDIATE_BINDER_RANK_MAX = 2.0


# ── Parameter setup ───────────────────────────────────────────────────────────

def _parse_peptide_lengths(raw_input: str, default_lengths: list) -> list:
    """
    Parses peptide length input supporting comma-separated and range notation.

    Examples:
        "9"      → [9]
        "9,10,11"→ [9, 10, 11]
        "9-11"   → [9, 10, 11]
    """
    raw_input = raw_input.strip()
    if not raw_input:
        return list(default_lengths)

    if '-' in raw_input and ',' not in raw_input:
        range_parts = raw_input.split('-')
        if len(range_parts) == 2:
            try:
                start_length = int(range_parts[0].strip())
                end_length   = int(range_parts[1].strip())
                if 1 <= start_length <= end_length <= 20:
                    return list(range(start_length, end_length + 1))
            except ValueError:
                pass

    try:
        parsed = [int(token.strip()) for token in raw_input.split(',') if token.strip()]
        return parsed if parsed else list(default_lengths)
    except ValueError:
        return list(default_lengths)


def _ask_binding_params(project_name: str, project_config: dict) -> tuple[list[str], list[int]]:
    """
    Returns HLA alleles and peptide lengths from project_config if already saved,
    otherwise asks the user once and saves them.

    Default alleles: 27 standard MHC-I alleles from config.DEFAULT_HLA_ALLELES.
    Default peptide length: 9 (configurable 8–12 for MHC-I).
    """
    from config import DEFAULT_HLA_ALLELES, DEFAULT_PEPTIDE_LENGTHS

    if 'hla_alleles' in project_config and project_config['hla_alleles'] \
            and 'peptide_lengths' in project_config and project_config['peptide_lengths']:
        hla_alleles     = project_config['hla_alleles']
        peptide_lengths = project_config['peptide_lengths']
        console.print(f"[dim]  Alleles (saved): {len(hla_alleles)} alleles[/dim]")
        console.print(f"[dim]  Peptide lengths (saved): {peptide_lengths}[/dim]")
        return hla_alleles, peptide_lengths

    console.print("\n[bold]Binding prediction setup[/bold]")

    # ── HLA alleles ───────────────────────────────────────────────────────────
    console.print(
        f"\n[dim]Default: {len(DEFAULT_HLA_ALLELES)} standard MHC-I alleles "
        f"(A*01:01, A*02:01 ... B*57:01, B*58:01)[/dim]"
    )
    try:
        use_defaults = Prompt.ask("Use default 27 alleles?", choices=['y', 'n'], default='y')
    except Exception:
        use_defaults = 'y'

    if use_defaults.lower() == 'y':
        hla_alleles = list(DEFAULT_HLA_ALLELES)
        console.print(f"[dim]  → {len(hla_alleles)} default alleles selected.[/dim]")
    else:
        console.print("[dim]Enter HLA alleles in IMGT format, comma-separated.[/dim]")
        console.print("[dim]Example: HLA-A*02:01,HLA-B*07:02[/dim]")
        try:
            allele_input = Prompt.ask("HLA alleles").strip()
        except EOFError:
            allele_input = ''
        hla_alleles = [a.strip() for a in allele_input.split(',') if a.strip()]
        if not hla_alleles:
            console.print("[yellow]  No alleles entered — using defaults.[/yellow]")
            hla_alleles = list(DEFAULT_HLA_ALLELES)

    # ── Peptide lengths ───────────────────────────────────────────────────────
    console.print(
        "\n[dim]Peptide lengths to predict (MHC-I binding groove: 8–12 aa).[/dim]\n"
        "[dim]Examples:  9  |  9,10,11  |  9-11[/dim]"
    )
    try:
        length_input = Prompt.ask("Peptide lengths", default="9").strip()
    except EOFError:
        length_input = '9'

    peptide_lengths = _parse_peptide_lengths(length_input, DEFAULT_PEPTIDE_LENGTHS)

    project_config['hla_alleles']     = hla_alleles
    project_config['peptide_lengths'] = peptide_lengths

    from utils.project_manager import save_project_config
    save_project_config(project_name, project_config)

    return hla_alleles, peptide_lengths


# ── Sequence loading ──────────────────────────────────────────────────────────

def _load_sequences(track_id: str, input_dir: Path) -> list:
    """
    Loads sequences from the canonical FASTA produced by fetch_sequences.
    Falls back to a user-provided file path or terminal paste if not found.
    Returns a list of validated SeqRecord objects.
    """
    canonical_fasta_path = input_dir / track_id / get_step_filename("SEQUENCES", track_id, ext="fasta")

    if canonical_fasta_path.exists():
        records = list(SeqIO.parse(str(canonical_fasta_path), 'fasta'))
        if records:
            console.print(f"[dim]  Loaded {len(records)} sequence(s) from {canonical_fasta_path.name}[/dim]")
            return records

    console.print(f"[yellow]  Sequence file not found: {canonical_fasta_path}[/yellow]")
    console.print("  Provide sequences via file path or paste FASTA content directly.")
    console.print("[dim]  (Enter a file path, or press Enter to paste sequences)[/dim]")

    user_path_input = Prompt.ask("  File path (or Enter to paste)").strip()

    if user_path_input:
        provided_path = Path(user_path_input)
        if not provided_path.exists():
            raise FileNotFoundError(f"File not found: {provided_path}")
        raw_records = list(SeqIO.parse(str(provided_path), 'fasta'))
    else:
        console.print("[dim]  Paste FASTA content below. Enter a blank line when done.[/dim]")
        pasted_lines = []
        while True:
            line = input()
            if line == '':
                break
            pasted_lines.append(line)
        raw_records = list(SeqIO.parse(io.StringIO('\n'.join(pasted_lines)), 'fasta'))

    if not raw_records:
        raise ValueError("No sequences found in the provided input.")

    valid_records = []
    for record in raw_records:
        sequence_is_valid, rejection_reason = is_valid_sequence(record)
        if sequence_is_valid:
            valid_records.append(record)
        else:
            console.print(f"[yellow]  Skipping {record.id}: {rejection_reason}[/yellow]")

    if not valid_records:
        raise ValueError("No valid sequences after validation.")

    console.print(f"[dim]  {len(valid_records)} valid sequence(s) accepted.[/dim]")
    return valid_records


# ── NetMHCpan 4.1 EL via IEDB classic API ────────────────────────────────────

def _run_netmhcpan_iedb_silent(
    sequence_records: list,
    hla_alleles:      list[str],
    peptide_lengths:  list[int],
) -> pd.DataFrame:
    """
    Pure HTTP worker — emits NO console output. Designed to run in a worker
    thread while the main thread drives the Rich Progress display. All
    user-facing messaging stays in the main thread.

    Returns DataFrame with columns:
      allele, peptide, netmhcpan_el_score, netmhcpan_el_percentile,
      and any additional columns returned by the API (seq_num, start, end, etc.)
    """
    fasta_text = ''.join(
        f'>{record.id}\n{str(record.seq).replace("*", "")}\n'
        for record in sequence_records
    )

    # The IEDB API requires an equal number of allele and length entries.
    # We expand to the cartesian product so each allele maps to each length.
    allele_length_pairs = [
        (allele, length)
        for allele in hla_alleles
        for length in peptide_lengths
    ]

    request_payload = {
        'method':        'netmhcpan_el',
        'sequence_text': fasta_text,
        'allele':        ','.join(allele for allele, _ in allele_length_pairs),
        'length':        ','.join(str(length) for _, length in allele_length_pairs),
    }

    def _post_to_iedb():
        try:
            return requests.post(IEDB_MHCI_URL, data=request_payload, timeout=300)
        except requests.exceptions.RequestException as request_error:
            raise ConnectionError(str(request_error)) from request_error

    api_response = retry_network_call(_post_to_iedb, 'IEDB NetMHCpan 4.1 EL')

    if api_response.status_code != 200:
        raise RuntimeError(
            f"IEDB API returned HTTP {api_response.status_code}: {api_response.text[:300]}"
        )

    non_comment_lines = [
        line for line in api_response.text.splitlines()
        if not line.startswith('#')
    ]
    net_dataframe = pd.read_csv(io.StringIO('\n'.join(non_comment_lines)), sep='\t')
    net_dataframe.columns = net_dataframe.columns.str.strip()

    net_dataframe = net_dataframe.rename(columns={
        'score':           'netmhcpan_el_score',
        'percentile_rank': 'netmhcpan_el_percentile',
    })

    columns_to_keep = [
        col for col in [
            'allele', 'peptide',
            'netmhcpan_el_score', 'netmhcpan_el_percentile',
            'seq_num', 'start', 'end', 'length', 'core', 'icore',
        ]
        if col in net_dataframe.columns
    ]
    return net_dataframe[columns_to_keep].copy()


# ── MHCFlurry 2.0 (per-allele, with output capture per call) ─────────────────

def _run_mhcflurry_with_progress(
    sequence_records: list,
    hla_alleles:      list[str],
    peptide_lengths:  list[int],
    progress_handle:  Progress,
    progress_task_id: int,
) -> tuple[pd.DataFrame, str]:
    """
    Loads MHCFlurry models and predicts presentation for every (peptide, allele)
    pair. Each model load and predict call is wrapped in `capture_fd_output` so
    TensorFlow / Keras INFO logs do not pollute the terminal. The captured text
    is concatenated and returned alongside the dataframe so the caller can show
    it if something goes wrong.

    Updates `progress_handle[progress_task_id]` after every allele finishes.
    """
    from mhcflurry import Class1PresentationPredictor

    captured_chunks: list[str] = []

    progress_handle.update(
        progress_task_id,
        description="loading TensorFlow models — please don't press keys",
    )

    def _load_mhcflurry_models():
        with capture_fd_output() as load_capture:
            loaded_predictor = Class1PresentationPredictor.load()
        captured_chunks.append(load_capture.stderr)
        captured_chunks.append(load_capture.stdout)
        return loaded_predictor

    predictor = retry_network_call(
        _load_mhcflurry_models,
        operation_description="MHCFlurry model load",
        max_attempts=2,
        initial_backoff_seconds=2.0,
    )

    all_unique_kmers = sorted(set(
        peptide
        for record in sequence_records
        for peptide in generate_peptides(str(record.seq).replace('*', ''), peptide_lengths)
    ))

    progress_handle.update(
        progress_task_id,
        description=f"predicting {len(all_unique_kmers)} k-mers",
        total=len(hla_alleles),
    )

    per_allele_dataframes: list[pd.DataFrame] = []
    for allele in hla_alleles:
        with capture_fd_output() as allele_capture:
            allele_batch = predictor.predict(
                peptides=all_unique_kmers,
                alleles=[allele],
                verbose=0,
            )
        captured_chunks.append(allele_capture.stderr)
        captured_chunks.append(allele_capture.stdout)

        allele_batch = allele_batch[['peptide', 'presentation_percentile']].copy()
        allele_batch['allele'] = allele
        allele_batch = allele_batch.rename(columns={
            'presentation_percentile': 'mhcflurry_presentation_percentile',
        })
        per_allele_dataframes.append(allele_batch)

        progress_handle.update(progress_task_id, advance=1)

    flurry_dataframe = pd.concat(per_allele_dataframes, ignore_index=True)
    flurry_dataframe = flurry_dataframe[['peptide', 'allele', 'mhcflurry_presentation_percentile']]

    captured_text = "".join(chunk for chunk in captured_chunks if chunk)
    return flurry_dataframe, captured_text


# ── Step class ────────────────────────────────────────────────────────────────

def _count_strong_and_intermediate_binders(net_dataframe: pd.DataFrame) -> tuple[int, int]:
    """Return (strong_binder_count, intermediate_binder_count) by EL %rank."""
    if net_dataframe is None or 'netmhcpan_el_percentile' not in net_dataframe.columns:
        return (0, 0)
    rank_series = pd.to_numeric(net_dataframe['netmhcpan_el_percentile'], errors='coerce')
    strong_count = int((rank_series <= SUMMARY_STRONG_BINDER_RANK_MAX).sum())
    intermediate_count = int(
        ((rank_series > SUMMARY_STRONG_BINDER_RANK_MAX) & (rank_series <= SUMMARY_INTERMEDIATE_BINDER_RANK_MAX)).sum()
    )
    return strong_count, intermediate_count


def _count_presentation_ready(flurry_dataframe: pd.DataFrame) -> int:
    """Count peptides MHCFlurry considered presentation-ready (percentile ≤ 2)."""
    if flurry_dataframe is None or 'mhcflurry_presentation_percentile' not in flurry_dataframe.columns:
        return 0
    percentile_series = pd.to_numeric(flurry_dataframe['mhcflurry_presentation_percentile'], errors='coerce')
    return int((percentile_series <= 2.0).sum())


def _peptides_in_common(net_dataframe: pd.DataFrame, flurry_dataframe: pd.DataFrame) -> int:
    """Return the number of unique peptides reported by both tools."""
    if net_dataframe is None or flurry_dataframe is None:
        return 0
    return len(set(net_dataframe['peptide']).intersection(set(flurry_dataframe['peptide'])))


class PredictBindingStep(BaseTrackStep):
    step_name   = 'predict_binding'
    description = (
        "Runs NetMHCpan 4.1 EL (via the IEDB HTTP API) and MHCFlurry 2.0 "
        "(local Python package) in parallel across every peptide × HLA "
        "combination, producing two separate prediction tables."
    )
    long_description = (
        "Generates every possible peptide of the configured lengths "
        "(default: 9-mers) from the reference protein, then asks two "
        "independent predictors to score each peptide × HLA pair:\n\n"
        "  • [bold]NetMHCpan-4.1 EL[/bold] (Eluted Ligand mode) — pan-allele "
        "neural network trained on mass-spec ligandome + binding affinity.\n"
        "  • [bold]MHCFlurry 2.0[/bold] (presentation predictor) — combines "
        "affinity + antigen-processing into a single presentation score.\n\n"
        "Both tools output a percentile rank per (peptide, allele); both are "
        "needed because the [bold]consensus_filter[/bold] step only keeps "
        "peptides flagged as binders by [italic]both[/italic]."
    )
    methodology = (
        "1. Reads the reference FASTA produced by fetch_sequences.\n"
        "2. Enumerates k-mers of each configured length, deduplicated.\n"
        "3. NetMHCpan: submits the FASTA + alleles to the IEDB classic API "
        "endpoint (mhci tools_api) and parses the EL %rank column.\n"
        "4. MHCFlurry: loads Class1PresentationPredictor (TF models) locally "
        "and calls predict() per allele, capturing TF logs to avoid spinner "
        "collisions.\n"
        "5. Both runs happen in parallel (ThreadPoolExecutor + main thread).\n"
        "6. Output tables share the column contract used by consensus_filter."
    )
    references = [
        {
            'authors': 'Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M.',
            'title':   'NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data',
            'journal': 'Nucleic Acids Research',
            'year':    2020,
            'doi':     '10.1093/nar/gkaa379',
        },
        {
            'authors': "O'Donnell TJ, Rubinsteyn A, Laserson U.",
            'title':   'MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing',
            'journal': 'Cell Systems',
            'year':    2020,
            'doi':     '10.1016/j.cels.2020.06.010',
        },
    ]
    data_format = (
        "Input is automatic — uses the FASTA from fetch_sequences. "
        "You will be asked once (at the first step run) for:\n"
        "  • [bold]HLA-I alleles[/bold] — IMGT format (e.g. HLA-A*02:01); "
        "default is the 27-allele panel covering global diversity.\n"
        "  • [bold]Peptide lengths[/bold] — typically 9 for MHC-I (8-10 are "
        "biologically valid; we default to 9 because it dominates the eluted "
        "ligandome literature)."
    )
    outputs_overview = (
        "[bold]PRED_NET_{track_id}.csv[/bold]      — NetMHCpan EL predictions (peptide, allele, percentile).\n"
        "[bold]PRED_FLURRY_{track_id}.csv[/bold]   — MHCFlurry presentation predictions (peptide, allele, percentile).\n"
        "[bold]PREDICT_AUDIT_{track_id}.json[/bold] — run metadata + timing per tool."
    )
    tips = [
        "First MHCFlurry call loads TF models — may take 20-40s on the first run; later steps reuse the cache.",
        "More alleles → linearly more API calls to IEDB. Keep the default 27-allele panel unless you have a reason to extend.",
        "Peptide length 9 is the gold standard for MHC-I; using 8-10 is valid but increases volume ~3×.",
        "If IEDB returns a transient error, the step retries automatically up to 5 times with backoff.",
    ]

    def run(self, input_data=None):
        hla_alleles, peptide_lengths = _ask_binding_params(
            self.project_name, self.project_config
        )
        sequence_records = _load_sequences(self.track_id, self.input_dir)

        predictions_dir = self.track_dir / 'predictions'
        predictions_dir.mkdir(parents=True, exist_ok=True)

        net_output_path    = predictions_dir / get_prediction_filename("NET_PRED",    self.track_id)
        flurry_output_path = predictions_dir / get_prediction_filename("FLURRY_PRED", self.track_id)
        audit_output_path  = predictions_dir / get_step_filename("PREDICT_AUDIT", self.track_id, ext="json")

        console.print(Panel(
            Text.from_markup(
                f"[bold]Track:[/bold] {self.track_id}\n"
                f"[bold]Sequences:[/bold] {len(sequence_records)}\n"
                f"[bold]Alleles:[/bold] {len(hla_alleles)}  "
                f"[dim]({', '.join(hla_alleles[:3])}{'…' if len(hla_alleles) > 3 else ''})[/dim]\n"
                f"[bold]Peptide lengths:[/bold] {peptide_lengths}"
            ),
            title="Binding prediction", border_style="cyan", box=box.ROUNDED,
        ))

        console.print()
        console.print(Panel(
            Text.from_markup(
                "[bold yellow]⚠ MHCFlurry first-time load: ~30-60 seconds.[/bold yellow]\n"
                "TensorFlow models are being read from disk and warmed up. The bar may "
                "briefly pause during native imports — that's expected.\n"
                "[dim]Please don't press any keys until the load finishes.[/dim]"
            ),
            border_style="yellow", box=box.ROUNDED, padding=(0, 1),
        ))

        run_start_time = time.time()

        net_dataframe: pd.DataFrame | None = None
        flurry_dataframe: pd.DataFrame | None = None
        net_error: str | None = None
        flurry_error: str | None = None
        flurry_captured_log: str = ""

        # Both tools run in parallel: NetMHCpan in a worker thread (HTTP, silent),
        # MHCFlurry in the main thread (loads TF models, captured per allele).
        with Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.fields[label]:<10}[/bold blue]"),
            BarColumn(
                bar_width=30,
                style="yellow",
                complete_style="green",
                finished_style="green",
                pulse_style="bold yellow",
            ),
            MofNCompleteColumn(),
            TextColumn("{task.description}"),
            console=console,
            transient=False,
        ) as progress_bar:
            netmhcpan_task_id = progress_bar.add_task(
                "sending request to IEDB",
                total=1,
                label="NetMHCpan",
            )
            mhcflurry_task_id = progress_bar.add_task(
                "preparing",
                total=None,
                label="MHCFlurry",
            )

            with ThreadPoolExecutor(max_workers=1) as netmhcpan_executor:
                netmhcpan_future = netmhcpan_executor.submit(
                    _run_netmhcpan_iedb_silent,
                    sequence_records, hla_alleles, peptide_lengths,
                )

                # Run MHCFlurry serially in the main thread so the progress
                # bar can be updated after each allele completes.
                try:
                    flurry_dataframe, flurry_captured_log = _run_mhcflurry_with_progress(
                        sequence_records=sequence_records,
                        hla_alleles=hla_alleles,
                        peptide_lengths=peptide_lengths,
                        progress_handle=progress_bar,
                        progress_task_id=mhcflurry_task_id,
                    )
                    progress_bar.update(
                        mhcflurry_task_id,
                        description=f"[green]done — {len(flurry_dataframe):,} rows[/green]",
                    )
                except Exception as caught_flurry_error:
                    flurry_error = str(caught_flurry_error)
                    progress_bar.update(
                        mhcflurry_task_id,
                        description=f"[red]failed: {caught_flurry_error}[/red]",
                    )

                # Wait for NetMHCpan to finish (it has been running concurrently).
                try:
                    net_dataframe = netmhcpan_future.result()
                    progress_bar.update(
                        netmhcpan_task_id,
                        completed=1,
                        description=f"[green]done — {len(net_dataframe):,} rows[/green]",
                    )
                except Exception as caught_netmhcpan_error:
                    net_error = str(caught_netmhcpan_error)
                    progress_bar.update(
                        netmhcpan_task_id,
                        description=f"[red]failed: {caught_netmhcpan_error}[/red]",
                    )

        elapsed_seconds = time.time() - run_start_time

        flush_stdin()

        if net_error or flurry_error:
            error_lines: list[str] = []
            if net_error:
                error_lines.append(f"NetMHCpan: {net_error}")
            if flurry_error:
                error_lines.append(f"MHCFlurry: {flurry_error}")
            if flurry_captured_log.strip():
                console.print(Panel(
                    Text.from_ansi(flurry_captured_log[-4000:]),
                    title="MHCFlurry stderr/stdout (captured)",
                    border_style="red", box=box.ROUNDED,
                ))
            raise RuntimeError(
                "Prediction failed (both outputs are required by consensus_filter):\n  "
                + "\n  ".join(error_lines)
            )

        if net_dataframe is not None:
            if 'netmhcpan_el_percentile' in net_dataframe.columns:
                net_dataframe['netmhcpan_el_percentile'] = pd.to_numeric(
                    net_dataframe['netmhcpan_el_percentile'], errors='coerce'
                ).round(2)
            if 'netmhcpan_el_score' in net_dataframe.columns:
                net_dataframe['netmhcpan_el_score'] = pd.to_numeric(
                    net_dataframe['netmhcpan_el_score'], errors='coerce'
                ).round(4)
            net_dataframe.to_csv(net_output_path, index=False)
        if flurry_dataframe is not None:
            if 'mhcflurry_presentation_percentile' in flurry_dataframe.columns:
                flurry_dataframe['mhcflurry_presentation_percentile'] = pd.to_numeric(
                    flurry_dataframe['mhcflurry_presentation_percentile'], errors='coerce'
                ).round(2)
            flurry_dataframe.to_csv(flurry_output_path, index=False)

        audit_data = {
            'track_id':        self.track_id,
            'timestamp':       datetime.datetime.now().isoformat(),
            'hla_alleles':     hla_alleles,
            'peptide_lengths': peptide_lengths,
            'sequence_count':  len(sequence_records),
            'elapsed_seconds': round(elapsed_seconds, 1),
            'netmhcpan': {
                'status':    'done'  if net_dataframe    is not None else 'error',
                'row_count': len(net_dataframe)    if net_dataframe    is not None else 0,
                'output':    str(net_output_path)  if net_dataframe    is not None else None,
                'error':     net_error,
            },
            'mhcflurry': {
                'status':    'done'  if flurry_dataframe is not None else 'error',
                'row_count': len(flurry_dataframe) if flurry_dataframe is not None else 0,
                'output':    str(flurry_output_path) if flurry_dataframe is not None else None,
                'error':     flurry_error,
            },
        }
        with open(audit_output_path, 'w') as audit_file:
            json.dump(audit_data, audit_file, indent=2)

        # Narrative summary — strong/intermediate (never "weak") + plain prose
        strong_binder_count, intermediate_binder_count = _count_strong_and_intermediate_binders(net_dataframe)
        presentation_ready_count = _count_presentation_ready(flurry_dataframe)
        peptides_in_common_count = _peptides_in_common(net_dataframe, flurry_dataframe)

        narrative_lines: list[str] = [
            f"[bold]NetMHCpan 4.1 EL[/bold] — {len(net_dataframe):,} rows",
            f"  {strong_binder_count:,} strong binders (rank ≤ {SUMMARY_STRONG_BINDER_RANK_MAX}%) + "
            f"{intermediate_binder_count:,} intermediate "
            f"({SUMMARY_STRONG_BINDER_RANK_MAX}% < rank ≤ {SUMMARY_INTERMEDIATE_BINDER_RANK_MAX}%)",
            f"[bold]MHCFlurry 2.0[/bold] — {len(flurry_dataframe):,} rows",
            f"  {presentation_ready_count:,} presentation-ready peptides (percentile ≤ 2)",
            f"[bold]Peptides reported by both tools:[/bold] {peptides_in_common_count:,} "
            f"[dim](these continue into consensus_filter)[/dim]",
        ]

        print_step_summary(
            step_title=f"Binding prediction complete for {self.track_id}",
            elapsed_seconds=elapsed_seconds,
            narrative_lines=narrative_lines,
            output_files=[net_output_path, flurry_output_path, audit_output_path],
        )

        return {
            'net_csv':    str(net_output_path)    if net_dataframe    is not None else None,
            'flurry_csv': str(flurry_output_path) if flurry_dataframe is not None else None,
            'audit_json': str(audit_output_path),
        }

    def describe_outputs(self) -> dict[Path, str]:
        predictions_dir = self.track_dir / 'predictions'
        return {
            predictions_dir / get_prediction_filename("NET_PRED", self.track_id):
                "NetMHCpan 4.1 EL predictions — one row per (peptide, allele) with EL score and %rank.",
            predictions_dir / get_prediction_filename("FLURRY_PRED", self.track_id):
                "MHCFlurry 2.0 predictions — one row per (peptide, allele) with presentation percentile.",
            predictions_dir / get_step_filename("PREDICT_AUDIT", self.track_id, ext="json"):
                "Run audit — alleles used, lengths, row counts per tool, elapsed time, and any errors.",
        }
