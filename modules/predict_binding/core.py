"""
Prediction engines for predict_binding: the NetMHCpan 4.1 EL (IEDB HTTP API)
and MHCFlurry 2.0 (local) runners, plus the binder-count summaries. The TF /
MHCFlurry warning suppression lives here because the lazy tensorflow import
happens inside the MHCFlurry runner.
"""

# Silence TF/Keras/MHCFlurry noise — env vars MUST be set before tensorflow's
# first import (which happens lazily inside _run_mhcflurry_with_progress).
import logging as _logging_module
import os as _os_module
import warnings as _warnings_module

_os_module.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
_os_module.environ.setdefault("TF_ENABLE_ONEDNN_OPTS", "0")
_os_module.environ.setdefault("CUDA_VISIBLE_DEVICES", _os_module.environ.get("CUDA_VISIBLE_DEVICES", ""))
_warnings_module.filterwarnings("ignore", category=DeprecationWarning)
_warnings_module.filterwarnings("ignore", category=FutureWarning)
# UserWarning broadened: MHCFlurry's import chain raises pkg_resources warnings that collide with the Progress spinner.
_warnings_module.filterwarnings("ignore", category=UserWarning)
_logging_module.getLogger("tensorflow").setLevel(_logging_module.ERROR)
_logging_module.getLogger("mhcflurry").setLevel(_logging_module.ERROR)
_logging_module.getLogger("absl").setLevel(_logging_module.ERROR)
_logging_module.getLogger("h5py").setLevel(_logging_module.ERROR)


import io

import pandas as pd
import requests
from rich.progress import Progress

from utils.fasta_utils import generate_peptides
from utils.output_capture import capture_fd_output
from utils.retry_helpers import retry_network_call

IEDB_MHCI_URL = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'

# Cutoffs used only for the narrative summary (the actual cutoffs that drive
# the consensus filter live in config.py and are not duplicated here).
SUMMARY_STRONG_BINDER_RANK_MAX = 0.5
SUMMARY_INTERMEDIATE_BINDER_RANK_MAX = 2.0
# ── NetMHCpan 4.1 EL via IEDB classic API ────────────────────────────────────

def _run_netmhcpan_iedb_silent(
    sequence_records: list,
    hla_alleles:      list[str],
    peptide_lengths:  list[int],
) -> pd.DataFrame:
    """Pure HTTP worker — emits no console output.
    Designed to run in a worker thread while the main thread drives the Progress display.

    Returns DataFrame with columns: allele, peptide, netmhcpan_el_score,
    netmhcpan_el_percentile, plus passthroughs from the API (seq_num, start, end, etc.)."""
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
    """Loads MHCFlurry models and predicts presentation for every (peptide, allele) pair.
    Each load/predict is wrapped in capture_fd_output to silence TF/Keras INFO logs;
    captured text is returned alongside the dataframe for error reporting.
    Updates progress_handle[progress_task_id] after each allele."""
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


