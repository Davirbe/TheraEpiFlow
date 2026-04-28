"""
Step — Predict Binding

Runs MHC-I binding predictions for all sequences in a track using two tools in parallel:
  1. NetMHCpan 4.1 EL — via IEDB classic API (http://tools-cluster-interface.iedb.org/tools_api/mhci/)
  2. MHCFlurry 2.0    — via Python API (mhcflurry.Class1PresentationPredictor)

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
"""

import datetime
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from rich.console import Console
from rich.prompt import Prompt

from modules.base_step import BaseTrackStep
from utils.fasta_utils import generate_peptides, is_valid_sequence
from utils.naming import get_prediction_filename, get_step_filename
from utils.retry_helpers import retry_network_call

console = Console()

IEDB_MHCI_URL = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'


# ── Parameter setup ───────────────────────────────────────────────────────────

def _ask_binding_params(project_name: str, project_config: dict) -> tuple[list[str], list[int]]:
    """
    Returns HLA alleles and peptide lengths from project_config if already saved,
    otherwise asks the user once and saves them.
    """
    if 'hla_alleles' in project_config and 'peptide_lengths' in project_config:
        hla_alleles     = project_config['hla_alleles']
        peptide_lengths = project_config['peptide_lengths']
        console.print(f"[dim]  Alleles (saved): {', '.join(hla_alleles)}[/dim]")
        console.print(f"[dim]  Peptide lengths (saved): {peptide_lengths}[/dim]")
        return hla_alleles, peptide_lengths

    console.print("\n[bold]Binding prediction setup[/bold]")
    console.print("[dim]Enter HLA alleles in IMGT format, comma-separated.[/dim]")
    console.print("[dim]Example: HLA-A*02:01,HLA-B*07:02,HLA-C*07:02[/dim]")

    allele_input = Prompt.ask("HLA alleles").strip()
    hla_alleles  = [allele.strip() for allele in allele_input.split(',') if allele.strip()]
    if not hla_alleles:
        raise ValueError("At least one HLA allele is required.")

    console.print("[dim]Peptide lengths to evaluate, comma-separated (default: 9).[/dim]")
    length_input    = Prompt.ask("Peptide lengths", default="9").strip()
    peptide_lengths = [int(token.strip()) for token in length_input.split(',') if token.strip()]
    if not peptide_lengths:
        peptide_lengths = [9]

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

def _run_netmhcpan_iedb(
    sequence_records: list,
    hla_alleles:      list[str],
    peptide_lengths:  list[int],
) -> pd.DataFrame:
    """
    Calls the IEDB classic API with method=netmhcpan_el.
    The API receives full FASTA text and generates k-mers internally.

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
            return requests.post(IEDB_MHCI_URL, data=request_payload, timeout=120)
        except requests.exceptions.RequestException as request_error:
            raise ConnectionError(str(request_error)) from request_error

    console.print("[dim]  NetMHCpan: sending request to IEDB API...[/dim]")
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


# ── MHCFlurry 2.0 ─────────────────────────────────────────────────────────────

def _run_mhcflurry(
    sequence_records: list,
    hla_alleles:      list[str],
    peptide_lengths:  list[int],
) -> pd.DataFrame:
    """
    Runs MHCFlurry 2.0 presentation predictions via the Python API.

    Calls predict() once per allele with all k-mers, so each output row
    corresponds to exactly one (peptide, allele) pair.

    Returns DataFrame with columns:
      peptide, allele, mhcflurry_presentation_percentile
    """
    from mhcflurry import Class1PresentationPredictor

    console.print("[dim]  MHCFlurry: loading models...[/dim]")
    predictor = Class1PresentationPredictor.load()

    all_unique_kmers = sorted(set(
        peptide
        for record in sequence_records
        for peptide in generate_peptides(str(record.seq).replace('*', ''), peptide_lengths)
    ))

    console.print(
        f"[dim]  MHCFlurry: {len(all_unique_kmers)} unique k-mers "
        f"× {len(hla_alleles)} alleles[/dim]"
    )

    per_allele_dataframes = []

    for allele in hla_alleles:
        # Pass a single-element genotype so predict() evaluates each peptide
        # against exactly this one allele, giving us per-allele rows.
        allele_batch = predictor.predict(
            peptides=all_unique_kmers,
            alleles=[allele],
            verbose=0,
        )
        allele_batch = allele_batch[['peptide', 'presentation_percentile']].copy()
        allele_batch['allele'] = allele
        allele_batch = allele_batch.rename(columns={
            'presentation_percentile': 'mhcflurry_presentation_percentile',
        })
        per_allele_dataframes.append(allele_batch)

    flurry_dataframe = pd.concat(per_allele_dataframes, ignore_index=True)
    return flurry_dataframe[['peptide', 'allele', 'mhcflurry_presentation_percentile']]


# ── Step class ────────────────────────────────────────────────────────────────

class PredictBindingStep(BaseTrackStep):
    step_name = 'predict_binding'

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

        console.print(f"\n[bold]Binding predictions for {self.track_id}[/bold]")
        console.print(f"  Alleles       : {', '.join(hla_alleles)}")
        console.print(f"  Peptide lengths: {peptide_lengths}")
        console.print(f"  Sequences      : {len(sequence_records)}")
        console.print("  Running NetMHCpan and MHCFlurry in parallel...\n")

        run_start_time = time.time()

        net_dataframe    = None
        flurry_dataframe = None
        net_error        = None
        flurry_error     = None

        with ThreadPoolExecutor(max_workers=2) as thread_executor:
            net_future    = thread_executor.submit(
                _run_netmhcpan_iedb, sequence_records, hla_alleles, peptide_lengths
            )
            flurry_future = thread_executor.submit(
                _run_mhcflurry, sequence_records, hla_alleles, peptide_lengths
            )

            for completed_future in as_completed([net_future, flurry_future]):
                if completed_future is net_future:
                    try:
                        net_dataframe = completed_future.result()
                        console.print(f"[green]  ✓ NetMHCpan done: {len(net_dataframe)} rows[/green]")
                    except Exception as caught_net_error:
                        net_error = str(caught_net_error)
                        console.print(f"[red]  ✗ NetMHCpan failed: {caught_net_error}[/red]")
                else:
                    try:
                        flurry_dataframe = completed_future.result()
                        console.print(f"[green]  ✓ MHCFlurry done: {len(flurry_dataframe)} rows[/green]")
                    except Exception as caught_flurry_error:
                        flurry_error = str(caught_flurry_error)
                        console.print(f"[red]  ✗ MHCFlurry failed: {caught_flurry_error}[/red]")

        elapsed_seconds = time.time() - run_start_time

        if net_error and flurry_error:
            raise RuntimeError(
                f"Both predictions failed.\n"
                f"  NetMHCpan: {net_error}\n"
                f"  MHCFlurry: {flurry_error}"
            )

        if net_dataframe is not None:
            net_dataframe.to_csv(net_output_path, index=False)
        if flurry_dataframe is not None:
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

        console.print(f"\n  Elapsed: {elapsed_seconds:.1f}s")
        console.print(f"  Audit  : {audit_output_path}")

        return {
            'net_csv':    str(net_output_path)    if net_dataframe    is not None else None,
            'flurry_csv': str(flurry_output_path) if flurry_dataframe is not None else None,
            'audit_json': str(audit_output_path),
        }
