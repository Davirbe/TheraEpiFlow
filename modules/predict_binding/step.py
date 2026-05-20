"""Orchestration for predict_binding: runs both predictors over a track and
writes PRED_NET / PRED_FLURRY CSVs + the audit JSON."""

import datetime
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from rich import box
from rich.panel import Panel
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn
from rich.text import Text

from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.naming import get_prediction_filename, get_step_filename
from utils.step_summary import print_step_summary

from .core import (
    SUMMARY_INTERMEDIATE_BINDER_RANK_MAX,
    SUMMARY_STRONG_BINDER_RANK_MAX,
    _count_presentation_ready,
    _count_strong_and_intermediate_binders,
    _peptides_in_common,
    _run_mhcflurry_with_progress,
    _run_netmhcpan_iedb_silent,
)
from .io import _load_sequences
from .prompts import _ask_binding_params

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
        "[bold]PRED_VIEW_{track_id}.csv[/bold]     — slim view: one row per (peptide, allele) with both tools' percentiles side-by-side.\n"
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

        # Slim VIEW: per (peptide, allele), one row with NetMHCpan + MHCFlurry percentiles side-by-side.
        view_path = predictions_dir / get_step_filename("PRED_VIEW", self.track_id)
        if net_dataframe is not None and flurry_dataframe is not None:
            net_slim    = net_dataframe[['peptide', 'allele', 'netmhcpan_el_percentile']]
            flurry_slim = flurry_dataframe[['peptide', 'allele', 'mhcflurry_presentation_percentile']]
            view_dataframe = net_slim.merge(flurry_slim, on=['peptide', 'allele'], how='outer')
            view_dataframe.to_csv(view_path, index=False)

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
            predictions_dir / get_step_filename("PRED_VIEW", self.track_id):
                "Slim view — one row per (peptide, allele) with NetMHCpan + MHCFlurry percentiles side-by-side.",
            predictions_dir / get_step_filename("PREDICT_AUDIT", self.track_id, ext="json"):
                "Run audit — alleles used, lengths, row counts per tool, elapsed time, and any errors.",
        }
