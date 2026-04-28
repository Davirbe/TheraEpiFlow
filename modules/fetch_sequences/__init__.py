"""
Step 01 — Fetch Sequences

Downloads or imports protein sequences for a single project track.
Runs once per track (e.g. once for HPV16_E6, once for CHIKV_E1, etc.).

GenBank mode uses the unified search (utils.genbank_utils.search_proteins_comprehensive)
which always runs two strategies in parallel:
  A) Direct search by protein name (for standalone protein records)
  B) Polyprotein search + mat_peptide extraction (for CHIKV, ZIKV, DENV, HCV...)

The user never has to decide in advance whether the virus uses a polyprotein —
the tool checks both and merges the results into a single selection table.

Output files (saved to data/input/{track_id}/):
  sequences_{track_id}.fasta            — VALIDATED sequences (canonical input for step03)
  sequence_registry_{track_id}.json     — metadata for ALL downloaded records
                                          (used by step08 to exclude reference IDs from variant search)
  validation_report_{track_id}.json     — what was rejected and why (ambiguous AAs, too short, etc.)

Validation phase (formerly step 02, merged in here on 27/04/2026):
After downloading, each record is checked by utils.fasta_utils.is_valid_sequence:
  - minimum length (≥ MIN_SEQUENCE_LENGTH)
  - no ambiguous amino acids (X, B, Z, U, O)
Records that fail are excluded from the canonical FASTA but preserved in the
registry (for traceability) and listed in the validation report (for transparency).
"""

import datetime
import json
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box

from modules.base_step import BaseTrackStep
from utils.genbank_utils import (
    search_proteins_comprehensive,
    suggest_reference_sequence,
    extract_source_qualifiers,
    record_is_refseq,
    save_sequences_as_fasta,
    build_sequence_registry,
)
from utils.fasta_utils import is_valid_sequence

console = Console(width=120)


class FetchSequencesStep(BaseTrackStep):
    step_name = 'fetch_sequences'

    def run(self, input_data=None):
        track_config  = self.project_config['tracks'][self.track_id]
        organism_name = track_config['organism_name']
        protein_name  = track_config.get('protein_name')
        target_host   = self.project_config.get('target_host', 'Homo sapiens')
        entrez_email  = self.project_config.get('entrez_email', '')
        input_source  = track_config.get('input_source', 'genbank')

        Entrez.email = entrez_email

        console.print(f'\n[bold cyan]━━━ Track: {self.track_id} ━━━[/bold cyan]')
        console.print(f'[dim]Organism : {organism_name}[/dim]')
        console.print(f'[dim]Protein  : {protein_name or "(not specified)"}[/dim]')
        console.print(f'[dim]Host     : {target_host}[/dim]')
        console.print(f'[dim]Source   : {input_source}[/dim]\n')

        if input_source == 'local':
            selected_records = _load_local_fasta(
                local_file_path=track_config.get('local_file_path')
            )
        else:
            selected_records = self._run_genbank_flow(
                organism_name=organism_name,
                protein_name=protein_name,
                target_host=target_host,
            )

        if not selected_records:
            raise ValueError(f'No sequences selected for {self.track_id}. Step aborted.')

        # ── Phase 1 of step 01: validate downloaded sequences ─────────────────
        # Triage: drop records that fail length/ambiguity checks. Keep a full
        # log so the rejection trail is auditable downstream.
        validated_records  = []
        rejected_records_log = []
        for downloaded_record in selected_records:
            record_is_valid, rejection_reason = is_valid_sequence(downloaded_record)
            if record_is_valid:
                validated_records.append(downloaded_record)
            else:
                rejected_records_log.append({
                    'id':          downloaded_record.id,
                    'description': downloaded_record.description[:80],
                    'length':      len(downloaded_record.seq),
                    'reason':      rejection_reason,
                })

        if not validated_records:
            raise ValueError(
                f'No sequences passed validation for {self.track_id}. '
                f'All {len(selected_records)} downloaded record(s) were rejected. '
                f'See validation_report_{self.track_id}.json for details.'
            )

        # ── Save outputs ──────────────────────────────────────────────────────
        track_input_dir          = self.input_dir / self.track_id
        fasta_output_path        = track_input_dir / f'sequences_{self.track_id}.fasta'
        registry_output_path     = track_input_dir / f'sequence_registry_{self.track_id}.json'
        validation_report_path   = track_input_dir / f'validation_report_{self.track_id}.json'

        # Canonical FASTA carries ONLY validated records — this is what step03 consumes.
        save_sequences_as_fasta(
            selected_records=validated_records,
            output_fasta_path=fasta_output_path,
        )
        # Registry covers ALL downloaded records (validated + rejected) so that
        # step08 (variant search) can exclude every original ID, not only the
        # ones that survived validation.
        build_sequence_registry(
            selected_records=selected_records,
            output_json_path=registry_output_path,
        )
        # Validation report — full audit trail of rejections.
        validation_report_payload = {
            'track_id':         self.track_id,
            'validated_at':     datetime.datetime.now().isoformat(),
            'total_downloaded': len(selected_records),
            'total_validated':  len(validated_records),
            'total_rejected':   len(rejected_records_log),
            'rejected_records': rejected_records_log,
        }
        with open(validation_report_path, 'w', encoding='utf-8') as report_file_handle:
            json.dump(validation_report_payload, report_file_handle, indent=2, ensure_ascii=False)

        # ── Console summary ───────────────────────────────────────────────────
        console.print(
            f'\n[bold green]Validated {len(validated_records)} '
            f'of {len(selected_records)} downloaded sequence(s):[/bold green]'
        )
        for surviving_record in validated_records:
            console.print(
                f'  [cyan]{surviving_record.id}[/cyan]  '
                f'{surviving_record.description[:65]}'
            )
        if rejected_records_log:
            console.print(
                f'\n[bold yellow]Rejected {len(rejected_records_log)} '
                f'sequence(s) during validation:[/bold yellow]'
            )
            for rejected_entry in rejected_records_log:
                console.print(
                    f'  [yellow]{rejected_entry["id"]}[/yellow]  '
                    f'({rejected_entry["length"]} aa) — {rejected_entry["reason"]}'
                )
        console.print(f'\n  FASTA    → {fasta_output_path}')
        console.print(f'  Registry → {registry_output_path}')
        console.print(f'  Report   → {validation_report_path}')

        return {
            'fasta_path':              str(fasta_output_path),
            'registry_path':           str(registry_output_path),
            'validation_report_path':  str(validation_report_path),
            'total_downloaded':        len(selected_records),
            'total_validated':         len(validated_records),
            'total_rejected':          len(rejected_records_log),
        }

    # ── GenBank flow (interactive) ────────────────────────────────────────────

    def _run_genbank_flow(self, organism_name, protein_name, target_host):
        """
        Runs the unified NCBI search and lets the user pick which sequences
        to download. On zero results, offers interactive fallbacks.
        """
        protein_name_being_searched = protein_name

        while True:
            console.print('\n[yellow]Searching NCBI (direct + polyprotein in parallel)...[/yellow]')
            comprehensive_search_result = search_proteins_comprehensive(
                organism_name=organism_name,
                protein_name=protein_name_being_searched,
                target_host=target_host,
            )

            _print_search_summary(
                search_result=comprehensive_search_result,
                organism_name=organism_name,
                protein_name=protein_name_being_searched,
            )

            candidate_records_from_search = comprehensive_search_result['records']

            if candidate_records_from_search:
                break

            # Zero results → interactive fallback
            zero_results_recovery_choice = _prompt_zero_results_recovery(
                organism_name=organism_name,
                protein_name=protein_name_being_searched,
                search_result=comprehensive_search_result,
            )
            if zero_results_recovery_choice == 'abort':
                raise ValueError(
                    f'No sequences found for organism="{organism_name}" '
                    f'protein="{protein_name_being_searched}". User aborted.'
                )
            if zero_results_recovery_choice == 'organism_only':
                protein_name_being_searched = None
                continue
            if isinstance(zero_results_recovery_choice, tuple) \
                    and zero_results_recovery_choice[0] == 'refine':
                protein_name_being_searched = zero_results_recovery_choice[1]
                continue

        suggested_reference_record, suggestion_reason_text = suggest_reference_sequence(
            candidate_records=candidate_records_from_search,
        )

        _display_sequences_table(
            records=candidate_records_from_search,
            suggested_record=suggested_reference_record,
            suggestion_reason=suggestion_reason_text,
            search_result=comprehensive_search_result,
        )

        user_selected_records = _prompt_user_selection(
            records=candidate_records_from_search,
            suggested_record=suggested_reference_record,
        )
        if not user_selected_records:
            raise ValueError('No sequences selected.')

        console.print(
            f'\n[green]Selected {len(user_selected_records)} sequence(s).[/green]'
        )
        return user_selected_records


# ── Local FASTA loader ────────────────────────────────────────────────────────

def _load_local_fasta(local_file_path: Optional[str] = None) -> list:
    """
    Loads sequences from a local FASTA file.
    Path is read from project_config (set at track setup).
    """
    if local_file_path:
        fasta_file_path = Path(local_file_path)
    else:
        console.print('\n[bold]Path to local FASTA file:[/bold]')
        fasta_file_path = Path(input('> ').strip())

    if not fasta_file_path.exists():
        raise FileNotFoundError(f'File not found: {fasta_file_path}')
    if fasta_file_path.suffix.lower() not in ('.fasta', '.fa', '.faa', '.fas'):
        console.print('[yellow]Warning: file extension not recognized as FASTA. '
                      'Attempting to parse anyway...[/yellow]')

    loaded_records = list(SeqIO.parse(str(fasta_file_path), 'fasta'))
    if not loaded_records:
        raise ValueError(f'No sequences found in file: {fasta_file_path}')

    console.print(f'[green]Loaded {len(loaded_records)} sequence(s) from {fasta_file_path}[/green]')
    return loaded_records


# ── Search summary and fallback prompt ────────────────────────────────────────

def _print_search_summary(search_result: dict, organism_name: str,
                          protein_name: Optional[str]):
    """
    Shows a short panel with the raw totals from each search strategy.
    Helps the user see why X results came back (or why none did).
    """
    total_direct_hits_in_ncbi     = search_result['total_direct_found']
    total_polyprotein_entries     = search_result['total_polyprotein_found']
    extracted_mat_peptide_count   = search_result['extracted_from_polyprotein_count']
    final_candidate_count         = len(search_result['records'])
    direct_search_tier_used       = search_result.get('direct_search_tier', 'none')

    tier_human_description = {
        'strict':        '[green]strict[/green] (product name match)',
        'loose':         '[yellow]loose[/yellow] (all fields — fallback)',
        'organism_only': '[cyan]organism only[/cyan]',
        'none':          '[red]no match found[/red]',
    }.get(direct_search_tier_used, direct_search_tier_used)

    summary_panel_lines = [
        f'[bold]Query:[/bold] {organism_name} / {protein_name or "(no protein)"}',
        f'  Direct search tier                     : {tier_human_description}',
        f'  Direct matches in NCBI                 : [cyan]{total_direct_hits_in_ncbi}[/cyan]',
        f'  Polyprotein entries in NCBI            : [cyan]{total_polyprotein_entries}[/cyan]',
        f'  Mat_peptide matches extracted          : [cyan]{extracted_mat_peptide_count}[/cyan]',
        f'  Final candidates (host-filtered, dedup): [bold cyan]{final_candidate_count}[/bold cyan]',
    ]
    console.print(Panel(
        '\n'.join(summary_panel_lines),
        box=box.ROUNDED, border_style='dim', title='[dim]Search summary[/dim]',
        title_align='left',
    ))


def _prompt_zero_results_recovery(
    organism_name: str,
    protein_name: Optional[str],
    search_result: dict,
):
    """
    Interactive fallback when no sequences were found.
    Returns one of:
      ('refine', new_protein_name)
      'organism_only'
      'abort'
    """
    console.print(Panel(
        f'[bold yellow]No sequences found[/bold yellow]\n'
        f'[dim]Direct query:     [/dim] {search_result["direct_query"]}\n'
        f'[dim]Polyprotein query:[/dim] {search_result["polyprotein_query"]}',
        box=box.ROUNDED, border_style='yellow',
    ))

    while True:
        console.print(
            '\n[bold]What to do?[/bold]\n'
            '  [cyan][r][/cyan] refine — try a different protein name (e.g. full name instead of abbreviation)\n'
            '  [cyan][o][/cyan] organism only — list every protein of this organism (manual filtering)\n'
            '  [cyan][a][/cyan] abort this step'
        )
        recovery_option_input = input('> ').strip().lower()

        if recovery_option_input in ('r', 'refine'):
            console.print(
                '[bold]New protein name[/bold] '
                '[dim](examples: "envelope glycoprotein E1", "spike", "surface glycoprotein")[/dim]'
            )
            refined_protein_name_input = input('> ').strip()
            if refined_protein_name_input:
                return ('refine', refined_protein_name_input)
            console.print('[dim]Empty input — try again.[/dim]')
            continue

        if recovery_option_input in ('o', 'organism_only'):
            return 'organism_only'

        if recovery_option_input in ('a', 'abort', 'q'):
            return 'abort'

        console.print('[dim]Unrecognized option — use r / o / a.[/dim]')


# ── Display table ─────────────────────────────────────────────────────────────

def _display_sequences_table(records, suggested_record, suggestion_reason,
                             search_result):
    """
    Renders the Rich table of candidate sequences.
    The Source column tells the user whether each row came from the direct
    search or was extracted from a polyprotein.
    """
    total_returned = len(records)
    direct_count    = sum(1 for r in records
                          if r.annotations.get('search_source') == 'direct')
    extracted_count = total_returned - direct_count

    table = Table(
        box=box.ROUNDED,
        show_header=True,
        header_style='bold white',
        title=(
            f'Sequences available — {total_returned} total  '
            f'(direct: {direct_count}  /  extracted from polyprotein: {extracted_count})'
        ),
        title_style='bold cyan',
    )

    table.add_column('#',                min_width=3,  no_wrap=True, justify='right')
    table.add_column('',                 min_width=1,  no_wrap=True)   # ★ marker
    table.add_column('Accession',        min_width=14, no_wrap=True, style='cyan')
    table.add_column('aa',               min_width=5,  no_wrap=True, justify='right')
    table.add_column('Strain / Isolate', min_width=18, no_wrap=True, max_width=24)
    table.add_column('Location',         min_width=14, no_wrap=True, max_width=22)
    table.add_column('Date',             min_width=10, no_wrap=True)
    table.add_column('Source',           min_width=10, no_wrap=True)

    for row_index, sequence_record in enumerate(records, start=1):
        qualifiers = extract_source_qualifiers(sequence_record)

        is_suggested = sequence_record.id == suggested_record.id
        row_style    = 'yellow' if is_suggested else ''
        star_marker  = '★' if is_suggested else ''

        strain_value  = qualifiers.get('strain',  'N/A')
        isolate_value = qualifiers.get('isolate', 'N/A')
        strain_or_isolate = strain_value if strain_value != 'N/A' else isolate_value

        # Build source label combining provenance + RefSeq flag
        search_source = sequence_record.annotations.get('search_source', 'direct')
        is_refseq     = record_is_refseq(sequence_record)
        if search_source == 'extracted':
            source_label = 'Extracted'
        elif is_refseq:
            source_label = 'RefSeq'
        else:
            source_label = 'GenBank'

        table.add_row(
            str(row_index),
            star_marker,
            sequence_record.id,
            str(len(sequence_record.seq)),
            strain_or_isolate[:22],
            qualifiers.get('geographic_location', 'N/A')[:20],
            qualifiers.get('collection_date', 'N/A'),
            source_label,
            style=row_style,
        )

    console.print(table)
    console.print(
        f'[dim]★ Suggestion: {suggested_record.id} — {suggestion_reason}[/dim]\n'
    )


# ── Selection prompt ──────────────────────────────────────────────────────────

def _prompt_user_selection(records, suggested_record):
    """
    Prompts user to select one or more sequences from the table.

    Accepted formats:
      1        → row 1 only
      1,3,5    → rows 1, 3 and 5
      1-4      → rows 1 through 4
      all      → all rows
      (blank)  → accept the suggested row (★)
    """
    suggested_row_number = next(
        (index for index, record in enumerate(records, start=1)
         if record.id == suggested_record.id),
        1,
    )

    console.print(
        '[bold]Select sequences to download[/bold] '
        '(e.g. [cyan]1[/cyan], [cyan]1,3[/cyan], [cyan]1-4[/cyan], [cyan]all[/cyan])'
    )
    console.print(f'[dim]Press Enter to accept suggestion (★ row {suggested_row_number})[/dim]')

    raw_selection_input = input('> ').strip()

    if raw_selection_input == '':
        return [suggested_record]

    if raw_selection_input.lower() == 'all':
        return list(records)

    selected_row_numbers = _parse_row_selection(
        raw_selection_input=raw_selection_input,
        total_rows=len(records),
    )

    if not selected_row_numbers:
        console.print('[red]Invalid selection.[/red]')
        return []

    return [records[row_number - 1] for row_number in selected_row_numbers]


def _parse_row_selection(raw_selection_input, total_rows):
    """
    Parses a selection string into a sorted list of valid row numbers.

    Handles: '1', '1,3,5', '1-4', '1,3-5,7'
    Silently ignores numbers outside [1, total_rows].
    """
    selected_row_numbers = set()
    try:
        for selection_token in raw_selection_input.split(','):
            selection_token = selection_token.strip()
            if '-' in selection_token:
                range_start, range_end = selection_token.split('-', 1)
                for row_number in range(int(range_start.strip()), int(range_end.strip()) + 1):
                    if 1 <= row_number <= total_rows:
                        selected_row_numbers.add(row_number)
            else:
                row_number = int(selection_token)
                if 1 <= row_number <= total_rows:
                    selected_row_numbers.add(row_number)
    except (ValueError, IndexError):
        return []

    return sorted(selected_row_numbers)
