"""
Step 01 — Fetch Sequences

Downloads or imports protein sequences for a single project track.
Runs once per track (e.g. once for HPV16_E6, once for HPV16_E7, etc.).

At the start of each track, asks the user three questions:
  1. Source: GenBank search or local FASTA file?
  2. If GenBank: single protein or polyprotein?
  3. If polyprotein: which mature peptide to extract?

Two-phase GenBank workflow:
  Phase 1 — Discovery: search NCBI, get IDs, fetch metadata → show table (no sequences yet)
  Phase 2 — Download:  fetch full sequences only for the IDs the user selected

Output files (saved to data/input/{track_id}/):
  sequences_{track_id}.fasta         — selected sequences
  sequence_registry_{track_id}.json  — metadata for each sequence (used by step08)
"""

from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from rich.console import Console
from rich.table import Table
from rich import box

from modules.base_step import BaseTrackStep
from utils.genbank_utils import (
    search_ncbi_protein_ids,
    fetch_records_by_accession_ids,
    suggest_reference_sequence,
    extract_source_qualifiers,
    extract_protein_from_polyprotein,
    record_is_refseq,
    record_is_polyprotein,
    save_sequences_as_fasta,
    build_sequence_registry,
)

console = Console(width=120)

MAX_RECORDS_TO_DISPLAY = 50


class FetchSequencesStep(BaseTrackStep):
    step_number = 1
    step_name   = 'fetch_sequences'

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

        # Input source read from project_config — defined at track setup
        if input_source == 'local':
            selected_records = _load_local_fasta(
                local_file_path=track_config.get('local_file_path')
            )
        else:
            selected_records = self._run_genbank_flow(
                organism_name=organism_name,
                protein_name=protein_name,
                target_host=target_host,
                include_polyproteins=track_config.get('search_include_polyproteins', False),
            )

        if not selected_records:
            raise ValueError(f'No sequences selected for {self.track_id}. Step aborted.')

        # ── Save outputs ──────────────────────────────────────────────────────
        track_input_dir      = self.input_dir / self.track_id
        fasta_output_path    = track_input_dir / f'sequences_{self.track_id}.fasta'
        registry_output_path = track_input_dir / f'sequence_registry_{self.track_id}.json'

        save_sequences_as_fasta(
            selected_records=selected_records,
            output_fasta_path=fasta_output_path,
        )
        build_sequence_registry(
            selected_records=selected_records,
            output_json_path=registry_output_path,
        )

        console.print(f'\n[bold green]Saved {len(selected_records)} sequence(s):[/bold green]')
        for saved_record in selected_records:
            console.print(f'  [cyan]{saved_record.id}[/cyan]  {saved_record.description[:65]}')
        console.print(f'\n  FASTA    → {fasta_output_path}')
        console.print(f'  Registry → {registry_output_path}')

        return {
            'fasta_path':     str(fasta_output_path),
            'registry_path':  str(registry_output_path),
            'sequence_count': len(selected_records),
        }

    def _run_genbank_flow(self, organism_name, protein_name, target_host,
                          include_polyproteins=False):
        """
        Handles the two-phase GenBank workflow:
          Phase 1 — search IDs + fetch metadata for display table
          Phase 2 — download full sequences only for selected IDs
        include_polyproteins is read from project_config (set at track setup).
        """
        # ── Phase 1: search IDs (no sequence download) ────────────────────────
        console.print('\n[yellow]Searching NCBI...[/yellow]')
        total_available, all_discovered_ids = search_ncbi_protein_ids(
            organism_name=organism_name,
            protein_name=protein_name,
            target_host=target_host,
            include_polyproteins=include_polyproteins,
            max_ids_to_return=200,
        )

        if not all_discovered_ids:
            console.print('[bold red]No sequences found. Check organism/protein names.[/bold red]')
            raise ValueError(f'No sequences found for organism="{organism_name}" protein="{protein_name}"')

        console.print(
            f'[green]{total_available} sequences found in NCBI.[/green] '
            f'Fetching metadata for the first {min(MAX_RECORDS_TO_DISPLAY, len(all_discovered_ids))} '
            f'to display...'
        )

        # Fetch only the first N IDs for display (metadata + sequence in GenBank format)
        ids_for_display = all_discovered_ids[:MAX_RECORDS_TO_DISPLAY]
        records_for_display = fetch_records_by_accession_ids(
            accession_ids=ids_for_display,
            target_host=target_host,
        )

        if not records_for_display:
            console.print('[bold red]No sequences remain after host filtering.[/bold red]')
            raise ValueError(f'No host-matching sequences found for {organism_name}')

        # ── Question 3 (conditional): polyprotein → which protein to extract? ─
        target_protein_to_extract = None
        if include_polyproteins and protein_name:
            polyprotein_records = [r for r in records_for_display
                                   if record_is_polyprotein(r)]
            if polyprotein_records:
                target_protein_to_extract = _ask_protein_to_extract_from_polyprotein(
                    protein_name_hint=protein_name,
                )

        # Separate and optionally extract from polyproteins
        display_records = _prepare_display_records(
            downloaded_records=records_for_display,
            include_polyproteins=include_polyproteins,
            target_protein_to_extract=target_protein_to_extract,
        )

        suggested_record, suggestion_reason = suggest_reference_sequence(display_records)

        # ── Display table ─────────────────────────────────────────────────────
        _display_sequences_table(
            records=display_records,
            suggested_record=suggested_record,
            suggestion_reason=suggestion_reason,
            total_available=total_available,
        )

        # ── User selection ────────────────────────────────────────────────────
        selected_rows = _prompt_user_selection(display_records, suggested_record)

        if not selected_rows:
            raise ValueError('No sequences selected.')

        # Phase 2: download only the selected sequences fully
        # (records are already in memory from the display fetch — reuse them)
        selected_records = selected_rows
        console.print(
            f'\n[green]Selected {len(selected_records)} sequence(s) — '
            f'no additional download needed.[/green]'
        )

        return selected_records


# ── Interactive prompts ────────────────────────────────────────────────────────

def _ask_protein_to_extract_from_polyprotein(protein_name_hint: str) -> str:
    """
    When the user opts to include polyproteins, asks which mature peptide
    to extract from the polyprotein's mat_peptide features.

    The protein_name_hint pre-fills the prompt with the track's protein name.
    Returns the protein name string to search in mat_peptide /product qualifiers.
    """
    console.print(
        f'\n[bold]Which protein to extract from polyprotein?[/bold]\n'
        f'[dim]This will search the mat_peptide features of each polyprotein '
        f'for a /product matching this name (case-insensitive, partial match).[/dim]\n'
        f'[dim]Press Enter to use: "{protein_name_hint}"[/dim]'
    )

    raw_input_value = input('> ').strip()
    return raw_input_value if raw_input_value else protein_name_hint


def _load_local_fasta(local_file_path: Optional[str] = None) -> list:
    """
    Loads sequences from a local FASTA file.
    The path is read from project_config (set at track setup).
    Returns a list of SeqRecord objects.
    """
    if local_file_path:
        fasta_file_path = Path(local_file_path)
    else:
        # Fallback: ask if path was not saved in config
        console.print('\n[bold]Path to local FASTA file:[/bold]')
        fasta_file_path = Path(input('> ').strip())

    if not fasta_file_path.exists():
        raise FileNotFoundError(f'File not found: {fasta_file_path}')
    if not fasta_file_path.suffix.lower() in ('.fasta', '.fa', '.faa', '.fas'):
        console.print('[yellow]Warning: file extension not recognized as FASTA. '
                      'Attempting to parse anyway...[/yellow]')

    loaded_records = list(SeqIO.parse(str(fasta_file_path), 'fasta'))

    if not loaded_records:
        raise ValueError(f'No sequences found in file: {fasta_file_path}')

    console.print(f'[green]Loaded {len(loaded_records)} sequence(s) from {fasta_file_path}[/green]')
    return loaded_records


# ── Record preparation ─────────────────────────────────────────────────────────

def _prepare_display_records(
    downloaded_records: list,
    include_polyproteins: bool,
    target_protein_to_extract: Optional[str],
) -> list:
    """
    Prepares the list of records to display in the table.

    If include_polyproteins is False: returns only non-polyprotein records.
    If include_polyproteins is True and target_protein_to_extract is set:
      extracts the target mature peptide from each polyprotein and adds it
      to the display list alongside any standalone records.
    """
    standalone_records  = [r for r in downloaded_records if not record_is_polyprotein(r)]
    polyprotein_records = [r for r in downloaded_records if record_is_polyprotein(r)]

    if not include_polyproteins:
        return standalone_records if standalone_records else downloaded_records

    # Include polyproteins: optionally extract specific mature peptide
    extracted_records = []
    if target_protein_to_extract:
        for polyprotein_record in polyprotein_records:
            extracted = extract_protein_from_polyprotein(
                polyprotein_record=polyprotein_record,
                target_protein_name=target_protein_to_extract,
            )
            if extracted is not None:
                extracted_records.append(extracted)

        if extracted_records:
            console.print(
                f'[green]Extracted {len(extracted_records)} "{target_protein_to_extract}" '
                f'sequence(s) from polyprotein records.[/green]'
            )

    all_display_records = standalone_records + extracted_records
    return all_display_records if all_display_records else downloaded_records


# ── Display ────────────────────────────────────────────────────────────────────

def _display_sequences_table(records, suggested_record, suggestion_reason, total_available):
    """
    Renders a Rich table listing all candidate sequences.
    Suggested sequence is highlighted in yellow with ★.
    RefSeq entries are labeled; extracted-from-polyprotein entries are marked.
    """
    table = Table(
        box=box.ROUNDED,
        show_header=True,
        header_style='bold white',
        title=f'Sequences available — showing {len(records)} of {total_available} total in NCBI',
        title_style='bold cyan',
    )

    table.add_column('#',                min_width=3,  no_wrap=True, justify='right')
    table.add_column('',                min_width=1,  no_wrap=True)   # ★ marker
    table.add_column('Accession',       min_width=14, no_wrap=True, style='cyan')
    table.add_column('aa',              min_width=5,  no_wrap=True, justify='right')
    table.add_column('Strain / Isolate', min_width=18, no_wrap=True, max_width=24)
    table.add_column('Location',        min_width=16, no_wrap=True, max_width=26)
    table.add_column('Date',            min_width=12, no_wrap=True)
    table.add_column('Source',          min_width=8,  no_wrap=True)

    for row_index, sequence_record in enumerate(records, start=1):
        qualifiers = extract_source_qualifiers(sequence_record)

        is_suggested = sequence_record.id == suggested_record.id
        row_style    = 'yellow' if is_suggested else ''
        star_marker  = '★' if is_suggested else ''

        strain_value  = qualifiers.get('strain',  'N/A')
        isolate_value = qualifiers.get('isolate', 'N/A')
        strain_or_isolate = strain_value if strain_value != 'N/A' else isolate_value

        extracted_from_poly = sequence_record.annotations.get(
            'extracted_from_polyprotein', False
        )
        if record_is_refseq(sequence_record):
            source_label = 'RefSeq'
        elif extracted_from_poly:
            source_label = 'Extracted'
        else:
            source_label = 'GenBank'

        table.add_row(
            str(row_index),
            star_marker,
            sequence_record.id,
            str(len(sequence_record.seq)),
            strain_or_isolate[:22],
            qualifiers.get('geographic_location', 'N/A')[:24],
            qualifiers.get('collection_date', 'N/A'),
            source_label,
            style=row_style,
        )

    console.print(table)
    console.print(
        f'[dim]★ Suggestion: {suggested_record.id} — {suggestion_reason}[/dim]\n'
    )


# ── Selection ──────────────────────────────────────────────────────────────────

def _prompt_user_selection(records, suggested_record):
    """
    Prompts user to select one or more sequences from the table.

    Accepted formats:
      1        → row 1 only
      1,3,5    → rows 1, 3 and 5
      1-4      → rows 1 through 4
      all      → all rows
      (blank)  → accept the suggested row (★)

    Returns a list of selected SeqRecord objects.
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
    Silently ignores row numbers outside [1, total_rows].
    Returns empty list if input cannot be parsed.
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


