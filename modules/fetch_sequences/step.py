"""Orchestration for fetch_sequences: resolves input, fetches sequences and
writes the per-track FASTA + registry."""

import csv
import datetime
import io
import json
from pathlib import Path
from typing import Optional

import requests
from Bio import SeqIO

from modules.base_step import BaseTrackStep
from utils.console import console, is_interactive_session
from utils.naming import get_step_filename
from utils.project_manager import save_project_config

from .core import (
    ORGANISM_ALIASES,
    _build_registry_entry,
    _normalize_organism,
    _search_uniprot,
    _sort_and_flag,
    _validate_records,
)
from .io import _download_fasta, _load_local_fasta
from .prompts import _prompt_selection
from .render import _display_uniprot_table

# ── Step class ─────────────────────────────────────────────────────────────────

class FetchSequencesStep(BaseTrackStep):
    step_name   = 'fetch_sequences'
    description = (
        "Pulls the reference protein sequence for each track from UniProt, "
        "validates it (minimum length, no ambiguous residues), and writes "
        "the FASTA the rest of the pipeline reads."
    )
    long_description = (
        "For each organism/protein pair, queries the UniProt REST API for "
        "the protein you named, ranks the candidate hits (Swiss-Prot "
        "preferred over TrEMBL, then by length), lets you pick the seed "
        "sequence, validates it for MHC-I-suitable amino acids (no ambiguous "
        "B, J, O, U, X, Z and a minimum length), and writes the reference "
        "FASTA used by every downstream step.\n\n"
        "You can also bypass the UniProt search and supply your own local "
        "FASTA — the validation rules are the same."
    )
    methodology = (
        "1. Query: UniProt REST API (/uniprotkb/search) with `protein_name AND organism_id:<taxid>` "
        "(taxid resolved automatically from common name / scientific name / accession).\n"
        "2. Ranking: Swiss-Prot reviewed entries first; ties broken by sequence length "
        "(closer to the median of all hits wins) and then by isoform canonicity.\n"
        "3. Validation: rejects sequences shorter than 30 aa or containing non-canonical "
        "residues (B, J, O, U, X, Z) — those would crash NetMHCpan / MHCFlurry later.\n"
        "4. Persistence: writes the FASTA, a registry of every candidate considered, and a "
        "validation report (accepted/rejected with reasons)."
    )
    references = [
        {
            'authors': 'The UniProt Consortium',
            'title':   'UniProt: the Universal Protein Knowledgebase in 2023',
            'journal': 'Nucleic Acids Research',
            'year':    2023,
            'doi':     '10.1093/nar/gkac1052',
        },
    ]
    data_format = (
        "Input is the organism + protein you defined when setting up the project:\n"
        "  • [bold]Organism[/bold] — alias (HPV16, ZIKV, CHIKV) or scientific name "
        "([italic]Human papillomavirus 16[/italic]).\n"
        "  • [bold]Protein[/bold] — name as it appears in UniProt "
        "(e.g. [italic]E6[/italic], [italic]envelope protein[/italic], [italic]nsP1[/italic]).\n\n"
        "If you choose [cyan]local FASTA[/cyan] instead, point to a `.fasta` / `.fa` file "
        "containing one or more protein sequences in standard FASTA format."
    )
    outputs_overview = (
        "[bold]SEQUENCES_{track_id}.fasta[/bold]      — reference protein FASTA (input for every downstream step).\n"
        "[bold]SEQUENCES_VIEW_{track_id}.csv[/bold]   — slim per-step view (accession, organism, length, source).\n"
        "[bold]REGISTRY_{track_id}.json[/bold]        — every UniProt candidate considered (full audit).\n"
        "[bold]VALIDATION_REPORT_{track_id}.json[/bold] — accepted vs rejected sequences + reasons."
    )
    tips = [
        "Use the same wording UniProt uses — \"E6\", \"envelope protein\", \"nsP1\" — to maximise hit precision.",
        "Lowercase / uppercase does not matter for the search, but [bold]avoid accents[/bold] (Vírus → Virus).",
        "Aliases (HPV16, ZIKV, CHIKV) are resolved automatically; full scientific names also work.",
        "Uncommon proteins may return few candidates — pick the one with the closest length and Swiss-Prot status.",
        "After this step, the seed FASTA can be inspected via [cyan]b[/cyan] (file browser) in the main menu.",
    ]

    def describe_outputs(self) -> dict:
        track_input_dir = self.input_dir / self.track_id
        return {
            track_input_dir / get_step_filename('SEQUENCES_VIEW', self.track_id):
                "Slim per-step view — one row with track_id, accession, organism, protein, length, source.",
            track_input_dir / get_step_filename('SEQUENCES', self.track_id, ext='fasta'):
                "Reference protein FASTA — used as query for every downstream step.",
            track_input_dir / get_step_filename('REGISTRY', self.track_id, ext='json'):
                "Registry of every candidate UniProt hit considered during the search.",
            track_input_dir / get_step_filename('VALIDATION_REPORT', self.track_id, ext='json'):
                "Validation report — counts of accepted/rejected sequences with reasons.",
        }

    def run(self, input_data=None):
        track_config  = self.project_config['tracks'][self.track_id]
        organism_name = track_config['organism_name']
        protein_name  = track_config.get('protein_name')
        input_source  = track_config.get('input_source', 'uniprot')

        console.print(f'\n[bold cyan]━━━ Track: {self.track_id} ━━━[/bold cyan]')
        console.print(f'[dim]Organism : {organism_name}[/dim]')
        console.print(f'[dim]Protein  : {protein_name or "(not specified)"}[/dim]')
        console.print(f'[dim]Source   : {input_source}[/dim]\n')

        if input_source == 'local':
            selected_records = _load_local_fasta(track_config.get('local_file_path'))
            selected_hit     = None
            validated_records, rejected_log = _validate_records(selected_records)
            if not validated_records:
                raise ValueError(
                    f'No local sequences passed validation for {self.track_id}.'
                )
        else:
            selected_records, selected_hit, validated_records, rejected_log = \
                self._run_uniprot_flow(organism_name, protein_name)

        # ── Save outputs ──────────────────────────────────────────────────────
        track_input_dir        = self.input_dir / self.track_id
        fasta_path             = track_input_dir / get_step_filename('SEQUENCES', self.track_id, ext='fasta')
        registry_path          = track_input_dir / get_step_filename('REGISTRY', self.track_id, ext='json')
        validation_report_path = track_input_dir / get_step_filename('VALIDATION_REPORT', self.track_id, ext='json')

        with open(fasta_path, 'w') as fh:
            SeqIO.write(validated_records, fh, 'fasta')

        registry_payload = {
            'total_sequences': len(selected_records),
            'sequences': [_build_registry_entry(r, selected_hit) for r in selected_records],
        }
        with open(registry_path, 'w', encoding='utf-8') as fh:
            json.dump(registry_payload, fh, indent=2, ensure_ascii=False)

        validation_payload = {
            'track_id':         self.track_id,
            'validated_at':     datetime.datetime.now().isoformat(),
            'total_downloaded': len(selected_records),
            'total_validated':  len(validated_records),
            'total_rejected':   len(rejected_log),
            'rejected_records': rejected_log,
        }
        with open(validation_report_path, 'w', encoding='utf-8') as fh:
            json.dump(validation_payload, fh, indent=2, ensure_ascii=False)

        view_path = track_input_dir / get_step_filename('SEQUENCES_VIEW', self.track_id)
        seed_accession = selected_hit['accession'] if selected_hit else (
            validated_records[0].id if validated_records else ''
        )
        seed_length = selected_hit['length'] if selected_hit else (
            len(validated_records[0].seq) if validated_records else 0
        )
        with open(view_path, 'w', newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(['track_id', 'accession', 'organism', 'protein', 'length', 'source'])
            writer.writerow([
                self.track_id,
                seed_accession,
                organism_name,
                protein_name or '',
                seed_length,
                input_source,
            ])

        # ── Save seed metadata to project_config ──────────────────────────────
        if selected_hit:
            cfg = self.project_config
            cfg['tracks'][self.track_id]['seed_accession'] = selected_hit['accession']
            cfg['tracks'][self.track_id]['seed_size']      = selected_hit['length']
            cfg['tracks'][self.track_id]['tax_id']         = selected_hit['tax_id']
            save_project_config(self.project_name, cfg)

        # ── Console summary ───────────────────────────────────────────────────
        console.print(
            f'\n[bold green]Validated {len(validated_records)} '
            f'of {len(selected_records)} sequence(s):[/bold green]'
        )
        for r in validated_records:
            console.print(f'  [cyan]{r.id}[/cyan]  {r.description[:65]}')
        if rejected_log:
            console.print(f'\n[bold yellow]Rejected {len(rejected_log)} sequence(s):[/bold yellow]')
            for entry in rejected_log:
                console.print(
                    f'  [yellow]{entry["id"]}[/yellow]  '
                    f'({entry["length"]} aa) — {entry["reason"]}'
                )
        console.print(f'\n  FASTA    → {fasta_path}')
        console.print(f'  Registry → {registry_path}')
        console.print(f'  Report   → {validation_report_path}')

        return {
            'fasta_path':              str(fasta_path),
            'registry_path':           str(registry_path),
            'validation_report_path':  str(validation_report_path),
            'total_downloaded':        len(selected_records),
            'total_validated':         len(validated_records),
            'total_rejected':          len(rejected_log),
        }

    def _run_uniprot_flow(
        self,
        organism_name: str,
        protein_name: Optional[str],
    ) -> tuple[list, dict, list, list]:
        """Full UniProt search + selection + download + validation flow.

        Auto-advances to the next valid hit if the downloaded sequence fails validation.
        Returns (selected_records, selected_hit, validated_records, rejected_log).

        TODO — Flavivirus polyprotein extraction:
        For organisms like DENV, ZIKV, HCV, the individual proteins (E, NS5, etc.) are
        stored as Chain features inside a single "Genome polyprotein" UniProt entry.
        Planned fix: if the selected hit is flagged (Polyprotein?) AND the search had a
        protein_name, fetch /{accession}.json, search features[type='Chain'] for a description
        matching protein_name (fuzzy), and slice seq[start-1:end] before saving to FASTA.
        Reference: UniProt JSON features[].type == "Chain", fields: location.{start,end}.value, description.
        """
        resolved_name, tax_id, was_aliased = _normalize_organism(organism_name)
        if was_aliased:
            console.print(
                f'[dim]→ Organism resolved: '
                f'[bold]{organism_name}[/bold] → [bold cyan]{resolved_name}[/bold cyan][/dim]'
            )

        with console.status(
            f'[yellow]Searching UniProt for "{protein_name or "all proteins"}" '
            f'in {resolved_name}…[/yellow]',
            spinner='dots',
        ):
            hits, query_used = _search_uniprot(
                resolved_name, protein_name, tax_id=tax_id,
            )

        if not hits:
            raise ValueError(
                f'No results found in UniProt for '
                f'"{resolved_name}" / "{protein_name}".'
            )

        hits = _sort_and_flag(hits)
        console.print(f'[dim]Query: {query_used}[/dim]')
        _display_uniprot_table(hits, resolved_name, protein_name)

        non_interactive = not is_interactive_session()

        tried: set[str] = set()

        while True:
            remaining = [h for h in hits if h['accession'] not in tried]
            if not remaining:
                raise ValueError(
                    f'All tested hits failed validation for '
                    f'{self.track_id}. Please review results manually.'
                )

            selected_hit = _prompt_selection(remaining, non_interactive=non_interactive)
            tried.add(selected_hit['accession'])

            with console.status(
                f'[yellow]Downloading FASTA: {selected_hit["accession"]} '
                f'({selected_hit["protein_name"][:40]})…[/yellow]',
                spinner='dots',
            ):
                records = list(SeqIO.parse(
                    io.StringIO(_download_fasta(selected_hit['accession'])), 'fasta',
                ))

            if not records:
                console.print(
                    f'[yellow]Empty FASTA for {selected_hit["accession"]}. '
                    f'Trying next...[/yellow]'
                )
                continue

            console.print(f'[green]Download complete: {len(records[0].seq)} aa[/green]')

            validated, rejected = _validate_records(records)
            if validated:
                return records, selected_hit, validated, rejected

            reason = rejected[0]['reason'] if rejected else 'unknown reason'
            next_candidates = [h for h in hits if h['accession'] not in tried]
            if next_candidates:
                console.print(
                    f'[yellow]→ {selected_hit["accession"]} rejected '
                    f'({reason}). Trying next: '
                    f'{next_candidates[0]["accession"]}...[/yellow]'
                )
            else:
                raise ValueError(
                    f'No valid sequence found for {self.track_id}. '
                    f'Last error: {reason}'
                )
