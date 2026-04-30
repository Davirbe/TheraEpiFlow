"""
Step fetch_sequences — UniProt-based sequence retrieval.

Replaces GenBank/Entrez with UniProt REST API:
  Phase 1: Normalize organism name (alias dict + difflib fuzzy match)
  Phase 2: Lightweight metadata search (no FASTA yet)
  Phase 3: Sort — Swiss-Prot (reviewed) first, TrEMBL second
           Flag polyproteins/fragments based on global length median
  Phase 4: Rich table display + user selection
           Enter / EOF (non-interactive) = first non-flagged hit
  Phase 5: Download + validate selected hit; auto-advance if invalid
  Phase 6: Save outputs (SEQUENCES, REGISTRY, VALIDATION_REPORT)

Saves seed_size and tax_id to project_config for use by analyze_conservation.
Local FASTA mode (input_source='local') is preserved unchanged.
"""

import datetime
import difflib
import io
import json
import time
from pathlib import Path
from statistics import median
from typing import Optional

import requests
from Bio import SeqIO
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box

from modules.base_step import BaseTrackStep
from utils.fasta_utils import is_valid_sequence
from utils.naming import get_step_filename
from utils.project_manager import save_project_config

console = Console(width=120)

UNIPROT_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
UNIPROT_FASTA_URL  = 'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
_REQUEST_TIMEOUT   = 15
_MAX_RESULTS       = 25

# (scientific_name, ncbi_taxonomy_id)
ORGANISM_ALIASES: dict[str, tuple[str, int]] = {
    'CHIKV':  ('Chikungunya virus',                                         37124),
    'ZIKV':   ('Zika virus',                                               64320),
    'HPV16':  ('Human papillomavirus 16',                                 333760),
    'HPV18':  ('Human papillomavirus 18',                                 333761),
    'HPV31':  ('Human papillomavirus 31',                                  10582),
    'HPV33':  ('Human papillomavirus 33',                                  10583),
    'HPV45':  ('Human papillomavirus 45',                                  10585),
    'HPV52':  ('Human papillomavirus 52',                                  10588),
    'HPV58':  ('Human papillomavirus 58',                                  10590),
    'DENV':   ('Dengue virus',                                             12637),
    'DENV1':  ('Dengue virus 1',                                           11053),
    'DENV2':  ('Dengue virus 2',                                           11060),
    'DENV3':  ('Dengue virus 3',                                           11069),
    'DENV4':  ('Dengue virus 4',                                           11070),
    'HCV':    ('Hepatitis C virus',                                        11103),
    'MPOX':   ('Monkeypox virus',                                          10244),
    'EBOV':   ('Zaire ebolavirus',                                        186538),
    'HIV':    ('Human immunodeficiency virus 1',                           11676),
    'SARS2':  ('Severe acute respiratory syndrome coronavirus 2',        2697049),
    'MERS':   ('Middle East respiratory syndrome-related coronavirus',    1335626),
    'INFA':   ('Influenza A virus',                                       11320),
    'INFB':   ('Influenza B virus',                                       11520),
    'RSV':    ('Human respiratory syncytial virus',                       11250),
    'RABV':   ('Rabies lyssavirus',                                       11292),
}

# Reverse lookup: canonical scientific name → tax_id (built from ORGANISM_ALIASES)
_NAME_TO_TAX_ID: dict[str, int] = {name: tid for name, tid in ORGANISM_ALIASES.values()}


# ── HTTP helper ────────────────────────────────────────────────────────────────

def _http_get(url: str, params: dict = None, max_attempts: int = 3) -> requests.Response:
    """GET with exponential-backoff retry for transient network errors."""
    last_err: Exception = RuntimeError('no attempts made')
    for attempt in range(max_attempts):
        try:
            response = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as err:
            last_err = err
            if attempt < max_attempts - 1:
                wait = 2 ** attempt
                console.print(f'[yellow]⚠ UniProt request failed (attempt {attempt + 1}/{max_attempts}): '
                               f'{err} — retrying in {wait}s[/yellow]')
                time.sleep(wait)
    raise last_err


# ── Input normalization ────────────────────────────────────────────────────────

def _normalize_organism(raw_name: str) -> tuple[str, Optional[int], bool]:
    """
    Resolves organism abbreviation/typo to canonical scientific name + tax_id.
    Returns (resolved_name, tax_id_or_None, was_aliased).
    Handles three cases:
      1. Alias shorthand (HPV16 → ...)
      2. Fuzzy match on alias keys
      3. Already a canonical scientific name stored in _NAME_TO_TAX_ID
    """
    upper = raw_name.strip().upper()
    if upper in ORGANISM_ALIASES:
        name, tax_id = ORGANISM_ALIASES[upper]
        return name, tax_id, True
    matches = difflib.get_close_matches(upper, ORGANISM_ALIASES.keys(), n=1, cutoff=0.75)
    if matches:
        name, tax_id = ORGANISM_ALIASES[matches[0]]
        return name, tax_id, True
    stripped = raw_name.strip()
    tax_id = _NAME_TO_TAX_ID.get(stripped)
    return stripped, tax_id, False


# ── UniProt search ─────────────────────────────────────────────────────────────

def _extract_protein_name(result: dict) -> str:
    """Extracts best protein name string from a UniProt JSON result entry."""
    desc = result.get('proteinDescription', {})
    recommended = desc.get('recommendedName', {})
    if recommended:
        return recommended.get('fullName', {}).get('value', 'N/A')
    submitted = desc.get('submittedName', [])
    if submitted:
        return submitted[0].get('fullName', {}).get('value', 'N/A')
    return 'N/A'


def _build_hits(results: list) -> list[dict]:
    """Converts raw UniProt JSON result list to normalized hit dicts."""
    hits = []
    for r in results:
        entry_type = r.get('entryType', '')
        reviewed   = 'Swiss-Prot' in entry_type
        organism   = r.get('organism', {})
        hits.append({
            'accession':    r.get('primaryAccession', ''),
            'reviewed':     reviewed,
            'protein_name': _extract_protein_name(r),
            'length':       r.get('sequence', {}).get('length', 0),
            'organism':     organism.get('scientificName', ''),
            'tax_id':       organism.get('taxonId', 0),
            'flag':         '',
        })
    return hits


def _search_uniprot(organism: str, protein_name: Optional[str], tax_id: Optional[int] = None) -> tuple[list[dict], str]:
    """
    Searches UniProt for metadata only (no FASTA download).
    Returns (hits, query_used).

    Uses taxonomy_id when available (exact match); falls back to organism_name substring.
    Strategy:
      1. organism + protein_name (protein field)
      2. organism + gene name
      3. organism only (last resort — warns user)
    """
    fields = 'accession,reviewed,protein_name,length,organism_name,organism_id'
    org_filter = f'(taxonomy_id:{tax_id})' if tax_id else f'(organism_name:"{organism}")'

    def _query(q: str) -> list:
        r = _http_get(UNIPROT_SEARCH_URL, params={
            'query': q, 'fields': fields, 'format': 'json', 'size': str(_MAX_RESULTS),
        })
        return r.json().get('results', [])

    if protein_name:
        q1 = f'{org_filter} AND (protein_name:"{protein_name}")'
        results = _query(q1)
        if results:
            return _build_hits(results), q1

        q2 = f'{org_filter} AND (gene:"{protein_name}")'
        results = _query(q2)
        if results:
            return _build_hits(results), q2

    q_org = org_filter
    if protein_name:
        console.print(
            f'[yellow]⚠ No results for protein "{protein_name}" — '
            f'falling back to organism-only search.[/yellow]'
        )
    results = _query(q_org)
    return _build_hits(results), q_org


# ── Sort & flag ────────────────────────────────────────────────────────────────

def _sort_and_flag(hits: list[dict]) -> list[dict]:
    """
    Sorts Swiss-Prot first, TrEMBL second.
    Flags polyproteins and fragments based on the global length median
    of ALL results — regardless of Swiss-Prot / TrEMBL status.
    This catches cases where even the reviewed entry is a polyprotein
    (e.g. alphavirus nsP1 stored as a feature of the non-structural polyprotein).
    """
    reviewed   = [h for h in hits if h['reviewed']]
    unreviewed = [h for h in hits if not h['reviewed']]
    sorted_hits = reviewed + unreviewed

    lengths = [h['length'] for h in sorted_hits if h['length'] > 0]
    if len(lengths) >= 2:
        med = median(lengths)
        for h in sorted_hits:
            if h['length'] > med * 2.0:
                h['flag'] = '(Poliproteína?)'
            elif h['length'] > 0 and h['length'] < med * 0.4:
                h['flag'] = '(Fragmento?)'

    return sorted_hits


# ── Display table ──────────────────────────────────────────────────────────────

def _display_uniprot_table(hits: list[dict], organism: str, protein_name: Optional[str]):
    """Renders Rich selection table with Swiss-Prot / TrEMBL status and size flags."""
    reviewed_count   = sum(1 for h in hits if h['reviewed'])
    unreviewed_count = len(hits) - reviewed_count

    table = Table(
        box=box.ROUNDED, show_header=True, header_style='bold white',
        title=(
            f'UniProt — {organism} / {protein_name or "all proteins"}  '
            f'([bold yellow]{reviewed_count}[/bold yellow] Swiss-Prot  '
            f'[dim]{unreviewed_count} TrEMBL[/dim])'
        ),
        title_style='bold cyan',
    )
    table.add_column('#',         no_wrap=True, justify='right', min_width=3)
    table.add_column('Status',    no_wrap=True, min_width=30)
    table.add_column('Accession', no_wrap=True, style='cyan', min_width=12)
    table.add_column('Proteína',  no_wrap=False, min_width=30, max_width=45)
    table.add_column('Tamanho',   no_wrap=True, justify='right', min_width=12)

    for i, hit in enumerate(hits, start=1):
        if hit['reviewed']:
            status_str = '[bold yellow]⭐[/bold yellow] Revisado (Swiss-Prot)'
            row_style  = 'yellow'
        else:
            status_str = '[dim]⚠️  Não Revisado (TrEMBL)[/dim]'
            row_style  = ''

        size_str = f'{hit["length"]} aa'
        if hit.get('flag'):
            size_str += f'  [dim red]{hit["flag"]}[/dim red]'
            row_style = 'dim'

        table.add_row(
            str(i),
            status_str,
            hit['accession'],
            hit['protein_name'][:45],
            size_str,
            style=row_style,
        )

    console.print(table)


# ── Selection prompt ───────────────────────────────────────────────────────────

def _prompt_selection(hits: list[dict], non_interactive: bool = False) -> dict:
    """
    Prompts user to select a hit by row number.
    Non-interactive priority: Swiss-Prot (any) > non-flagged TrEMBL > first hit.
    Interactive (Enter): first non-flagged hit, or first hit if all flagged.
    """
    if non_interactive:
        # Swiss-Prot is human-curated — always prefer it, even if flagged
        best = (
            next((h for h in hits if h['reviewed']), None)
            or next((h for h in hits if not h.get('flag')), hits[0])
        )
        console.print(
            f'[dim]→ Auto-selecionado: {best["accession"]} '
            f'— {best["protein_name"]}[/dim]'
        )
        return best

    best_default = next((h for h in hits if not h.get('flag')), hits[0])

    console.print(
        f'\n[bold]Selecionar sequência[/bold] '
        f'[dim](1–{len(hits)}, ou Enter para aceitar melhor candidato)[/dim]'
    )
    try:
        raw = input('> ').strip()
    except EOFError:
        raw = ''

    if not raw:
        console.print(
            f'[dim]→ Selecionado: {best_default["accession"]} '
            f'— {best_default["protein_name"]}[/dim]'
        )
        return best_default

    try:
        idx = int(raw) - 1
        if 0 <= idx < len(hits):
            return hits[idx]
    except ValueError:
        pass

    console.print('[yellow]Seleção inválida, usando melhor candidato.[/yellow]')
    return best_default


# ── Validation helper ──────────────────────────────────────────────────────────

def _validate_records(records: list) -> tuple[list, list]:
    """Returns (validated_records, rejected_log)."""
    validated = []
    rejected  = []
    for record in records:
        is_valid, reason = is_valid_sequence(record)
        if is_valid:
            validated.append(record)
        else:
            rejected.append({
                'id':          record.id,
                'description': record.description[:80],
                'length':      len(record.seq),
                'reason':      reason,
            })
    return validated, rejected


# ── FASTA download ─────────────────────────────────────────────────────────────

def _download_fasta(accession: str) -> str:
    """Downloads FASTA text for a single UniProt accession."""
    url = UNIPROT_FASTA_URL.format(accession=accession)
    response = _http_get(url)
    return response.text


# ── Local FASTA loader (unchanged) ────────────────────────────────────────────

def _load_local_fasta(local_file_path: Optional[str] = None) -> list:
    """Loads sequences from a local FASTA file."""
    if local_file_path:
        fasta_path = Path(local_file_path)
    else:
        console.print('\n[bold]Path to local FASTA file:[/bold]')
        fasta_path = Path(input('> ').strip())

    if not fasta_path.exists():
        raise FileNotFoundError(f'File not found: {fasta_path}')
    if fasta_path.suffix.lower() not in ('.fasta', '.fa', '.faa', '.fas'):
        console.print('[yellow]Warning: file extension not recognized as FASTA. '
                      'Attempting to parse anyway...[/yellow]')

    records = list(SeqIO.parse(str(fasta_path), 'fasta'))
    if not records:
        raise ValueError(f'No sequences found in file: {fasta_path}')

    console.print(f'[green]Loaded {len(records)} sequence(s) from {fasta_path}[/green]')
    return records


# ── Registry entry builder ─────────────────────────────────────────────────────

def _build_registry_entry(record, hit: Optional[dict]) -> dict:
    """Builds a registry JSON entry from a BioPython record + optional UniProt hit."""
    if hit:
        return {
            'accession_id':       hit['accession'],
            'description':        hit['protein_name'],
            'organism':           hit['organism'],
            'sequence_length_aa': hit['length'],
            'search_source':      'uniprot',
            'reviewed':           hit['reviewed'],
            'tax_id':             hit['tax_id'],
        }
    return {
        'accession_id':       record.id,
        'description':        record.description[:120],
        'organism':           'local',
        'sequence_length_aa': len(record.seq),
        'search_source':      'local',
        'reviewed':           None,
        'tax_id':             None,
    }


# ── Step class ─────────────────────────────────────────────────────────────────

class FetchSequencesStep(BaseTrackStep):
    step_name = 'fetch_sequences'

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
        """
        Full UniProt search + selection + download + validation flow.

        Shows the table once. Lets the user pick a hit (or auto-picks in
        non-interactive mode). If the downloaded sequence fails validation,
        automatically advances to the next valid candidate until one passes
        or all hits are exhausted.

        Returns (selected_records, selected_hit, validated_records, rejected_log).
        """
        resolved_name, tax_id, was_aliased = _normalize_organism(organism_name)
        if was_aliased:
            console.print(
                f'[dim]→ Organismo resolvido: '
                f'[bold]{organism_name}[/bold] → [bold cyan]{resolved_name}[/bold cyan][/dim]'
            )

        console.print('[yellow]Buscando no UniProt...[/yellow]')
        hits, query_used = _search_uniprot(resolved_name, protein_name, tax_id=tax_id)

        if not hits:
            raise ValueError(
                f'Nenhum resultado encontrado no UniProt para '
                f'"{resolved_name}" / "{protein_name}".'
            )

        hits = _sort_and_flag(hits)
        console.print(f'[dim]Query: {query_used}[/dim]')
        _display_uniprot_table(hits, resolved_name, protein_name)

        # Detect non-interactive mode upfront
        non_interactive = _is_non_interactive()

        tried: set[str] = set()

        while True:
            remaining = [h for h in hits if h['accession'] not in tried]
            if not remaining:
                raise ValueError(
                    f'Todos os hits testados falharam na validação para '
                    f'{self.track_id}. Revise os resultados manualmente.'
                )

            selected_hit = _prompt_selection(remaining, non_interactive=non_interactive)
            tried.add(selected_hit['accession'])

            console.print(
                f'\n[yellow]Baixando FASTA: {selected_hit["accession"]} '
                f'({selected_hit["protein_name"][:40]})...[/yellow]'
            )
            fasta_text = _download_fasta(selected_hit['accession'])
            records    = list(SeqIO.parse(io.StringIO(fasta_text), 'fasta'))

            if not records:
                console.print(
                    f'[yellow]FASTA vazio para {selected_hit["accession"]}. '
                    f'Tentando próximo...[/yellow]'
                )
                continue

            console.print(f'[green]Download concluído: {len(records[0].seq)} aa[/green]')

            validated, rejected = _validate_records(records)
            if validated:
                return records, selected_hit, validated, rejected

            reason = rejected[0]['reason'] if rejected else 'motivo desconhecido'
            next_candidates = [h for h in hits if h['accession'] not in tried]
            if next_candidates:
                console.print(
                    f'[yellow]→ {selected_hit["accession"]} rejeitado '
                    f'({reason}). Tentando próximo: '
                    f'{next_candidates[0]["accession"]}...[/yellow]'
                )
            else:
                raise ValueError(
                    f'Nenhuma sequência válida encontrada para {self.track_id}. '
                    f'Último erro: {reason}'
                )


# ── Non-interactive detection ──────────────────────────────────────────────────

def _is_non_interactive() -> bool:
    """Returns True if stdin is not a TTY (piped, redirected, or CI environment)."""
    import sys
    return not sys.stdin.isatty()
