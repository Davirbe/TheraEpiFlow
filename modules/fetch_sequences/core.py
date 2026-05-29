"""
UniProt search and sequence-fetch logic for fetch_sequences: organism
normalization, the search/scoring helpers and record validation. Pure-ish
domain functions (network GET included) — no Rich tables, no prompts.
"""

import difflib
from statistics import median
from typing import Optional

from utils.console import console
from utils.fasta_utils import is_valid_sequence
from utils.http import http_get

UNIPROT_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
UNIPROT_FASTA_URL  = 'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
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


# ── Input normalization ────────────────────────────────────────────────────────

def _normalize_organism(raw_name: str) -> tuple[str, Optional[int], bool]:
    """Resolves organism abbreviation/typo to canonical scientific name + tax_id.
    Returns (resolved_name, tax_id_or_None, was_aliased).
    Order: exact alias → fuzzy alias → canonical scientific name in _NAME_TO_TAX_ID."""
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
    """Searches UniProt for metadata only (no FASTA download); returns (hits, query_used).
    Uses taxonomy_id when available; falls back to organism_name substring. Fallback chain:
      1. organism + protein_name (protein field)
      2. organism + gene name
      3. organism only (last resort — warns user)
    """
    fields = 'accession,reviewed,protein_name,length,organism_name,organism_id'
    org_filter = f'(taxonomy_id:{tax_id})' if tax_id else f'(organism_name:"{organism}")'

    def _query(q: str) -> list:
        r = http_get(UNIPROT_SEARCH_URL, params={
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
            f'[yellow]⚠ No results for protein "{protein_name}"; '
            f'falling back to organism-only search.[/yellow]'
        )
    results = _query(q_org)
    return _build_hits(results), q_org


# ── Sort & flag ────────────────────────────────────────────────────────────────

def _sort_and_flag(hits: list[dict]) -> list[dict]:
    """Sorts Swiss-Prot first, TrEMBL second.
    Flags polyproteins (>2× median length) and fragments (<0.4× median) against the
    global length median — catches reviewed-but-polyprotein entries like alphavirus nsP1."""
    reviewed   = [h for h in hits if h['reviewed']]
    unreviewed = [h for h in hits if not h['reviewed']]
    sorted_hits = reviewed + unreviewed

    lengths = [h['length'] for h in sorted_hits if h['length'] > 0]
    if len(lengths) >= 2:
        med = median(lengths)
        for h in sorted_hits:
            if h['length'] > med * 2.0:
                h['flag'] = '(Polyprotein?)'
            elif h['length'] > 0 and h['length'] < med * 0.4:
                h['flag'] = '(Fragment?)'

    return sorted_hits


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


# ── Registry entry builder ─────────────────────────────────────────────────────

def _build_registry_entry(record, hit: Optional[dict]) -> dict:
    """Builds a registry JSON entry from a BioPython record + optional UniProt hit."""
    if hit:
        registry_entry = {
            'accession_id':       hit['accession'],
            'description':        hit['protein_name'],
            'organism':           hit['organism'],
            'sequence_length_aa': hit['length'],
            'search_source':      'uniprot',
            'reviewed':           hit['reviewed'],
            'tax_id':             hit['tax_id'],
        }
        if hit.get('chain_slice'):
            registry_entry['chain_slice'] = hit['chain_slice']
        return registry_entry
    return {
        'accession_id':       record.id,
        'description':        record.description[:120],
        'organism':           'local',
        'sequence_length_aa': len(record.seq),
        'search_source':      'local',
        'reviewed':           None,
        'tax_id':             None,
    }


