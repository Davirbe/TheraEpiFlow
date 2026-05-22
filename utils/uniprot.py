"""
Shared UniProt helpers — Chain-feature matching and entry fetch.

Used by both fetch_sequences (slice the reference out of a polyprotein) and
search_variants (slice variant candidates the same way). The matcher works on any
dict shaped like a UniProt JSON record — a full entry (`/{accession}.json`) or a
single search result returned with `fields=...,ft_chain` — since both expose
`features[]` and `sequence.length`.
"""

import difflib
import re
from typing import Optional

from utils.http import http_get

UNIPROT_ENTRY_JSON_URL = 'https://rest.uniprot.org/uniprotkb/{accession}.json'

# Minimum score to accept a Chain as the target protein. High on purpose: trusts
# exact / whole-word-token matches and rejects fuzzy near-misses like NS1↔NS3.
CHAIN_MATCH_CUTOFF = 0.85

# Flaviviruses store mature proteins as Chain features whose description is the long
# form "Non-structural protein N" (e.g. NS1, NS2A, NS4B), while users type the standard
# abbreviation "NSN". This pattern recovers the abbreviation so the two forms match.
_NONSTRUCTURAL_PROTEIN_PATTERN = re.compile(r'non[\s-]*structural protein\s*([0-9]+[a-z]?)')


def _nonstructural_protein_alias(description_lower: str) -> Optional[str]:
    """Returns the 'nsN' abbreviation for a 'Non-structural protein N' description, else None."""
    alias_match = _NONSTRUCTURAL_PROTEIN_PATTERN.search(description_lower)
    return f"ns{alias_match.group(1)}" if alias_match else None


def fetch_uniprot_entry_json(accession: str) -> dict:
    """Downloads the full UniProt entry JSON for a single accession (for Chain features)."""
    return http_get(UNIPROT_ENTRY_JSON_URL.format(accession=accession)).json()


def score_chain_match(protein_name: str, description: str) -> float:
    """Fuzzy score (0-1) between a requested protein_name and a Chain description.
    Whole-word token match (e.g. 'E' in 'Envelope protein E') scores high; otherwise
    falls back to difflib ratio over the full string and over each token."""
    name = protein_name.strip().lower()
    desc = description.strip().lower()
    if not name or not desc:
        return 0.0
    if name == desc:
        return 1.0
    desc_tokens = desc.replace('-', ' ').split()
    if name in desc_tokens:
        return 0.95
    # "Non-structural protein 1" ≡ "NS1" (and 2A, 4B, …) — flavivirus convention.
    nonstructural_alias = _nonstructural_protein_alias(desc)
    if nonstructural_alias and name.replace(' ', '') == nonstructural_alias:
        return 0.95
    best_token_ratio = max(
        (difflib.SequenceMatcher(None, name, token).ratio() for token in desc_tokens),
        default=0.0,
    )
    full_ratio = difflib.SequenceMatcher(None, name, desc).ratio()
    return max(full_ratio, best_token_ratio)


def find_chain_for_protein(uniprot_record: dict, protein_name: str) -> Optional[dict]:
    """Finds the mature-chain region matching protein_name in a UniProt record's features.

    Accepts a full entry JSON or a search result (both expose `features[]` and
    `sequence.length`). Scans features[type=='Chain'], skips the full 'Genome polyprotein'
    (and any chain spanning the whole sequence), fuzzy-matches protein_name to each Chain
    description, and returns {'start', 'end', 'description', 'score'} for the best match
    above CHAIN_MATCH_CUTOFF — or None when nothing matches (caller keeps the full sequence).
    Coordinates are 1-based inclusive, as UniProt reports them."""
    if not protein_name:
        return None

    total_length = uniprot_record.get('sequence', {}).get('length', 0)
    best_match: Optional[dict] = None

    for feature in uniprot_record.get('features', []):
        if feature.get('type') != 'Chain':
            continue
        description = feature.get('description', '') or ''
        if 'polyprotein' in description.lower():
            continue

        location = feature.get('location', {})
        start = location.get('start', {}).get('value')
        end   = location.get('end', {}).get('value')
        if not isinstance(start, int) or not isinstance(end, int) or start > end:
            continue
        # Skip a chain spanning the whole sequence — that is the polyprotein itself.
        if total_length and (end - start + 1) >= total_length:
            continue

        score = score_chain_match(protein_name, description)
        if best_match is None or score > best_match['score']:
            best_match = {'start': start, 'end': end, 'description': description, 'score': score}

    if best_match and best_match['score'] >= CHAIN_MATCH_CUTOFF:
        return best_match
    return None
