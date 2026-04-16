"""
Utility functions for searching and downloading protein sequences from NCBI GenBank.

Used by:
  step01_fetch_sequences — download reference sequences for each project track
  step08_search_variants — download variant sequences for conservation analysis

Entrez.email must be set by the caller before using any search function:
  from Bio import Entrez
  Entrez.email = project_config["entrez_email"]
"""

import json
from collections import Counter
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ── Constants ─────────────────────────────────────────────────────────────────

REFSEQ_ACCESSION_PREFIXES = ('YP_', 'NP_', 'XP_', 'WP_', 'AP_')

# Polyproteins contain multiple viral proteins joined in one entry.
# We exclude them because they don't represent a single protein of interest.
MINIMUM_POLYPROTEIN_LENGTH_AA = 1000


# ── Record classification ──────────────────────────────────────────────────────

def record_is_polyprotein(genbank_record: SeqRecord) -> bool:
    """
    Returns True if the record is a polyprotein — a concatenated multi-protein
    entry common in flaviviruses (e.g. Zika, Dengue).

    Detection: description contains 'polyprotein' AND sequence length > 1000 aa.
    Both conditions required to avoid false positives on short partial entries.
    """
    description_lowercase = genbank_record.description.lower()
    sequence_exceeds_polyprotein_threshold = len(genbank_record.seq) > MINIMUM_POLYPROTEIN_LENGTH_AA
    return 'polyprotein' in description_lowercase and sequence_exceeds_polyprotein_threshold


def record_is_refseq(genbank_record: SeqRecord) -> bool:
    """
    Returns True if the record is a curated NCBI RefSeq entry.

    RefSeq accessions are manually reviewed and use standard prefixes:
      YP_ = viral/prokaryotic proteins
      NP_ = RefSeq protein (well-characterized)
      XP_ = predicted RefSeq protein
      WP_ = non-redundant prokaryotic protein
      AP_ = protein from alternative assembly
    """
    return genbank_record.id.startswith(REFSEQ_ACCESSION_PREFIXES)


# ── Metadata extraction ────────────────────────────────────────────────────────

def extract_source_qualifiers(genbank_record: SeqRecord) -> dict:
    """
    Extracts biological metadata from the 'source' feature of a GenBank record.

    Returns a flat dictionary with the following keys:
      strain             — viral strain name (e.g. 'Natal RGN', 'MR 766')
      isolate            — clinical or environmental isolate identifier
      host               — host organism (e.g. 'Homo sapiens', 'Aedes aegypti')
      geographic_location — country or region of collection
      collection_date    — date the sample was collected
      isolation_source   — biological material (e.g. 'serum', 'brain tissue')

    All values default to 'N/A' when the qualifier is absent.
    Note: NCBI is migrating 'country' → 'geo_loc_name'; both are checked.
    """
    source_features = [feature for feature in genbank_record.features
                       if feature.type == 'source']
    if not source_features:
        return {
            'strain': 'N/A', 'isolate': 'N/A', 'host': 'N/A',
            'geographic_location': 'N/A', 'collection_date': 'N/A',
            'isolation_source': 'N/A',
        }

    raw_qualifiers = source_features[0].qualifiers

    # geo_loc_name is the new NCBI standard; country is the legacy field name
    geographic_location_value = raw_qualifiers.get(
        'geo_loc_name', raw_qualifiers.get('country', ['N/A'])
    )[0]

    return {
        'strain':              raw_qualifiers.get('strain',           ['N/A'])[0],
        'isolate':             raw_qualifiers.get('isolate',          ['N/A'])[0],
        'host':                raw_qualifiers.get('host',             ['N/A'])[0],
        'geographic_location': geographic_location_value,
        'collection_date':     raw_qualifiers.get('collection_date',  ['N/A'])[0],
        'isolation_source':    raw_qualifiers.get('isolation_source', ['N/A'])[0],
    }


# ── Polyprotein extraction ────────────────────────────────────────────────────

def extract_protein_from_polyprotein(
    polyprotein_record: SeqRecord,
    target_protein_name: str,
) -> Optional[SeqRecord]:
    """
    Extracts a specific mature peptide from a polyprotein GenBank record.

    Polyproteins (common in flaviviruses like Zika, Dengue, Yellow Fever) encode
    multiple proteins in a single continuous sequence. GenBank annotates each
    protein boundary using 'mat_peptide' features with a /product qualifier.

    Example for ZIKV polyprotein:
      mat_peptide  291..790  /product="envelope protein E"
      → extracts amino acids 291 to 790 as a standalone SeqRecord

    Matching is case-insensitive and partial:
      target_protein_name='envelope' matches /product="envelope protein E"

    Returns a new SeqRecord with the extracted sequence, or None if no matching
    mat_peptide feature is found in the polyprotein record.
    """
    target_name_lowercase = target_protein_name.lower()

    # Feature types that can annotate individual proteins within a polyprotein
    protein_feature_types = ('mat_peptide', 'Protein', 'Region')

    for genbank_feature in polyprotein_record.features:
        if genbank_feature.type not in protein_feature_types:
            continue

        # Check /product, /name, and /note qualifiers for a match
        product_value = genbank_feature.qualifiers.get('product', [''])[0].lower()
        name_value    = genbank_feature.qualifiers.get('name',    [''])[0].lower()
        note_value    = genbank_feature.qualifiers.get('note',    [''])[0].lower()

        name_matches_target = (
            target_name_lowercase in product_value
            or target_name_lowercase in name_value
            or target_name_lowercase in note_value
        )

        if not name_matches_target:
            continue

        extracted_sequence = genbank_feature.extract(polyprotein_record.seq)
        source_qualifiers  = extract_source_qualifiers(polyprotein_record)

        extracted_record = SeqRecord(
            seq=extracted_sequence,
            id=polyprotein_record.id,
            description=(
                f'{genbank_feature.qualifiers.get("product", ["extracted protein"])[0]}'
                f' [extracted from polyprotein] [{polyprotein_record.annotations.get("organism", "")}]'
            ),
            annotations={
                'organism': polyprotein_record.annotations.get('organism', 'N/A'),
                'extracted_from_polyprotein': True,
                'source_accession': polyprotein_record.id,
            },
        )

        # Carry over source qualifiers so metadata is preserved after extraction
        extracted_record.features = [f for f in polyprotein_record.features
                                     if f.type == 'source']

        return extracted_record

    return None


# ── Reference suggestion ───────────────────────────────────────────────────────

def suggest_reference_sequence(candidate_records: list) -> tuple:
    """
    Suggests the most appropriate reference sequence from a list of candidates.

    This function does NOT make the final selection — it returns a suggestion
    that step01 highlights in the display table. The user always decides.

    Selection priority:
      1. RefSeq non-polyprotein entries → longest sequence (most complete)
      2. No RefSeq available → modal length across all non-polyprotein records
         (the most common length = the canonical reference length for that protein),
         then the first record at that length

    Returns:
      (suggested_record, suggestion_reason_text)
    """
    non_polyprotein_records = [record for record in candidate_records
                               if not record_is_polyprotein(record)]
    if not non_polyprotein_records:
        non_polyprotein_records = candidate_records  # fallback: all are polyproteins

    refseq_non_polyprotein = [record for record in non_polyprotein_records
                              if record_is_refseq(record)]

    if refseq_non_polyprotein:
        suggested_record = max(refseq_non_polyprotein, key=lambda record: len(record.seq))
        suggestion_reason = 'RefSeq — longest curated entry'
        return suggested_record, suggestion_reason

    # No RefSeq: use modal length as proxy for canonical reference length
    sequence_lengths = [len(record.seq) for record in non_polyprotein_records]
    canonical_length = Counter(sequence_lengths).most_common(1)[0][0]
    records_at_canonical_length = [record for record in non_polyprotein_records
                                   if len(record.seq) == canonical_length]
    suggested_record = records_at_canonical_length[0]
    suggestion_reason = f'GenBank — canonical length ({canonical_length} aa, most common)'
    return suggested_record, suggestion_reason


# ── NCBI search — phase 1: get IDs (no sequence download) ────────────────────

def search_ncbi_protein_ids(
    organism_name: str,
    protein_name: Optional[str],
    target_host: str,
    include_polyproteins: bool = False,
    max_ids_to_return: int = 200,
) -> tuple:
    """
    Searches NCBI protein database and returns accession IDs only — no sequences
    are downloaded. This is the fast discovery phase.

    Strategy:
      1. RefSeq IDs first (curated, preferred)
      2. General GenBank IDs fill remaining slots
      3. Optionally excludes polyproteins from the query itself

    Args:
      organism_name        — e.g. 'Human papillomavirus 16', 'Zika virus'
      protein_name         — e.g. 'E6', 'envelope protein E'. None = organism only
      target_host          — e.g. 'Homo sapiens'
      include_polyproteins — if False, adds NOT polyprotein[Title] to query
      max_ids_to_return    — maximum number of IDs to return

    Returns:
      (total_available_in_ncbi, ordered_accession_ids)
      IDs are ordered: RefSeq first, then general GenBank.
    """
    base_search_term = f'"{organism_name}"[Organism]'
    if protein_name:
        base_search_term += f' AND "{protein_name}"[Protein Name]'
    if not include_polyproteins:
        base_search_term += ' NOT polyprotein[Title]'

    # Count total available (no download)
    count_handle = Entrez.esearch(db='protein', term=base_search_term, retmax=0)
    total_available_in_ncbi = int(Entrez.read(count_handle)['Count'])
    count_handle.close()

    if total_available_in_ncbi == 0:
        return 0, []

    # Fetch RefSeq IDs first
    refseq_handle = Entrez.esearch(
        db='protein',
        term=base_search_term + ' AND RefSeq[Filter]',
        retmax=max_ids_to_return,
    )
    refseq_ids = Entrez.read(refseq_handle)['IdList']
    refseq_handle.close()

    # Fill remaining slots with general GenBank IDs
    remaining_slots = max_ids_to_return - len(refseq_ids)
    general_ids = []
    if remaining_slots > 0:
        general_handle = Entrez.esearch(
            db='protein', term=base_search_term, retmax=max_ids_to_return,
        )
        all_general_ids = Entrez.read(general_handle)['IdList']
        general_handle.close()
        refseq_ids_set = set(refseq_ids)
        general_ids = [
            accession_id for accession_id in all_general_ids
            if accession_id not in refseq_ids_set
        ][:remaining_slots]

    ordered_accession_ids = refseq_ids + general_ids
    return total_available_in_ncbi, ordered_accession_ids


# ── NCBI fetch — phase 2: download records by ID ─────────────────────────────

def fetch_records_by_accession_ids(
    accession_ids: list,
    target_host: str,
) -> list:
    """
    Downloads full GenBank records for a list of accession IDs.
    Filters by target host: keeps records that match or have no host declared.

    This is the download phase — call only for IDs the user actually selected,
    not for the full discovery set.

    Returns list of SeqRecord with full sequence and metadata.
    """
    if not accession_ids:
        return []

    records_handle = Entrez.efetch(
        db='protein',
        id=','.join(accession_ids),
        rettype='gb',
        retmode='text',
    )
    downloaded_records = list(SeqIO.parse(records_handle, 'genbank'))
    records_handle.close()

    target_host_lowercase = target_host.lower()
    host_filtered_records = []
    for downloaded_record in downloaded_records:
        record_qualifiers  = extract_source_qualifiers(downloaded_record)
        declared_host      = record_qualifiers.get('host', 'N/A')
        host_matches       = target_host_lowercase in declared_host.lower()
        host_not_declared  = declared_host == 'N/A'
        if host_matches or host_not_declared:
            host_filtered_records.append(downloaded_record)

    return host_filtered_records


# ── Output ─────────────────────────────────────────────────────────────────────

def save_sequences_as_fasta(
    selected_records: list,
    output_fasta_path: Path,
) -> Path:
    """
    Saves a list of selected GenBank records as a FASTA file.
    Creates parent directories if they do not exist.
    Returns the path to the saved file.
    """
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(selected_records, str(output_fasta_path), 'fasta')
    return output_fasta_path


def build_sequence_registry(
    selected_records: list,
    output_json_path: Path,
) -> Path:
    """
    Saves a JSON registry of downloaded sequences and their metadata.

    Purpose: step08 (variant search) reads this file to exclude the original
    reference sequences from variant results — preventing the reference from
    appearing as its own variant during conservation analysis.

    Saved fields per sequence:
      accession_id, description, organism, sequence_length_aa,
      strain, isolate, host, geographic_location, collection_date.
    """
    registry_entries = []
    for selected_record in selected_records:
        source_qualifiers = extract_source_qualifiers(selected_record)
        registry_entries.append({
            'accession_id':        selected_record.id,
            'description':         selected_record.description,
            'organism':            selected_record.annotations.get('organism', 'N/A'),
            'sequence_length_aa':  len(selected_record.seq),
            'strain':              source_qualifiers.get('strain',              'N/A'),
            'isolate':             source_qualifiers.get('isolate',             'N/A'),
            'host':                source_qualifiers.get('host',                'N/A'),
            'geographic_location': source_qualifiers.get('geographic_location', 'N/A'),
            'collection_date':     source_qualifiers.get('collection_date',     'N/A'),
        })

    output_json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json_path, 'w', encoding='utf-8') as registry_file:
        json.dump(
            {
                'total_sequences': len(registry_entries),
                'sequences':       registry_entries,
            },
            registry_file,
            indent=2,
            ensure_ascii=False,
        )
    return output_json_path
