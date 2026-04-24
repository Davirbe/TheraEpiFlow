"""
Utility functions for searching and downloading protein sequences from NCBI GenBank.

Used by:
  step01_fetch_sequences — download reference sequences for each project track
  step08_search_variants — download variant sequences for conservation analysis

Entrez.email must be set by the caller before using any search function:
  from Bio import Entrez
  Entrez.email = project_config["entrez_email"]

Design note — unified search:
  Many viruses (CHIKV, ZIKV, DENV, HCV, alphaviruses, flaviviruses) deposit
  their individual proteins ONLY as mat_peptide features inside a polyprotein
  record. Others (HPV, SARS-CoV-2, HIV-1) deposit each protein as its own
  standalone entry.

  Instead of asking the user up front which case applies, search_proteins_comprehensive()
  runs BOTH strategies in parallel and returns a merged list. Each record is
  annotated with annotations['search_source'] = 'direct' | 'extracted' so the
  caller can show the provenance in the selection table.
"""

import json
from collections import Counter
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from utils.retry_helpers import retry_network_call


# ── Constants ─────────────────────────────────────────────────────────────────

REFSEQ_ACCESSION_PREFIXES = ('YP_', 'NP_', 'XP_', 'WP_', 'AP_')

# Polyproteins contain multiple viral proteins joined in one entry.
# Detection is done after fetching — see record_is_polyprotein() below.
MINIMUM_POLYPROTEIN_LENGTH_AA = 1000

# Default ceilings for the comprehensive search — keeps NCBI traffic reasonable.
DEFAULT_MAX_DIRECT_IDS_TO_FETCH     = 50
DEFAULT_MAX_POLYPROTEINS_TO_FETCH   = 10   # we only need a few polyproteins
                                            # to extract mat_peptide candidates
DEFAULT_MAX_RECORDS_TO_RETURN       = 50


# ── Record classification ──────────────────────────────────────────────────────

def record_is_polyprotein(genbank_record: SeqRecord) -> bool:
    """
    Returns True if the record looks like a polyprotein — a concatenated
    multi-protein entry common in flaviviruses (Zika, Dengue) and
    alphaviruses (Chikungunya, Sindbis).

    Detection uses two independent signals:
      1. Description contains 'polyprotein' AND sequence length > 1000 aa
      2. Record has at least one mat_peptide feature
    Either signal is enough — NCBI is inconsistent about labeling.
    """
    description_lowercase = genbank_record.description.lower()
    sequence_is_long = len(genbank_record.seq) > MINIMUM_POLYPROTEIN_LENGTH_AA
    has_polyprotein_label = (
        'polyprotein' in description_lowercase and sequence_is_long
    )
    has_mat_peptide_features = any(
        feature.type == 'mat_peptide'
        for feature in genbank_record.features
    )
    return has_polyprotein_label or has_mat_peptide_features


def record_is_refseq(genbank_record: SeqRecord) -> bool:
    """
    Returns True if the record is a curated NCBI RefSeq entry.

    RefSeq accessions use standard prefixes:
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

    Polyproteins (flaviviruses, alphaviruses, HCV, HIV gag-pol) encode multiple
    proteins in a single continuous sequence. GenBank annotates each protein
    boundary with 'mat_peptide' features that carry a /product qualifier.

    Example for Chikungunya structural polyprotein:
      mat_peptide  782..1249   /product="envelope glycoprotein E1"
      → extracts residues 782..1249 as a standalone SeqRecord

    Matching is case-insensitive and partial — 'E1' matches 'glycoprotein E1'.

    Returns a new SeqRecord or None if no matching mat_peptide is found.
    """
    target_protein_name_lowercase = target_protein_name.lower().strip()

    feature_types_that_can_label_proteins = ('mat_peptide', 'Protein', 'Region')

    for candidate_feature in polyprotein_record.features:
        if candidate_feature.type not in feature_types_that_can_label_proteins:
            continue

        product_qualifier_text = candidate_feature.qualifiers.get('product', [''])[0].lower()
        name_qualifier_text    = candidate_feature.qualifiers.get('name',    [''])[0].lower()
        note_qualifier_text    = candidate_feature.qualifiers.get('note',    [''])[0].lower()

        feature_labels_this_protein = (
            target_protein_name_lowercase in product_qualifier_text
            or target_protein_name_lowercase in name_qualifier_text
            or target_protein_name_lowercase in note_qualifier_text
        )
        if not feature_labels_this_protein:
            continue

        extracted_peptide_sequence = candidate_feature.extract(polyprotein_record.seq)

        annotated_product_name = candidate_feature.qualifiers.get(
            'product', ['extracted protein']
        )[0]

        extracted_mat_peptide_record = SeqRecord(
            seq=extracted_peptide_sequence,
            id=polyprotein_record.id,
            description=(
                f'{annotated_product_name} [extracted from polyprotein] '
                f'[{polyprotein_record.annotations.get("organism", "")}]'
            ),
            annotations={
                'organism': polyprotein_record.annotations.get('organism', 'N/A'),
                'extracted_from_polyprotein': True,
                'source_accession':           polyprotein_record.id,
                'search_source':              'extracted',
            },
        )

        # Preserve source-feature metadata (strain, host, location, etc.)
        extracted_mat_peptide_record.features = [
            polyprotein_feature for polyprotein_feature in polyprotein_record.features
            if polyprotein_feature.type == 'source'
        ]

        return extracted_mat_peptide_record

    return None


# ── Reference suggestion ───────────────────────────────────────────────────────

def suggest_reference_sequence(candidate_records: list) -> tuple:
    """
    Suggests the most appropriate reference from a list of candidate records.
    The user always decides — this is only a hint highlighted in the table.

    Priority:
      1. RefSeq + 'direct' (standalone RefSeq is usually best)
      2. RefSeq + 'extracted' (extracted from a curated polyprotein)
      3. GenBank + modal length (the most common length = canonical)

    Returns (suggested_record, human_readable_reason).
    """
    def _search_source(record: SeqRecord) -> str:
        return record.annotations.get('search_source', 'direct')

    # Tier 1 — RefSeq direct standalone
    refseq_direct = [
        record for record in candidate_records
        if record_is_refseq(record)
        and _search_source(record) == 'direct'
        and not record_is_polyprotein(record)
    ]
    if refseq_direct:
        suggested_record = max(refseq_direct, key=lambda record: len(record.seq))
        return suggested_record, 'RefSeq standalone — curated direct entry'

    # Tier 2 — RefSeq extracted from polyprotein
    refseq_extracted = [
        record for record in candidate_records
        if record_is_refseq(record) and _search_source(record) == 'extracted'
    ]
    if refseq_extracted:
        suggested_record = max(refseq_extracted, key=lambda record: len(record.seq))
        return suggested_record, 'RefSeq — extracted from curated polyprotein'

    # Tier 3 — modal length among non-polyproteins
    non_polyprotein_records = [record for record in candidate_records
                               if not record_is_polyprotein(record)]
    if not non_polyprotein_records:
        non_polyprotein_records = candidate_records

    sequence_lengths = [len(record.seq) for record in non_polyprotein_records]
    if not sequence_lengths:
        # Nothing to suggest — fall back to first record
        return candidate_records[0], 'first candidate'

    canonical_length = Counter(sequence_lengths).most_common(1)[0][0]
    records_at_canonical_length = [record for record in non_polyprotein_records
                                   if len(record.seq) == canonical_length]
    return (
        records_at_canonical_length[0],
        f'GenBank — canonical length ({canonical_length} aa, most common)',
    )


# ── Low-level Entrez helpers (private) ────────────────────────────────────────

def _esearch_protein_database_count_and_ids(
    entrez_search_query: str,
    maximum_ids_to_return: int,
) -> tuple:
    """
    Runs Entrez.esearch against the NCBI Protein database once and returns
    a (total_count_matching_query, list_of_accession_ids) tuple.
    Wrapped in retry for transient network failures.
    """
    def _perform_esearch_request():
        esearch_response_handle = Entrez.esearch(
            db='protein', term=entrez_search_query, retmax=maximum_ids_to_return,
        )
        esearch_parsed_response = Entrez.read(esearch_response_handle)
        esearch_response_handle.close()
        return esearch_parsed_response

    esearch_parsed_response = retry_network_call(
        callable_to_execute=_perform_esearch_request,
        operation_description=f'NCBI esearch ({entrez_search_query[:60]}...)',
    )
    total_count_matching_query  = int(esearch_parsed_response['Count'])
    matching_accession_ids      = list(esearch_parsed_response['IdList'])
    return total_count_matching_query, matching_accession_ids


def _efetch_genbank_records_by_accession_ids(accession_ids_to_fetch: list) -> list:
    """
    Downloads full GenBank records for a list of accession IDs.
    Wrapped in retry. Returns [] for an empty id list.
    """
    if not accession_ids_to_fetch:
        return []

    def _perform_efetch_request():
        efetch_response_handle = Entrez.efetch(
            db='protein', id=','.join(accession_ids_to_fetch),
            rettype='gb', retmode='text',
        )
        downloaded_genbank_records = list(SeqIO.parse(efetch_response_handle, 'genbank'))
        efetch_response_handle.close()
        return downloaded_genbank_records

    return retry_network_call(
        callable_to_execute=_perform_efetch_request,
        operation_description=f'NCBI efetch ({len(accession_ids_to_fetch)} records)',
    )


def _filter_records_by_declared_host(
    records_to_filter: list, target_host: str,
) -> list:
    """
    Keeps records where the 'source' feature's /host qualifier matches target_host
    (case-insensitive, partial match), OR where no host is declared at all.
    """
    target_host_lowercase = target_host.lower()
    host_matching_records = []
    for downloaded_genbank_record in records_to_filter:
        source_qualifiers           = extract_source_qualifiers(downloaded_genbank_record)
        declared_host_text          = source_qualifiers.get('host', 'N/A')
        host_text_matches_target    = target_host_lowercase in declared_host_text.lower()
        host_qualifier_is_absent    = declared_host_text == 'N/A'
        if host_text_matches_target or host_qualifier_is_absent:
            host_matching_records.append(downloaded_genbank_record)
    return host_matching_records


def _build_strict_direct_term(organism_name: str, protein_name: str) -> str:
    """
    Strict direct search: '{protein}[Protein Name]' WITHOUT quotes.

    NCBI's [Protein Name] field indexes the /product qualifier from CDS/mat_peptide
    features. Without quotes, NCBI does token-level match, so 'E1[Protein Name]'
    matches records with product='E1' OR product='envelope glycoprotein E1' —
    but NOT random occurrences of 'E1' in strain codes, notes, or titles.

    This is the sweet spot: flexible enough to catch compound protein names,
    tight enough to reject noise.
    """
    return f'"{organism_name}"[Organism] AND {protein_name}[Protein Name]'


def _build_loose_direct_term(organism_name: str, protein_name: str) -> str:
    """
    Fallback direct search: '{protein}[All Fields]'. Used only when the
    strict [Protein Name] search returns zero results — typically when the
    user typed a synonym that isn't the annotated product name.
    """
    return f'"{organism_name}"[Organism] AND {protein_name}[All Fields]'


def _build_polyprotein_search_term(organism_name: str) -> str:
    """
    Polyproteins are reliably indexed under [Protein Name]='polyprotein' in NCBI.
    This query hits the canonical polyprotein records for the organism.
    """
    return f'"{organism_name}"[Organism] AND polyprotein[Protein Name]'


def _reorder_accession_ids_with_refseq_first(
    all_matching_accession_ids: list,
    refseq_only_accession_ids: list,
) -> list:
    """
    Returns all_matching_accession_ids reordered so that RefSeq accessions
    appear first, followed by non-RefSeq accessions, preserving each group's
    original relative order. Duplicates between the two inputs are collapsed.
    """
    refseq_only_accession_ids_set = set(refseq_only_accession_ids)
    refseq_ids_in_matches = [
        accession_id for accession_id in all_matching_accession_ids
        if accession_id in refseq_only_accession_ids_set
    ]
    non_refseq_ids_in_matches = [
        accession_id for accession_id in all_matching_accession_ids
        if accession_id not in refseq_only_accession_ids_set
    ]
    return refseq_ids_in_matches + non_refseq_ids_in_matches


# ── Public API: unified comprehensive search ──────────────────────────────────

def search_proteins_comprehensive(
    organism_name: str,
    protein_name: Optional[str],
    target_host: str,
    max_records_to_return: int = DEFAULT_MAX_RECORDS_TO_RETURN,
) -> dict:
    """
    Runs two parallel search strategies against NCBI Protein and returns
    a merged, deduplicated, host-filtered list of SeqRecord ready to display.

    Strategy A — Direct search (for proteins deposited as standalone entries):
        '"{organism}"[Organism] AND {protein_name}[All Fields]'

    Strategy B — Polyprotein extraction (for proteins annotated as mat_peptide):
        '"{organism}"[Organism] AND polyprotein[Protein Name]'
        Fetches top polyprotein records, then extract_protein_from_polyprotein()
        pulls out the mat_peptide whose /product matches protein_name.

    Each returned record carries annotations['search_source']:
      'direct'    — downloaded as an independent protein entry
      'extracted' — carved out of a polyprotein's mat_peptide feature

    Args:
      organism_name          — e.g. 'Chikungunya virus', 'Zika virus'
      protein_name           — e.g. 'E1', 'envelope protein E'. None = organism only
      target_host            — e.g. 'Homo sapiens'. Used to filter records
      max_records_to_return  — cap on the final merged list size (default 50)

    Returns:
      {
        'records': list[SeqRecord],
        'total_direct_found':       int,  # how many IDs exist in NCBI (strategy A)
        'total_polyprotein_found':  int,  # how many polyprotein IDs exist
        'extracted_from_polyprotein_count': int,
        'direct_search_tier': str,  # 'strict' | 'loose' | 'organism_only' | 'none'
        'direct_query':      str,   # the actual Entrez term used
        'polyprotein_query': str,
      }
    """
    # ── Strategy A — Direct search (two-tier) ─────────────────────────────────
    direct_search_tier_used:             str  = 'none'
    direct_search_query_used:            str  = f'"{organism_name}"[Organism]'
    total_direct_hits_in_ncbi:           int  = 0
    direct_search_accession_ids:         list = []

    if protein_name:
        # Tier 1 — strict match on /product via [Protein Name]
        strict_protein_name_query = _build_strict_direct_term(organism_name, protein_name)
        total_strict_hits, strict_matching_accession_ids = (
            _esearch_protein_database_count_and_ids(
                entrez_search_query=strict_protein_name_query,
                maximum_ids_to_return=DEFAULT_MAX_DIRECT_IDS_TO_FETCH,
            )
        )

        if total_strict_hits > 0:
            direct_search_tier_used     = 'strict'
            direct_search_query_used    = strict_protein_name_query
            total_direct_hits_in_ncbi   = total_strict_hits
            direct_search_accession_ids = strict_matching_accession_ids
        else:
            # Tier 2 — fall back to [All Fields] (broader, noisier)
            loose_all_fields_query = _build_loose_direct_term(organism_name, protein_name)
            total_loose_hits, loose_matching_accession_ids = (
                _esearch_protein_database_count_and_ids(
                    entrez_search_query=loose_all_fields_query,
                    maximum_ids_to_return=DEFAULT_MAX_DIRECT_IDS_TO_FETCH,
                )
            )
            if total_loose_hits > 0:
                direct_search_tier_used     = 'loose'
                direct_search_query_used    = loose_all_fields_query
                total_direct_hits_in_ncbi   = total_loose_hits
                direct_search_accession_ids = loose_matching_accession_ids
    else:
        # Organism-only search (user chose 'organism_only' fallback in step01)
        organism_only_query = f'"{organism_name}"[Organism]'
        total_organism_only_hits, organism_only_matching_accession_ids = (
            _esearch_protein_database_count_and_ids(
                entrez_search_query=organism_only_query,
                maximum_ids_to_return=DEFAULT_MAX_DIRECT_IDS_TO_FETCH,
            )
        )
        if total_organism_only_hits > 0:
            direct_search_tier_used     = 'organism_only'
            direct_search_query_used    = organism_only_query
            total_direct_hits_in_ncbi   = total_organism_only_hits
            direct_search_accession_ids = organism_only_matching_accession_ids

    # Reorder direct hits RefSeq first (second esearch limited to RefSeq)
    if direct_search_accession_ids:
        _, direct_search_refseq_only_accession_ids = (
            _esearch_protein_database_count_and_ids(
                entrez_search_query=direct_search_query_used + ' AND RefSeq[Filter]',
                maximum_ids_to_return=DEFAULT_MAX_DIRECT_IDS_TO_FETCH,
            )
        )
        direct_search_accession_ids = _reorder_accession_ids_with_refseq_first(
            all_matching_accession_ids=direct_search_accession_ids,
            refseq_only_accession_ids=direct_search_refseq_only_accession_ids,
        )

    direct_hit_records_raw = _efetch_genbank_records_by_accession_ids(
        accession_ids_to_fetch=direct_search_accession_ids,
    )

    # Tag each direct hit with provenance and drop records that are actually
    # polyproteins. Polyproteins are handled via the extraction path below —
    # keeping them here as 'direct' hits would bury the real protein-of-interest.
    direct_hit_records_without_polyproteins = []
    for direct_hit_record in direct_hit_records_raw:
        if record_is_polyprotein(direct_hit_record):
            continue
        direct_hit_record.annotations['search_source'] = 'direct'
        direct_hit_records_without_polyproteins.append(direct_hit_record)

    # ── Strategy B — Polyprotein extraction (only if protein_name set) ────────
    polyprotein_search_query             = _build_polyprotein_search_term(organism_name)
    total_polyprotein_hits_in_ncbi       = 0
    downloaded_polyprotein_records       = []
    extracted_mat_peptide_records        = []

    if protein_name:
        total_polyprotein_hits_in_ncbi, polyprotein_accession_ids = (
            _esearch_protein_database_count_and_ids(
                entrez_search_query=polyprotein_search_query,
                maximum_ids_to_return=DEFAULT_MAX_POLYPROTEINS_TO_FETCH,
            )
        )

        if polyprotein_accession_ids:
            _, polyprotein_refseq_only_accession_ids = (
                _esearch_protein_database_count_and_ids(
                    entrez_search_query=polyprotein_search_query + ' AND RefSeq[Filter]',
                    maximum_ids_to_return=DEFAULT_MAX_POLYPROTEINS_TO_FETCH,
                )
            )
            polyprotein_accession_ids = _reorder_accession_ids_with_refseq_first(
                all_matching_accession_ids=polyprotein_accession_ids,
                refseq_only_accession_ids=polyprotein_refseq_only_accession_ids,
            )

        downloaded_polyprotein_records = _efetch_genbank_records_by_accession_ids(
            accession_ids_to_fetch=polyprotein_accession_ids,
        )

        for downloaded_polyprotein_record in downloaded_polyprotein_records:
            extracted_mat_peptide_record = extract_protein_from_polyprotein(
                polyprotein_record=downloaded_polyprotein_record,
                target_protein_name=protein_name,
            )
            if extracted_mat_peptide_record is not None:
                extracted_mat_peptide_records.append(extracted_mat_peptide_record)

    # ── Merge, deduplicate, host-filter ───────────────────────────────────────
    # Dedup by accession ID. Direct hits win if the same accession appears in
    # both lists (unlikely, but possible if a polyprotein's accession was also
    # returned by the direct search and managed to pass the polyprotein drop).
    already_added_accession_ids = set()
    merged_candidate_records    = []

    for direct_hit_record in direct_hit_records_without_polyproteins:
        if direct_hit_record.id in already_added_accession_ids:
            continue
        already_added_accession_ids.add(direct_hit_record.id)
        merged_candidate_records.append(direct_hit_record)

    for extracted_mat_peptide_record in extracted_mat_peptide_records:
        if extracted_mat_peptide_record.id in already_added_accession_ids:
            continue
        already_added_accession_ids.add(extracted_mat_peptide_record.id)
        merged_candidate_records.append(extracted_mat_peptide_record)

    host_filtered_candidate_records = _filter_records_by_declared_host(
        records_to_filter=merged_candidate_records,
        target_host=target_host,
    )

    # Cap to max_records_to_return (keeps RefSeq-first order from the search)
    final_capped_candidate_records = host_filtered_candidate_records[:max_records_to_return]

    return {
        'records':                           final_capped_candidate_records,
        'total_direct_found':                total_direct_hits_in_ncbi,
        'total_polyprotein_found':           total_polyprotein_hits_in_ncbi,
        'extracted_from_polyprotein_count':  len(extracted_mat_peptide_records),
        'direct_search_tier':                direct_search_tier_used,
        'direct_query':                      direct_search_query_used,
        'polyprotein_query':                 polyprotein_search_query,
    }


# ── Output writers ────────────────────────────────────────────────────────────

def save_sequences_as_fasta(
    selected_records: list,
    output_fasta_path: Path,
) -> Path:
    """
    Saves a list of selected GenBank records as a FASTA file.
    Creates parent directories if they do not exist.
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

    Purpose: step08 (variant search) reads this to exclude the original
    reference sequences from variant results — preventing the reference
    from appearing as its own variant during conservation analysis.
    """
    registry_entries = []
    for selected_record in selected_records:
        source_qualifiers = extract_source_qualifiers(selected_record)
        registry_entries.append({
            'accession_id':        selected_record.id,
            'description':         selected_record.description,
            'organism':            selected_record.annotations.get('organism', 'N/A'),
            'sequence_length_aa':  len(selected_record.seq),
            'search_source':       selected_record.annotations.get('search_source', 'direct'),
            'source_accession':    selected_record.annotations.get('source_accession'),
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
