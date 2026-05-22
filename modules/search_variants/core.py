"""
Variant-search domain logic for search_variants: the HTTP helper, taxonomy
lineage, UniProt variant search/scoring, identity computation and the
validate-and-build-SeqRecords step. Network + compute, no Rich tables.
"""

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from config import VARIANTS_MAX_RESULTS
from utils.console import console
from utils.http import http_get
from utils.uniprot import find_chain_for_protein

UNIPROT_SEARCH_URL   = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy/{tax_id}"
_MAX_VARIANTS        = VARIANTS_MAX_RESULTS
_REQUEST_TIMEOUT     = 20
_AMBIGUOUS_AAS       = set("XBZUO")
_INFORMATIVE_RANKS   = frozenset({
    "genus", "subgenus", "subfamily", "family", "superfamily",
    "suborder", "order",
})


# ── Taxonomy lineage ─────────────────────────────────────────────────────────

def _fetch_taxonomy_lineage(tax_id: int) -> list[dict]:
    """Returns genus/family/order ancestors for tax_id from UniProt, closest-first.
    Returns [] on any failure — non-critical, just skips the family prompt."""
    try:
        resp = http_get(
            UNIPROT_TAXONOMY_URL.format(tax_id=tax_id), max_attempts=2,
            timeout=_REQUEST_TIMEOUT,
        )
        data    = resp.json()
        lineage = data.get("lineage", [])
        result  = []
        for entry in reversed(lineage):   # root→species → reverse to closest-first
            rank = (entry.get("rank") or "").lower().replace(" ", "")
            if rank in _INFORMATIVE_RANKS:
                result.append({
                    "name":  entry.get("scientificName", ""),
                    "rank":  rank,
                    "taxid": entry.get("taxonId", 0),
                })
        return result
    except Exception:
        return []


def _validate_variant(record: SeqRecord) -> tuple[bool, str]:
    """Lenient validation for variants (conservation, not binding) — rejects only empty seqs
    or all-ambiguous residues. Short sequences (< 50 aa) pass; caller emits a warning."""
    seq = str(record.seq).upper().strip()
    if not seq:
        return False, "empty sequence"
    if all(aa in _AMBIGUOUS_AAS for aa in seq):
        return False, "sequence consists entirely of ambiguous residues (X/B/Z/U/O)"
    return True, ""


def _extract_protein_name(result: dict) -> str:
    desc = result.get("proteinDescription", {})
    recommended = desc.get("recommendedName", {})
    if recommended:
        return recommended.get("fullName", {}).get("value", "N/A")
    submitted = desc.get("submittedName", [])
    if submitted:
        return submitted[0].get("fullName", {}).get("value", "N/A")
    return "N/A"


def _build_candidates(
    results: list,
    protein_name: str = "",
    should_slice_variants: bool = False,
) -> list[dict]:
    """Builds candidate dicts from raw UniProt results.

    When should_slice_variants is set, each candidate that is a polyprotein (has a Chain
    feature matching protein_name) is sliced down to that mature chain — the same
    treatment fetch_sequences applies to the reference — so identity is computed
    chain-vs-chain rather than chain-vs-polyprotein. Candidates with no matching Chain
    keep their full sequence."""
    candidates = []
    for raw_result in results:
        is_reviewed = "Swiss-Prot" in raw_result.get("entryType", "")
        organism    = raw_result.get("organism", {})
        sequence_data = raw_result.get("sequence", {})
        candidate = {
            "accession":    raw_result.get("primaryAccession", ""),
            "reviewed":     is_reviewed,
            "protein_name": _extract_protein_name(raw_result),
            "length":       sequence_data.get("length", 0),
            "organism":     organism.get("scientificName", ""),
            "tax_id":       organism.get("taxonId", 0),
            "sequence":     sequence_data.get("value", ""),
            "identity":     None,
        }
        if should_slice_variants and protein_name:
            _slice_candidate_to_chain(candidate, raw_result, protein_name)
        candidates.append(candidate)
    return candidates


def _slice_candidate_to_chain(candidate: dict, raw_result: dict, protein_name: str) -> None:
    """Slices a polyprotein candidate down to its matching mature chain, in place.

    Reads Chain features from the search result (requested via fields=...,ft_chain) and,
    when one matches protein_name, trims candidate['sequence'] to seq[chain_start-1:chain_end],
    updating 'length' and recording 'chain_slice'. No matching chain → left untouched."""
    matched_chain = find_chain_for_protein(raw_result, protein_name)
    if matched_chain is None:
        return
    full_sequence = candidate["sequence"]
    if not full_sequence:
        return
    chain_start = matched_chain["start"]
    chain_end   = matched_chain["end"]
    sliced_sequence = full_sequence[chain_start - 1:chain_end]
    candidate["sequence"]    = sliced_sequence
    candidate["length"]      = len(sliced_sequence)
    candidate["chain_slice"] = {
        "start":       chain_start,
        "end":         chain_end,
        "description": matched_chain["description"],
        "match_score": round(matched_chain["score"], 3),
    }


# ── Genotype grouping (interspecific presentation) ────────────────────────────

def _genotype_label(candidate: dict) -> str:
    """Human-readable genotype label for a candidate (organism name, tax_id fallback)."""
    return candidate.get("organism") or f"tax_id {candidate.get('tax_id', '?')}"


def _group_candidates_by_genotype(candidates: list[dict]) -> dict[int, list[dict]]:
    """Groups candidates by genotype (organism tax_id), preserving input order.
    Returns {tax_id: [candidates]} — each bucket is one genotype (e.g. one HPV type)."""
    groups: dict[int, list[dict]] = {}
    for candidate in candidates:
        groups.setdefault(candidate.get("tax_id", 0), []).append(candidate)
    return groups


def _pick_best_per_genotype(candidates: list[dict]) -> list[dict]:
    """One best representative per genotype: highest identity, Swiss-Prot as tiebreak.
    Returned list is sorted by identity descending (best genotype first)."""
    best_by_genotype: dict[int, dict] = {}
    for tax_id, bucket in _group_candidates_by_genotype(candidates).items():
        best_by_genotype[tax_id] = max(
            bucket,
            key=lambda c: (c.get("identity") or 0.0, 1 if c.get("reviewed") else 0),
        )
    return sorted(
        best_by_genotype.values(),
        key=lambda c: c.get("identity") or 0.0,
        reverse=True,
    )


def _search_uniprot_variants(
    protein_name: str,
    tax_id: int | None,
    scope: str,
    host_filter: str | None,
    family_taxid: int | None = None,
    should_slice_variants: bool = False,
) -> list[dict]:
    """Searches UniProt for protein variants, including sequence data.

    When should_slice_variants is set, also requests Chain features (ft_chain) so
    polyprotein candidates can be sliced down to their mature chain — see
    _build_candidates / _slice_candidate_to_chain."""
    fields = "accession,reviewed,protein_name,organism_name,organism_id,length,sequence"
    if should_slice_variants:
        fields += ",ft_chain"

    if scope == "intraspecific" and tax_id:
        query = f'(taxonomy_id:{tax_id}) AND (protein_name:"{protein_name}")'
    elif scope == "interspecific" and family_taxid:
        query = f'(taxonomy_id:{family_taxid}) AND (protein_name:"{protein_name}")'
        if host_filter:
            query += f' AND (virus_host_name:"{host_filter}")'
    else:
        query = f'(protein_name:"{protein_name}")'
        if host_filter:
            query += f' AND (virus_host_name:"{host_filter}")'

    console.print(f"[dim]Query: {query}[/dim]")

    with console.status(
        f"[yellow]Searching UniProt variants ({scope})…[/yellow]",
        spinner="dots",
    ):
        uniprot_response = http_get(UNIPROT_SEARCH_URL, params={
            "query":  query,
            "fields": fields,
            "format": "json",
            "size":   str(_MAX_VARIANTS),
        }, timeout=_REQUEST_TIMEOUT)
        raw_uniprot_results = uniprot_response.json().get("results", [])

    console.print(f"[dim]→ {len(raw_uniprot_results)} raw results returned.[/dim]")
    return _build_candidates(raw_uniprot_results, protein_name, should_slice_variants)


# ── Identity computation ──────────────────────────────────────────────────────

def _compute_identity(seq_a: str, seq_b: str) -> float:
    """Returns percent identity (0–100) via global alignment (match=1, gaps=0).
    Denominator is min(len_a, len_b) — a short fragment that perfectly matches its region
    in a longer sequence scores 100% instead of being penalised by max(). The right metric
    when comparing full proteins to fragments or to proteins of different lengths."""
    aligner = PairwiseAligner()
    aligner.mode             = "global"
    aligner.match_score      = 1
    aligner.mismatch_score   = 0
    aligner.open_gap_score   = 0
    aligner.extend_gap_score = 0
    score   = aligner.score(seq_a, seq_b)
    min_len = min(len(seq_a), len(seq_b))
    if min_len == 0:
        return 0.0
    return float(score) / min_len * 100.0


def _build_and_validate(candidates: list[dict]) -> tuple[list[SeqRecord], list[dict]]:
    """Builds SeqRecords from candidates and applies lenient validation; returns (valid_records, rejected_log).
    Short sequences (< 50 aa) pass with a console warning."""
    valid_records:  list[SeqRecord] = []
    rejected_log:   list[dict]      = []
    short_warnings: list[str]       = []

    for c in candidates:
        seq_str = c.get("sequence", "")
        record  = SeqRecord(
            seq=Seq(seq_str),
            id=c.get("accession", "unknown"),
            description=(
                f"{c.get('organism', '')} | {c.get('protein_name', '')[:60]}"
                + (f" | identity={c['identity']:.1f}%" if c.get("identity") is not None else "")
            ).strip(" |"),
        )
        ok, reason = _validate_variant(record)
        if not ok:
            rejected_log.append({"id": record.id, "reason": reason})
            continue

        if len(seq_str) < 50:
            short_warnings.append(f"{record.id} ({len(seq_str)} aa)")

        valid_records.append(record)

    if short_warnings:
        console.print(
            f"[yellow]⚠ {len(short_warnings)} short sequence(s) included "
            f"(< 50 aa — valid for conservation, not for binding prediction):[/yellow]"
        )
        for w in short_warnings:
            console.print(f"  [yellow]{w}[/yellow]")

    return valid_records, rejected_log


