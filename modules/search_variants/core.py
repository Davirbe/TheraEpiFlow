"""
Variant-search domain logic for search_variants: the HTTP helper, taxonomy
lineage, UniProt variant search/scoring, identity computation and the
validate-and-build-SeqRecords step. Network + compute, no Rich tables.
"""

import time

import requests
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils.console import console

UNIPROT_SEARCH_URL   = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy/{tax_id}"
_MAX_VARIANTS        = 100
_REQUEST_TIMEOUT     = 20
_AMBIGUOUS_AAS       = set("XBZUO")
_INFORMATIVE_RANKS   = frozenset({
    "genus", "subgenus", "subfamily", "family", "superfamily",
    "suborder", "order",
})


def _http_get(url: str, params: dict = None, max_attempts: int = 3) -> requests.Response:
    """GET with exponential-backoff retry for transient network errors."""
    last_err: Exception = RuntimeError("no attempts made")
    for attempt in range(max_attempts):
        try:
            response = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as err:
            last_err = err
            if attempt < max_attempts - 1:
                wait = 2 ** attempt
                console.print(
                    f"[yellow]⚠ UniProt request failed (attempt {attempt + 1}/{max_attempts}): "
                    f"{err} — retrying in {wait}s[/yellow]"
                )
                time.sleep(wait)
    raise last_err


# ── Taxonomy lineage ─────────────────────────────────────────────────────────

def _fetch_taxonomy_lineage(tax_id: int) -> list[dict]:
    """Returns genus/family/order ancestors for tax_id from UniProt, closest-first.
    Returns [] on any failure — non-critical, just skips the family prompt."""
    try:
        resp = _http_get(
            UNIPROT_TAXONOMY_URL.format(tax_id=tax_id), max_attempts=2
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


def _build_candidates(results: list) -> list[dict]:
    candidates = []
    for r in results:
        reviewed = "Swiss-Prot" in r.get("entryType", "")
        organism = r.get("organism", {})
        seq_data = r.get("sequence", {})
        candidates.append({
            "accession":    r.get("primaryAccession", ""),
            "reviewed":     reviewed,
            "protein_name": _extract_protein_name(r),
            "length":       seq_data.get("length", 0),
            "organism":     organism.get("scientificName", ""),
            "tax_id":       organism.get("taxonId", 0),
            "sequence":     seq_data.get("value", ""),
            "identity":     None,
        })
    return candidates


def _search_uniprot_variants(
    protein_name: str,
    tax_id: int | None,
    scope: str,
    host_filter: str | None,
    family_taxid: int | None = None,
) -> list[dict]:
    """Searches UniProt for protein variants, including sequence data."""
    fields = "accession,reviewed,protein_name,organism_name,organism_id,length,sequence"

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
        uniprot_response = _http_get(UNIPROT_SEARCH_URL, params={
            "query":  query,
            "fields": fields,
            "format": "json",
            "size":   str(_MAX_VARIANTS),
        })
        raw_uniprot_results = uniprot_response.json().get("results", [])

    console.print(f"[dim]→ {len(raw_uniprot_results)} raw results returned.[/dim]")
    return _build_candidates(raw_uniprot_results)


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


