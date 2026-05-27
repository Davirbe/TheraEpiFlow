"""IEDB IQ-API client for the clinical-validation experiment.

Port (English + library-shaped) of the original prototype at
`existing_scripts/teste.py`. For each peptide, queries two PostgREST
endpoints exposed by the IEDB Query API:

  - `tcell_search` — T-cell assays, MHC Class I only, host = Human
  - `mhc_search`   — MHC binding/elution assays, MHC Class I only, host = Human

B-cell assays are excluded by simply not querying that endpoint (there is no
bcell endpoint here). MHC Class I is enforced server-side via the same
`mhc_class=eq.I` filter on both endpoints; if the column name turns out to be
different on `tcell_search` (it has happened historically) the client retries
without that filter and flags the result as `class_filter_dropped`.

Returns one record per peptide containing the assay counts and a deduplicated
list of supporting publications (PMIDs + titles). Callers decide how to
aggregate (per organism×protein, per ★ vs negative-control, etc.).

Usage as a library:
    from tests.validation.iedb_query import query_peptides
    records = query_peptides(["SIINFEKL", "GILGFVFTL"], host="Human")

Usage on the command line:
    python -m tests.validation.iedb_query --peptides SIINFEKL,GILGFVFTL --out hits.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import requests


IEDB_BASE_URL = "https://query-api.iedb.org"
DEFAULT_HOST = "Human"
DEFAULT_TIMEOUT_SECONDS = 30
DEFAULT_RETRY_COUNT = 2
DEFAULT_RETRY_BACKOFF_SECONDS = 2.0

# Server-side class filter applied to both endpoints. If `tcell_search` rejects
# this column we retry without it (see `_search_with_class_fallback`).
DEFAULT_CLASS_FILTER_VALUE = "I"
CLASS_FILTER_COLUMN = "mhc_class"

# IEDB returns publication identifiers under several alternative keys
# depending on the endpoint, so we look at each candidate in order.
PMID_FIELD_CANDIDATES   = ("pubmed_id", "reference_pubmed_id", "reference_id")
TITLE_FIELD_CANDIDATES  = ("reference_titles", "reference_title", "title", "reference_name")


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class PeptideHits:
    """Single peptide's IEDB hit summary."""
    peptide:                str
    tcell_assay_count:      int                          = 0
    mhc_assay_count:        int                          = 0
    publications:           list[tuple[str, str]]        = field(default_factory=list)
    class_filter_dropped:   bool                         = False
    errors:                 list[str]                    = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "peptide":              self.peptide,
            "tcell_assay_count":    self.tcell_assay_count,
            "mhc_assay_count":      self.mhc_assay_count,
            "publication_count":    len(self.publications),
            "publications":         [
                {"pmid": pmid, "title": title} for pmid, title in self.publications
            ],
            "class_filter_dropped": self.class_filter_dropped,
            "errors":               list(self.errors),
        }


# ── Low-level HTTP helper ─────────────────────────────────────────────────────

def _http_get_json(
    url: str,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
    retry_count: int = DEFAULT_RETRY_COUNT,
) -> tuple[Optional[list], Optional[str]]:
    """GET `url`, returning (json_list, error_message).

    Treats anything other than a JSON list as an error — the PostgREST layer
    answers errors as JSON OBJECTS (`{"code": ..., "message": ...}`) so the
    is-list check distinguishes data from error responses without parsing
    HTTP status codes.
    """
    last_error: Optional[str] = None
    for attempt_index in range(retry_count + 1):
        try:
            response = requests.get(url, timeout=timeout_seconds)
        except requests.RequestException as connection_exception:
            last_error = f"connection error: {connection_exception}"
        else:
            try:
                decoded_payload = response.json()
            except ValueError:
                last_error = (
                    f"non-json response (HTTP {response.status_code}): "
                    f"{response.text[:200]}"
                )
            else:
                if isinstance(decoded_payload, list):
                    return decoded_payload, None
                last_error = (
                    f"API error (HTTP {response.status_code}): {decoded_payload}"
                )
        if attempt_index < retry_count:
            time.sleep(DEFAULT_RETRY_BACKOFF_SECONDS * (attempt_index + 1))
    return None, last_error


def _build_url(endpoint: str, filters: list[str]) -> str:
    return f"{IEDB_BASE_URL}/{endpoint}?" + "&".join(filters)


# ── Field extraction ──────────────────────────────────────────────────────────

def _first_present(record: dict, candidate_keys: Iterable[str]):
    for key in candidate_keys:
        value = record.get(key)
        if value not in (None, "", [], {}):
            return value
    return None


def _normalize_to_text(value) -> str:
    if isinstance(value, list):
        return "; ".join(str(item) for item in value if item not in (None, ""))
    return str(value) if value is not None else ""


def _extract_publications(records: list[dict]) -> list[tuple[str, str]]:
    """Deduplicate publications across an endpoint's records (by PMID + title)."""
    seen: dict[tuple[str, str], None] = {}
    for record in records:
        pmid_raw  = _first_present(record, PMID_FIELD_CANDIDATES)
        title_raw = _first_present(record, TITLE_FIELD_CANDIDATES)
        pmid_text  = _normalize_to_text(pmid_raw)  if pmid_raw  is not None else "N/A"
        title_text = _normalize_to_text(title_raw) if title_raw is not None else "Title not available"
        seen.setdefault((pmid_text, title_text), None)
    return list(seen.keys())


# ── Per-endpoint search with class-filter fallback ────────────────────────────

def _search_with_class_fallback(
    endpoint: str,
    sequence_filter: str,
    host_filter: str,
) -> tuple[Optional[list], bool, Optional[str]]:
    """Try the class-restricted search; if the column is missing, retry without it.

    Returns (records_or_None, class_filter_dropped_flag, error_or_None).
    """
    class_filter = f"{CLASS_FILTER_COLUMN}=eq.{DEFAULT_CLASS_FILTER_VALUE}"
    primary_url  = _build_url(endpoint, [sequence_filter, host_filter, class_filter])
    records, error_message = _http_get_json(primary_url)
    if error_message is None:
        return records, False, None

    # Heuristic: column-not-found errors mention the column name; retry without it.
    fallback_needed = CLASS_FILTER_COLUMN in (error_message or "")
    if not fallback_needed:
        return None, False, error_message

    fallback_url = _build_url(endpoint, [sequence_filter, host_filter])
    fallback_records, fallback_error = _http_get_json(fallback_url)
    if fallback_error is not None:
        return None, False, fallback_error
    return fallback_records, True, None


# ── Public API ────────────────────────────────────────────────────────────────

def query_peptide(
    peptide:    str,
    host:       str = DEFAULT_HOST,
) -> PeptideHits:
    """Query both endpoints for a single peptide and return a `PeptideHits` summary."""
    peptide_clean = (peptide or "").strip().upper()
    if not peptide_clean.isalpha():
        hits = PeptideHits(peptide=peptide)
        hits.errors.append("invalid peptide: only A-Z characters allowed")
        return hits

    sequence_filter = f"linear_sequence=eq.{peptide_clean}"
    host_filter     = f"host_organism_name=ilike.%25{host}%25"

    hits = PeptideHits(peptide=peptide_clean)

    tcell_records, tcell_dropped, tcell_error = _search_with_class_fallback(
        "tcell_search", sequence_filter, host_filter,
    )
    if tcell_error:
        hits.errors.append(f"tcell_search: {tcell_error}")
    if tcell_records is not None:
        hits.tcell_assay_count = len(tcell_records)
        hits.class_filter_dropped = hits.class_filter_dropped or tcell_dropped
        for publication in _extract_publications(tcell_records):
            hits.publications.append(publication)

    mhc_records, mhc_dropped, mhc_error = _search_with_class_fallback(
        "mhc_search", sequence_filter, host_filter,
    )
    if mhc_error:
        hits.errors.append(f"mhc_search: {mhc_error}")
    if mhc_records is not None:
        hits.mhc_assay_count = len(mhc_records)
        hits.class_filter_dropped = hits.class_filter_dropped or mhc_dropped
        for publication in _extract_publications(mhc_records):
            hits.publications.append(publication)

    # Deduplicate publications union (T-cell + MHC may cite the same paper).
    deduplicated: dict[tuple[str, str], None] = {pub: None for pub in hits.publications}
    hits.publications = list(deduplicated.keys())
    return hits


def query_peptides(
    peptides:           Iterable[str],
    host:               str = DEFAULT_HOST,
    sleep_between_secs: float = 0.25,
) -> list[PeptideHits]:
    """Query every peptide; small inter-request sleep to be polite to the IEDB API."""
    results: list[PeptideHits] = []
    for peptide in peptides:
        results.append(query_peptide(peptide=peptide, host=host))
        if sleep_between_secs > 0:
            time.sleep(sleep_between_secs)
    return results


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--peptides",
        help="Comma-separated peptide list (e.g. SIINFEKL,GILGFVFTL).",
    )
    group.add_argument(
        "--peptides-file", type=Path,
        help="Plain-text file with one peptide per line.",
    )
    parser.add_argument(
        "--host", default=DEFAULT_HOST,
        help='Host organism filter (default: "Human").',
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Destination JSON path. Default: prints to stdout.",
    )
    args = parser.parse_args()

    if args.peptides:
        peptide_list = [item.strip() for item in args.peptides.split(",") if item.strip()]
    else:
        peptide_list = [
            line.strip()
            for line in args.peptides_file.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]

    if not peptide_list:
        parser.error("no peptides provided.")

    results = query_peptides(peptide_list, host=args.host)
    payload = {
        "host":            args.host,
        "peptide_count":   len(results),
        "records":         [r.to_dict() for r in results],
    }
    serialized = json.dumps(payload, indent=2, ensure_ascii=False)

    if args.out is None:
        sys.stdout.write(serialized + "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(serialized)
        print(f"Wrote {len(results)} record(s) to {args.out}")


if __name__ == "__main__":
    main()
