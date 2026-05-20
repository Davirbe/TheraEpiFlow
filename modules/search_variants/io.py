"""I/O for search_variants: reference-sequence/accession loading and the
empty-output writer used when no variants are found."""

import datetime
import json
import time
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from utils.fasta_utils import write_fasta
from utils.naming import get_step_filename

def _load_reference_sequence(input_dir: Path, track_id: str) -> SeqRecord:
    fasta_path = input_dir / track_id / get_step_filename("SEQUENCES", track_id, ext="fasta")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in reference FASTA: {fasta_path}")
    return records[0]


def _load_reference_accessions(input_dir: Path, track_id: str) -> set[str]:
    registry_path = input_dir / track_id / get_step_filename("REGISTRY", track_id, ext="json")
    if not registry_path.exists():
        return set()
    with open(registry_path, encoding="utf-8") as fh:
        data = json.load(fh)
    return {s["accession_id"] for s in data.get("sequences", [])}


def _write_empty_outputs(
    fasta_path: Path,
    audit_path: Path,
    track_id: str,
    scope: str,
    host_filter: str | None,
    family_taxid: int | None,
    ref_accessions: set[str],
    note: str,
):
    fasta_path.touch()
    audit = {
        "track_id":           track_id,
        "generated_at":       datetime.datetime.now().isoformat(),
        "scope":              scope,
        "family_taxid":       family_taxid,
        "host_filter":        host_filter,
        "reference_excluded": list(ref_accessions),
        "total_valid":        0,
        "note":               note,
    }
    with open(audit_path, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, ensure_ascii=False)


