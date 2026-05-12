"""
Shared FASTA parsing and validation utilities.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


AMBIGUOUS_AMINO_ACIDS = set("XBZUO")
MIN_SEQUENCE_LENGTH   = 50


def write_fasta(records: list[SeqRecord], output_path: str | Path):
    """Writes a list of SeqRecord objects to a FASTA file."""
    SeqIO.write(records, str(output_path), "fasta")


def has_ambiguous_residues(sequence: str) -> bool:
    """Returns True if the sequence contains ambiguous amino acid characters."""
    return bool(set(sequence.upper()) & AMBIGUOUS_AMINO_ACIDS)


def is_valid_sequence(record: SeqRecord) -> tuple[bool, str]:
    """
    Validates a SeqRecord.
    Returns (is_valid, reason_if_invalid).
    """
    seq = str(record.seq).upper()

    if len(seq) < MIN_SEQUENCE_LENGTH:
        return False, f"Sequence too short ({len(seq)} aa, minimum {MIN_SEQUENCE_LENGTH})"

    if has_ambiguous_residues(seq):
        return False, f"Contains ambiguous residues: {set(seq) & AMBIGUOUS_AMINO_ACIDS}"

    return True, ""


def generate_peptides(sequence: str, lengths: list[int]) -> list[str]:
    """Generates all possible peptides of given lengths from a sequence."""
    peptides = []
    for length in lengths:
        peptides.extend(
            sequence[i: i + length]
            for i in range(len(sequence) - length + 1)
        )
    return peptides
