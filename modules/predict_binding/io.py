"""Input loading for predict_binding: reads the per-track FASTA (or prompts for
a path / pasted sequences) and expands it into peptides."""

import io
from pathlib import Path

from Bio import SeqIO
from rich.prompt import Prompt

from utils.console import console
from utils.fasta_utils import is_valid_sequence
from utils.naming import get_step_filename

def _load_sequences(track_id: str, input_dir: Path) -> list:
    """Loads validated SeqRecords from the canonical FASTA produced by fetch_sequences.
    Falls back to a user-provided file path or terminal paste if not found."""
    canonical_fasta_path = input_dir / track_id / get_step_filename("SEQUENCES", track_id, ext="fasta")

    if canonical_fasta_path.exists():
        records = list(SeqIO.parse(str(canonical_fasta_path), 'fasta'))
        if records:
            console.print(f"[dim]  Loaded {len(records)} sequence(s) from {canonical_fasta_path.name}[/dim]")
            return records

    console.print(f"[yellow]  Sequence file not found: {canonical_fasta_path}[/yellow]")
    console.print("  Provide sequences via file path or paste FASTA content directly.")
    console.print("[dim]  (Enter a file path, or press Enter to paste sequences)[/dim]")

    user_path_input = Prompt.ask("  File path (or Enter to paste)").strip()

    if user_path_input:
        provided_path = Path(user_path_input)
        if not provided_path.exists():
            raise FileNotFoundError(f"File not found: {provided_path}")
        raw_records = list(SeqIO.parse(str(provided_path), 'fasta'))
    else:
        console.print("[dim]  Paste FASTA content below. Enter a blank line when done.[/dim]")
        pasted_lines = []
        while True:
            line = input()
            if line == '':
                break
            pasted_lines.append(line)
        raw_records = list(SeqIO.parse(io.StringIO('\n'.join(pasted_lines)), 'fasta'))

    if not raw_records:
        raise ValueError("No sequences found in the provided input.")

    valid_records = []
    for record in raw_records:
        sequence_is_valid, rejection_reason = is_valid_sequence(record)
        if sequence_is_valid:
            valid_records.append(record)
        else:
            console.print(f"[yellow]  Skipping {record.id}: {rejection_reason}[/yellow]")

    if not valid_records:
        raise ValueError("No valid sequences after validation.")

    console.print(f"[dim]  {len(valid_records)} valid sequence(s) accepted.[/dim]")
    return valid_records


