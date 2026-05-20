"""I/O for fetch_sequences: UniProt FASTA download and local-FASTA loading."""

import io
from typing import Optional

from Bio import SeqIO

from utils.console import console

from .core import UNIPROT_FASTA_URL, _http_get

# ── FASTA download ─────────────────────────────────────────────────────────────

def _download_fasta(accession: str) -> str:
    """Downloads FASTA text for a single UniProt accession."""
    url = UNIPROT_FASTA_URL.format(accession=accession)
    response = _http_get(url)
    return response.text


# ── Local FASTA loader (unchanged) ────────────────────────────────────────────

def _load_local_fasta(local_file_path: Optional[str] = None) -> list:
    """Loads sequences from a local FASTA file."""
    if local_file_path:
        fasta_path = Path(local_file_path)
    else:
        console.print('\n[bold]Path to local FASTA file:[/bold]')
        fasta_path = Path(input('> ').strip())

    if not fasta_path.exists():
        raise FileNotFoundError(f'File not found: {fasta_path}')
    if fasta_path.suffix.lower() not in ('.fasta', '.fa', '.faa', '.fas'):
        console.print('[yellow]Warning: file extension not recognized as FASTA. '
                      'Attempting to parse anyway...[/yellow]')

    records = list(SeqIO.parse(str(fasta_path), 'fasta'))
    if not records:
        raise ValueError(f'No sequences found in file: {fasta_path}')

    console.print(f'[green]Loaded {len(records)} sequence(s) from {fasta_path}[/green]')
    return records


