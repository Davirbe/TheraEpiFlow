# fetch_sequences

Downloads the reference protein sequence for each track from the UniProt REST API. This is the entry point of the pipeline. Every later step assumes the FASTA produced here is the canonical sequence for the track.

## What it does

For each track (organism plus protein), the step:

1. Resolves the organism. Short aliases like `HPV16`, `MPOX`, `CHIKV` are mapped to a scientific name and an NCBI tax ID through the `ORGANISM_ALIASES` dictionary. Unknown organisms fall back to a free-text search.
2. Builds a UniProt query. Known tax IDs use `(taxonomy_id:{id})`, which is exact. Unknown organisms use `(organism_name:"{name}")`, which is fuzzier.
3. Tries three search strategies in order until something hits: protein name, gene name, organism alone (the last one prints a warning so the user knows the search loosened).
4. Flags candidates relative to the global median length of the result set. Anything longer than median * 2 gets a `(Polyprotein?)` mark, anything shorter than median * 0.4 gets `(Fragment?)`.
5. In interactive mode, prompts the user to pick one of the candidates. In non-interactive mode, prefers Swiss-Prot, then non-flagged TrEMBL, then the first hit.
6. If the chosen sequence fails validation (presence of `X`, length too short, etc.), retries with the next candidate automatically.
7. Saves the sequence, a JSON registry with metadata, and a validation report.

## Why UniProt instead of NCBI

The earlier version used Entrez/GenBank. UniProt was chosen because Swiss-Prot entries are curated, the REST API is more predictable than Entrez, and the metadata schema is consistent across organisms. The trade-off: flaviviruses (DENV, ZIKV, HCV) only appear in UniProt as Chain features inside a single "Genome polyprotein" record, so the pipeline currently uses the full polyprotein for those organisms. A planned fix slices out the mature protein region using the Chain feature coordinates.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `FetchSequencesStep` orchestration — `run` / `describe_outputs` |
| `core.py` | Organism normalization, UniProt search/scoring, record validation, constants (`ORGANISM_ALIASES`) |
| `io.py` | UniProt FASTA download + local-FASTA loading |
| `prompts.py` | Interactive candidate selection |
| `render.py` | Rich candidates table |
| `__init__.py` | Facade — re-exports `FetchSequencesStep` and `ORGANISM_ALIASES` |

## Input

`fetch_sequences` reads from `project_config.json`:

```json
"tracks": {
  "HPV16_E6": {
    "organism_name":  "Human papillomavirus 16",
    "organism_label": "HPV16",
    "protein_name":   "E6",
    "protein_label":  "E6",
    "input_source":   "uniprot",
    "local_file_path": null
  }
}
```

When `input_source` is `local_file`, the step skips the API and reads from `local_file_path` instead.

## Output

Saved to `data/input/{track_id}/`:

| File | Contents |
|---|---|
| `SEQUENCES_{track_id}.fasta` | Selected sequence in FASTA format |
| `SEQUENCES_VIEW_{track_id}.csv` | Slim per-step view — one row with `track_id, accession, organism, protein, length, source` for quick inspection |
| `REGISTRY_{track_id}.json` | UniProt accession, tax ID, organism, length, source DB |
| `VALIDATION_REPORT_{track_id}.json` | Length checks, ambiguous residue counts, retry trail |

After the run, `project_config['tracks'][track_id]` gains `seed_accession`, `seed_size`, `tax_id`, and `input_source='uniprot'`. These are used downstream and survive across reruns.

## Selection table

When more than one candidate is shown:

| Column | Description |
|---|---|
| `#` | Row number for selection |
| `★` | Suggested sequence (Swiss-Prot wins, then closest to median length) |
| Accession | UniProt accession |
| aa | Sequence length in amino acids |
| Source | `Swiss-Prot` (curated) or `TrEMBL` |
| Flag | `(Polyprotein?)`, `(Fragment?)`, or empty |
| Description | UniProt protein description |

Pressing Enter picks the suggested row. Other accepted inputs: `1`, `1,3,5`, `1-4`, `all`.

## Known organisms

`ORGANISM_ALIASES` covers the targets used so far: HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58, ZIKV, DENV1-4, CHIKV, MPOX, HCV, SARS-CoV-2. Adding a new organism only requires appending one entry with its scientific name and tax ID.
