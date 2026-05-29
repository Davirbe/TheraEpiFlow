# fetch_sequences

Resolves each track to its reference protein sequence and writes the FASTA that every later step reads. This is the entry point of the pipeline: searches the UniProt REST API (or loads a local FASTA), lets the user pick among candidates, validates the chosen sequence, and persists it together with a full search audit.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `FetchSequencesStep` orchestration: `run`, `describe_outputs`, the UniProt flow and polyprotein-chain slicing |
| `core.py` | Organism normalization, UniProt search + ranking/flagging, record validation, constants (`ORGANISM_ALIASES`, UniProt URLs) |
| `io.py` | UniProt FASTA download (`_download_fasta`) + local-FASTA loading (`_load_local_fasta`) |
| `prompts.py` | Interactive candidate selection + FASTA-first prompt (validates local paths via `utils.input_validation`) |
| `render.py` | Rich candidates table |
| `__init__.py` | Facade that re-exports `FetchSequencesStep` and `ORGANISM_ALIASES` |

`core.py` performs the network GET (through `utils.http.http_get`) but holds no Rich tables or prompts, keeping the domain logic isolated from the UI.

## What it does

For each track (organism + protein):

1. **Resolve the organism** (`core._normalize_organism`). Exact alias → fuzzy alias (`difflib`, cutoff 0.75) → canonical scientific name. Aliases map to a scientific name + NCBI tax ID via `ORGANISM_ALIASES`; unknown organisms fall back to a free-text name search.
2. **Build the UniProt query** (`core._search_uniprot`). Known tax IDs use the exact `(taxonomy_id:{id})`; unknown organisms use `(organism_name:"{name}")`. Three strategies are tried in order until one hits: protein name → gene name → organism only (the last prints a warning so the user knows the search loosened).
3. **Rank & flag** (`core._sort_and_flag`). Swiss-Prot before TrEMBL; relative to the median length of the result set, entries longer than `median × 2` are marked `(Polyprotein?)` and entries shorter than `median × 0.4` are marked `(Fragment?)`.
4. **Select.** Interactive: the user picks a row (or Enter for the suggested one). Non-interactive: Swiss-Prot first, then the first non-flagged TrEMBL hit, then the first hit.
5. **Polyprotein slicing** (`step._slice_polyprotein_chain`). Flaviviruses (ZIKV/DENV/HCV) exist in UniProt only as a single "Genome polyprotein"; the target protein is a Chain feature inside it. When the entry looks like a polyprotein, the mature chain is sliced out using `utils.uniprot.find_chain_for_protein`. On no match it warns and keeps the full polyprotein.
6. **Validate & retry** (`core._validate_records` → `utils.fasta_utils.is_valid_sequence`). Rejects sequences shorter than 50 aa or containing ambiguous residues (B, J, O, U, X, Z). If the chosen candidate fails, the next candidate is tried automatically.
7. **Persist** the FASTA, a registry of every candidate, and a validation report.

## Why UniProt

Swiss-Prot entries are curated, the REST API is more predictable than Entrez, and the metadata schema is consistent across organisms. The trade-off is polyprotein viruses, handled by the Chain-slicing logic above.

## Inputs

Read from `project_config.json` under `tracks[track_id]`:

| Field | Meaning |
|---|---|
| `organism_name` | Alias (`HPV16`, `ZIKV`) or scientific name (`Human papillomavirus 16`) |
| `protein_name` | Protein as named in UniProt (`E6`, `envelope protein`, `nsP1`) |
| `input_source` | `uniprot` (default) or `local` |
| `local_file_path` | Path to a local `.fasta`/`.fa`/`.faa`/`.fas` when `input_source == 'local'` |

Even when `input_source` is `uniprot`, the run offers a local FASTA up front (FASTA-first prompt); accepting it switches the source to `local` for that run.

## Outputs

Saved to `data/input/{track_id}/`:

| File | Contents |
|---|---|
| `SEQUENCES_{track_id}.fasta` | Selected, validated sequence(s), the query for every downstream step |
| `SEQUENCES_VIEW_{track_id}.csv` | Slim per-step view: one row `track_id, accession, organism, protein, length, source` |
| `REGISTRY_{track_id}.json` | Every candidate considered (accession, organism, length, reviewed flag, tax ID, chain slice if any) |
| `VALIDATION_REPORT_{track_id}.json` | Counts of downloaded / validated / rejected, with the rejection reason per record |

After a UniProt run, `project_config['tracks'][track_id]` gains `seed_accession`, `seed_size` and `tax_id`, which persist across reruns.

## Selection table

| Column | Description |
|---|---|
| `#` | Row number for selection |
| `★` | Suggested sequence (Swiss-Prot first, else first non-flagged hit) |
| Accession | UniProt accession |
| aa | Sequence length |
| Source | `Swiss-Prot` (reviewed) or `TrEMBL` |
| Flag | `(Polyprotein?)`, `(Fragment?)` or empty |
| Description | UniProt protein description |

Press Enter for the suggested row, or type a single row number. (Multi-selection like `1,3,5`, `1-4`, or `all` belongs to `search_variants`, not here: a track resolves to exactly one seed sequence.)

## Known organisms

`ORGANISM_ALIASES` (in `core.py`) ships 24 aliases: HPV16/18/31/33/45/52/58, ZIKV, DENV + DENV1–4, CHIKV, HCV, MPOX, EBOV, HIV, SARS2, MERS, INFA, INFB, RSV, RABV. Adding one is a single dict entry of `(scientific_name, tax_id)`.

## References

UniProt is the sequence source for this step:

- The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2023.* Nucleic Acids Research, 2023. doi:10.1093/nar/gkac1052
