# Utils

Shared utility modules used across multiple pipeline steps.

---

## genbank_utils.py

**Used by:** Step 01 (fetch sequences), Step 08 (search variants)

Functions for searching, downloading and processing sequences from NCBI GenBank.

| Function | Purpose |
|---|---|
| `search_ncbi_protein_ids()` | Phase 1: returns accession IDs only (no download). RefSeq entries returned first. |
| `fetch_records_by_accession_ids()` | Phase 2: downloads full GenBank records for a list of IDs. Filters by target host. |
| `suggest_reference_sequence()` | Suggests the best reference from a list of records. RefSeq-first, then modal length. |
| `extract_source_qualifiers()` | Extracts metadata (strain, isolate, host, location, date) from a GenBank record. |
| `record_is_polyprotein()` | Returns True if record is a polyprotein (description + length heuristic). |
| `record_is_refseq()` | Returns True if record is a curated NCBI RefSeq entry. |
| `extract_protein_from_polyprotein()` | Extracts a specific mature peptide from a polyprotein using `mat_peptide` features. |
| `save_sequences_as_fasta()` | Saves selected records as a FASTA file. |
| `build_sequence_registry()` | Saves a JSON registry of downloaded sequences with full metadata. |

**Note:** `Entrez.email` must be set by the caller before using search functions.

---

## project_manager.py

**Used by:** `main.py`

Handles the full project lifecycle: creation, track setup, listing, loading, and deletion.

| Function | Purpose |
|---|---|
| `create_project_interactive()` | Interactive wizard: 4 questions → creates minimal project structure |
| `setup_project_tracks_interactive()` | Asks organisms, proteins, labels and input sources → creates tracks |
| `create_project()` | Non-interactive project creation (called by `create_project_interactive`) |
| `load_project_config()` | Returns the `project_config.json` dict for a project |
| `save_project_config()` | Saves updated config (each step adds its own configuration here) |
| `list_projects()` | Returns all projects with metadata for display |
| `tracks_are_defined()` | Returns True if the project has at least one track defined |
| `get_track_intermediate_dir()` | Returns the intermediate data path for a track |
| `get_tracks_for_protein()` | Returns all tracks sharing a protein name (e.g. all E6 tracks) |
| `delete_project()` | Permanently deletes a project folder and registry entry |
| `_suggest_organism_label()` | Suggests standard abbreviation (HPV16, ZIKV, DENV2, SARS2...) |
| `_suggest_protein_label()` | Suggests protein abbreviation (E6, E7, E, NS1, NS5, S, HA...) |

---

## pipeline_state.py

**Used by:** `modules/base_step.py`

Low-level helpers for reading and writing `pipeline.json`.

| Function | Purpose |
|---|---|
| `load_pipeline_state()` | Loads `pipeline.json` for a project |
| `save_pipeline_state()` | Saves updated pipeline state |
| `get_track_step_status()` | Returns status of a step for a specific track |
| `set_track_step_status()` | Marks a step as done/error with timestamp |
| `get_global_step_status()` | Returns status of a global step (steps 13-14) |
| `set_global_step_status()` | Marks a global step as done/error |

---

## naming.py

**Used by:** `project_manager.py`, all step modules

Standardized naming conventions for files and columns.

| Item | Pattern | Example |
|---|---|---|
| Track ID | `{ORGANISM}_{PROTEIN}` | `HPV16_E6` |
| Prediction file | `PRED_{TOOL}_{TRACK_ID}.csv` | `PRED_NET_HPV16_E6.csv` |
| Consensus file | `CONSENSUS_{TRACK_ID}.csv` | `CONSENSUS_HPV16_E6.csv` |
| HLA allele (NetMHCpan) | `HLA-A*02:01` → `HLA-A0201` | Format conversion |
| Column suffix | `_{tool}_{type}` | `_net_pred`, `_flurry_pred` |

Key function: `build_track_id(organism_label, protein_label)` → `"HPV16_E6"`

---

## fasta_utils.py

**Used by:** Step 02 (validate input), Step 03 (predict binding)

Utilities for handling FASTA sequences.

| Function | Purpose |
|---|---|
| `parse_fasta()` | Parses a FASTA file into a list of SeqRecord objects |
| `validate_sequence()` | Checks for ambiguous amino acids and minimum length |
| `remove_duplicate_sequences()` | Removes records with identical amino acid sequences |
| `generate_peptides()` | Generates all overlapping peptides of given lengths from a sequence |

---

## column_finder.py

**Used by:** Step 04 (consensus filter), Step 06 (select representatives)

Flexible column lookup for DataFrames where column names may vary between runs.

| Function | Purpose |
|---|---|
| `find_column()` | Finds a column by partial name match, raising an error if not found |
| `find_columns_by_suffix()` | Returns all columns ending with a given suffix |

Adapted from the original HPV pipeline scripts.
