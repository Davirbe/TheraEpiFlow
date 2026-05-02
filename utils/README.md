# utils

Shared helpers used across pipeline steps. Each file groups functions by concern, with no thin wrappers and no helpers used by a single caller.

## naming.py

Standard naming for files, columns, and aliases. Used everywhere.

| Item | Pattern or behavior |
|---|---|
| `build_track_id(organism_label, protein_label)` | Returns `"HPV16_E6"` |
| `get_prediction_filename(tool, track_id)` | Returns `"PRED_NET_HPV16_E6.csv"` for `tool="NET"` etc. |
| `get_step_filename(step, track_id, tool="", ext="csv")` | Returns `"CONSENSUS_HPV16_E6.csv"`, `"TOXICITY_AUDIT_HPV16_E6.json"`, etc. |
| `allele_to_netmhcpan_format(allele)` | `"HLA-A*02:01"` becomes `"HLA-A0201"` (strips `*` and `:`) |
| `find_column_name(df, candidates)` | Returns the first matching column name, or `None` |
| `require_column(df, candidates, context="")` | Same, but raises `ValueError` if no match |

Constants:

| Name | Value |
|---|---|
| `SUFFIX_NET_PRED`, `SUFFIX_FLURRY_PRED` | `"_net_pred"`, `"_flurry_pred"` |
| `SUFFIX_NET_PROC`, `SUFFIX_FLURRY_PROC` | `"_net_proc"`, `"_flurry_proc"` |
| `COLUMN_PEPTIDE` | `"peptide"` |
| `COLUMN_NETMHC_EL_PERCENTILE` | `"netmhcpan_el_percentile_net_pred"` |
| `COLUMN_FLURRY_PERCENTILE` | `"mhcflurry_presentation_percentile_flurry_pred"` |

`column_finder.py` was absorbed into this file in April 2026, so older references to `from utils.column_finder import ...` should now point at `utils.naming`.

## project_manager.py

Project lifecycle. Used by `main.py` and most steps when they need to read or save shared config.

| Function | Purpose |
|---|---|
| `create_project_interactive()` | 4-question wizard, creates the minimal project structure |
| `setup_project_tracks_interactive()` | Asks for organisms, proteins, labels, and input source on the first run of `fetch_sequences` |
| `create_project()` | Non-interactive project creation |
| `load_project_config()` / `save_project_config()` | Read or write `project_config.json` |
| `list_projects()` | Returns project metadata for the menu and `--list` |
| `tracks_are_defined()` | True if at least one track is configured |
| `get_track_intermediate_dir()` | Path to `data/intermediate/{track_id}` |
| `get_tracks_for_protein()` | All tracks that share a protein name (for example all E6 tracks) |
| `delete_project()` | Permanently removes the project directory |
| `_suggest_organism_label()` | Default abbreviation for known organisms (HPV16, ZIKV, DENV2, SARS2, MPOX, CHIKV, etc.) |
| `_suggest_protein_label()` | Default protein abbreviation (E6, E7, E, NS1, NS5, S, HA, etc.) |

## pipeline_state.py

Low-level helpers around `pipeline.json`. Used by `BaseTrackStep.execute()` and the REPL.

| Function | Purpose |
|---|---|
| `load_pipeline_state()` / `save_pipeline_state()` | Read or write the per-track step status file |
| `get_track_step_status()` / `set_track_step_status()` | Status of a step for one track |
| `get_global_step_status()` / `set_global_step_status()` | Status of a global step (`integrate_data`, `generate_report`) |
| `reset_track_step()` | Clears a step's cached status to force a rerun |

## retry_helpers.py

Wraps any callable in an exponential-backoff retry loop. Used for IEDB and UniProt requests.

```python
from utils.retry_helpers import retry_network_call

result = retry_network_call(
    fn=lambda: requests.post(url, data=payload, timeout=120),
    max_attempts=5,
    initial_delay=2.0,
    backoff=2.0,
    description="IEDB NetMHCpan",
)
```

The wrapper retries on connection errors, timeouts, and HTTP 5xx responses. HTTP 4xx errors are raised immediately.

## fasta_utils.py

FASTA-specific helpers. Used by `fetch_sequences` and `predict_binding`.

| Function | Purpose |
|---|---|
| `parse_fasta()` | Parses a FASTA file into a list of `SeqRecord` objects |
| `validate_sequence()` | Checks for ambiguous amino acids and minimum length |
| `remove_duplicate_sequences()` | Drops records with identical amino acid sequences |
| `generate_peptides()` | Generates all overlapping peptides of given lengths from a sequence |

## genbank_utils.py

GenBank-specific helpers, kept around for `search_variants` (planned). The earlier `fetch_sequences` step also used these but has since moved to UniProt.

| Function | Purpose |
|---|---|
| `search_ncbi_protein_ids()` | Returns accession IDs without downloading full records |
| `fetch_records_by_accession_ids()` | Pulls full GenBank records for a list of IDs |
| `extract_source_qualifiers()` | Strain, isolate, host, location, collection date |
| `record_is_polyprotein()` / `record_is_refseq()` | Heuristic flags |
| `extract_protein_from_polyprotein()` | Slices a mature peptide using `mat_peptide` features |
| `save_sequences_as_fasta()` | Saves selected records as FASTA |
| `build_sequence_registry()` | JSON registry with full metadata |

`Entrez.email` must be set before calling any search function, otherwise NCBI will throttle or reject the request.
