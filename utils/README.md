# utils

Shared helpers used across pipeline steps. Each file groups functions by concern, with no thin wrappers and no helpers used by a single caller.

## console.py

Centralised Rich `Console` plus interactive-prompt helpers. Every module imports the shared `console` from here instead of constructing its own `Console()`, so width and styling stay consistent across the whole CLI.

| Item | Purpose |
|---|---|
| `console` | Module-level `Console(width=120)` instance |
| `is_interactive_session()` | True when `stdin` is a TTY (gates interactive prompts) |
| `ask(text, default=None, choices=None, required=False)` | Rich-styled prompt with TTY fallback; `required=True` rejects empty input |
| `confirm(text, default=False)` | y/n prompt with TTY fallback |
| `press_enter_to_continue(text=...)` | Pauses execution until Enter is pressed (no-op when non-interactive) |
| `confirm_value(label, value, indent='')` | Echoes a value the user just typed + asks `[y] confirm / [n] re-enter` (returns True/False so the caller can loop) |
| `show_recap_and_confirm(title, fields, proceed_label='...')` | Renders a recap Panel of `(label, value)` pairs and asks a single final y/n. Used right before any wizard saves to disk |

## naming.py

Standard naming for files, columns, and aliases. Used everywhere.

| Item | Pattern or behavior |
|---|---|
| `build_track_id(organism_label, protein_label)` | Returns `"HPV16_E6"` |
| `get_prediction_filename(tool, track_id)` | Returns `"PRED_NET_HPV16_E6.csv"` for `tool="NET_PRED"` etc. |
| `get_step_filename(step, track_id, tool="", ext="csv")` | Returns `"CONSENSUS_HPV16_E6.csv"`, `"TOXICITY_AUDIT_HPV16_E6.json"`, etc. |
| `allele_to_netmhcpan_format(allele)` | `"HLA-A*02:01"` becomes `"HLA-A0201"` (strips `*` and `:`) |
| `find_column_name(df, candidates)` | Returns the first matching column name, or `None` |

Constants (English, prefix-only — the project-wide column-naming convention):

| Name | Value |
|---|---|
| `COLUMN_PEPTIDE` | `"peptide"` |
| `COLUMN_NETMHC_BEST_ALLELE` | `"netmhcpan_best_allele"` |
| `COLUMN_NETMHC_EL_PERCENTILE` | `"netmhcpan_el_percentile"` |
| `COLUMN_NETMHC_ALLELES` | `"netmhcpan_alleles"` |
| `COLUMN_NETMHC_NUM_ALLELES` | `"netmhcpan_num_alleles"` |
| `COLUMN_NETMHC_EL_PERCENTILES_ALL` | `"netmhcpan_el_percentiles_all"` |
| `COLUMN_FLURRY_BEST_ALLELE` | `"mhcflurry_best_allele"` |
| `COLUMN_FLURRY_PERCENTILE` | `"mhcflurry_presentation_percentile"` |
| `COLUMN_FLURRY_ALLELES` | `"mhcflurry_alleles"` |
| `COLUMN_FLURRY_NUM_ALLELES` | `"mhcflurry_num_alleles"` |
| `COLUMN_FLURRY_PERCENTILES_ALL` | `"mhcflurry_presentation_percentiles_all"` |
| `COLUMN_ALLELES_UNITED` | `"alleles_united"` |
| `COLUMN_NUM_ALLELES_UNITED` | `"num_alleles_united"` |
| `COLUMN_BEST_REPRESENTATIVE` | `"BEST_REPRESENTATIVE"` |
| `STAR_MARKER` | `"★"` — the Unicode glyph written into the `BEST_REPRESENTATIVE` column by `select_representatives` and read by every downstream step (`analyze_conservation`, `population_coverage`, `predict_murine`, …). |

`column_finder.py` was absorbed into this file in April 2026, so older references to `from utils.column_finder import ...` should now point at `utils.naming`.

## project_manager.py

Project lifecycle. Used by `main.py` and most steps when they need to read or save shared config.

| Function | Purpose |
|---|---|
| `create_project_interactive()` | 3-question wizard (project name, description, Entrez email); target host is fixed at `config.TARGET_HOST = 'Homo sapiens'` since HLA-I is human-only |
| `setup_project_tracks_interactive()` | Outer save-wrapper that retries `_collect_tracks_interactive` until the final recap is confirmed, then persists to `project_config.json` |
| `_collect_tracks_interactive()` | Walks the user through organisms × proteins with per-field y/n confirmations + per-pair recap + final recap. Returns the tracks dict; the caller saves it |
| `create_project()` | Non-interactive project creation |
| `edit_track_interactive()` | Lets the user edit one existing track and rename / re-fetch as needed |
| `load_project_config()` / `save_project_config()` | Read or write `project_config.json` |
| `list_projects()` | Returns project metadata for the menu and `--list` |
| `tracks_are_defined()` | True if at least one organism/protein pair is configured |
| `delete_project()` | Permanently removes the project directory |
| `update_last_used()` | Bumps the `last_used` timestamp in the registry |
| `_suggest_organism_label()` | Default abbreviation for known organisms (HPV16, ZIKV, DENV2, SARS2, MPOX, CHIKV, etc.) |
| `_suggest_protein_label()` | Default protein abbreviation (E6, E7, E, NS1, NS5, S, HA, etc.) |
| `_render_prompt_with_example()` | Renders a Rich two-column layout (explanation left, filled-in example right) above a prompt |
| `_prompt_required_nonempty()` | `input()` that re-asks until the user types something non-empty |

The wizard never saves to disk on a `[n]` answer at the final recap — collection restarts from scratch.

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
