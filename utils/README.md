# utils

Shared helpers used across pipeline steps. Each file groups functions by concern, with no thin wrappers and no helpers used by a single caller. The list below is intentionally short — many helpers (private `_foo`) only matter to the file that owns them; only the public surface is documented here.

## archive.py

Packages project outputs into archives. Provides two format pairs — tar.gz (Linux / macOS) and zip (Windows / WSL). The `download_ui` selects the format automatically: `.zip` when running under WSL (natively openable in Windows Explorer), `.tar.gz` otherwise. Stdlib `tarfile` / `zipfile`, no extra dependency.

| Function | Purpose |
|---|---|
| `archive_project(project_name, destination_dir, include_predictions=False)` | Bundles `projects/{project}/` into `{project}_full_{stamp}.tar.gz`. Default excludes `data/intermediate/*/predictions/` (large, re-derivable). |
| `archive_project_zip(project_name, destination_dir, include_predictions=False)` | Same scope, same exclusion logic — writes `{project}_full_{stamp}.zip` (Windows-friendly). |
| `archive_step(project_name, step_name, destination_dir)` | Bundles a single step's outputs across every track into `{project}_{step}_{stamp}.tar.gz`. |
| `archive_step_zip(project_name, step_name, destination_dir)` | Same scope — writes `{project}_{step}_{stamp}.zip` (Windows-friendly). |

## console.py

Centralised Rich `Console` plus interactive-prompt helpers. Every module imports the shared `console` from here instead of constructing its own `Console()`, so width and styling stay consistent across the whole CLI.

| Item | Purpose |
|---|---|
| `console` | Module-level `Console(width=120)` instance |
| `DEFAULT_TABLE_ROW_STYLES` | Zebra striping for long data tables (`["", "on grey7"]`) |
| `is_interactive_session()` | True when `stdin` is a TTY (gates interactive prompts) |
| `ask(text, default=None, choices=None, required=False)` | Rich-styled prompt with TTY fallback; `required=True` rejects empty input |
| `confirm(text, default=False)` | y/n prompt with TTY fallback |
| `press_enter_to_continue(text=...)` | Pauses execution until Enter is pressed (no-op when non-interactive) |
| `flush_stdin()` | Discards pending keystrokes (POSIX-only) so impatient Enter presses don't skip the next prompt |
| `confirm_value(label, value, indent='')` | Echoes a value the user just typed + asks `[y] confirm / [n] re-enter` (returns True/False so the caller can loop) |
| `show_recap_and_confirm(title, fields, proceed_label='...')` | Renders a recap Panel of `(label, value)` pairs and asks a single final y/n. Used right before any wizard saves to disk |
| `format_file_size_human(num_bytes)` | Short human-readable size string (`12.3 MB`) |

## csv_write.py

One-line wrapper around `pandas.DataFrame.to_csv` that all user-facing CSVs go through. Guarantees consistent UTF-8 BOM + numeric quoting so the same files open cleanly in Excel, LibreOffice and Python without quoting surprises. Per project convention, every step writes its CSVs via this helper instead of calling `df.to_csv` directly.

| Function | Purpose |
|---|---|
| `write_user_facing_csv(df, output_path, *, sep=',', quoting=csv.QUOTE_NONNUMERIC, encoding='utf-8-sig')` | Save a DataFrame with portable defaults (UTF-8 BOM for Excel; `QUOTE_NONNUMERIC` so cells with `;` or `,` don't break the schema) |

## fasta_utils.py

FASTA-specific helpers. Used by `fetch_sequences`, `predict_binding`, `search_variants`, and `analyze_conservation`.

| Function | Purpose |
|---|---|
| `write_fasta(records, output_path)` | Writes a list of `SeqRecord` to a FASTA file |
| `has_ambiguous_residues(sequence)` | True if the sequence contains any of `X B Z J U O` |
| `is_valid_sequence(record)` | Returns `(bool, reason)` — checks min length + ambiguous-residue rules |
| `generate_peptides(sequence, lengths)` | Yields every overlapping peptide of the requested lengths (deduplicated) |

## file_browser.py

Rich-styled inspector for step outputs (CSVs, JSONs, FASTAs, XLSXs, HTMLs, tar.gz). Used by `BaseTrackStep.execute()` to offer a post-step "peek at what just got written" menu, and by the top-level `[b]` REPL key for project-wide browsing.

| Function | Purpose |
|---|---|
| `browse_step_outputs(output_descriptions: dict[Path, str])` | Per-step pop-up: lists artefacts with size + description; user picks a file to preview (with format-specific renderer) |
| `run_project_browser(project_name)` | Top-level chooser: Project outputs / Downloads / Tracks, with HTML opener + tar.gz inspector |

## http.py

Centralised HTTP wrapper around `requests.get` with the project's standard timeout + retry policy. Anything that hits an external HTTP endpoint should call this rather than `requests.get` directly.

| Function | Purpose |
|---|---|
| `http_get(url, params=None, max_attempts=3, ...)` | GET with exponential backoff on connection errors / timeouts / HTTP 5xx; HTTP 4xx raises immediately |

## input_validation.py

Per-field-type input validation for interactive prompts. Identifier fields block shell-, glob- and quoting-breaking characters (so values can flow safely into file paths and command lines) and suggest an accent-free NFKD form instead of mangling accents; descriptions allow prose; paths keep separators and drive colons; peptides accept only the 20 standard amino acids. All `validate_*` functions are pure (no I/O) and return a `ValidationResult`; `prompt_validated()` wires a validator to `input()`.

| Function | Purpose |
|---|---|
| `validate_organism_name(raw)` / `validate_protein_name(raw)` | Identifier fields — block `/ \ ( ) < > \| & ; " ' * ?` + control chars; suggest accent-free form |
| `validate_description(raw)` | Free prose — blocks only control chars; empty allowed |
| `validate_local_path(raw)` | Strips wrapping quotes; blocks redirection/glob/quote chars; keeps separators and `:` |
| `validate_peptide(raw)` | Accepts only the 20 standard amino acids (uppercased) |
| `strip_accents(text)` | NFKD diacritic removal helper |
| `prompt_validated(validator, indent='')` | `input()` loop enforcing a validator; re-asks on a hard block, offers the accent-free suggestion; returns `''` non-interactively |

`ValidationResult` carries `ok`, `value`, `suggestion` (accent-free alternative) and `error`. Wired into the project wizard (`project_manager.py`) and the local-FASTA path prompts in `fetch_sequences` / `analyze_conservation`.

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
| `STAR_MARKER` | `"★"` — the Unicode glyph written into the `BEST_REPRESENTATIVE` column by `select_representatives` and read by every downstream step (`analyze_conservation`, `population_coverage`, `predict_murine`, …) |

## output_capture.py

Captures stdout/stderr at the file-descriptor level (not the Python level) so we can swallow noisy logs from third-party libraries — MHCFlurry / TensorFlow in particular — without losing the panel layout of the Rich console. Used inside `predict_binding` and `predict_murine`.

| Item | Purpose |
|---|---|
| `CapturedOutput` | Holds the captured bytes after the `with` block exits |
| `capture_fd_output()` | Context manager: returns a `CapturedOutput` you can read after the block |

## pipeline_state.py

Low-level helpers around `pipeline.json`. Used by `BaseTrackStep.execute()` and the REPL to track per-step status.

| Function | Purpose |
|---|---|
| `load_pipeline_state(project_name)` / `save_pipeline_state(project_name, state)` | Read or write the per-track step status file |
| `get_track_step_status(project_name, track_id, step_key)` / `set_track_step_status(...)` | Status of a step for one track (`'pending'` / `'done'` / `'error'`) |
| `get_global_step_status(project_name, step_key)` / `set_global_step_status(...)` | Same for global steps |
| `reset_track_step(project_name, track_id, step_key)` | Clears a step's cached status to force a rerun |
| `_migrate_legacy_step_keys()` | Auto-strips `stepNN_` prefixes from old `pipeline.json` keys on load |

## project_manager.py

Project lifecycle. Used by `main.py` and most steps when they need to read or save shared config.

| Function | Purpose |
|---|---|
| `create_project_interactive()` | 2-question wizard (project name, optional description); target host is fixed at `config.TARGET_HOST = 'Homo sapiens'` since HLA-I is human-only |
| `setup_project_tracks_interactive()` | Outer save-wrapper that retries `_collect_tracks_interactive` until the final recap is confirmed, then persists to `project_config.json` |
| `_collect_tracks_interactive()` | Walks the user through organisms × proteins with per-field y/n confirmations + per-pair recap + final recap |
| `create_project(project_name, description='')` | Non-interactive project creation (used by `tests/validation/seed_project.py`) |
| `edit_track_interactive(project_name, track_id)` | Lets the user edit one existing track; renames folders and resets step statuses if the track ID changes |
| `load_project_config(project_name)` / `save_project_config(project_name, config)` | Read or write `project_config.json` |
| `list_projects(expected_track_step_names=None)` | Returns project metadata for the menu and `--list` |
| `tracks_are_defined(project_name)` | True if at least one organism/protein pair is configured |
| `delete_project(project_name)` | Permanently removes the project directory + registry entry |
| `update_last_used(project_name)` | Bumps the `last_used` timestamp in the registry |

The wizard never saves to disk on a `[n]` answer at the final recap — collection restarts from scratch.

## retry_helpers.py

Wraps any callable in an exponential-backoff retry loop. Used for IEDB and UniProt requests when the calling code wants fine-grained control (most callers should prefer `http.py` instead).

```python
from utils.retry_helpers import retry_network_call

result = retry_network_call(
    fn=lambda: requests.post(url, data=payload, timeout=120),
    description="IEDB NetMHCpan",
)
```

The wrapper retries on connection errors, timeouts, and HTTP 5xx responses. HTTP 4xx errors are raised immediately.

## step_summary.py

Standardised end-of-step block printed after each step's `run()` returns. Renders a Panel with the step's headline + a small "files written" table sourced from the step's `describe_outputs()` dict.

| Function | Purpose |
|---|---|
| `print_step_summary(title, headline, file_paths=...)` | One-shot Rich Panel for the post-step recap |

## text_format.py

Tiny string helpers for table cells and log lines.

| Function | Purpose |
|---|---|
| `compact_num(value)` | `1234` → `'1.2k'`; `1_500_000` → `'1.5M'`. Returns `''` for non-numeric input |
| `truncate_with_ellipsis(text, max_len=60)` | Trims long strings to fit a table column with `…` suffix |

## uniprot.py

UniProt-specific helpers used by `fetch_sequences` to handle viral polyproteins (DENV, ZIKV, HCV — where the mature peptide is stored as a Chain feature inside one big record).

| Function | Purpose |
|---|---|
| `fetch_uniprot_entry_json(accession)` | Pulls the full UniProt JSON for an accession |
| `score_chain_match(protein_name, chain_description)` | Heuristic 0–1 score for how well a chain feature matches the user-requested protein name |
| `find_chain_for_protein(uniprot_record, protein_name)` | Picks the best chain (or `None`) for slicing a mature peptide out of a polyprotein |

---

`column_finder.py` was absorbed into `naming.py` in April 2026, so older references to `from utils.column_finder import ...` should now point at `utils.naming`.
