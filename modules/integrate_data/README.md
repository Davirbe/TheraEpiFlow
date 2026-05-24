# integrate_data

Global step — runs once after every track finished `curate_murine`. Stacks the per-track master tables into a project-wide table and projects a user-configurable VIEW that feeds `generate_report`.

## Why a global step

`curate_murine` already does the heavy per-track JOIN (★ epitopes + conservation + coverage + murine). What was still missing was a project-level union. `integrate_data` is that union — plus the user choice of which columns matter for the construct calculator.

## SRP layout

| File | Role |
|---|---|
| `__init__.py` | Facade — re-exports `IntegrateDataStep` for `STEP_REGISTRY`. |
| `step.py` | `IntegrateDataStep(BaseGlobalStep)` — orchestrates load → prompt → project → write → audit. Owns the pre-step page ClassVars. |
| `core.py` | Pure pandas — track loading, project-wide stacking, anchor-mutation aggregation, VIEW projection. No Rich, no openpyxl, no `input()`. |
| `io.py` | XLSX writers (FULL + VIEW + CSV sidecar). Owns the canonical palette mirrored from other modules. |
| `prompts.py` | Rich-based numbered-list checkbox UI for VIEW column customization. |

## Outputs (in `data/output/`)

| File | What |
|---|---|
| `MASTER_TABLE_FULL_{project}.xlsx` | Every column from every track, with `track_id` + `organism` + `protein` prepended. Audit-grade dump, no cell coloring. |
| `MASTER_TABLE_VIEW_{project}.xlsx` | User-selected columns, English display headers, canonical palette. No `track_id` — identity is `organism` + `protein`. |
| `MASTER_TABLE_VIEW_{project}.csv` | Same VIEW with display headers, in CSV. |
| `MASTER_TABLE_AUDIT_{project}.json` | Track counts, skipped tracks, populations, chosen columns, generation timestamp. |

## Customization flow

On the first run for a project, the user is asked:

```
Use default VIEW columns? [Y/n]
```

- `Y` → default selection (Organism, Protein, Peptide, Best percentile, HLA alleles bound, Coverage (per population), HLA alleles list, Conservation 100% / @ threshold, Conservation label, Murine binding + tooltip cols).
- `n` → opens a numbered checkbox UI with the full catalog grouped by topic (binding / conservation / coverage / murine). Toggle by typing numbers; Enter confirms.

The choice is persisted to `project_config.json` under:

```json
"step_overrides": {
  "integrate_data": {
    "view_columns": "default"            // or a list of column names
  }
}
```

Reruns reuse the persisted choice silently. To re-open the prompt, run the step interactively in the REPL with `retry` (which sets `reconfigure=True`).

## Optional columns

Off by default; opt in via the checkbox UI:

- `netmhcpan_el_percentile`, `mhcflurry_presentation_percentile` — per-method percentiles for the best allele
- `calis_score` — immunogenicity index (legacy, kept for completeness)
- `pct_identity_90`, `pct_identity_80` — additional conservation thresholds
- `n_mutations_safe`, `n_mutations_risky` — anchor-mutation aggregation read on demand from `CONSERVATION_MUTATIONS_{track_id}.xlsx`

## Reuse map

- `BaseGlobalStep` contract → `modules/base_step.py`
- Track → organism/protein lookup → `utils.naming.parse_track_id`
- XLSX palette → mirrors `modules/select_representatives/__init__.py`, `modules/analyze_conservation/io.py`, `modules/population_coverage/io.py`
- Persistence helper → `utils.project_manager.save_project_config`
- End-of-step summary → `utils.step_summary.print_step_summary`

## Verification

```bash
python main.py --project hpv16     --step integrate_data    # 5-track real case
python main.py --project scer_test --step integrate_data    # 2-track small case
pyflakes modules/integrate_data utils/naming.py
```
