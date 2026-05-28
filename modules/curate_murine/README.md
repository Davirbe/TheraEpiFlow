# curate_murine

Assembles the per-track master table by joining each ★ representative's human
qualification with its conservation, population coverage and (when available)
murine prediction.

> **Code layout:** single-file module (`__init__.py`) — small enough that the per-role split (see the top-level README) would only fragment it. Splits along the standard seams (core/io/step) if it grows.

## Where it sits in the pipeline

```
analyze_conservation + population_coverage (+ predict_murine)  →  curate_murine
```

It is the last per-track step: everything downstream (`integrate_data`,
`generate_report`) is global.

## What it does

A **JOIN-only** step — it does not re-rank, re-score or relabel anything. It
left-joins, on `peptide`, the ★ representatives with:

- **conservation** (`CONSERVATION_{track_id}.csv`) — **required**
- **population coverage** (`COVERAGE_{track_id}.csv`, pivoted long→wide) — **required**
- **murine** (`MURINE_AGG_{track_id}.csv`) — **optional** (only if `predict_murine` ran)

The `★` marker and all upstream columns are preserved verbatim.

## Inputs

- `clusters/CLUSTER_REPR_{track_id}.csv` — ★ rows (human qualification).
- `conservation/CONSERVATION_{track_id}.csv` — conservation label + identity.
- `coverage/COVERAGE_{track_id}.csv` — long-format coverage (pivoted per population).
- `murine/MURINE_AGG_{track_id}.csv` — optional aggregated murine row per peptide.

## Outputs (`murine/`)

| File | Description |
|---|---|
| `CURATE_MURINE_{track_id}.csv` | One row per ★ epitope with every annotation joined. |
| `CURATE_MURINE_VIEW_{track_id}.csv` | Slim view — key columns only. |
| `CURATE_MURINE_AUDIT_{track_id}.json` | Join provenance + label histogram. |

## Notes

- Qualitative — never removes epitopes.
- Conservation and coverage are hard requirements; the step errors clearly if
  either is missing. Murine is optional and simply omitted when absent.
- This per-track table is what `integrate_data` (global) stacks across all tracks
  into `output/master_table.xlsx`.
