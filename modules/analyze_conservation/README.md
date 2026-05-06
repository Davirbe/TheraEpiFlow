# analyze_conservation

Measures how faithfully each ★ representative epitope appears across the variant sequences produced by `search_variants`. The step is qualitative — no epitopes are removed. All representatives are kept; the output annotates them with conservation metrics.

## Method: sliding window identity

For each `(peptide, variant_sequence)` pair, a window of `len(peptide)` residues is slid across the full variant sequence. The window with the highest per-position identity is recorded:

```
identity = matches / peptide_length    (0.0 – 1.0)
```

The best-matching window and its identity score are stored for every variant.

## Analysis threshold

The threshold (default `1.0` = exact match) is configurable per project. It determines:

- `n_variants_passed_threshold` — count of variants with identity ≥ threshold
- `variants_passed_threshold` — list of passing accessions
- `variants_failed_threshold` — list of failing accessions with their actual best window
- `position_mutation_profile` — aggregated per-position substitutions from failed variants (e.g. `pos3:L>I(2x); pos6:Q>R(1x)`)
- `conservation_summary` — tier summary including the chosen threshold plus fixed 90% and 80% reference tiers

**What does NOT change with threshold:**

- `conservation_label` — always based on `mean_max_identity`
- Row colours in the XLSX — always based on `conservation_label`
- `n_variants_exact_match`, `n_variants_identity_90pct`, `n_variants_identity_80pct` — always computed at fixed reference tiers

The threshold is saved to `project_config["conservation_threshold"]`. It can be changed on re-run by resetting the step.

## Conservation labels (fixed)

| Label | Condition |
|---|---|
| `perfect` | `mean_max_identity == 1.0` |
| `high` | `mean_max_identity >= 0.90` |
| `moderate` | `mean_max_identity >= 0.80` |
| `low` | `mean_max_identity < 0.80` |
| `conservation_unknown` | No variant FASTA available |

## Input

- `clusters/CLUSTER_REPR_{track_id}.csv` — ★ representatives from `select_representatives`
- `variants/VARIANTS_{track_id}.fasta` — variant sequences from `search_variants` (or user-supplied)

When no FASTA is found, interactive mode offers to provide a local path. Non-interactive mode labels all epitopes `conservation_unknown` and continues.

## Output

| File | Contents |
|---|---|
| `conservation/CONSERVATION_{track_id}.csv` | One row per ★ representative with all conservation columns |
| `conservation/CONSERVATION_{track_id}.xlsx` | Same table, rows colour-coded by `conservation_label` |
| `conservation/CONSERVATION_VISUAL_{track_id}.xlsx` | Position heatmap — one row per epitope, one column per position |
| `conservation/CONSERVATION_AUDIT_{track_id}.json` | Run metadata: threshold, FASTA source, label counts, mean identity |

## Visual heatmap XLSX

Each row is a ★ representative, each column is a peptide position (P1–P9 for 9-mers). Cell content:

- Line 1: conservation percentage at that position (e.g. `67%`)
- Line 2: most common substitution if not 100% (e.g. `Q→R`)

**Cell colours by position conservation:**

| Conservation | Colour |
|---|---|
| 100% | Dark green `#00B050` |
| ≥ 90% | Light green `#92D050` |
| ≥ 70% | Yellow `#FFFF99` |
| ≥ 50% | Orange `#FFC000` |
| < 50% | Red `#FF9999` |

Anchor positions P2 and Pc (C-terminus) are highlighted with a blue header — these positions are critical for MHC-I binding. Rows are sorted by anchor score descending (min of P2 and Pc conservation), placing the best vaccine candidates at the top.

## Summary XLSX colour coding

**Header cells:**

- Numeric columns (`n_variants_*`, `mean_max_identity`, `analysis_threshold`) → orange `#FFD966`
- Text columns (`conservation_summary`, `conservation_label`, `variants_*`, `position_mutation_profile`) → gray `#D9D9D9`

**Row background by `conservation_label`:**

| Label | Colour |
|---|---|
| perfect | `#00B050` |
| high | `#92D050` |
| moderate | `#FFFF99` |
| low | `#FF9999` |
| conservation_unknown | `#D9D9D9` |
