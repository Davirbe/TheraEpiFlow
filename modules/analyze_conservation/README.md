# analyze_conservation

Measures how faithfully each ★ representative epitope appears across the variant sequences produced by `search_variants`. The step is qualitative — no epitopes are removed. All representatives are kept; the output annotates them with conservation metrics and (where applicable) per-variant mutation verdicts.

## Method: sliding window identity

For each `(peptide, variant_sequence)` pair, a window of `len(peptide)` residues is slid across the full variant sequence. The window with the highest per-position identity is recorded:

```
identity = matches / peptide_length    (0.0 – 1.0)
```

The best-matching window and its identity score are stored for every variant.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `AnalyzeConservationStep` orchestration — preflight (FASTA inspection) / run / postflight |
| `core.py` | BLOSUM62 scoring, MHC-I anchor verdicts, sliding-window identity, position stats, summary builder |
| `io.py` | FASTA loading + conservation/mutations XLSX writers (colour palette) |
| `charts.py` | Dual-panel conservation PNG |
| `prompts.py` | Identity threshold + local-FASTA override prompts |
| `render.py` | Rich conservation table + preflight FASTA-status table |

## Analysis threshold

The threshold (default `1.0` = exact match) is configurable per project and saved to `project_config["conservation_threshold"]`. It feeds the `pct_identity_threshold` column in the summary and the `≥threshold` column in the heatmap PNG. The fixed tiers (100%, 90%, 80%) are always present regardless of threshold.

## Conservation labels (fixed)

| Label | Condition |
|---|---|
| `perfect` | `mean_max_identity == 1.0` |
| `high` | `mean_max_identity >= 0.90` |
| `moderate` | `mean_max_identity >= 0.80` |
| `low` | `mean_max_identity < 0.80` |
| `conservation_unknown` | No variant FASTA available |

## Mutation tolerance verdict (MHC-I heuristic)

For every `(epitope, variant)` pair where the variant's best window has **1 or 2 substitutions**, the step emits a verdict combining:

- **Anchor flag**: does any substitution land at P2 or PΩ (the universal MHC-I anchor positions, B and F pockets)?
- **BLOSUM62 chemistry**: minimum BLOSUM62 score across the substitutions. `≥ 0` = conservative (preserves polarity, size, hydrophobicity).

| Verdict | Condition | Colour |
|---|---|---|
| `excellent_match` | 1 mutation, non-anchor, BLOSUM62 ≥ 0 | Green |
| `tolerated` | 1–2 non-anchor mutations, otherwise | Yellow |
| `likely_lost` | Any mutation in P2 or PΩ | Red |

This is a **screening heuristic**, not re-prediction. Pairs with 0 or > 2 mutations are excluded from the mutation file.

## Input

- `clusters/CLUSTER_REPR_{track_id}.csv` — ★ representatives from `select_representatives`
- `variants/VARIANTS_{track_id}.fasta` — variant sequences from `search_variants` (or user-supplied)

When no FASTA is found, interactive mode offers to provide a local path. Non-interactive mode labels all epitopes `conservation_unknown` and continues.

## Output

| File | Contents |
|---|---|
| `conservation/CONSERVATION_{track_id}.csv` | IEDB-style summary: one row per ★ rep with tier % + fraction, min/max/avg identity, label |
| `conservation/CONSERVATION_{track_id}.xlsx` | Same table, header gray, only `conservation_label` cell coloured |
| `conservation/CONSERVATION_VIEW_{track_id}.csv` | Slim per-step view — `peptide, length, min_identity, max_identity, avg_identity, conservation_label` only |
| `conservation/CONSERVATION_HEATMAP_{track_id}.png` | Dual-panel heatmap: position conservation + identity tiers |
| `conservation/CONSERVATION_MUTATIONS_{track_id}.xlsx` | Per (epitope, variant) breakdown for ≤ 2-mut variants, row coloured by verdict |
| `conservation/CONSERVATION_AUDIT_{track_id}.json` | Run metadata, label counts, verdict counts |

## Summary CSV/XLSX columns

| Column | Type | Example |
|---|---|---|
| `#` | int | `1` |
| `peptide` | str | `YLQPRTFLL` |
| `length` | int | `9` |
| `pct_identity_100` | str | `39.34% (24/61)` |
| `pct_identity_90` | str | `95.08% (58/61)` |
| `pct_identity_80` | str | `100.00% (61/61)` |
| `pct_identity_threshold` | str | matches user-chosen threshold |
| `min_identity` | str | `26.67%` |
| `max_identity` | str | `100.00%` |
| `avg_identity` | str | `96.62%` |
| `conservation_label` | str | `perfect` / `high` / `moderate` / `low` |
| `n_excellent_match` | int | `2` — pairs with 1 mut, non-anchor, BLOSUM62 ≥ 0 |
| `n_tolerated` | int | `0` — pairs 1-2 mut, non-anchor, otherwise |
| `variants_exact_match` | str | `A0A8B1JN83(SARS-CoV-2); ...` — variants at identity 100% |
| `variants_tolerable` | str | `A0A8B1JJ74[P5:A→S]; ...` — variants with verdict ∈ {excellent, tolerated} |
| `alleles_united` | str | `HLA-A*02:01;HLA-A*68:02` |
| `num_alleles_united` | int | `5` |

The `variants_exact_match` and `variants_tolerable` columns are key for **multi-target vaccine design**: filter the CSV for "contains `HPV31`" to find epitopes that cross-protect against multiple strains.

## Heatmap PNG layout

Single figure, three axes sharing the Y axis (peptides):

```
[label sidebar] [position conservation P1..PΩ] [identity tiers]
```

- **Filter**: epitopes with at least one variant having ≤ 2 substitutions. Perfect-conserved epitopes are included (0 mutations qualifies). Highly divergent epitopes (no variant within 2 mutations) are excluded.
- **Order**: anchor conservation score descending (best vaccine candidates on top).
- **Position panel**: cell colour = % conservation at that position; P2 and PΩ marked with thicker blue borders.
- **Tier panel**: 4 columns (`≥100%`, `≥90%`, `≥80%`, `≥threshold`). Cell text shows `n/N`. Same colour scale as the position panel.

## Mutation XLSX columns

| Column | Example |
|---|---|
| `peptide` | `YLQPRTFLL` |
| `length` | `9` |
| `variant_accession` | `QHD43416.1` |
| `variant_window` | `YLQPRTFLM` |
| `identity` | `88.89% (8/9)` |
| `n_mutations` | `1` |
| `mutations` | `P9:L→M` |
| `mutation_positions` | `9` |
| `anchor_hit` | `Y` (P2 or PΩ?) |
| `blosum62_min` | `2` |
| `mhc_verdict` | `excellent_match` / `tolerated` / `likely_lost` |
| `alleles_united` | `HLA-A*02:01;HLA-A*68:02` |

Rows sorted: verdict (excellent → tolerated → likely_lost), then peptide, then `n_mutations`.
