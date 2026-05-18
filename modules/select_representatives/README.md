# select_representatives

Selects the single best epitope representative from each cluster produced by `cluster_epitopes`.

## Where it sits in the pipeline

```
cluster_epitopes  →  select_representatives  →  search_variants
```

Every cluster — whether it contains one peptide (singleton) or many — contributes exactly one representative to the next step. Singletons pass through automatically; for multi-member clusters the step scores each candidate and marks the winner with ★.

## Scoring

Two independent metrics are computed for every peptide across the whole track, then combined:

### 1. `best_combined_percentile`
The minimum value across **all** allele-level percentile values from both tools:

```
best_combined_percentile = min(
    netmhcpan_el_percentiles_all   (semicolon-separated)
    + mhcflurry_presentation_percentiles_all
)
```

This captures the absolute peak binding quality: the strongest single MHC-I binding event the peptide can achieve, regardless of which allele or tool produced it. A peptide that is an exceptional binder for even one allele will score well here.

### 2. `num_alleles_united`
The number of unique HLA alleles in the union of both tools:

```
alleles_united = set(netmhcpan_alleles) ∪ set(mhcflurry_alleles)
num_alleles_united = |alleles_united|
```

This captures promiscuity: a peptide that binds many alleles is more likely to cover a broader population.

### 3. Final score

Both metrics are min-max normalised over the whole track (percentile is inverted, since lower is better):

```
norm_best_percentile = (max_pct - pct) / (max_pct - min_pct)   → 1.0 if all equal
norm_alleles         = (n - min_n)     / (max_n - min_n)        → 1.0 if all equal

final_score = (norm_best_percentile + norm_alleles) / 2
```

The peptide with the highest `final_score` in each cluster receives `BEST_REPRESENTATIVE = ★`. Ties are broken by row order (first row wins).

## Why minimum percentile, not mean?

Using the mean across all tested alleles penalises a peptide that binds HLA-A\*02:01 (present in ~45 % of the global population) at 0.1 % but performs poorly on less frequent alleles. The minimum captures the single strongest binding event without that penalty. Population-weighted coverage analysis is handled downstream in `population_coverage`.

## Input

`{track_dir}/clusters/CLUSTER_{track_id}.csv` produced by `cluster_epitopes`.

Required columns: `peptide`, `cluster_id`, `netmhcpan_el_percentiles_all`, `mhcflurry_presentation_percentiles_all`, `netmhcpan_alleles`, `mhcflurry_alleles`.

## Output

Saved to `{track_dir}/clusters/`:

| File | Contents |
|---|---|
| `CLUSTER_REPR_{track_id}.csv` | All input rows with scoring columns and ★ added — feeds every downstream step. |
| `CLUSTER_REPR_{track_id}.xlsx` | Same data, colour-coded for review |
| `REPRESENTATIVES_VIEW_{track_id}.csv` | Slim per-step view — `peptide, cluster_id, BEST_REPRESENTATIVE (★), final_score` only. |
| `CLUSTER_REPR_AUDIT_{track_id}.json` | Run summary and parameters |

### Columns added

| Column | Description |
|---|---|
| `best_combined_percentile` | Minimum EL% across all allele–tool pairs |
| `alleles_united` | Semicolon-separated union of NetMHCpan + MHCFlurry alleles |
| `num_alleles_united` | Count of unique alleles in the union |
| `norm_best_percentile` | Min-max normalised, inverted (0–1, higher = better binder) |
| `norm_alleles` | Min-max normalised allele count (0–1, higher = more promiscuous) |
| `final_score` | Average of `norm_best_percentile` and `norm_alleles` (0–1) |
| `BEST_REPRESENTATIVE` | `★` for the highest-scoring peptide per cluster; empty otherwise |

### Excel colour scheme

| Colour | Columns |
|---|---|
| Orange header | Percentile-related columns |
| Pink header | HLA/allele-related columns |
| Grey header | All other columns |
| Yellow row | Best representative (★) |

## Audit JSON

```json
{
  "timestamp": "2026-05-04T14:00:00.000000",
  "track_id": "HPV16_E5",
  "input_file": ".../CLUSTER_HPV16_E5.csv",
  "n_epitopes": 9,
  "n_clusters": 4,
  "n_multi_member_clusters": 3,
  "n_singletons": 1,
  "n_representatives": 4,
  "scoring": {
    "best_combined_percentile": "min of all values in netmhcpan_el_percentiles_all + mhcflurry_presentation_percentiles_all",
    "norm_best_percentile": "min-max normalised, inverted (lower percentile → higher score)",
    "alleles_united": "set-union of netmhcpan_alleles and mhcflurry_alleles",
    "norm_alleles": "min-max normalised num_alleles_united",
    "final_score": "(norm_best_percentile + norm_alleles) / 2"
  },
  "output_csv": ".../CLUSTER_REPR_HPV16_E5.csv",
  "output_xlsx": ".../CLUSTER_REPR_HPV16_E5.xlsx"
}
```

## Downstream dependency

The `★` symbol in `BEST_REPRESENTATIVE` is the key that downstream steps (`population_coverage`, `curate_murine`, `integrate_data`) use to filter for selected candidates. Do not rename or remove this column.
