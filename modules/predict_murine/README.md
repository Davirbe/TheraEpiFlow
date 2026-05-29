# predict_murine

Re-runs the two MHC-I binding predictors with H-2 (murine) alleles against the ★ representatives so the human-selected epitopes can be ranked by how well they translate to mouse models. **Qualitative — no peptides are removed.**

## Where it sits in the pipeline

```
select_representatives  →  predict_murine  →  curate_murine
```

Used by `curate_murine` to flag candidates that are both strong human binders and translatable to mouse strains — the "MAXIMUM priority" set for in vivo validation.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `PredictMurineStep` orchestration — runs NetMHCpan + MHCFlurry, writes tables |
| `core.py` | ★ peptide loading, synthetic records, tier labels, per-peptide aggregation |
| `prompts.py` | Interactive murine-strain selection |
| `__init__.py` | Facade — re-exports `PredictMurineStep` |

## What it does

For every ★ representative epitope:

1. Predicts MHC-I presentation against the configured H-2 allele set using **both** NetMHCpan EL (via the IEDB classic HTTP API) and **MHCFlurry 2.0** (local Python package).
2. Operates in **peptide-direct mode** — the input is the ★ peptide list, not a FASTA. No k-mer enumeration here; the lengths are inherited from `predict_binding`.
3. Aggregates per peptide: best percentile across both tools, the list of H-2 alleles bound (best first), the count, and a 4-tier binder label.

## Strain groups

| Group | Alleles | Use case |
|---|---|---|
| `C57BL/6` | H-2Kb, H-2Db | Most common immunology lab strain |
| `BALB/c` | H-2Kd, H-2Dd, H-2Ld | Th2-biased strain, very common for vaccine work |
| `complete` | All H-2 alleles available in both predictors | Maximum coverage, slower |

The chosen group is asked once and saved to `project_config["murine_strain_group"]`.

## Binder tiers

Per-peptide label based on the best percentile rank across all H-2 alleles tested (lower = stronger):

| Tier (`best_percentile_label`) | Condition | Meaning |
|---|---|---|
| `optimal` | min %rank ≤ 0.5 | Strong binder, top priority for in vivo |
| `good` | min %rank ≤ 2.0 | Established binder |
| `borderline` | min %rank ≤ 2.5 | Marginal binder |
| `non_binder` | min %rank > 2.5 | Not translatable to this strain |

The label is the best (lowest) percentile rank a peptide achieves across all H-2 alleles tested, taking the better of the two predictors per (peptide, allele) pair. The number of H-2 alleles bound at tier `borderline` or better is reported as `num_murine_alleles_bound`.

Thresholds live in `config.py`: `MURINE_OPTIMAL_BINDER_RANK_MAX` (0.5), `MURINE_STRONG_BINDER_EL_RANK` (2.0), `MURINE_BORDERLINE_BINDER_RANK_MAX` (2.5).

## Inputs

- `clusters/CLUSTER_REPR_{track_id}.csv` (from `select_representatives`) — only rows with `BEST_REPRESENTATIVE == "★"` are processed.

Required columns on the input: `peptide`, `alleles_united` (used for cross-checks), `BEST_REPRESENTATIVE`.

## Outputs

Under `data/intermediate/{track_id}/murine/`:

| File | Contents |
|---|---|
| `MURINE_{track_id}.csv` | Long format — one row per `(peptide, allele, tool)` with raw percentile. |
| `MURINE_AGG_{track_id}.csv` | Slim, one row per ★ peptide: `peptide`, `best_percentile_label` (the tier — `optimal`/`good`/`borderline`/`non_binder`), `best_percentile_value` (the numeric rank), `murine_alleles_bound` (semicolon-joined, best first), `num_murine_alleles_bound`. **This is the file `curate_murine` reads.** |
| `MURINE_VIEW_{track_id}.csv` | Slim per-step view — `peptide`, `best_percentile_label`, `best_percentile_value`, `num_murine_alleles_bound`. |
| `MURINE_AUDIT_{track_id}.json` | Strain group used, allele list, peptide lengths, and `label_counts` per tier (optimal/good/borderline/non_binder). |

## External tools — provenance and version

Same dependencies as `predict_binding`:

| Tool | Source | Notes |
|---|---|---|
| **NetMHCpan EL** | IEDB classic API | One HTTP POST per H-2 allele, `method=netmhcpan_el`. Subject to IEDB availability. |
| **MHCFlurry** | `mhcflurry` package (≥ 2.0) | Reuses the same `Class1PresentationPredictor` model load as `predict_binding`; the first call in a session pays the ~30s TF model load. |

The IEDB endpoint accepts H-2 alleles in NetMHCpan compact form (e.g. `H-2-Kb`); the conversion happens in `utils.naming.allele_to_netmhcpan_format`.

## References

Same predictors as `predict_binding`:

- Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. *NetMHCpan-4.1 and NetMHCIIpan-4.0.* Nucleic Acids Research. 2020;48(W1):W449–W454. doi:10.1093/nar/gkaa379
- O'Donnell TJ, Rubinsteyn A, Laserson U. *MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing.* Cell Systems. 2020;11(1):42–48.e7. doi:10.1016/j.cels.2020.06.010

## Operational notes

- Pick the strain group your wet lab actually uses — running `complete` is heavier (more API calls + more MHCFlurry batches) and rarely yields extra useful candidates.
- Peptides labelled `optimal` or `strong` are immediate candidates for in vivo validation. `borderline` warrants an experimental check, not skipped outright.
- Because the predictors are the same as `predict_binding`, network failures and the first-load MHCFlurry latency apply here too. Both are retried automatically (`utils.retry_helpers.retry_network_call`).
- Input is peptide-direct, not from FASTA — no length re-enumeration happens. If you want to test additional lengths, rerun `predict_binding` with a wider `peptide_lengths` first.
