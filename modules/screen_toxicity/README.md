# screen_toxicity

Removes toxic peptides from the immunogenic candidate list using ToxinPred3.

> **Code layout:** single-file module (`__init__.py`) — small enough that the per-role split (see the top-level README) would only fragment it. Splits along the standard seams (core/io/prompts/step) if it grows.

## Where it sits in the pipeline

```
consensus_filter  →  screen_toxicity  →  cluster_epitopes
```

Toxicity is screened before clustering on purpose. Clustering toxic peptides only to throw them away later wastes effort and can pull the cluster representative toward a peptide that would never make it to the vaccine anyway.

## How prediction works

Instead of calling the ToxinPred3 standalone tool through subprocess (which writes temp files for AAC, DPC, motifs, predictions, then reads them back), this step uses the underlying scikit-learn model directly:

1. `importlib.util.find_spec('toxinpred3')` locates the installed package.
2. `joblib.load()` reads `toxinpred3.0_model.pkl` (Extra Trees classifier, trained by the ToxinPred3 authors on AAC + DPC features).
3. For each peptide we build the feature vector inline: 20 amino-acid composition values plus 400 dipeptide composition values, in percentage units.
4. `clf.predict_proba()` returns the toxicity probability. The PPV column is computed with the same linear calibration the upstream tool applies (`PPV = score * 1.2341 - 0.1182`, clipped to [0, 1]). Coefficients come from `toxinpred3/python_scripts/toxinpred3.py:346` (Model 1, no MERCI).

The tool is faster this way (no I/O, no perl invocation for the MERCI motif step we do not use), and the output matches what ToxinPred3 standalone would produce for Model 1.

## Threshold

Default is 0.38, asked the first time the step runs:

| Option | Value | Notes |
|---|---|---|
| `[1]` | 0.38 | ToxinPred3 default for short peptides |
| `[2]` | 0.50 | More restrictive, fewer false toxics, more false safes |
| `[3]` | custom | Any value in (0, 1] |

The chosen threshold is saved to `project_config['toxicity_threshold']`. Reruns reuse it without prompting again. Use `r` (rerun) in the REPL or remove the key from `project_config.json` if you want a fresh prompt.

## Defensive cleanup

Before prediction the step drops:

- Rows with `NaN` peptide
- Peptides that are not strings or are empty
- Peptides containing residues outside the canonical 20 (`X`, `B`, `Z`, etc.)

These would either crash the AAC vector (division by zero on empty strings) or produce silently truncated feature vectors. The number of dropped rows is reported in the audit JSON as `n_dropped_invalid`.

## Input

`{track_dir}/consensus/CONSENSUS_IMMUNOGENIC_{track_id}.csv` produced by `consensus_filter`. The only column the step touches by name is `peptide`; everything else is carried through to the output.

## Output

Saved to `{track_dir}/toxicity/`:

| File | Contents |
|---|---|
| `TOXICITY_ALL_{track_id}.csv` | Every input row, with the three new columns added |
| `TOXICITY_SAFE_{track_id}.csv` | Only the `Non-Toxin` rows (used by the next step) |
| `TOXICITY_VIEW_{track_id}.csv` | Slim per-step view — `peptide, toxinpred3_score, toxinpred3_ppv, toxinpred3_label` only |
| `TOXICITY_AUDIT_{track_id}.json` | Run timestamp, threshold, counts, output paths |

Three columns are added:

| Column | Range | Meaning |
|---|---|---|
| `toxinpred3_score` | 0 to 1 | Raw probability from the Extra Trees classifier |
| `toxinpred3_ppv` | 0 to 1 | Calibrated PPV, comparable to the ToxinPred3 web tool output |
| `toxinpred3_label` | `Toxin` or `Non-Toxin` | Verdict at the chosen threshold |

## Audit JSON

```json
{
  "timestamp": "2026-05-01T10:32:11.452103",
  "track_id": "HPV16_E6",
  "tool": "ToxinPred3",
  "model": "AAC+DPC Extra Trees (Model 1)",
  "threshold": 0.38,
  "n_input": 218,
  "n_dropped_invalid": 0,
  "n_toxic": 14,
  "n_safe": 204,
  "pct_removed": 6.4,
  "output_all": ".../TOXICITY_ALL_HPV16_E6.csv",
  "output_safe": ".../TOXICITY_SAFE_HPV16_E6.csv"
}
```

## Reproducibility notes

- `scikit-learn` is pinned to `1.2.2` in `environment.yml`. ToxinPred3's `model.pkl` was serialized with this version, and newer scikit-learn releases change internal Cython structures that break unpickling.
- The PPV linear coefficients are not magic numbers, they are reproduced verbatim from upstream (commented at the call site).
- The chosen threshold lives in `project_config.json`. Reruns are deterministic.
