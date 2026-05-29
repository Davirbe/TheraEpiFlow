# predict_binding

Runs two MHC-I binding predictors in parallel against every peptide × HLA allele combination for a track. The two predictors capture complementary signal:

- **NetMHCpan 4.1 EL** — eluted-ligand likelihood scoring, accessed via the IEDB classic HTTP API (`http://tools-cluster-interface.iedb.org/tools_api/mhci/`). No local install required; the API does the work and returns a tab-separated table per allele.
- **MHCFlurry 2.0** — presentation predictor running entirely locally through the `mhcflurry` Python package (`Class1PresentationPredictor`). Requires `mhcflurry >= 2.0` and a one-time `mhcflurry-downloads fetch` to pull the trained models.

The two are kept separate (no consensus yet) so the downstream `consensus_filter` step can decide what to do with disagreement.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `PredictBindingStep` orchestration — `run` / `describe_outputs` |
| `core.py` | NetMHCpan + MHCFlurry runners, binder-count summaries, TF-warning suppression |
| `io.py` | FASTA loading / peptide expansion |
| `prompts.py` | Interactive allele + peptide-length selection |
| `__init__.py` | Facade — re-exports `PredictBindingStep` and the two runners |

The two runners (`_run_netmhcpan_iedb_silent`, `_run_mhcflurry_with_progress`) are re-exported from `__init__.py` because `predict_murine` reuses them.

## Inputs

- `data/input/{track_id}/SEQUENCES_{track_id}.fasta` — the reference FASTA written by `fetch_sequences`. If missing, the step falls back to an interactive prompt for a local path or pasted sequence.

## Configuration

Asked once and saved to `project_config.json`:

| Key | Default | Description |
|---|---|---|
| `hla_alleles` | `config.DEFAULT_HLA_ALLELES` (27 HLA-A / HLA-B alleles) | IMGT format, e.g. `["HLA-A*02:01", "HLA-B*07:02"]`. |
| `peptide_lengths` | `[9]` | Any subset of `[8, 9, 10, 11, 12]` for MHC-I. |

Both predictors require the IMGT format (`HLA-A*02:01`). The wizard validates each user-entered allele via `utils.naming.parse_hla_allele`: missing-asterisk and lowercase inputs are auto-corrected with user confirmation; unparseable tokens re-prompt the whole list.

## Outputs

Under `data/intermediate/{track_id}/predictions/`:

| File | Schema |
|---|---|
| `PRED_NET_{track_id}.csv` | `allele, peptide, netmhcpan_el_score, netmhcpan_el_percentile, seq_num, start, end, length, core, icore` |
| `PRED_FLURRY_{track_id}.csv` | `peptide, allele, mhcflurry_presentation_percentile` |
| `PREDICT_AUDIT_{track_id}.json` | Run metadata, allele list, peptide-length list, row counts, strong/intermediate-binder counts, peptides found by both tools. |

Both CSVs have a multi-row schema (one row per peptide × allele) — they are consolidated to one row per peptide by `consensus_filter`.

## Thresholds (informational)

The thresholds are not applied here — they are reported in the audit JSON and printed in the run summary so the user has a sense of the data before the next step.

| Constant (in `config.py`) | Value | Meaning |
|---|---|---|
| `NETMHC_STRONG_BINDER_EL_RANK` | `0.5` | NetMHCpan EL %Rank ≤ 0.5 is reported as "strong binder". |
| `NETMHC_WEAK_BINDER_EL_RANK` | `2.0` | NetMHCpan EL %Rank ≤ 2.0 is reported as "intermediate binder" (renamed from NetMHC's "weak" to avoid implying poor T-cell response). |
| `CONSENSUS_NETMHCPAN_EL_RANK_MAX_PERCENT` | `2.0` | Hard cut applied later by `consensus_filter`. |
| `CONSENSUS_MHCFLURRY_PRESENTATION_PERCENTILE_MAX` | `2.0` | Hard cut applied later by `consensus_filter`. |

## External tools — provenance and version

| Tool | Source | Version | How invoked |
|---|---|---|---|
| **NetMHCpan EL** | IEDB HTTP API `http://tools-cluster-interface.iedb.org/tools_api/mhci/` | 4.1 (whatever the IEDB endpoint currently serves) | One POST per allele, parameters `method=netmhcpan_el`. |
| **MHCFlurry** | `mhcflurry` package on PyPI | `>= 2.0` (requires `Class1PresentationPredictor`) | In-process Python API; one batch call per track. |

**Note.** Because NetMHCpan is a third-party API, the run depends on IEDB availability. Network errors are caught, recorded in the audit JSON, and surfaced to the user via `retry_network_call` (in `utils.retry_helpers`).

## First-load latency and warning suppression

The first call to `Class1PresentationPredictor.load()` in a session deserialises TensorFlow weights and may take 20–40 s; the progress description reflects this with `"loading models (first run may take ~30s)"`. The load is wrapped in `retry_network_call` (`max_attempts=2`) so a transient file-system or model-cache error does not abort the step.

`UserWarning` is suppressed at module load (broader than the previous `module="tensorflow"` filter) because MHCFlurry's import chain triggers a `pkg_resources is deprecated` notice from setuptools that otherwise collides with the Rich spinner output and looks like a freeze.

## References

The two predictors used by this step:

- Reynisson B et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data.* Nucleic Acids Research. 2020;48(W1):W449–W454. doi:10.1093/nar/gkaa379
- O'Donnell TJ, Rubinsteyn A, Laserson U. *MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating Antigen Processing.* Cell Systems. 2020;11(1):42–48.e7. doi:10.1016/j.cels.2020.06.010
