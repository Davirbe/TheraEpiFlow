# consensus_filter

Takes the two prediction tables from `predict_binding` (NetMHCpan EL and MHCFlurry presentation, both at the peptide × allele level) and produces a single consensus table at the peptide level, then trims it by IEDB Calis 2013 immunogenicity scoring.

The step has two stages:

1. **Presentation filter** — keep only peptides predicted as binders by *both* tools (`%Rank ≤ 2.0` for each), consolidate to one row per peptide carrying the joined HLA list.
2. **Immunogenicity filter** — score every surviving peptide with the IEDB Class I immunogenicity tool (Calis 2013) and drop anything with a score at or below zero.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `ConsensusFilterStep` orchestration — `run` / `describe_outputs` |
| `core.py` | Column resolution, Stage-1 presentation filter, Stage-2 Calis scoring |
| `prompts.py` | Interactive threshold selection |
| `render.py` | Rich progressive per-stage tables |
| `predict_immunogenicity.py` | Vendored Calis 2013 implementation (**DO NOT MODIFY**) |

## Inputs

Under `data/intermediate/{track_id}/predictions/`:

- `PRED_NET_{track_id}.csv` — NetMHCpan EL multi-row schema (one row per peptide × allele).
- `PRED_FLURRY_{track_id}.csv` — MHCFlurry presentation multi-row schema.

## Configuration

Asked once and saved to `project_config["consensus_threshold"]`. Default: the `CONSENSUS_NETMHCPAN_EL_RANK_MAX_PERCENT` and `CONSENSUS_MHCFLURRY_PRESENTATION_PERCENTILE_MAX` constants in `config.py` (both `2.0`).

## Algorithm — five-phase trace

Audit files are written for each phase under `data/intermediate/{track_id}/consensus/` so the run is fully reproducible.

| Phase | Action | Output |
|---|---|---|
| `0a` | Drop rows with NaN in the percentile column. | `0a_PRED_NET_no_nan.csv`, `0a_PRED_FLURRY_no_nan.csv` |
| `0b` | Apply the percentile threshold (default `≤ 2.0%`). | `0b_PRED_NET_thresholded.csv`, `0b_PRED_FLURRY_thresholded.csv` |
| `1` | Consolidate to one row per peptide: best allele, `{tool}_alleles` (semicolon-joined), `{tool}_num_alleles`. | `1_PRED_NET_consolidated.csv`, `1_PRED_FLURRY_consolidated.csv` |
| `2` | Inner-join on `peptide`: keep only peptides present in *both* tools. | `2_intersection.csv` |
| `3` | Merge per-tool columns into a single table. | `CONSENSUS_{track_id}.csv` |
| Calis | Score with the Calis 2013 immunogenicity tool; keep `calis_score > 0`. | `CONSENSUS_IMMUNOGENIC_{track_id}.csv` |

Final summary metadata is written to `consensus_audit_summary.json` (per-phase counts, threshold used, Calis allele picked per peptide).

## Outputs

| File | Role |
|---|---|
| `CONSENSUS_IMMUNOGENIC_{track_id}.csv` | Final output: cumulative table feeding `screen_toxicity` / `cluster_epitopes`. |
| `CONSENSUS_VIEW_{track_id}.csv` | Slim per-step view — `peptide`, best allele/percentile per tool, and `calis_score`. For quick inspection only. |
| `CONSENSUS_{track_id}.csv` | Same peptide set as `CONSENSUS_IMMUNOGENIC` **before** the Calis filter (Stage 1 only). |
| `consensus_audit_summary.json` | Funnel audit — counts at each stage. |

## Output schema (`CONSENSUS_IMMUNOGENIC_{track_id}.csv`)

Column names follow the project-wide convention (prefix-only, English) defined in `utils/naming.py`:

| Column | Meaning |
|---|---|
| `peptide` | Peptide sequence (primary key). |
| `netmhcpan_best_allele` | Best NetMHCpan allele for this peptide (lowest %Rank). |
| `netmhcpan_alleles` | Semicolon-joined NetMHCpan alleles for this peptide. |
| `netmhcpan_num_alleles` | Count of NetMHCpan alleles. |
| `netmhcpan_el_percentile` | Best NetMHCpan EL %Rank. |
| `netmhcpan_el_percentiles_all` | All per-allele NetMHCpan percentiles (semicolon-joined, same order as `netmhcpan_alleles`). |
| `mhcflurry_best_allele` | Best MHCFlurry allele for this peptide. |
| `mhcflurry_alleles` | Semicolon-joined MHCFlurry alleles. |
| `mhcflurry_num_alleles` | Count of MHCFlurry alleles. |
| `mhcflurry_presentation_percentile` | Best MHCFlurry presentation percentile. |
| `mhcflurry_presentation_percentiles_all` | All per-allele MHCFlurry percentiles (semicolon-joined). |
| `calis_score` | Calis 2013 immunogenicity score (> 0 for survivors). |
| `calis_allele_used` | HLA allele the Calis tool was queried with. |

## Calis 2013 — provenance

The file `modules/consensus_filter/predict_immunogenicity.py` is **vendored verbatim** from the IEDB Class I Immunogenicity tool. Its algorithm is the one described in:

> Calis JJA, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Keşmir C, Peters B. *Properties of MHC class I presented peptides that enhance immunogenicity.* PLOS Computational Biology. 2013;9(10):e1003266.

The score combines per-position amino-acid contributions weighted by an anchor mask. **Do not modify the body of `predict_immunogenicity.py`** — it is the reference implementation and any deviation breaks comparability with IEDB-published scores. Bug fixes (e.g. Python compatibility) are allowed; algorithm changes are not.

Citation for this step's choice to use Calis 2013 as a sequential filter after presentation:

> Bui H-H, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A. *Predicting population coverage of T-cell epitope-based diagnostics and vaccines.* BMC Bioinformatics. 2006;7:153. (For the broader IEDB workflow rationale.)

## Thresholds

Defined in `config.py`:

| Constant | Value | Used for |
|---|---|---|
| `CONSENSUS_NETMHCPAN_EL_RANK_MAX_PERCENT` | `2.0` | NetMHCpan threshold in phase `0b`. |
| `CONSENSUS_MHCFLURRY_PRESENTATION_PERCENTILE_MAX` | `2.0` | MHCFlurry threshold in phase `0b`. |
| `CONSENSUS_CALIS_IMMUNOGENICITY_SCORE_MIN` | `0.0` | Calis cut-off (peptides with score ≤ 0 are dropped). |
| `IEDB_IMMUNOGENICITY_API_URL` | `http://tools-cluster-interface.iedb.org/tools_api/immunogenicity/` | Endpoint hit if the local Calis script falls back to the API. |
