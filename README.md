# TheraEPIflow

Pipeline for identifying and selecting MHC class I (CTL/CD8+) epitopes for therapeutic vaccine design. The tool runs prediction, filtering, clustering, toxicity screening, variant search, conservation analysis, population coverage, and optional murine validation as a single reproducible workflow driven by an interactive CLI.

## Pipeline at a glance

```
fetch_sequences        Pull reference protein sequences from UniProt
predict_binding        NetMHCpan 4.1 EL (IEDB API) + MHCFlurry 2.0 (local)
consensus_filter       Stage 1: NET ∩ FLURRY presentation, Stage 2: Calis 2013 immunogenicity
screen_toxicity        ToxinPred3 in-memory predictor, default cutoff 0.38
cluster_epitopes       Pairwise identity + NetworkX clustering (cluster_break, default 80%)
select_representatives Best peptide per cluster (min percentile + allele breadth)
search_variants        UniProt variant lookup — intraspecific or interspecific scope
analyze_conservation   Sliding window identity across variants, visual heatmap XLSX
population_coverage    IEDB allele-frequency pickle, diploid per-epitope coverage
predict_murine         NetMHCpan + MHCFlurry with H-2 (murine) alleles
curate_murine          Per-track master: joins human ★ + conservation + coverage + murine
integrate_data         (planned, global) Merge all tracks into a master table
generate_report        (planned, global) Self-contained HTML report
```

Steps 1–11 (per-track) are implemented and validated end-to-end — most recently a clean two-track run (hantavirus Nucleoprotein + Zika NS1, **11/11 steps**, ~6 min with a 4-allele panel). The two global steps — `integrate_data` and `generate_report` — remain.

## Why this design

A peptide only earns a spot in the candidate list if three independent pieces of evidence agree.

| Biological event | Tool | Field | Cutoff |
|---|---|---|---|
| MHC presentation (EL) | NetMHCpan 4.1 via IEDB classic API | `Rnk_EL` | ≤ 2% |
| MHC presentation (BA + AP) | MHCFlurry 2.0 local | `presentation_percentile` | ≤ 2% |
| TCR recognition | Calis et al. 2013 via IEDB API | `score` | > 0 |

Immunogenicity (Calis) is only computed for peptides that already passed the binding intersection. Scoring TCR contact on something that failed presentation would just waste API calls.

Toxicity is screened right after the consensus filter, before clustering. Removing toxic peptides first keeps clusters meaningful and avoids work on candidates that would be discarded anyway.

### Per-step module layout

Every step is a Python package under `modules/<step>/`. Larger steps are split by
**responsibility** (the Single Responsibility Principle / separation of concerns),
so each file answers one question — *where is the maths? where is the file I/O?
where is the prompt?* The file names form a consistent vocabulary across steps:

| File | Responsibility |
|---|---|
| `step.py` | Orchestration — the `BaseTrackStep`/`BaseGlobalStep` subclass (`run`, `preflight`, `postflight`, `describe_outputs`) |
| `core.py` | Domain logic — pure computation (scoring, classification, the maths) and input loaders |
| `io.py` | Data writers — CSV / XLSX / JSON |
| `charts.py` | Matplotlib PNG generation |
| `prompts.py` | Interactive terminal prompts |
| `render.py` | Rich console output (tables, panels, summaries) |
| `__init__.py` | Facade — re-exports the step class (and any symbol other modules import) |

The guiding rule: **the domain logic in `core.py` never imports Rich, openpyxl
or `input()`** — only `step.py` knows about all the layers and wires them
together. A file is only created when the step actually has that responsibility
(no empty `charts.py` for a step that emits no PNGs), and small steps that would
fragment into trivial files stay as a single `__init__.py`.

## Installation

You need Linux or WSL2, Conda (Miniconda or Anaconda), and Git.

```bash
git clone <repository-url>
cd TheraEPIflow
bash setup.sh
conda activate TheraEPIflow
```

The setup script creates the `TheraEPIflow` environment from `environment.yml` with Python 3.10, Biopython, MHCFlurry 2.0, ToxinPred3, scikit-learn 1.2.2, NetworkX, Rich, Pandas, matplotlib + seaborn (for heatmap/hit-chart PNGs), Jinja2, and the rest of the stack. MHCFlurry presentation models are downloaded automatically.

### Known constraints (do not bump these pins lightly)

Three packages are pinned to exact versions because newer releases break the rest of the stack:

| Pin | Reason |
|---|---|
| `tensorflow==2.15.0` | MHCFlurry 2.0 internals depend on the TF 2.15 API; TF 2.16+ removes / renames functions and MHCFlurry stops loading its presentation model. |
| `scikit-learn==1.2.2` | ToxinPred3's shipped model `.pkl` was trained on sklearn 1.2.2; loading it under 1.3+ raises an `AttributeError` because of internal class changes. |
| `mhcflurry==2.0.6` | MHCFlurry 2.2 changed the `Class1PresentationPredictor.predict_to_dataframe` signature; the predict_binding step's batching logic assumes the 2.0.x API. |

If you upgrade any of these, the affected step will fail at import or at first call. A future dependency-management pass will revisit this triplet once the upstream tools converge again.

## Quick start

```bash
conda activate TheraEPIflow

# List existing projects and open the menu
python main.py

# Create a new project (asks name and description)
python main.py --new-project

# Open an interactive REPL session for an existing project
python main.py --project hpv_study

# Run a specific step in non-interactive mode and exit
python main.py --project hpv_study --step cluster_epitopes

# Show progress without launching the REPL
python main.py --project hpv_study --status

# List every project on disk
python main.py --list
```

Inside the REPL: `Enter` runs the next pending step, `a` runs all remaining steps, `r` reruns the last step, `j NAME` jumps to a specific step, `s` shows the status table, `q` quits.

## How a project is organized

```
projects/
  {project_name}/
    project_config.json     Shared settings — grows as each step adds its own config
    pipeline.json           Per-track step states (done, error, pending)
    data/
      input/{track_id}/                       Reference FASTA + sequence registry
      intermediate/{track_id}/
        predictions/  PRED_NET_*.csv, PRED_FLURRY_*.csv
        consensus/    CONSENSUS_FILTERED_*.csv, CONSENSUS_IMMUNOGENIC_*.csv, audit JSON
        toxicity/     TOXICITY_ALL_*.csv, TOXICITY_SAFE_*.csv, audit JSON
        clusters/     CLUSTER_*.csv, CLUSTER_REPR_*.csv (with ★ column), audit JSON
        variants/     VARIANTS_*.fasta (permanent cache), audit JSON
        conservation/ CONSERVATION_*.csv, CONSERVATION_*.xlsx, CONSERVATION_VISUAL_*.xlsx
        coverage/     (created by population_coverage)
        murine/       (created by predict_murine and curate_murine)
      output/
        master_table.xlsx   (created by integrate_data)
        report.html         (created by generate_report)
```

A track is one organism plus one protein. All tracks in a project share the same HLA alleles and pipeline parameters. Track IDs follow `{ORGANISM_LABEL}_{PROTEIN_LABEL}`, for example `HPV16_E6` or `SARS2_S`. Labels are suggested automatically based on standard abbreviations and can be overridden when the project is created.

## Configuration as you go

The CLI never asks for parameters you do not need yet. Each step collects its own settings the first time it runs and saves them to `project_config.json`, so reruns stay deterministic.

| Stage | Asked at first run |
|---|---|
| `--new-project` | Project name, description (target host is fixed at Homo sapiens — HLA-I is human) |
| `fetch_sequences` | Organism aliases, proteins, labels, input source |
| `predict_binding` | HLA alleles (default: 27 globally diverse), peptide lengths |
| `consensus_filter` | Stage 1 percentile threshold |
| `screen_toxicity` | ToxinPred3 cutoff (default 0.38) |
| `cluster_epitopes` | Identity cutoff (default 80%), clustering method |
| `search_variants` | Scope (intraspecific / interspecific) and optional host filter |
| `analyze_conservation` | Analysis threshold (default 1.0 = exact match) |
| `population_coverage` | Populations to evaluate (human only) |
| `predict_murine` | Murine strains (optional) |

## Implemented steps

### fetch_sequences

Pulls the reference protein sequence for each track from the UniProt REST API. Maps short aliases (`HPV16`, `MPOX`, `CHIKV`, `SARS2`) to scientific names and taxonomy IDs, runs three search strategies (protein name, gene name, organism only), prefers Swiss-Prot entries, and flags polyproteins or fragments based on the median length of the result set. Saves a cleaned FASTA, a JSON registry with metadata, and a validation report. See `modules/fetch_sequences/README.md`.

### predict_binding

Runs NetMHCpan 4.1 EL via IEDB's classic API and MHCFlurry 2.0 locally. Both predictors run in parallel via `ThreadPoolExecutor`. The IEDB call iterates over the cartesian product of alleles and peptide lengths. MHCFlurry uses `Class1PresentationPredictor.predict_to_dataframe`, called once per allele to keep the per-(peptide, allele) granularity that consensus_filter expects. Output: `PRED_NET_{track_id}.csv`, `PRED_FLURRY_{track_id}.csv`, and an audit JSON.

### consensus_filter

Two stages.

Stage 1 keeps only peptides predicted as presented by both NetMHCpan and MHCFlurry. Per-peptide aggregation collapses N rows (one per allele) into one row per peptide, recording the best allele and the joined HLA list.

Stage 2 sends the surviving peptides to the IEDB Calis 2013 immunogenicity endpoint and keeps those with `score > 0`.

Audit files preserve every intermediate state for reproducibility.

### screen_toxicity

Loads the ToxinPred3 Model 1 (AAC + DPC features, Extra Trees) once with joblib and predicts in memory, no temporary files or subprocess. Computes the ML score and the calibrated PPV. Default cutoff is 0.38. Drops rows with non-canonical residues (`X`, `B`, `Z`) before predicting. Outputs a full table and a safe-only table for the next step. See `modules/screen_toxicity/README.md`.

### cluster_epitopes

Groups safe epitopes by pairwise sequence identity using Biopython global alignment and NetworkX. Three methods available; `cluster_break` is the default.

`cluster_break` starts from connected components and iteratively removes the minimum-weight edge from any component whose average intra-cluster similarity falls below the threshold. This avoids "bridge" groupings where peptides with distinct TCR recognition land in the same cluster — relevant because pMHC-TCR recognition is highly sensitive to anchor position differences.

Default identity threshold: **80%**. See `modules/cluster_epitopes/README.md`.

### select_representatives

Selects one representative peptide per cluster using a combined score: minimum percentile across all allele-level values from both tools (`best_combined_percentile`) plus allele breadth (`num_alleles_united`). Marks the winner with `★` in the `best_representative` column. Produces a colour-coded XLSX deliverable alongside the CSV. See `modules/select_representatives/README.md`.

### search_variants

Queries UniProt REST API for protein variants. Two scopes:

- **Intraspecific** — variants within the same taxonomy ID (e.g. different HPV16 isolates, SARS-CoV-2 strains)
- **Interspecific** — same protein name across related species, with optional host filter (e.g. E5 from HPV16, HPV18, HPV31 infecting Homo sapiens)

The multi-FASTA is written once and cached permanently. On re-run the user can choose to keep or regenerate the file. Near-identical sequences (≥ 99% vs. reference) are filtered out before writing. See `modules/search_variants/README.md`.

### analyze_conservation

Measures how faithfully each ★ representative appears across the variant sequences from `search_variants`. Uses a sliding window: for each variant the best-matching window of the same length as the peptide is found and its per-position identity is recorded.

The analysis threshold is configurable (default 1.0 = exact match). It controls which variants are labelled "passed" vs "failed" and what mutations are surfaced in `position_mutation_profile`. The `conservation_label` (perfect / high / moderate / low / conservation_unknown) and all row colours are based on `mean_max_identity` and never change regardless of the chosen threshold.

Also produces a `CONSERVATION_VISUAL_{track_id}.xlsx` heatmap with one row per ★ epitope, one column per position, cells colour-coded from green (100% conserved) to red (< 50%). See `modules/analyze_conservation/README.md`.

### population_coverage

For each ★ epitope, computes the fraction of one or more human populations that carries at least one of the epitope's HLA alleles, using a vendored copy of the IEDB Population Coverage tool's allele-frequency pickle (IEDB v3.0.2, pickle v1.1.2, sourced from AlleleFrequencies.net). Diploid model: `p_locus = 1 − (1 − q)²` per locus, combined as `1 − Π(1 − p_locus)`. Math is bit-for-bit reproducible against `iedb.org/population_coverage`.

Supports multiple populations in a single run; per-population output includes an IEDB-style Hit Chart PNG plus a long-format CSV. With two or more populations, a comparative heatmap PNG is also generated. See `modules/population_coverage/README.md` and `modules/population_coverage/data/SOURCE.md` for the pickle's provenance.

### predict_murine

Re-runs NetMHCpan EL + MHCFlurry on the ★ representatives against murine **H-2** alleles instead of human HLA-I, to flag which human-selected epitopes would also bind in a mouse model. Strain groups (`C57BL/6`, `BALB/c`, `complete`) live in `config.MURINE_ALLELES`. Each ★ epitope gets per-allele records plus an aggregated row with best percentile, H-2 alleles bound, and a four-tier binder label (`optimal` / `good` / `borderline` / `non_binder`). Qualitative — never removes epitopes. See `modules/predict_murine/README.md`.

### curate_murine

Assembles the per-track master table by joining each ★ representative's human qualification with its conservation, population coverage and (when present) murine prediction. A JOIN-only step — no ranking or priority relabelling; conservation and coverage are required inputs, murine is optional. Output: `CURATE_MURINE_{track_id}.csv` (one row per ★ epitope with all annotations) + a slim view + audit JSON.

## Master-table strategy

`curate_murine` already produces a per-track joined table. The upcoming `integrate_data` global step stacks those across all tracks into a single `output/master_table.xlsx`. The pipeline emits its files in the right shape: every implemented step except `search_variants` (FASTA, by design — consumed by `analyze_conservation`) and `population_coverage` (long format, by design — pivoted at merge time) produces a one-row-per-peptide CSV with `peptide` as the primary key. Per-track assembly order:

```
CONSENSUS_IMMUNOGENIC_{track}.csv     (1 row / peptide — base)
  + TOXICITY_SAFE_{track}.csv         (appends toxicity columns)
  + CLUSTER_{track}.csv               (appends cluster_id, etc.)
  + CLUSTER_REPR_{track}.csv          (filters to ★ rows, appends norm scores)
  + CONSERVATION_{track}.csv          (appends conservation_label etc.)
  + COVERAGE_{track}.csv (pivoted)    (long → wide: coverage_World, coverage_Brazil, …)
```

Track context (`track_id`, `organism_label`, `protein_label`) comes from `project_config.json`. The master table is what feeds the future HTML report's interactive selection calculator.

## Tools and licenses

| Tool | License | Used for |
|---|---|---|
| MHCFlurry 2.0 | Apache 2.0 | Local MHC-I presentation prediction |
| NetMHCpan 4.1 | Academic, accessed through IEDB | EL rank prediction |
| IEDB classic API | Free for academic use | NetMHCpan and Calis immunogenicity |
| ToxinPred3 | Academic | Toxicity screening |
| Biopython | BSD | FASTA handling, pairwise alignment |
| NetworkX | BSD | Clustering graph algorithms |
| Rich | MIT | Terminal UI |
| openpyxl | MIT | XLSX generation |

## Development status

| Step | Status |
|---|---|
| fetch_sequences | implemented, validated on HPV, MPOX, CHIKV, DENV, SARS-CoV-2 |
| predict_binding | implemented, validated end-to-end |
| consensus_filter | implemented, two stages plus Calis 2013 |
| screen_toxicity | implemented, ToxinPred3 in memory |
| cluster_epitopes | implemented, default 80% threshold, cluster_break method |
| select_representatives | implemented, min-percentile + allele-breadth scoring |
| search_variants | implemented, intraspecific + interspecific with host filter |
| analyze_conservation | implemented, sliding window, heatmap XLSX |
| population_coverage | implemented, IEDB pickle (v1.1.2) vendored, multi-population |
| predict_murine | implemented, NetMHCpan + MHCFlurry with H-2 alleles |
| curate_murine | implemented, per-track join of human ★ + conservation + coverage + murine |
| integrate_data | planned (global step) |
| generate_report | planned (global step) |

All step modules were refactored into the per-step layout described above (Single Responsibility split), with no behaviour change verified module-by-module, and the pipeline was re-validated end-to-end on a two-track run (hantavirus Nucleoprotein + Zika NS1) — all 11 implemented steps completed cleanly through the interactive CLI.

Earlier validation projects: `hpv_study` (HPV16/18/31 × E5/E6/E7), `chikungunya_study` (CHIKV M, NSP1), `dengue_study` (DENV E, NS5), `monkeypox_study` (MPOX M), `oropouche2` (OROV N), `sars2_test` (SARS2 S, N).
