# TheraEPIflow

Pipeline for identifying and selecting MHC class I (CTL/CD8+) epitopes for therapeutic vaccine design. The tool runs prediction, filtering, clustering, toxicity screening, conservation analysis, population coverage, and optional murine validation as a single reproducible workflow driven by an interactive CLI.

## Pipeline at a glance

```
fetch_sequences        Pull reference protein sequences from UniProt
predict_binding        NetMHCpan 4.1 EL (IEDB API) + MHCFlurry 2.0 (local)
consensus_filter       Stage 1: NET ∩ FLURRY presentation, Stage 2: Calis 2013 immunogenicity
screen_toxicity        ToxinPred3 in-memory predictor, default cutoff 0.38
cluster_epitopes       (planned) Pairwise alignment + NetworkX-based clustering
select_representatives (planned) Best peptide per cluster
search_variants        (planned) GenBank variant lookup
analyze_conservation   (planned) Pairwise conservation across variants
population_coverage    (planned) IEDB allele frequency database (human only)
predict_murine         (planned) NetMHCpan + MHCFlurry with H-2 alleles
curate_murine          (planned) Cross-species priority labelling
integrate_data         (planned, global) Merge all tracks into a master table
generate_report        (planned, global) Self-contained HTML report
```

The first four steps are implemented and validated end-to-end. The remaining ones are scheduled for May 2026.

## Why this design

A peptide only earns a spot in the candidate list if three independent pieces of evidence agree.

| Biological event | Tool | Field | Cutoff |
|---|---|---|---|
| MHC presentation (EL) | NetMHCpan 4.1 via IEDB classic API | `Rnk_EL` | <= 2 percent |
| MHC presentation (BA + AP) | MHCFlurry 2.0 local | `presentation_percentile` | <= 2 percent |
| TCR recognition | Calis et al. 2013 via IEDB API | `score` | > 0 |

Immunogenicity (Calis) is only computed for peptides that already passed the binding intersection. Scoring TCR contact on something that failed presentation would just waste API calls.

Toxicity is screened right after the consensus filter, before clustering. Removing toxic peptides first keeps clusters meaningful and avoids work on candidates that would be discarded anyway.

## Installation

You need Linux or WSL2, Conda (Miniconda or Anaconda), and Git.

```bash
git clone <repository-url>
cd TheraEPIflow
bash setup.sh
conda activate theraEPIflow
```

The setup script creates the `theraEPIflow` environment from `environment.yml` with Python 3.10, Biopython, MHCFlurry 2.0, ToxinPred3, scikit-learn 1.2.2 (pinned for ToxinPred3 model compatibility), Rich, Pandas, Jinja2, and the rest of the stack.

## Quick start

```bash
conda activate theraEPIflow

# List existing projects and open the menu
python main.py

# Create a new project (asks 4 questions: name, description, email, host)
python main.py --new-project

# Open an interactive REPL session for an existing project
python main.py --project hpv_study

# Run a specific step in non-interactive mode and exit
python main.py --project hpv_study --step screen_toxicity

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
    project_config.json     Shared settings, grows as each step adds its own config
    pipeline.json           Per-track step states (done, error, pending)
    data/
      input/{track_id}/                       Reference FASTA + sequence registry
      intermediate/{track_id}/
        predictions/  PRED_NET_*.csv, PRED_FLURRY_*.csv
        consensus/    CONSENSUS_FILTERED_*.csv, CONSENSUS_IMMUNOGENIC_*.csv, audit JSON
        toxicity/     TOXICITY_ALL_*.csv, TOXICITY_SAFE_*.csv, audit JSON
        clusters/     (created by cluster_epitopes)
        variants/     (created by search_variants)
        conservation/ (created by analyze_conservation)
        coverage/     (created by population_coverage)
        murine/       (created by predict_murine and curate_murine)
      output/
        master_table.xlsx   (created by integrate_data)
        report.html         (created by generate_report)
```

A track is one organism plus one protein. All tracks in a project share the same HLA alleles and pipeline parameters. Track IDs follow `{ORGANISM_LABEL}_{PROTEIN_LABEL}`, for example `HPV16_E6` or `ZIKV_E`. Labels are suggested automatically based on standard abbreviations and can be overridden when the project is created.

## Configuration as you go

The CLI never asks for parameters you do not need yet. Each step collects its own settings the first time it runs and saves them to `project_config.json`, so reruns stay deterministic.

| Stage | Asked at this point |
|---|---|
| `--new-project` | Project name, description, Entrez email, target host |
| `fetch_sequences` first run | Organism aliases, proteins, labels, input source |
| `predict_binding` first run | HLA alleles (defaults to 27 globally diverse alleles), peptide lengths |
| `consensus_filter` first run | Stage 1 percentile threshold |
| `screen_toxicity` first run | ToxinPred3 cutoff (default 0.38) |
| `cluster_epitopes` first run | Identity cutoff (planned) |
| `population_coverage` first run | Populations to evaluate (planned, human only) |
| `predict_murine` first run | Murine strains, optional (planned) |

## Implemented steps

### fetch_sequences

Pulls the reference protein sequence for each track from the UniProt REST API. Maps short aliases (`HPV16`, `MPOX`, `CHIKV`) to scientific names and tax IDs, runs three search strategies (protein name, gene name, organism only), prefers Swiss-Prot entries, and flags polyproteins or fragments based on the median length of the result set. Saves a cleaned FASTA, a JSON registry with metadata, and a validation report. See `modules/fetch_sequences/README.md`.

### predict_binding

Runs NetMHCpan 4.1 EL via IEDB's classic API and MHCFlurry 2.0 locally. Both predictors run in parallel via `ThreadPoolExecutor`. The IEDB call iterates over the cartesian product of alleles and peptide lengths. MHCFlurry uses `Class1PresentationPredictor.predict_to_dataframe`, called once per allele to keep the per-(peptide, allele) granularity that consensus_filter expects. Output: `PRED_NET_{track_id}.csv`, `PRED_FLURRY_{track_id}.csv`, and an audit JSON.

### consensus_filter

Two stages.

Stage 1 keeps only peptides predicted as presented by both NetMHCpan and MHCFlurry. Per-peptide aggregation collapses N rows (one per allele) into one row per peptide, recording the best allele and the joined HLA list under suffixes `_net_pred` and `_flurry_pred`.

Stage 2 sends the surviving peptides to the IEDB Calis 2013 immunogenicity endpoint and keeps those with `score > 0`.

Audit files (`0a_`, `0b_`, `1_`, `2_`, `3_` plus a JSON summary) preserve every intermediate state for reproducibility.

### screen_toxicity

Loads the ToxinPred3 Model 1 (AAC + DPC features, Extra Trees) once with joblib and predicts in memory, no temporary files or subprocess. Computes the ML score and the calibrated PPV with the exact linear coefficients used upstream (`PPV = score * 1.2341 - 0.1182`, clipped to [0, 1]). Default cutoff is 0.38, which is the value the ToxinPred3 paper uses for short peptides. Drops NaN rows and any peptide containing non-canonical residues (`X`, `B`, `Z`) before predicting, reporting the count when that happens. Outputs:

| File | Contents |
|---|---|
| `TOXICITY_ALL_{track_id}.csv` | Every input row plus `toxinpred3_score`, `toxinpred3_ppv`, `toxinpred3_label` |
| `TOXICITY_SAFE_{track_id}.csv` | Only the `Non-Toxin` rows, used as input for the next step |
| `TOXICITY_AUDIT_{track_id}.json` | Run metadata, counts, and chosen threshold |

## Tools and licenses

| Tool | License | Used for |
|---|---|---|
| MHCFlurry 2.0 | Apache 2.0 | Local MHC-I presentation prediction |
| NetMHCpan 4.1 | Academic, accessed through IEDB | EL rank prediction |
| IEDB classic API | Free for academic use | NetMHCpan and Calis immunogenicity |
| ToxinPred3 | Academic | Toxicity screening |
| Biopython | BSD | FASTA handling and GenBank lookups |
| Rich | MIT | Terminal UI |

## Development status

| Step | Status |
|---|---|
| fetch_sequences | implemented, validated on HPV, MPOX, CHIKV, DENV |
| predict_binding | implemented, validated end-to-end |
| consensus_filter | implemented, two stages plus Calis 2013 |
| screen_toxicity | implemented, ToxinPred3 in memory |
| cluster_epitopes onwards | planned, target delivery May 2026 |

The pipeline is under active development. See `CLAUDE.md` for ongoing project notes (internal, not committed).
