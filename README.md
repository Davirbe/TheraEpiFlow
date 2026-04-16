# TheraEPIflow

**MHC-I Epitope Pipeline for Therapeutic Vaccine Design**

Automated pipeline for identifying and selecting MHC class I (CTL/CD8+) epitopes. Integrates binding prediction, quality filtering, variant analysis, conservation scoring, population coverage, and murine validation into a single reproducible workflow.

---

## Pipeline Overview

```
Step 01 — Fetch Sequences        Download protein sequences from NCBI GenBank
Step 02 — Validate Input         Clean and validate FASTA sequences
Step 03 — Predict Binding        NetMHCpan 4.1 (EL) + MHCFlurry 2.2 predictions
Step 04 — Consensus Filter       Intersection of both predictors + immunogenicity
Step 05 — Cluster Epitopes       Group by sequence identity (local, core-based)
Step 06 — Select Representatives Pick best per cluster (HLA count + mean rank)
Step 07 — Screen Toxicity        ToxinPred2 toxicity screening
Step 08 — Search Variants        GenBank variant/isolate sequences
Step 09 — Analyze Conservation   Pairwise conservation across variants
Step 10 — Population Coverage    IEDB allele frequency database (human only)
Step 11 — Predict Murine         NetMHCpan + MHCFlurry with H-2 alleles (optional)
Step 12 — Curate Murine          Filter promiscuous H-2 binders
Step 13 — Integrate Data         Merge all tracks into master table
Step 14 — Generate Report        Interactive HTML report with construct builder
```

### Consensus Filter Design

Epitopes must pass three independent criteria simultaneously:

| Biological event | Tool | Parameter | Cutoff |
|---|---|---|---|
| MHC binding + processing (EL) | NetMHCpan 4.1 via IEDB API | `Rnk_EL` | ≤ 2% |
| MHC binding + processing (BA+AP) | MHCFlurry 2.2 (local) | `presentation_percentile` | ≤ 2% |
| TCR recognition | Calis et al. 2013 via IEDB API | `score` | > 0 |

Immunogenicity is scored only on epitopes that pass the binding intersection — not in parallel.

---

## Installation

### Requirements

- Linux or WSL2 (Windows Subsystem for Linux)
- Conda (Miniconda or Anaconda)
- Git

### Setup

```bash
git clone <repository-url>
cd TheraEPIflow

# Install Miniconda (if not already installed) and create the environment
bash setup.sh

# Activate the environment
conda activate theraEPIflow
```

The setup script installs all Python dependencies including Biopython, MHCFlurry 2.2, Rich, Pandas, and Jinja2.

---

## Quick Start

```bash
conda activate theraEPIflow
cd TheraEPIflow

# Show all projects and quick-start menu
python main.py

# Create a new project
python main.py --new-project

# Run all pending steps for a project
python main.py --project PROJECT_NAME --run

# Run a specific step
python main.py --project PROJECT_NAME --step 1

# Check step-by-step progress
python main.py --project PROJECT_NAME --status

# List all projects
python main.py --list
```

---

## Project Structure

```
projects/
  {project_name}/
    project_config.json       Shared settings (grows as each step adds config)
    pipeline.json             Per-track step states (done/pending/error)
    data/
      input/{track_id}/       Reference FASTA + sequence registry
      intermediate/{track_id}/
        predictions/          NetMHCpan + MHCFlurry output CSVs
        consensus/            Filtered epitope tables
        clusters/             Cluster assignments + representative selection
        toxicity/             ToxinPred2 results
        variants/             GenBank variant sequences (FASTA)
        conservation/         Per-epitope conservation scores
        coverage/             Population coverage tables
        murine/               H-2 prediction results (if enabled)
      output/
        master_table.xlsx     All epitopes with all scores
        report.html           Interactive HTML report (self-contained)
```

A **track** is one organism + one protein pair (e.g. `HPV16_E6`, `ZIKV_E`).
All tracks in a project share the same HLA alleles and pipeline parameters.

---

## Track ID Convention

Track IDs follow the pattern `{ORGANISM_LABEL}_{PROTEIN_LABEL}`:

| Organism | Protein | Track ID |
|---|---|---|
| Human papillomavirus 16 | E6 | `HPV16_E6` |
| Human papillomavirus 18 | E7 | `HPV18_E7` |
| Zika virus | envelope protein E | `ZIKV_E` |
| Dengue virus 2 | NS5 | `DENV2_NS5` |
| SARS-CoV-2 | spike protein | `SARS2_S` |

Labels are suggested automatically and can be overridden at project creation.

---

## Configuration Philosophy

Parameters are asked contextually — only when they are needed:

| Stage | Parameters collected |
|---|---|
| `--new-project` | Project name, description, Entrez email, target host |
| Step 01 | Organisms, proteins, input source (GenBank or local file) |
| Step 03 | HLA alleles, peptide lengths |
| Step 05 | Clustering identity cutoff |
| Step 07 | Toxicity threshold |
| Step 10 | Populations for coverage analysis |
| Step 11 | Murine validation strains (optional) |

---

## Tools and Licenses

| Tool | License | Usage |
|---|---|---|
| MHCFlurry 2.2 | Apache 2.0 | MHC-I presentation prediction (local) |
| NetMHCpan 4.1 | Academic (via IEDB API) | EL rank prediction |
| IEDB API | Free for academic use | NetMHCpan, Calis immunogenicity |
| ToxinPred2 | Academic | Toxicity screening |
| Biopython | BSD | GenBank access, sequence handling |
| Rich | MIT | Terminal UI |

---

## Development Status

| Step | Status |
|---|---|
| Step 01 — Fetch Sequences | ✅ Complete |
| Step 02 — Validate Input | 🔄 In progress |
| Steps 03-14 | ⏳ Planned |

*Pipeline under active development — April/May 2026*
