# TheraEpiFlow

[![Tutorial](https://img.shields.io/badge/%F0%9F%93%96_Tutorial-online-2563eb)](https://davirbe.github.io/TheraEpiFlow/)

> **New here? Start with the [interactive tutorial](https://davirbe.github.io/TheraEpiFlow/)**, a self-contained walkthrough served from the `docs/` folder via GitHub Pages.

Pipeline for identifying and selecting MHC class I (CTL/CD8+) epitopes for therapeutic vaccine design. The tool runs prediction, filtering, clustering, toxicity screening, variant search, conservation analysis, population coverage, and optional murine validation as a single reproducible workflow driven by an interactive CLI.

## Pipeline at a glance

```
fetch_sequences        Pull reference protein sequences from UniProt
predict_binding        NetMHCpan 4.1 EL (IEDB API) + MHCFlurry 2.0 (local)
consensus_filter       Stage 1: NET ∩ FLURRY presentation, Stage 2: Calis 2013 immunogenicity
screen_toxicity        ToxinPred3 in-memory predictor, default cutoff 0.38
cluster_epitopes       Pairwise identity + NetworkX clustering (cluster_break, default 80%)
select_representatives Best peptide per cluster (min percentile + allele breadth)
search_variants        UniProt variant lookup (intraspecific or interspecific scope)
analyze_conservation   Sliding window identity across variants, visual heatmap XLSX
population_coverage    IEDB allele-frequency pickle, diploid per-epitope coverage
predict_murine         NetMHCpan + MHCFlurry with H-2 (murine) alleles
curate_murine          Per-track master: joins human ★ + conservation + coverage + murine
integrate_data         (global) Stack every track into MASTER_TABLE_FULL/VIEW + audit JSON
generate_report        (global) Self-contained interactive HTML calculator (REPORT_*.html)
```

Plus a non-step utility, reachable from the project menu (`[z]`):

```
download menu          .zip (WSL) or .tar.gz (Linux), to in-project / ~/Downloads / WSL Windows
```

All **13 pipeline steps** (11 per-track + 2 global) are implemented and validated end-to-end. The most recent multi-track validation: `hpv16` (5 tracks, 86 ★ peptides) and `scer_test` (2 tracks, 83 ★ peptides), both producing the master tables and the HTML calculator without intervention. Once the pipeline completes, the `[z]` menu key in the REPL packages the whole project into an archive. See [`utils/download_ui.py`](utils/download_ui.py).

## Why this design

A peptide only earns a spot in the candidate list if three independent pieces of evidence agree.

| Biological event | Tool | Field | Cutoff |
|---|---|---|---|
| MHC presentation (EL) | NetMHCpan 4.1 via IEDB classic API | `Rnk_EL` | ≤ 2% |
| MHC presentation (BA + AP) | MHCFlurry 2.0 local | `presentation_percentile` | ≤ 2% |
| TCR recognition | Calis et al. 2013 via IEDB API | `score` | > 0 |

Immunogenicity (Calis) is only computed for peptides that already passed the binding intersection. Scoring TCR contact on something that failed presentation would just waste API calls.

Toxicity is screened right after the consensus filter, before clustering. Removing toxic peptides first keeps clusters meaningful and avoids work on candidates that would be discarded anyway.

Past the binder filter, the rest of the pipeline turns a raw hit list into a shortlist a researcher can defend:

- **Cluster, then keep one representative.** Binders that differ by a single anchor residue usually share a TCR footprint, so reporting all of them just pads the list with near-duplicates. `cluster_epitopes` groups by sequence identity and `select_representatives` marks one `★` per cluster (best percentile plus the widest allele coverage). Every analysis downstream runs on that non-redundant `★` set.
- **Check the variants the candidate must survive.** An epitope is only useful if it is conserved in the strains being targeted. `search_variants` pulls the relevant set (isolates of the same species, or the same protein across related species) and `analyze_conservation` scores how well each `★` holds up. These steps annotate and never silently drop a peptide: conservation comes back as a label and a heatmap, so the judgement call stays with the researcher.
- **Coverage and model translatability.** `population_coverage` computes, from the IEDB allele-frequency data, what fraction of a population carries an HLA that presents the epitope; `predict_murine` flags which human-selected epitopes would also bind mouse H-2, so a candidate can be tested in a mouse model. Both are qualitative add-ons, not filters.
- **One table, one report.** `curate_murine` joins each `★` peptide's evidence per track, `integrate_data` stacks every track into a single master table, and `generate_report` turns it into a self-contained offline HTML calculator: the researcher filters, sorts, selects epitopes, and assembles a vaccine construct (linkers, adjuvant, His-tag) with population coverage recomputed live. That HTML plus the master tables are the deliverable.

The order follows one rule throughout: the cheapest and most decisive filters run first (binding intersection, then immunogenicity, then toxicity), and the expensive or purely informational analyses run only on the small `★` set that survives. The whole workflow is validated end to end (see [Validation](#validation)).

### Per-step module layout

Every step is a Python package under `modules/<step>/`. Larger steps are split by
**responsibility** (the Single Responsibility Principle / separation of concerns),
so each file answers one clear question: where is the maths, where is the file
I/O, where is the prompt? The file names form a consistent vocabulary across steps:

| File | Responsibility |
|---|---|
| `step.py` | Orchestration. The `BaseTrackStep`/`BaseGlobalStep` subclass (`run`, `preflight`, `postflight`, `describe_outputs`) |
| `core.py` | Domain logic: pure computation (scoring, classification, the maths) and input loaders |
| `io.py` | Data writers: CSV, XLSX, JSON |
| `charts.py` | Matplotlib PNG generation |
| `prompts.py` | Interactive terminal prompts |
| `render.py` | Rich console output (tables, panels, summaries) |
| `__init__.py` | Facade. Re-exports the step class (and any symbol other modules import) |

The guiding rule: **the domain logic in `core.py` never imports Rich, openpyxl
or `input()`**. Only `step.py` knows about all the layers and wires them
together. A file is only created when the step actually has that responsibility
(no empty `charts.py` for a step that emits no PNGs), and small steps that would
fragment into trivial files stay as a single `__init__.py`.

## Installation

You need Linux or WSL2 and Git. Conda is **not required** beforehand. `setup.sh` installs Miniconda3 into `~/miniconda3` if it doesn't find a pre-existing conda (Miniconda, Anaconda, Miniforge, or system-wide).

```bash
git clone https://github.com/Davirbe/TheraEpiFlow.git
cd TheraEpiFlow
bash setup.sh
conda activate TheraEpiFlow
```

`setup.sh` is self-contained: locates (or installs) conda, creates the `TheraEpiFlow` environment from `environment.yml` with Python 3.10, Biopython, MHCFlurry 2.0, ToxinPred3, scikit-learn 1.2.2, NetworkX, Rich, Pandas, matplotlib + seaborn (for heatmap/hit-chart PNGs), Jinja2, and the rest of the stack, then downloads MHCFlurry presentation models. Re-running `setup.sh` on an existing checkout updates the env in place (`conda env update --prune`).

### Known constraints (do not bump these pins lightly)

Three packages are pinned to exact versions because newer releases break the rest of the stack:

| Pin | Reason |
|---|---|
| `tensorflow==2.15.0` | MHCFlurry 2.0 internals depend on the TF 2.15 API; TF 2.16+ removes / renames functions and MHCFlurry stops loading its presentation model. |
| `scikit-learn==1.2.2` | ToxinPred3's shipped model `.pkl` was trained on sklearn 1.2.2; loading it under 1.3+ raises an `AttributeError` because of internal class changes. |
| `mhcflurry==2.0.6` | MHCFlurry 2.2 changed the `Class1PresentationPredictor.predict_to_dataframe` signature; the predict_binding step's batching logic assumes the 2.0.x API. |

If you upgrade any of these, the affected step will fail at import or at first call. A future dependency-management pass will revisit this triplet once the upstream tools converge again.

### Troubleshooting

> A more comprehensive list of installation pitfalls observed across user
> machines (WSL env activation, `conda` already installed, `~/Downloads` vs
> `/mnt/c/Users/…/Downloads`, RAM allocation under WSL2, etc.) lives in
> [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md). Start there if anything below
> doesn't match your symptom.

**`ImportError: No module named 'pkg_resources'` (or MHCFlurry crashes on import)**

Your conda env has `setuptools >= 81`, which removed `pkg_resources`. MHCFlurry 2.0.6 still imports it. The env file pins `setuptools<81`, but `conda env update --prune` does not force-downgrade an already-installed `setuptools`. Recreate the env:

```bash
conda env remove -n TheraEpiFlow
bash setup.sh
```

Or downgrade in place:

```bash
conda activate TheraEpiFlow && pip install 'setuptools<81'
```

`python main.py` checks for `pkg_resources` on startup and exits with this message before any pipeline step runs.

**NetMHCpan times out / `predict_binding` hangs**

The pipeline talks to IEDB's classic API (`tools-cluster-interface.iedb.org`) over plain HTTP (port 80), the only outbound network dependency. If your network blocks it (university and corporate firewalls are common offenders), `predict_binding` will time out at 120 s. Diagnose with:

```bash
timeout 10 bash -c 'cat < /dev/tcp/tools-cluster-interface.iedb.org/80' && echo OK || echo BLOCKED
```

If `BLOCKED`, try a phone hotspot to confirm it's the network, then ask your network admin to whitelist `tools-cluster-interface.iedb.org`. MHCFlurry, Calis 2013, clustering, conservation, and coverage all run locally; only NetMHCpan needs IEDB.

## Quick start

```bash
conda activate TheraEpiFlow
python main.py
```

That opens the menu, and everything runs from there: one key per action, so you rarely type a full command. To watch the whole pipeline run once without configuring anything, type `demo` — it builds a small HPV16 E7 project and runs every step end to end in a couple of minutes. If a step fails it stops there and names the fix, so `demo` doubles as an install check; for anything unexpected see [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md). If you would rather read first, the [interactive tutorial](https://davirbe.github.io/TheraEpiFlow/) walks through the same flow.

**Project menu** (shown by `python main.py` with no project open):

| Key | Action |
|---|---|
| `1`..`N` | open a project by its number |
| `n` | create a new project |
| `e` | edit a project's track configuration |
| `d` | delete a project |
| `demo` | run a guided demo (HPV16 E7) through the whole pipeline |
| `q` | quit |

If you have just installed and want to see the pipeline run without configuring a project, type `demo`. It seeds a small HPV16 E7 project and runs every step end to end in a couple of minutes. If something in the environment is missing it stops at that step and tells you what to fix, so it is also the quickest way to confirm the install works.

**Step menu** (shown once a project is open):

| Key | Action |
|---|---|
| `Enter` | run the next pending step |
| `a` | run all remaining steps |
| `r` | redo from a chosen step (re-asks its config and wipes every output after it) |
| `h` | show the full intro of the next step |
| `b` | browse intermediate files |
| `t` | edit a track's configuration |
| `s` | show the status table |
| `z` | open the download menu (`.zip` on WSL, `.tar.gz` on Linux) |
| `q` | quit |

### Scripted / non-interactive use

The same actions are available as flags, for CI or scripting (the interactive menu is the normal path):

```bash
python main.py --new-project                                # create a project
python main.py --project hpv_study                          # open the REPL for a project
python main.py --project hpv_study --step cluster_epitopes  # run one step and exit
python main.py --project hpv_study --status                 # show progress
python main.py --list                                       # list every project
```

## Quick self-test after install

If everything in `setup.sh` ran cleanly, validate the full pipeline end-to-end against a tiny reference protein (HPV16 E7, 98 aa, 1 track):

```bash
conda activate TheraEpiFlow
python tests/selftest.py
```

The script seeds a disposable project (`bench_selftest`), runs every step in non-interactive mode, and prints a pass/fail summary. It requires network access for the IEDB API calls (`predict_binding`, `consensus_filter`). On a mid-range laptop expect roughly 80–120 seconds wall time.

If a step fails, its error message points at the underlying issue. Start there, then check [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) for the most common install pitfalls (IEDB API connectivity, WSL conda activation, Windows path quirks).

The project remains on disk afterward at `projects/bench_selftest/` so you can inspect intermediate files and open the generated HTML report in a browser.

The heavy validation suite (Experiments 1 and 2, about 3 hours wall on a Ryzen 5 mobile, used in the master's thesis to characterize the pipeline) is at `tests/validation/run_experiment_1.py` and `run_experiment_2.py`. End users do not need to run it; published figures live under `tests/validation/figures/`.

## How a project is organized

```
projects/
  {project_name}/
    project_config.json     Shared settings, extended as each step adds its own config
    pipeline.json           Per-track step states (done, error, pending)
    data/
      input/{track_id}/                       Reference FASTA + sequence registry
      intermediate/{track_id}/
        predictions/  PRED_NET_*.csv, PRED_FLURRY_*.csv, audit JSON
        consensus/    CONSENSUS_*.csv, CONSENSUS_IMMUNOGENIC_*.csv, audit JSON + stage CSVs
        toxicity/     TOXICITY_ALL_*.csv, TOXICITY_SAFE_*.csv, audit JSON
        clusters/     CLUSTER_*.csv, CLUSTER_REPR_*.csv (with ★ column), audit JSON + XLSX
        variants/     VARIANTS_*.fasta (permanent cache), audit JSON
        conservation/ CONSERVATION_*.csv/.xlsx, CONSERVATION_HEATMAP_*.png, CONSERVATION_MUTATIONS_*.xlsx, audit JSON
        coverage/     COVERAGE_*.csv/.xlsx, COVERAGE_HIT_CHART_*.png, audit JSON (population_coverage)
        murine/       MURINE_*.csv, MURINE_AGG_*.csv, CURATE_MURINE_*.csv, audit JSON (predict_murine + curate_murine)
      output/
        MASTER_TABLE_FULL_{project}.xlsx    (created by integrate_data, audit-grade)
        MASTER_TABLE_VIEW_{project}.xlsx    (created by integrate_data, styled VIEW)
        MASTER_TABLE_VIEW_{project}.csv     (created by integrate_data, portable CSV)
        MASTER_TABLE_AUDIT_{project}.json   (created by integrate_data)
        REPORT_{project}.html               (created by generate_report, offline calculator)
      downloads/
        {project}_full_{stamp}.tar.gz       (created by the [z] download menu, optional destination)
```

A track is one organism plus one protein. All tracks in a project share the same HLA alleles and pipeline parameters. Track IDs follow `{ORGANISM_LABEL}_{PROTEIN_LABEL}`, for example `HPV16_E6` or `SARS2_S`. Labels are suggested automatically based on standard abbreviations and can be overridden when the project is created.

## Configuration as you go

The CLI never asks for parameters you do not need yet. Each step collects its own settings the first time it runs and saves them to `project_config.json`, so reruns stay deterministic.

| Stage | Asked at first run |
|---|---|
| `--new-project` | Project name and optional description (target host is fixed at Homo sapiens, since HLA-I is human) |
| `fetch_sequences` | Organism aliases, proteins, labels, input source (UniProt search or local FASTA) |
| `predict_binding` | HLA alleles (default: 27 globally diverse) and peptide lengths |
| `consensus_filter` | Presentation percentile threshold (default 2.0%) |
| `screen_toxicity` | ToxinPred3 cutoff (default 0.38) |
| `cluster_epitopes` | Identity cutoff (default 80%) and clustering method |
| `search_variants` | Scope (intraspecific or interspecific), optional host filter, optional family taxid |
| `analyze_conservation` | Analysis threshold (default 1.0, meaning exact match) |
| `population_coverage` | Populations to evaluate, plus an optional informational coverage cutoff |
| `predict_murine` | Murine strain group (`default`, `C57BL/6`, `BALB/c`, `CBA/C3H/AKR`, or `all`) |
| `integrate_data` | Which columns the report VIEW should contain (default selection, or pick from the catalog) |

## Implemented steps

### fetch_sequences

Pulls the reference protein sequence for each track from the UniProt REST API. Maps short aliases (`HPV16`, `MPOX`, `CHIKV`, `SARS2`) to scientific names and taxonomy IDs, runs three search strategies (protein name, gene name, organism only), prefers Swiss-Prot entries, and flags polyproteins or fragments based on the median length of the result set. Saves a cleaned FASTA, a JSON registry with metadata, and a validation report. See `modules/fetch_sequences/README.md`.

### predict_binding

Runs NetMHCpan 4.1 EL via IEDB's classic API and MHCFlurry 2.0 locally. Both predictors run in parallel (NetMHCpan in a `ThreadPoolExecutor` worker, MHCFlurry on the main thread). The IEDB call iterates over the cartesian product of alleles and peptide lengths. MHCFlurry loads `Class1PresentationPredictor` once per process and calls `.predict()` once per allele, keeping the per-(peptide, allele) granularity that consensus_filter expects. Output: `PRED_NET_{track_id}.csv`, `PRED_FLURRY_{track_id}.csv`, and an audit JSON.

### consensus_filter

Two stages.

Stage 1 keeps only peptides predicted as presented by both NetMHCpan and MHCFlurry. Per-peptide aggregation collapses N rows (one per allele) into one row per peptide, recording the best allele and the joined HLA list.

Stage 2 sends the surviving peptides to the IEDB Calis 2013 immunogenicity endpoint and keeps those with `score > 0`.

Audit files preserve every intermediate state for reproducibility.

### screen_toxicity

Loads the ToxinPred3 Model 1 (AAC + DPC features, Extra Trees) once with joblib and predicts in memory, no temporary files or subprocess. Computes the ML score and the calibrated PPV. Default cutoff is 0.38. Drops rows with non-canonical residues (`X`, `B`, `Z`) before predicting. Outputs a full table and a safe-only table for the next step. See `modules/screen_toxicity/README.md`.

### cluster_epitopes

Groups safe epitopes by pairwise sequence identity using Biopython global alignment and NetworkX. Three methods available; `cluster_break` is the default.

`cluster_break` starts from connected components and iteratively removes the minimum-weight edge from any component whose average intra-cluster similarity falls below the threshold. This avoids "bridge" groupings where peptides with distinct TCR recognition land in the same cluster. It matters because pMHC-TCR recognition is highly sensitive to anchor position differences.

Default identity threshold: **80%**. See `modules/cluster_epitopes/README.md`.

### select_representatives

Selects one representative peptide per cluster using a combined score: minimum percentile across all allele-level values from both tools (`best_combined_percentile`) plus allele breadth (`num_alleles_united`). Marks the winner with `★` in the `best_representative` column. Produces a colour-coded XLSX deliverable alongside the CSV. See `modules/select_representatives/README.md`.

### search_variants

Queries UniProt REST API for protein variants. Two scopes:

- **Intraspecific**: variants within the same taxonomy ID (e.g. different HPV16 isolates, SARS-CoV-2 strains)
- **Interspecific**: same protein name across related species, with optional host filter (e.g. E5 from HPV16, HPV18, HPV31 infecting Homo sapiens)

The multi-FASTA is written once and cached permanently. On re-run the user can choose to keep or regenerate the file. Identity to the reference is computed and shown, but it never excludes: near-identical (≥ 99%), possibly-unrelated (< 30%) and uncut-polyprotein candidates are only **flagged**, and the user decides which to keep during selection (intraspecific conservation usually *wants* the near-identical isolates). Exact-duplicate sequences are collapsed to one row. See `modules/search_variants/README.md`.

### analyze_conservation

Measures how faithfully each ★ representative appears across the variant sequences from `search_variants`. Uses a sliding window: for each variant the best-matching window of the same length as the peptide is found and its per-position identity is recorded.

The analysis threshold is configurable (default 1.0 = exact match) and feeds the `pct_identity_threshold` summary column. The `conservation_label` (perfect / high / moderate / low / conservation_unknown) and all row colours are based on `mean_max_identity` and never change regardless of the chosen threshold. It is qualitative and removes no epitopes.

For every (epitope, variant) pair with 1–2 substitutions it also emits a mutation-tolerance verdict (MHC-I anchors P2/PΩ + BLOSUM62 chemistry) in `CONSERVATION_MUTATIONS_{track_id}.xlsx`, and a dual-panel `CONSERVATION_HEATMAP_{track_id}.png` (per-position conservation + identity tiers). See `modules/analyze_conservation/README.md`.

### population_coverage

For each ★ epitope, computes the fraction of one or more human populations that carries at least one of the epitope's HLA alleles, using a vendored copy of the IEDB Population Coverage tool's allele-frequency pickle (IEDB v3.0.2, pickle v1.1.2, sourced from AlleleFrequencies.net). Diploid model: `p_locus = 1 − (1 − q)²` per locus, combined as `1 − Π(1 − p_locus)`. Math is bit-for-bit reproducible against `iedb.org/population_coverage`.

Supports multiple populations in a single run; per-population output includes an IEDB-style Hit Chart PNG plus a long-format CSV. With two or more populations, a comparative heatmap PNG is also generated. See `modules/population_coverage/README.md` and `modules/population_coverage/data/SOURCE.md` for the pickle's provenance.

### predict_murine

Re-runs NetMHCpan EL + MHCFlurry on the ★ representatives against murine **H-2** alleles instead of human HLA-I, to flag which human-selected epitopes would also bind in a mouse model. Five strain groups live in `config.MURINE_ALLELES`: `default` (C57BL/6 + BALB/c, 5 alleles), `C57BL/6` (H-2-Db, H-2-Kb), `BALB/c` (H-2-Dd, H-2-Kd, H-2-Ld), `CBA/C3H/AKR` (H-2-Dk, H-2-Kk, the H-2k haplotype), and `all` (the 7-allele consolidated set). Each ★ epitope gets per-allele records plus an aggregated row with best percentile, H-2 alleles bound, and a four-tier binder label (`optimal`, `good`, `borderline`, `non_binder`). It is qualitative and never removes epitopes. See `modules/predict_murine/README.md`.

### curate_murine

Assembles the per-track master table by joining each ★ representative's human qualification with its conservation, population coverage and (when present) murine prediction. It is a JOIN-only step with no ranking or priority relabelling. Conservation and coverage are required inputs, murine is optional. Output: `CURATE_MURINE_{track_id}.csv` (one row per ★ epitope with all annotations), a slim view, and an audit JSON.

### integrate_data (global)

Runs once after every track has finished `curate_murine`. Stacks the per-track tables into a project-wide master table (`MASTER_TABLE_FULL_{project}.xlsx`, every column) plus a user-configurable VIEW (`MASTER_TABLE_VIEW_{project}.xlsx` and `.csv`) with display-ready headers. The user picks once which columns belong to the VIEW, and the choice is persisted to `project_config.step_overrides.integrate_data.view_columns` and re-used on reruns. To re-open the column prompt, redo the step with the `r` (redo) command in the REPL. It also writes `MASTER_TABLE_AUDIT_{project}.json` capturing tracks integrated/skipped, coverage populations, and the chosen columns. See `modules/integrate_data/README.md`.

### generate_report (global)

Renders an offline interactive HTML calculator (`REPORT_{project}.html`) from the master tables and the IEDB allele-frequency pickle. Every dataset is inlined as JSON and the in-browser exports use the native `Blob` API, so the report works in any browser with no network and no bundled JS library. The user filters/sorts epitopes and clicks rows to build the selection; the sidebar updates live with the organism × protein coverage map, cumulative population coverage, and the construct builder (editable linker `GPGPG`, optional adjuvant and His-tag with N/C/off placement). Three inline exports are offered: `.tsv` (selected epitopes), `.fasta` (construct + epitopes), and `matrix` (allele × population frequencies). Header chips, an SVG progress ring and per-allele tooltips mirror the original vaxbuilder prototype. See `modules/generate_report/README.md`.

### Download menu (REPL key `[z]`)

This is not a pipeline step. It packages the project **on demand** from the REPL. Interactive prompts ask scope (full project or one earlier step's outputs across all tracks), opt-in for the heavy `predictions/` folder, and destination, offering the in-project `downloads/` folder always, `~/Downloads` when it exists, and the Windows-side Downloads folder when running under WSL (auto-detected via `/proc/version`, `WSLENV`, and `/mnt/c/Users/{user}`).

Format is chosen automatically: **`.zip`** when running under WSL (natively openable in Windows Explorer with a double-click), **`.tar.gz`** on pure Linux or macOS. Lives at [`utils/download_ui.py`](utils/download_ui.py); the archive is written by `utils/archive.py` (stdlib `tarfile` / `zipfile`, no extra dependency).

## Master-table strategy

`curate_murine` produces a per-track joined table. `integrate_data` stacks those across all tracks into `output/MASTER_TABLE_FULL_{project}.xlsx` (every column) plus a user-configurable VIEW (`MASTER_TABLE_VIEW_{project}.xlsx` / `.csv`) that feeds `generate_report`. The pipeline emits its files in the right shape: every implemented step except `search_variants` (FASTA, by design, consumed by `analyze_conservation`) and `population_coverage` (long format, by design, pivoted at merge time) produces a one-row-per-peptide CSV with `peptide` as the primary key. Per-track assembly order:

```
CONSENSUS_IMMUNOGENIC_{track}.csv     (1 row / peptide, the base)
  + TOXICITY_SAFE_{track}.csv         (appends toxicity columns)
  + CLUSTER_{track}.csv               (appends cluster_id, etc.)
  + CLUSTER_REPR_{track}.csv          (filters to ★ rows, appends norm scores)
  + CONSERVATION_{track}.csv          (appends conservation_label etc.)
  + COVERAGE_{track}.csv (pivoted)    (long → wide: coverage_World, coverage_Brazil, …)
```

Track context (`track_id`, `organism_label`, `protein_label`) comes from `project_config.json`. The master table is what feeds the HTML report's interactive selection calculator (`generate_report` step).

## Validation

The pipeline is validated at three levels:

- **Automated end-to-end self-test.** `python tests/selftest.py` seeds a disposable project (HPV16 E7, 98 aa, 1 track) and runs all 13 steps in non-interactive mode, asserting each produces its expected outputs (about 80 to 120 s on a mid-range laptop; needs the IEDB API for `predict_binding` and `consensus_filter`).
- **Multi-track end-to-end runs.** Most recently `hpv16` (5 tracks, 86 ★ peptides) and `scer_test` (2 tracks, 83 ★ peptides), both producing the per-track master tables, the project-wide `MASTER_TABLE_*`, and the self-contained HTML report without manual intervention.
- **Manual verification.** The generated `REPORT_{project}.html` was opened in a browser and exercised end-to-end: filtering and sorting, epitope selection, the live organism × protein coverage map, cumulative population-coverage updates, the construct builder (linker, adjuvant, and tag with N/C placement), and the `.tsv` / `.fasta` / `matrix` exports.
- **Characterization suite** (thesis). Experiments 1 and 2 under `tests/validation/` quantify determinism and IEDB API jitter; published figures live in `tests/validation/figures/`.

## Citation

If you use TheraEpiFlow, please cite the archived release and the underlying method papers (see each step's README for the relevant tool citations).

- **Software archive:** TheraEpiFlow, Zenodo. DOI: _to be assigned on first release_ (`https://github.com/Davirbe/TheraEpiFlow`).

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
| integrate_data | implemented, project-wide master table (FULL + customizable VIEW + AUDIT) |
| generate_report | implemented, offline HTML calculator with vaxbuilder-style header chips + sidebar heatmap + per-allele tooltips |

Download is a menu utility (`[z]` in the REPL → `utils/download_ui.py`), not a pipeline step.

All step modules follow the per-step layout described above (Single Responsibility split), with no behaviour change verified module-by-module. The pipeline was re-validated end-to-end on `hpv16` (5 tracks, 86 ★ peptides) and `scer_test` (2 tracks, 83 ★ peptides), producing the master tables and the interactive HTML calculator without intervention.

Earlier validation projects: `hpv_study` (HPV16/18/31 × E5/E6/E7), `chikungunya_study` (CHIKV M, NSP1), `dengue_study` (DENV E, NS5), `monkeypox_study` (MPOX M), `oropouche2` (OROV N), `sars2_test` (SARS2 S, N).
