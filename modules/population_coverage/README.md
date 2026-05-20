# population_coverage

Computes, for every ★ representative epitope, the fraction of one or more human populations that carries at least one HLA-I allele the epitope binds. The step is **qualitative** — no epitopes are removed; coverage is recorded so the final HTML report can filter or rank them.

## Method: per-epitope diploid coverage

Replicates the IEDB Population Coverage tool's formula at the per-epitope level, using a vendored copy of the IEDB allele-frequency database (`data/population_genotype_map.p`, ~3.6 MB, snapshot of IEDB SVN `population-coverage-pickle v1.1.2`).

For each population × epitope pair:

```
For each locus the epitope binds (HLA-A, HLA-B, HLA-C, …):
    q       = Σ frequency(allele)   for each allele in the epitope's HLA list
    p_locus = 1 - (1 - q)²                         (probability of cover, diploid)

Combine independent loci:
    p_total      = 1 - Π (1 - p_locus)
    coverage_pct = p_total · 100
```

Frequencies are normalised per locus when the raw locus total exceeds 1.0 (matches the source data). Alleles that are not present in the IEDB database are silently skipped and counted in the audit JSON.

## Code layout

The step is split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `PopulationCoverageStep` orchestration — `preflight` / `run` / `postflight` / `describe_outputs` |
| `core.py` | Database loading + diploid coverage math (pure functions) |
| `prompts.py` | Interactive population / cutoff selection |
| `io.py` | CSV / XLSX writers |
| `charts.py` | Hit-chart + comparative-matrix PNGs |
| `render.py` | Rich console summary |
| `__init__.py` | Facade — re-exports `PopulationCoverageStep` |

## Inputs

- `clusters/CLUSTER_REPR_{track_id}.csv` (from `select_representatives`) — only rows with `BEST_REPRESENTATIVE == "★"` are used. Required columns: `peptide`, `alleles_united` (semicolon-joined IMGT, e.g. `HLA-A*02:01;HLA-B*07:02`), `num_alleles_united`.

## Configuration (per project)

Set during the interactive preflight, saved into `project_config.json`:

| Key | Default | Description |
|---|---|---|
| `coverage_populations` | `["World"]` | One or more population names from the IEDB database (239 available for class I). |
| `coverage_minimum_pct` | `null` | Informational cutoff. Used only to highlight low-coverage epitopes in the postflight summary and downstream HTML report — never filters here. |

Curated short list shown in the preflight prompt:
`World`, `Europe`, `East Asia`, `South Asia`, `Africa`, `North America (Native American)`, `South America`, `Brazil`. Free-form names (any of the 239 populations) accepted via the `other` option.

## Outputs

Under `data/intermediate/{track_id}/coverage/`:

| File | Contents |
|---|---|
| `COVERAGE_{track_id}.csv` | Long-format summary, one row per `(peptide, population)`. Columns: `peptide, population, mhc_class, coverage_pct, n_hlas_used, n_hlas_in_db, alleles_united`. |
| `COVERAGE_{track_id}.xlsx` | Same data; `coverage_pct` cell coloured by band (green ≥ 50%, yellow ≥ cutoff or 10%, red < cutoff, gray = no matching HLA). |
| `COVERAGE_DETAIL_{population}_{track_id}.csv` | IEDB-style per-population CSV: row per epitope, columns = each HLA seen in any ★ epitope (sorted by frequency desc), cells `+`/`-`, plus an `Epitope set` totals row. |
| `COVERAGE_HIT_CHART_{population}_{track_id}.png` | IEDB-style hit chart (PNG) — title + locus-grouped HLA headers + `+`/`-` markers + coverage column + totals row. |
| `COVERAGE_MATRIX_{track_id}.png` | Comparative heatmap (rows = epitopes, columns = populations), produced only when ≥ 2 populations are selected. |
| `COVERAGE_AUDIT_{track_id}.json` | Run metadata: timestamp, populations, MHC class, allele-match stats, mean coverage per population, output paths. |

## Algorithm cross-check

The math has been validated bit-for-bit against the user's prototype (`existing_scripts/cobertura_por_epitopo.py`), which itself was validated against the IEDB website. 117 coverage values (39 ★ epitopes × 3 populations) on `sars2_test/SARS2_N` matched with zero deviation.

## Notes

- This pipeline is MHC-I only; `mhc_class` is hard-coded to `"I"`.
- The pickle has 3 sequential objects; we load the first (`population_coverage`) and discard the other two (`country_ethnicity`, `ethnicity`).
- Alleles in `alleles_united` and in the pickle are both IMGT format (`HLA-A*02:01`) — no conversion needed.
- The vendored `.p` file is required at runtime; it lives inside the module to keep distribution self-contained.

## Data source and provenance

The allele-frequency database at `data/population_genotype_map.p` is a vendored copy of the IEDB Population Coverage tool's bundled pickle.

| Attribute | Value |
|---|---|
| IEDB tool version | 3.0.2 (LATEST, 2023-03-07) |
| Internal pickle version | `population-coverage-pickle v1.1.2` |
| Upstream raw data | AlleleFrequencies.net (AFND) |
| Vendored on | 2026-05-12 |
| Coverage | 2 MHC classes, 239 populations (class I), 5 loci, ~3,245 alleles |
| Source URL | `https://downloads.iedb.org/tools/population/3.0.2/` |

This version is kept because every `coverage_pct` we produce is bit-for-bit reproducible on `https://iedb.org/population_coverage` for the same input. Switching to AFND (which has finer-grained data, ~163,695 frequencies, and the latest IMGT/HLA nomenclature update on 2026-01-01) would require building our own pickle: AFND has no public bulk download as of 2026-05, only the per-study web interface. That migration would also break parity with the IEDB website, which is the standard reference in the field.

### Updating the pickle

1. Check `https://downloads.iedb.org/tools/population/` for a release above 3.0.2.
2. Download the tarball and extract `population_coverage_pickle/population_genotype_map.p`.
3. Replace `data/population_genotype_map.p` with the new file.
4. Re-run the math cross-check on a known epitope (see "Algorithm cross-check" above). If outputs shift, decide whether to keep the new version and note the change in the audit JSON of subsequent runs.

See `data/SOURCE.md` for the citable per-file record.
