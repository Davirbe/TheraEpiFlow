# generate_report

The final global step. Reads the master tables produced by `integrate_data` and the vendored IEDB allele-frequency pickle, and renders a single self-contained interactive HTML calculator under `data/output/REPORT_{project}.html`.

## Why a self-contained HTML

Researchers and collaborators don't always have Python set up — they want to open a file and start picking epitopes. The report is the pipeline's user-facing deliverable. Every dataset (epitopes, project metadata, allele frequencies for the chosen populations) is inlined as JSON; JSZip is vendored inside the template; the only remote dependency is Google Fonts, which degrades gracefully.

## SRP layout

| File | Role |
|---|---|
| `__init__.py` | Facade — re-exports `GenerateReportStep`. |
| `step.py` | `GenerateReportStep(BaseGlobalStep)` — orchestrates load → meta → payload → render. Owns the pre-step page ClassVars. |
| `core.py` | Pure data prep — reads VIEW/FULL/AUDIT, builds the per-peptide JSON list and project_meta. No Jinja, no Rich. |
| `coverage_db.py` | Reduces the IEDB pickle to `{population: {allele: freq}}` restricted to alleles actually used by the project. Reuses `modules.population_coverage.core`'s loader. |
| `io.py` | Jinja2 render to `data/output/REPORT_{project}.html`. |
| `templates/calculator.html.j2` | The HTML/CSS/JS app — generalized from `existing_scripts/vaxbuilder_updated.html` to any organism × protein matrix. |
| `templates/jszip.min.js` | Vendored JSZip 3.10.1 (~96 KB) — included inline so the ZIP export works offline. |

## In-browser flow

1. **Filter bar** — Organism / Protein / Conservation / free-text peptide search.
2. **Table** — 13 default columns, click-to-sort, hover tooltips on Conservation / HLA / Murine cells, pill colors for label columns and coverage bands.
3. **Construct builder (right panel)** — three configurable inputs: TAG (default `HHHHHH`, position selector Off / N-term / C-term), Adjuvant (default `PADRE` = `AKFVAAWTLKAAA`, same selector), Linker (default `AAY`). Live chip-rendered preview.
4. **Finalize construct → modal**:
   - Organism × protein coverage heatmap (cell value = n selected epitopes that hit that pair).
   - Cumulative population coverage per population in `project_config`, computed as `1 - Π(1 - f_i)` over the union of `alleles_united`.
   - Construct stats: length aa/nt, molecular weight (monoisotopic table in JS), mean conservation 100%, HLA union count, linker / TAG / adjuvant in use.
   - FASTA preview.
5. **Download ZIP** — packages five files: `vaccine_construct.fasta`, `selected_epitopes.csv` (full 46 columns for chosen peptides, with utf-8 BOM + quoted strings), `construction_stats.json`, `coverage_heatmap.png` (canvas-rendered), `selection_summary.txt`.

## Reuse map

- `BaseGlobalStep` → `modules/base_step.py`
- IEDB pickle loader + per-locus normalization → `modules/population_coverage/core.py`
- Master-table file naming → `data/output/MASTER_TABLE_*_{project}.{csv,xlsx,json}` produced by `modules/integrate_data/`
- Conservation label palette → mirrors `modules/analyze_conservation/io.py` (`perfect / high / moderate / low / conservation_unknown`)
- Coverage palette → mirrors `modules/population_coverage/io.py` (`high ≥ 80% / moderate ≥ 50% / low`)

## Verification

```bash
python main.py --project hpv16     --step generate_report
python main.py --project scer_test --step generate_report
pyflakes modules/generate_report
```

Open the generated `REPORT_*.html` in a browser (offline-tested). Select a few epitopes, click `Finalize construct`, download the ZIP, verify all 5 files are present.
