# generate_report

The final global step. Reads the master tables produced by `integrate_data` and the vendored IEDB allele-frequency pickle, and renders a single self-contained interactive HTML calculator under `data/output/REPORT_{project}.html`.

## Why a self-contained HTML

Researchers and collaborators don't always have Python set up; they want to open a file and start picking epitopes. The report is the pipeline's user-facing deliverable. Every dataset (epitopes, project metadata, allele frequencies for the chosen populations) is inlined as JSON, and the in-browser exports are produced with the native `Blob` + `URL.createObjectURL` API (no JSZip or any other bundled JS library). The only remote dependency is Google Fonts, which degrades gracefully to a system font.

## SRP layout

| File | Role |
|---|---|
| `__init__.py` | Facade that re-exports `GenerateReportStep`. |
| `step.py` | `GenerateReportStep(BaseGlobalStep)`: orchestrates load → meta → payload → render. Owns the pre-step page ClassVars. |
| `core.py` | Pure data prep: reads VIEW/FULL/AUDIT, builds the per-peptide JSON list and project_meta. No Jinja, no Rich. |
| `coverage_db.py` | Reduces the IEDB pickle to `{population: {allele: freq}}` restricted to alleles actually used by the project. Reuses `modules.population_coverage.core`'s loader. |
| `io.py` | Jinja2 render to `data/output/REPORT_{project}.html`. |
| `templates/calculator.html.j2` | The HTML/CSS/JS app, generalized from `existing_scripts/vaxbuilder_updated.html` to any organism × protein matrix. Exports use the native `Blob` API; no bundled JS library. |

## In-browser flow

1. **Filter bar:** Organism, Protein, Conservation, and free-text peptide search.
2. **Table:** sortable, with hover tooltips on Conservation / HLA / Murine cells and pill colors for label columns and coverage bands. Click a row to add or remove it from the construct selection; selected rows are highlighted and the sidebar updates live.
3. **Coverage map:** organism × protein distribution heatmap, always visible in the sidebar. Each cell shows how many selected epitopes cover that organism-protein pair. Sub-labels (`missing` / `1 ep.` / `N ep.`) and a per-organism totals column (`covered / total proteins`) update on every selection change.
4. **Population coverage:** scrollable table below the coverage map; cumulative diploid coverage per population (from `project_config`) computed as `1 − Π(1 − f_i)` over the union of `alleles_united`. Updates live with the selection.
5. **Construct builder:** three configurable inputs:
   - **Linker** (editable, default `GPGPG`): inserted between every adjacent epitope chip.
   - **Adjuvant** (editable, no default): prepended (N) or appended (C) to the whole construct; Off by default. PADRE (`AKFVAAWTLKAAA`) is a common choice.
   - **Tag** (editable, default `HHHHHH`): his-tag or other purification tag, same N/C/Off toggle.
   Live chip-rendered preview shows every part of the construct colour-coded (epitope / linker / tag / adjuvant). Size bar shows total length in aa and nt plus a breakdown by component.
6. **Export:** three buttons.
   - `.tsv`: selected epitopes with all annotation columns, UTF-8 BOM, tab-separated.
   - `.fasta`: one FASTA record per selected epitope plus, when a non-empty construct exists, the full vaccine construct sequence.
   - `.matrix`: allele × population HLA frequency matrix (sparse) for the alleles bound by the selected epitopes.

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

Open the generated `REPORT_*.html` in a browser (offline). Select a few epitopes from different tracks, and the coverage map and population coverage table should update live. Fill in an adjuvant (e.g. `AKFVAAWTLKAAA`), set its position to N or C, adjust the tag, and confirm the chip bar and size breakdown reflect the change. Test each export button (`.tsv`, `.fasta`, `.matrix`) and verify the downloaded files are non-empty and parseable.
