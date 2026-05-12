# population_genotype_map.p — source

**File:** `population_genotype_map.p` (~3.6 MB binary pickle)

**Origin.** IEDB Population Coverage tool **v3.0.2** (LATEST as of
2026-05-12), released **2023-03-07**, downloaded from
`https://downloads.iedb.org/tools/population/3.0.2/IEDB_Population_Coverage-3.0.2.tar.gz`.

**Internal pickle version.** `population-coverage-pickle v1.1.2`. The
historical SVN tag is `https://trac.liai.org/svn/iedbtools/djangotools-deps/population-coverage-pickle/tags/v1.1.2`
(server now retired).

**Upstream data source.** AlleleFrequencies.net (AFND). The IEDB
consolidates published HLA-typing studies from AFND into this single
pickle.

**Contents used by this module.**
- 2 MHC classes (`I`, `II`); only class `I` is used.
- 239 populations under class `I` (e.g. `World`, `Brazil`, `East Asia`).
- 5 loci: `HLA-A`, `HLA-B`, `HLA-C`, `HLA-E`, `HLA-G`.
- ~3,245 allele entries across all populations.

The pickle stream contains three serialised objects in order:
`population_coverage`, `country_ethnicity`, `ethnicity`. We use only the
first.

**Vendored on.** 2026-05-12.

**License.** IEDB tools are free for academic and non-commercial use.
See `https://www.iedb.org/` for current terms.

**Citation.** Bui H-H, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A.
*Predicting population coverage of T-cell epitope-based diagnostics and
vaccines.* BMC Bioinformatics. 2006;7:153.

## How to update

1. Check `https://downloads.iedb.org/tools/population/` for a version
   newer than 3.0.2.
2. Download the tarball and extract
   `population_coverage_pickle/population_genotype_map.p`.
3. Replace this file in place.
4. Re-run the math regression described in the module README and
   confirm coverage values for the test project still match the IEDB
   website. If they shift, decide whether the upgrade is worth the
   reproducibility break and document the change in the audit JSON of
   subsequent runs.
