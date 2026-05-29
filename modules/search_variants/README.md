# search_variants

Retrieves protein sequence variants from UniProt REST API and writes a permanent multi-FASTA that `analyze_conservation` uses to measure epitope conservation across natural sequence diversity.

## What it does

Queries UniProt for proteins related to each track's reference sequence. Two search scopes:

| Scope | Query strategy | Typical use |
|---|---|---|
| Intraspecific | `taxonomy_id:{tax_id} AND protein_name:"{name}"` | Strains and isolates of the same species, e.g. different HPV16 submissions or SARS-CoV-2 variants |
| Interspecific | `protein_name:"{name}"` ± optional host filter | Same protein across related species, e.g. E5 from HPV16, HPV18, HPV31 infecting Homo sapiens |

The scope and optional host filter are asked interactively once per track and saved to `project_config.json["tracks"][track_id]` for reproducibility.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `SearchVariantsStep` orchestration: `run` and `describe_outputs` |
| `core.py` | HTTP helper, taxonomy lineage, UniProt search, identity, record building |
| `io.py` | Reference loading (+ chain-slice detection) + empty-output writer |
| `prompts.py` | Cache-redo, scope/host filter, tax-id guard, candidate multi-selection, grouped/flat view choice |
| `render.py` | Rich candidates table + per-genotype grouped table |
| `__init__.py` | Facade that re-exports `SearchVariantsStep` |

## Filtering pipeline

After the raw UniProt results come back:

1. **Polyprotein slicing (gated):** only when the track needs it (the reference was itself sliced out of a polyprotein, or the user supplied a local FASTA). Each candidate that carries a Chain feature matching the protein name is sliced down to that mature chain, so identity is computed chain-vs-chain. Clean single-protein UniProt tracks skip this entirely (no `ft_chain` request).
2. **Reference exclusion:** removes any accession that `fetch_sequences` already fetched (reads the registry JSON).
3. **Identity + advisory flags (no exclusion):** pairwise identity is computed for every candidate and shown in the table, but it never excludes. Two advisory flags are set so the user decides during selection (intraspecific conservation *wants* the near-identical isolates):
   - `possibly_unrelated`: identity `< VARIANTS_UNRELATED_IDENTITY_PERCENT` (30%, Rost 1999 twilight zone).
   - `near_identical`: identity `≥ VARIANTS_NEAR_IDENTICAL_PERCENT` (99%).
   - `uncut_polyprotein`: slicing was expected but no Chain matched and the sequence is > 1.5× the reference length (likely an unmatched polyprotein name).
4. **Exact-duplicate collapse:** byte-identical sequences (the same sliced chain re-deposited under different accessions) collapse to one row, preferring the reviewed (Swiss-Prot) copy.
5. **Identity sort:** remaining candidates are sorted descending by identity to the reference.
6. **Interactive selection:** table displayed; user selects by index range (e.g. `1,3,5-8`, `all`, `none`). For interspecific scope the user can switch to a per-genotype grouped view (one best representative per genotype, so other types aren't buried under the reference's own isolates). Non-interactive mode selects all.
7. **Lenient revalidation:** only rejects empty sequences or those composed entirely of ambiguous residues (X/B/Z/U/O). Short sequences (< 50 aa) are kept with a warning (valid for conservation, not for binding prediction).

## Cache behaviour

The FASTA is **permanent**. When the file already exists the step asks whether to keep it or redo the search. In non-interactive mode the existing file is always kept and the step completes instantly. This avoids re-running expensive API searches on every pipeline re-run.

## Input

- `input/{track_id}/SEQUENCES_{track_id}.fasta`: reference sequence (from `fetch_sequences`)
- `input/{track_id}/REGISTRY_{track_id}.json`: reference accession list (for exclusion)

## Output

| File | Contents |
|---|---|
| `variants/VARIANTS_{track_id}.fasta` | Multi-FASTA of selected variant sequences, the input for `analyze_conservation`. |
| `variants/VARIANTS_VIEW_{track_id}.csv` | Slim per-step view, one row per selected variant: `accession, organism, protein, length, identity_to_seed`. |
| `variants/VARIANTS_AUDIT_{track_id}.json` | Scope, host filter, counts, selected accessions. |

An empty FASTA is written (with a note in the audit) when no variants are found after filtering, so downstream steps can always expect the file to exist.

## Scope, family restriction and host filter

Three fields are saved in `project_config["tracks"][track_id]` after the first interactive run:

| Field | Applies to | Purpose |
|---|---|---|
| `variants_scope` | all | `"intraspecific"` or `"interspecific"` |
| `variants_family_taxid` | interspecific only | Restrict to a virus family/genus (e.g. Orthopoxvirus 10242, Coronaviridae 11118) |
| `variants_host_filter` | interspecific only | Filter by host organism (e.g. `"Homo sapiens"`) |

**Query built per configuration:**

| Scope | `family_taxid` | Query |
|---|---|---|
| `intraspecific` | (ignored) | `(taxonomy_id:{species_tax_id}) AND (protein_name:"{name}")` |
| `interspecific` | set | `(taxonomy_id:{family_taxid}) AND (protein_name:"{name}")` ± host filter |
| `interspecific` | null | `(protein_name:"{name}")` ± host filter |

When scope is `intraspecific`, both `variants_family_taxid` and `variants_host_filter` should be `null`; the species taxonomy already restricts results.

## Taxonomic restriction (family_taxid): interactive auto-suggestion

When selecting interspecific scope interactively, the step fetches the full taxonomic lineage of the track's `tax_id` from the UniProt taxonomy API and presents genus/family/order options to choose from:

```
Lineage options (closest ancestor first):

  1  Orthopoxvirus      [genus]   tax_id: 10242
  2  Chordopoxviridae   [family]  tax_id: 10240
  3  Poxviridae         [family]  tax_id: 10240
  4  No restriction (all organisms)

Choice (1-4 or custom tax_id, Enter=4):
```

The chosen tax_id is saved as `variants_family_taxid` and reused on subsequent runs without prompting.

## Biological interpretation of scope

| Scope | Question answered | Typical outcome |
|---|---|---|
| Intraspecific | Are these epitopes conserved across strains/isolates of this virus? | High conservation expected for stable vaccine targets |
| Interspecific + family | Are these epitopes conserved across the virus family? | Lower conservation expected (a different question) |

**These answer different research questions.** Low conservation in an interspecific family analysis does not indicate a poor vaccine target; it means the epitope is virus-specific rather than family-conserved. For strain coverage (the primary vaccine concern), use intraspecific.

## Choosing the right scope

| Situation | Recommended |
|---|---|
| Strain/variant coverage for vaccine design | `intraspecific` |
| Cross-species protection potential (e.g. pan-poxvirus) | `interspecific` + genus/family `family_taxid` |
| Generic protein name ("membrane protein", "nucleoprotein") | `intraspecific` or `interspecific` + `family_taxid`, never unrestricted interspecific |
| Short reference (< 200 aa) | `intraspecific` preferred, since min-length inflation affects interspecific identity scores |

**Identity denominator note:** `_compute_identity` uses `min(len_a, len_b)` as denominator. A 70 aa reference vs. a 700 aa unrelated protein can score 70%+ because the global aligner finds ~50 incidental matches across 70 positions. Intraspecific scope or a tight `family_taxid` restriction eliminates this problem.

## Identity thresholds (flags, not filters)

Identity never excludes a candidate; both cutoffs only set advisory flags (see the filtering pipeline above). Values live in `config.py`:

| Constant | Value | Flag set |
|---|---|---|
| `VARIANTS_UNRELATED_IDENTITY_PERCENT` | `30.0` | `possibly_unrelated` (twilight zone, Rost 1999) |
| `VARIANTS_NEAR_IDENTICAL_PERCENT` | `99.0` | `near_identical` |
| `VARIANTS_MAX_RESULTS` | `100` | UniProt page size cap for the search |

The min-length-inflation caveat still applies: `_compute_identity` uses `min(len_a, len_b)` as denominator, so a short reference vs. a long unrelated protein can score high. Use intraspecific scope or a tight `family_taxid` as the primary fix for contaminated searches; the flags are only an advisory aid.

## UniProt coverage note

UniProt intraspecific coverage can be sparse, and many viral isolates are deposited in NCBI/GenBank but not in UniProt. For richer strain diversity, a pre-built FASTA (e.g. from NCBI) can be supplied directly in the `analyze_conservation` step.

If an intraspecific search returns 0 results due to a mismatched `protein_name` (e.g. the config says `"membrane protein"` but UniProt annotates the protein as `"structural polyprotein"`), conservation will be labelled `conservation_unknown` for all epitopes. Fix by correcting `protein_name` in `project_config.json` and resetting the step.

## References

Data source and method:

- UniProt (variant search): The UniProt Consortium, Nucleic Acids Research, 2023. doi:10.1093/nar/gkac1052
- Biopython `PairwiseAligner` (identity): Cock PJA et al., Bioinformatics, 2009. doi:10.1093/bioinformatics/btp163
- 30% "twilight zone" threshold: Rost B. *Twilight zone of protein sequence alignments.* Protein Engineering. 1999;12(2):85–94.
