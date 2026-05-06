# search_variants

Retrieves protein sequence variants from UniProt REST API and writes a permanent multi-FASTA that `analyze_conservation` uses to measure epitope conservation across natural sequence diversity.

## What it does

Queries UniProt for proteins related to each track's reference sequence. Two search scopes:

| Scope | Query strategy | Typical use |
|---|---|---|
| Intraspecific | `taxonomy_id:{tax_id} AND protein_name:"{name}"` | Strains / isolates of the same species — e.g. different HPV16 submissions, SARS-CoV-2 variants |
| Interspecific | `protein_name:"{name}"` ± optional host filter | Same protein across related species — e.g. E5 from HPV16, HPV18, HPV31 infecting Homo sapiens |

The scope and optional host filter are asked interactively once per track and saved to `project_config.json["tracks"][track_id]` for reproducibility.

## Filtering pipeline

After the raw UniProt results come back:

1. **Reference exclusion** — removes any accession that was fetched by `fetch_sequences` (reads the registry JSON).
2. **Near-identical filter** — removes full-length candidates (≥ 80% of reference length) with ≥ 99% identity vs. reference. Fragments are always kept.
3. **Identity sort** — remaining candidates are sorted descending by identity to the reference.
4. **Interactive selection** — table displayed; user selects by index range (e.g. `1,3,5-8`, `all`, `none`). Non-interactive mode selects all.
5. **Lenient validation** — only rejects empty sequences or those composed entirely of ambiguous residues (X/B/Z/U/O). Short sequences (< 50 aa) are kept with a warning.

## Cache behaviour

The FASTA is **permanent**. When the file already exists the step asks whether to keep it or redo the search. In non-interactive mode the existing file is always kept and the step completes instantly. This avoids re-running expensive API searches on every pipeline re-run.

## Input

- `input/{track_id}/SEQUENCES_{track_id}.fasta` — reference sequence (from `fetch_sequences`)
- `input/{track_id}/REGISTRY_{track_id}.json` — reference accession list (for exclusion)

## Output

| File | Contents |
|---|---|
| `variants/VARIANTS_{track_id}.fasta` | Multi-FASTA of selected variant sequences |
| `variants/VARIANTS_AUDIT_{track_id}.json` | Scope, host filter, counts, selected accessions |

An empty FASTA is written (with a note in the audit) when no variants are found after filtering, so downstream steps can always expect the file to exist.

## UniProt coverage note

UniProt intraspecific coverage can be sparse — many viral isolates are deposited in NCBI/GenBank but not in UniProt. For richer strain diversity, a pre-built FASTA (e.g. from NCBI) can be supplied directly in the `analyze_conservation` step.
