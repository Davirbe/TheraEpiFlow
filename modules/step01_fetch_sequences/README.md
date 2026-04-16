# Step 01 — Fetch Sequences

Downloads protein sequences from NCBI GenBank for each project track, or loads sequences from a local FASTA file.

---

## What this step does

For each track (organism + protein combination), step01:

1. Reads the search configuration from `project_config.json` (defined at project setup)
2. Searches NCBI GenBank protein database — RefSeq entries first
3. Displays a table of candidate sequences for the user to review
4. The user selects one or more sequences
5. Saves the selected sequences and their metadata

---

## Two-phase GenBank workflow

**Phase 1 — Discovery (no download):**
Searches NCBI and retrieves accession IDs. RefSeq entries (curated, reviewed) are prioritized over general GenBank entries. This phase is fast and downloads no sequences.

**Phase 2 — Display and selection:**
Fetches the first 50 records (full GenBank format) to populate the table. The table shows key metadata so the user can make an informed selection. Only the selected sequences are kept — the rest are discarded.

---

## Reference sequence suggestion

A `★` marker highlights the suggested reference sequence:

- **If RefSeq entries exist:** suggests the longest non-polyprotein RefSeq entry
- **If no RefSeq:** suggests the entry at the most common sequence length (canonical reference length)

The suggestion is a hint, not a decision. The user always makes the final choice.

---

## Polyprotein handling

Some viruses (e.g. Zika, Dengue) store multiple proteins in a single polyprotein GenBank record. If the user selects "include polyproteins" at project setup, this step:

1. Detects polyprotein records (description contains "polyprotein" AND length > 1000 aa)
2. Extracts the target protein using `mat_peptide` feature annotations
3. Presents the extracted sequence as a selectable option in the table

---

## Input

Reads from `project_config.json`, specifically:

```json
"tracks": {
  "HPV16_E6": {
    "organism_name": "Human papillomavirus 16",
    "protein_name": "E6",
    "input_source": "genbank",
    "search_include_polyproteins": false,
    "local_file_path": null
  }
}
```

---

## Output

Saved to `data/input/{track_id}/`:

| File | Contents |
|---|---|
| `sequences_{track_id}.fasta` | Selected sequences in FASTA format |
| `sequence_registry_{track_id}.json` | Metadata for each sequence (accession, organism, strain, isolate, host, location, date) |

The `sequence_registry` file is read by **Step 08** to exclude original reference sequences from variant searches.

---

## Table columns

| Column | Description |
|---|---|
| `#` | Row number for selection |
| `★` | Suggested reference sequence |
| Accession | GenBank accession ID |
| aa | Sequence length in amino acids |
| Strain / Isolate | Viral strain or clinical isolate identifier |
| Location | Geographic collection location |
| Date | Sample collection date |
| Source | RefSeq (curated) or GenBank |

---

## Selection formats

| Input | Effect |
|---|---|
| `1` | Select row 1 only |
| `1,3,5` | Select rows 1, 3 and 5 |
| `1-4` | Select rows 1 through 4 |
| `all` | Select all displayed rows |
| *(Enter)* | Accept the suggested sequence (★) |
