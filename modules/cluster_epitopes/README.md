# cluster_epitopes

Groups safe (non-toxic) epitopes into sequence-similarity clusters using pairwise global alignment and NetworkX graph algorithms.

## What it does

For every pair of peptides the step computes a sequence identity score via Biopython global alignment (`match=1, mismatch=0, gap=0`). Identity is normalized by `max(len_a, len_b)`, matching the definition used by the IEDB Cluster Analysis Tool.

Three clustering methods are available:

| Method | How it works | When to use |
|---|---|---|
| `cluster_break` **(default)** | Starts from connected components (same as single_linkage), then iteratively removes the minimum-weight edge from any component whose average intra-cluster similarity falls below the threshold. | Best biological relevance — avoids bridge groupings where peptides with distinct TCR contact patterns end up in the same cluster. |
| `single_linkage` | Connected components of the threshold graph. Peptide A and C join the same cluster if B is similar to both, even if A and C are not directly similar. | IEDB-compatible results; most permissive. |
| `clique` | Each peptide is assigned to the largest maximal clique it belongs to. Every member must be pairwise similar above the threshold. | Strictest — guarantees all members are mutually similar; may over-fragment large groups. |

**Default identity threshold: 80% (0.8).** This value is locked as the project default after first setup.

## Code layout

Split by responsibility (one role per file):

| File | Responsibility |
|---|---|
| `step.py` | `ClusterEpitopesStep` orchestration — `run` / `describe_outputs` |
| `core.py` | Input resolution, pairwise similarity matrix, the three clustering algorithms |
| `io.py` | CLUSTER XLSX writer |
| `prompts.py` | Interactive clustering-parameter selection |

## Input

First file found in this order:

1. `intermediate/{track_id}/toxicity/TOXICITY_SAFE_{track_id}.csv` — preferred (post-toxicity)
2. `intermediate/{track_id}/consensus/CONSENSUS_IMMUNOGENIC_{track_id}.csv` — fallback

## Output

| File | Contents |
|---|---|
| `clusters/CLUSTER_{track_id}.csv` | Input rows plus `cluster_id`, `cluster_method`, `cluster_size` — feeds `select_representatives`. |
| `clusters/CLUSTER_VIEW_{track_id}.csv` | Slim per-step view — `peptide, cluster_id, cluster_size, cluster_method` only. |
| `clusters/CLUSTER_AUDIT_{track_id}.json` | Run metadata: input, threshold, method, n_clusters, n_singletons. |

The cumulative CSV is the input for `select_representatives`.

## Parameters

Asked once per project and saved to `project_config.json`:

- `cluster_threshold` — identity cutoff (0.0–1.0, default 0.8)
- `cluster_method` — `cluster_break` / `single_linkage` / `clique` (default `cluster_break`)

In non-interactive mode (stdin not a tty) both defaults apply automatically.

## Why cluster_break is the default

pMHC-TCR recognition is highly sensitive to anchor position differences. A "bridge" peptide — one that shares one end with group A and the other end with group B — could link two functionally distinct epitopes in single_linkage clustering. `cluster_break` prevents this by enforcing cohesion at the cluster level, not just at the edge level.

## Console output

A Rich table shows:
- Multi-member clusters (2+ members)
- Singletons (1-member clusters)
- Total clusters formed
