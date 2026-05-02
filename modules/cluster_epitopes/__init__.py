"""
cluster_epitopes step.

Groups safe epitopes into sequence-similarity clusters using NetworkX.

Three clustering methods are available — all use the same pairwise similarity
matrix (Biopython global alignment, identity score), differing only in how
they decide cluster membership:

  cluster_break  — starts from connected components (like single_linkage), then
                   iteratively removes the minimum-weight edge from any component
                   whose average intra-cluster similarity falls below the threshold.
                   Intermediate: clusters are cohesive but not as strict as cliques.
                   DEFAULT — chosen for biological relevance: avoids "bridge" groupings
                   where peptides with distinct TCR recognition land in the same cluster.

  single_linkage — connected components of the threshold graph (most permissive).
                   Peptide A and C end up in the same cluster if B is similar to both,
                   even if A and C are not directly similar.

  clique         — each peptide is assigned to the largest maximal clique it belongs to
                   (most restrictive). Every member must be pairwise similar to every
                   other member above the threshold.

Input (first found):
    track_dir/toxicity/TOXICITY_SAFE_{track_id}.csv
    track_dir/consensus/CONSENSUS_IMMUNOGENIC_{track_id}.csv

Output:
    track_dir/clusters/CLUSTER_{track_id}.csv
    track_dir/clusters/CLUSTER_AUDIT_{track_id}.json
"""

import json
import datetime
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, BarColumn, MofNCompleteColumn, TextColumn
from rich.table import Table

import config
from modules.base_step import BaseTrackStep
from utils.naming import get_step_filename, COLUMN_PEPTIDE
from utils.project_manager import save_project_config

console = Console(width=120)

_CANONICAL_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# ── Input resolution ──────────────────────────────────────────────────────────

def _find_input_file(track_dir: Path, track_id: str) -> Path:
    """Returns the best available input file for clustering, in priority order."""
    candidate_paths = [
        track_dir / "toxicity"  / get_step_filename("TOXICITY_SAFE",        track_id),
        track_dir / "consensus" / get_step_filename("CONSENSUS_IMMUNOGENIC", track_id),
    ]
    for candidate_path in candidate_paths:
        if candidate_path.exists():
            return candidate_path
    raise FileNotFoundError(
        "No input file found for cluster_epitopes. Tried:\n"
        + "\n".join(f"  {candidate_path}" for candidate_path in candidate_paths)
    )


# ── Pairwise similarity ───────────────────────────────────────────────────────

def _build_similarity_matrix(epitope_sequences: list) -> np.ndarray:
    """
    Computes the pairwise sequence identity matrix for all epitope pairs.

    Identity = alignment score / max(len_a, len_b), using global alignment with
    match=1, mismatch=0, no gap penalties. This matches the IEDB clustering
    tool's percent-identity definition and produces a value in [0.0, 1.0].
    """
    sequence_aligner = PairwiseAligner()
    sequence_aligner.mode          = "global"
    sequence_aligner.match_score   = 1
    sequence_aligner.mismatch_score = 0
    sequence_aligner.open_gap_score = 0
    sequence_aligner.extend_gap_score = 0

    number_of_sequences = len(epitope_sequences)
    similarity_matrix = np.zeros((number_of_sequences, number_of_sequences), dtype=np.float32)
    np.fill_diagonal(similarity_matrix, 1.0)

    total_pairs = number_of_sequences * (number_of_sequences - 1) // 2

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        console=console,
        transient=True,
    ) as progress_bar:
        similarity_task = progress_bar.add_task(
            "  Computing pairwise sequence identity...", total=total_pairs
        )
        for row_index in range(number_of_sequences):
            for col_index in range(row_index + 1, number_of_sequences):
                alignment_score = sequence_aligner.score(
                    epitope_sequences[row_index], epitope_sequences[col_index]
                )
                max_sequence_length = max(
                    len(epitope_sequences[row_index]),
                    len(epitope_sequences[col_index])
                )
                identity = float(alignment_score) / max_sequence_length if max_sequence_length else 0.0
                similarity_matrix[row_index, col_index] = identity
                similarity_matrix[col_index, row_index] = identity
                progress_bar.advance(similarity_task)

    return similarity_matrix


# ── Graph construction ────────────────────────────────────────────────────────

def _build_similarity_graph(similarity_matrix: np.ndarray, identity_threshold: float) -> nx.Graph:
    """
    Builds an undirected weighted graph where an edge exists between peptides
    i and j iff their sequence identity >= identity_threshold.
    """
    number_of_nodes = similarity_matrix.shape[0]
    similarity_graph = nx.Graph()
    similarity_graph.add_nodes_from(range(number_of_nodes))
    for row_index in range(number_of_nodes):
        for col_index in range(row_index + 1, number_of_nodes):
            if similarity_matrix[row_index, col_index] >= identity_threshold:
                similarity_graph.add_edge(
                    row_index, col_index,
                    weight=float(similarity_matrix[row_index, col_index])
                )
    return similarity_graph


# ── Clustering algorithms ─────────────────────────────────────────────────────

def _run_single_linkage(similarity_matrix: np.ndarray, identity_threshold: float) -> np.ndarray:
    """
    Single-linkage clustering: connected components of the threshold graph.
    Most permissive — a "bridge" peptide can chain together peptides that are
    not directly similar to each other.
    """
    similarity_graph = _build_similarity_graph(similarity_matrix, identity_threshold)
    cluster_labels = np.empty(similarity_matrix.shape[0], dtype=int)
    for cluster_id, connected_component in enumerate(nx.connected_components(similarity_graph)):
        for node_index in connected_component:
            cluster_labels[node_index] = cluster_id
    return cluster_labels


def _run_clique_clustering(similarity_matrix: np.ndarray, identity_threshold: float) -> np.ndarray:
    """
    Clique-based clustering: each peptide is assigned to the largest maximal clique
    it belongs to (greedy, sorted by clique size descending). All members within a
    cluster are pairwise similar above the threshold.

    Most restrictive — no "bridge" groupings allowed. Peptides in overlapping cliques
    get assigned to the larger (or earlier) one; remaining singletons get their own cluster.
    """
    similarity_graph = _build_similarity_graph(similarity_matrix, identity_threshold)
    maximal_cliques  = sorted(nx.find_cliques(similarity_graph), key=len, reverse=True)

    number_of_nodes = similarity_matrix.shape[0]
    cluster_labels  = np.full(number_of_nodes, -1, dtype=int)
    next_cluster_id = 0

    for clique_members in maximal_cliques:
        unassigned_members = [node for node in clique_members if cluster_labels[node] == -1]
        if unassigned_members:
            for node in unassigned_members:
                cluster_labels[node] = next_cluster_id
            next_cluster_id += 1

    # Remaining unassigned nodes (isolated — not in any edge) become singletons
    for node_index in range(number_of_nodes):
        if cluster_labels[node_index] == -1:
            cluster_labels[node_index] = next_cluster_id
            next_cluster_id += 1

    return cluster_labels


def _compute_average_intra_cluster_similarity(member_indices: list, similarity_matrix: np.ndarray) -> float:
    """
    Returns the average of ALL pairwise similarities within the cluster,
    including pairs that may not have an edge in the graph.
    """
    number_of_members = len(member_indices)
    if number_of_members <= 1:
        return 1.0
    total_similarity = sum(
        similarity_matrix[member_indices[i], member_indices[j]]
        for i in range(number_of_members)
        for j in range(i + 1, number_of_members)
    )
    number_of_pairs = number_of_members * (number_of_members - 1) / 2
    return float(total_similarity) / number_of_pairs


def _run_cluster_break(similarity_matrix: np.ndarray, identity_threshold: float) -> np.ndarray:
    """
    Cluster-break clustering: starts from connected components (same as single_linkage),
    then iteratively removes the minimum-weight edge from any component whose
    average intra-cluster similarity falls below identity_threshold.

    Chosen as default for biological relevance: prevents "bridge" groupings where
    peptides with distinct TCR recognition are lumped together, while being less
    aggressive than cliques in fragmenting genuinely similar groups.
    """
    similarity_graph = _build_similarity_graph(similarity_matrix, identity_threshold)

    graph_modified = True
    while graph_modified:
        graph_modified = False
        for connected_component in list(nx.connected_components(similarity_graph)):
            if len(connected_component) <= 1:
                continue
            component_member_indices = list(connected_component)
            average_similarity = _compute_average_intra_cluster_similarity(
                component_member_indices, similarity_matrix
            )
            if average_similarity < identity_threshold:
                component_edges = [
                    (node_a, node_b, similarity_graph[node_a][node_b]["weight"])
                    for node_a, node_b in similarity_graph.edges(connected_component)
                ]
                if component_edges:
                    weakest_node_a, weakest_node_b, _ = min(
                        component_edges, key=lambda edge_tuple: edge_tuple[2]
                    )
                    similarity_graph.remove_edge(weakest_node_a, weakest_node_b)
                    graph_modified = True
                    break

    cluster_labels = np.empty(similarity_matrix.shape[0], dtype=int)
    for cluster_id, connected_component in enumerate(nx.connected_components(similarity_graph)):
        for node_index in connected_component:
            cluster_labels[node_index] = cluster_id
    return cluster_labels


# ── Parameter prompt ──────────────────────────────────────────────────────────

def _ask_clustering_params(project_name: str, project_config: dict) -> tuple:
    """
    Returns (identity_threshold, clustering_method) from project_config if already
    saved, otherwise asks the user once and saves both values.
    """
    saved_threshold = project_config.get("cluster_threshold")
    saved_method    = project_config.get("cluster_method")
    if saved_threshold is not None and saved_method is not None:
        console.print(
            f"[dim]  Clustering params (saved): "
            f"method={saved_method}, threshold={saved_threshold}[/dim]"
        )
        return float(saved_threshold), str(saved_method)

    default_threshold = config.CLUSTER_IDENTITY_CUTOFF

    console.print(Panel(
        "[bold]Clustering parameters[/bold]\n\n"
        "[dim]Epitopes with pairwise sequence identity ≥ threshold are grouped "
        "into the same cluster. One representative is selected per cluster.[/dim]\n\n"
        "  [cyan][1][/cyan] cluster_break   — cohesive clusters, no bridge groupings [bold](recommended)[/bold]\n"
        "  [cyan][2][/cyan] single_linkage  — most permissive (IEDB-compatible)\n"
        "  [cyan][3][/cyan] clique          — all members pairwise similar (strictest)",
        box=box.ROUNDED, title="Setup: cluster_epitopes", title_align="left",
    ))

    while True:
        try:
            raw_threshold = input(
                f"  Identity threshold (0.0–1.0) [{default_threshold}]: "
            ).strip()
            identity_threshold = float(raw_threshold) if raw_threshold else default_threshold
            if 0.0 < identity_threshold <= 1.0:
                break
            console.print("[yellow]  Must be > 0 and ≤ 1.[/yellow]")
        except (ValueError, EOFError):
            identity_threshold = default_threshold
            break

    method_choices = {"1": "cluster_break", "2": "single_linkage", "3": "clique"}
    try:
        raw_method = input("  Method [1]: ").strip()
        clustering_method = method_choices.get(raw_method, "cluster_break")
    except EOFError:
        clustering_method = "cluster_break"

    project_config["cluster_threshold"] = identity_threshold
    project_config["cluster_method"]    = clustering_method
    save_project_config(project_name, project_config)

    return identity_threshold, clustering_method


# ── Step ──────────────────────────────────────────────────────────────────────

class ClusterEpitopesStep(BaseTrackStep):
    step_name = "cluster_epitopes"

    def run(self, input_data=None):
        input_file_path = _find_input_file(self.track_dir, self.track_id)
        epitopes_dataframe = pd.read_csv(input_file_path)

        if COLUMN_PEPTIDE not in epitopes_dataframe.columns:
            raise ValueError(
                f"Column '{COLUMN_PEPTIDE}' not found in {input_file_path.name}."
            )

        epitopes_dataframe[COLUMN_PEPTIDE] = epitopes_dataframe[COLUMN_PEPTIDE].astype(str)
        valid_peptide_mask = epitopes_dataframe[COLUMN_PEPTIDE].apply(
            lambda peptide_sequence: (
                len(peptide_sequence) > 0
                and set(peptide_sequence.upper()) <= _CANONICAL_AMINO_ACIDS
            )
        )
        number_of_dropped_rows = int((~valid_peptide_mask).sum())
        if number_of_dropped_rows:
            console.print(
                f"[yellow]  Dropped {number_of_dropped_rows} row(s) with "
                f"empty or non-canonical peptide residues.[/yellow]"
            )
        epitopes_dataframe = epitopes_dataframe[valid_peptide_mask].reset_index(drop=True)

        epitope_sequences = epitopes_dataframe[COLUMN_PEPTIDE].tolist()
        if not epitope_sequences:
            raise ValueError("No valid peptide sequences found in input file.")

        identity_threshold, clustering_method = _ask_clustering_params(
            self.project_name, self.project_config
        )
        console.print(
            f"  [dim]Clustering {len(epitope_sequences)} epitopes — "
            f"method=[bold]{clustering_method}[/bold]  "
            f"threshold={identity_threshold}[/dim]"
        )

        similarity_matrix = _build_similarity_matrix(epitope_sequences)

        if clustering_method == "clique":
            cluster_labels = _run_clique_clustering(similarity_matrix, identity_threshold)
        elif clustering_method == "single_linkage":
            cluster_labels = _run_single_linkage(similarity_matrix, identity_threshold)
        else:
            cluster_labels = _run_cluster_break(similarity_matrix, identity_threshold)

        epitopes_dataframe["cluster_id"]     = cluster_labels
        epitopes_dataframe["cluster_method"] = clustering_method
        epitopes_dataframe["cluster_size"]   = (
            epitopes_dataframe
            .groupby("cluster_id")["cluster_id"]
            .transform("count")
        )

        clusters_output_dir = self.track_dir / "clusters"
        clusters_output_dir.mkdir(parents=True, exist_ok=True)
        output_csv_path = clusters_output_dir / get_step_filename("CLUSTER", self.track_id)
        epitopes_dataframe.to_csv(output_csv_path, index=False)

        number_of_clusters  = int(epitopes_dataframe["cluster_id"].nunique())
        number_of_singletons = int((epitopes_dataframe["cluster_size"] == 1).sum())

        audit_data = {
            "timestamp":          datetime.datetime.now().isoformat(),
            "track_id":           self.track_id,
            "input_file":         str(input_file_path),
            "n_epitopes":         len(epitope_sequences),
            "n_dropped_invalid":  number_of_dropped_rows,
            "identity_threshold": identity_threshold,
            "clustering_method":  clustering_method,
            "n_clusters":         number_of_clusters,
            "n_singletons":       number_of_singletons,
        }
        audit_json_path = clusters_output_dir / get_step_filename(
            "CLUSTER_AUDIT", self.track_id, ext="json"
        )
        audit_json_path.write_text(json.dumps(audit_data, indent=2))

        summary_table = Table(box=box.SIMPLE, show_header=False)
        summary_table.add_column(style="cyan")
        summary_table.add_column(justify="right")
        summary_table.add_row("Input file",          input_file_path.name)
        summary_table.add_row("Epitopes clustered",  str(len(epitope_sequences)))
        summary_table.add_row("Method",              clustering_method)
        summary_table.add_row("Identity threshold",  str(identity_threshold))
        summary_table.add_row("Clusters formed",     str(number_of_clusters))
        summary_table.add_row("Singletons",          str(number_of_singletons))
        console.print(summary_table)

        return {
            "cluster_csv":    str(output_csv_path),
            "n_clusters":     number_of_clusters,
            "n_singletons":   number_of_singletons,
        }
