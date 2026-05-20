"""
Clustering domain logic for cluster_epitopes: input resolution, the pairwise
identity matrix and the three graph-based clustering algorithms. Pure compute
(NetworkX / NumPy / Biopython aligner) — no XLSX, no prompts.
"""

from pathlib import Path

import networkx as nx
import numpy as np
from Bio.Align import PairwiseAligner
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

from utils.console import console, flush_stdin
from utils.naming import get_step_filename

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
        BarColumn(
            style="yellow",
            complete_style="green",
            finished_style="green",
            pulse_style="bold yellow",
        ),
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

    flush_stdin()

    return similarity_matrix


# ── Graph construction ────────────────────────────────────────────────────────

def _build_similarity_graph(similarity_matrix: np.ndarray, identity_threshold: float) -> nx.Graph:
    """Threshold graph: weighted edge i→j iff similarity[i,j] >= identity_threshold."""
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
    """Single-linkage: connected components of the threshold graph."""
    similarity_graph = _build_similarity_graph(similarity_matrix, identity_threshold)
    cluster_labels = np.empty(similarity_matrix.shape[0], dtype=int)
    for cluster_id, connected_component in enumerate(nx.connected_components(similarity_graph)):
        for node_index in connected_component:
            cluster_labels[node_index] = cluster_id
    return cluster_labels


def _run_clique_clustering(similarity_matrix: np.ndarray, identity_threshold: float) -> np.ndarray:
    """Clique clustering: assign each peptide to its largest maximal clique (greedy, size-desc).
    Overlaps go to the larger/earlier clique; isolated nodes become singletons."""
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
    """Mean of every pairwise similarity inside the cluster (including pairs without an edge)."""
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
    """Cluster-break: connected components, then iteratively drop the weakest edge from any
    component whose mean intra-cluster similarity falls below identity_threshold."""
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


