"""Orchestration for cluster_epitopes: builds the similarity matrix, clusters
the epitopes and writes the CLUSTER CSV + XLSX."""

import datetime
import json

import pandas as pd
from rich import box
from rich.table import Table

from modules.base_step import BaseTrackStep
from utils.console import console
from utils.csv_write import write_user_facing_csv
from utils.naming import COLUMN_PEPTIDE, get_step_filename

from .core import (
    _CANONICAL_AMINO_ACIDS,
    _build_similarity_matrix,
    _find_input_file,
    _run_clique_clustering,
    _run_cluster_break,
    _run_single_linkage,
)
from .io import _write_cluster_xlsx
from .prompts import _ask_clustering_params

# ── Step ──────────────────────────────────────────────────────────────────────

class ClusterEpitopesStep(BaseTrackStep):
    step_name   = "cluster_epitopes"
    description = (
        "Builds a pairwise similarity graph over the remaining peptides "
        "(NetworkX) and partitions them into clusters at the configured "
        "identity cutoff so each cluster can be represented by a single ★ "
        "epitope downstream."
    )
    long_description = (
        "Two peptides that share ≥ 80% identity bind very similar TCRs and "
        "are immunologically redundant. This step groups them so the next "
        "step (select_representatives) can pick one ★ per cluster instead "
        "of carrying duplicates through the rest of the pipeline.\n\n"
        "Three clustering algorithms are available — all build the same "
        "Biopython pairwise-alignment similarity graph, but differ in how "
        "they decide who belongs to which cluster."
    )
    methodology = (
        "Similarity matrix: Biopython PairwiseAligner global alignment, "
        "identity score (match=1, mismatch=0, gap=0).\n\n"
        "Then three method choices:\n"
        "  • [bold]cluster_break[/bold] (DEFAULT) — connected components of the "
        "thresholded graph, then iteratively prune the weakest edge from any "
        "component whose mean intra-cluster similarity drops below the "
        "cutoff. Biologically motivated: avoids 'bridge' groupings.\n"
        "  • [bold]single_linkage[/bold] — pure connected components. Most "
        "permissive (A and C cluster together if both link to B).\n"
        "  • [bold]clique[/bold] — every peptide assigned to its largest "
        "maximal clique. Most restrictive (members are pairwise similar)."
    )
    references = [
        {
            'authors': 'Cock PJ, Antao T, Chang JT, et al.',
            'title':   'Biopython: freely available Python tools for computational molecular biology and bioinformatics',
            'journal': 'Bioinformatics',
            'year':    2009,
            'doi':     '10.1093/bioinformatics/btp163',
        },
        {
            'authors': 'Bron C, Kerbosch J.',
            'title':   'Algorithm 457: Finding all cliques of an undirected graph',
            'journal': 'Communications of the ACM',
            'year':    1973,
            'doi':     '10.1145/362342.362367',
        },
    ]
    data_format = (
        "Input is automatic — uses TOXICITY_SAFE_{track_id}.csv from the "
        "previous step. You will be asked once for:\n"
        "  • [bold]Identity cutoff[/bold] (0.0–1.0, default 0.8 = 80%).\n"
        "  • [bold]Method[/bold] — cluster_break / single_linkage / clique."
    )
    outputs_overview = (
        "[bold]CLUSTER_VIEW_{track_id}.csv[/bold]   — slim per-step view (peptide + cluster_id + cluster_size + method).\n"
        "[bold]CLUSTER_{track_id}.csv[/bold]        — full per-peptide row with all upstream columns + cluster info.\n"
        "[bold]CLUSTER_{track_id}.xlsx[/bold]       — same data, sorted by cluster_id with alternating band colours per cluster (orange = cluster_id, pink = cluster_size).\n"
        "[bold]CLUSTER_AUDIT_{track_id}.json[/bold] — threshold, method, n_clusters, n_singletons."
    )
    tips = [
        "80% identity is the default in IEDB tooling and a sensible biological cutoff for MHC-I.",
        "single_linkage is permissive — pick it only if you want to maximise diversity capture.",
        "clique is strict and gives smaller clusters — useful when downstream analysis tolerates redundancy.",
        "Singletons (size=1) are kept and treated as their own clusters — they always become ★ representatives.",
    ]

    def describe_outputs(self) -> dict:
        clusters_dir = self.track_dir / "clusters"
        return {
            clusters_dir / get_step_filename("CLUSTER_VIEW", self.track_id):
                "Slim per-step view — peptide + cluster_id + cluster_size + cluster_method only.",
            clusters_dir / get_step_filename("CLUSTER", self.track_id):
                "Per-peptide clustering result — adds cluster_id, cluster_size, cluster_method columns.",
            clusters_dir / get_step_filename("CLUSTER", self.track_id, ext="xlsx"):
                "Same as the CSV, sorted by cluster_id with alternating band colours so groupings are visible.",
            clusters_dir / get_step_filename("CLUSTER_AUDIT", self.track_id, ext="json"):
                "Run audit — threshold, method, n_clusters, n_singletons.",
        }

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
            raise ValueError(
                f"No valid peptide sequences in {input_file_path.name} — nothing to cluster. "
                "'consensus_filter' kept no binders; re-run it (try a looser threshold) "
                "before 'cluster_epitopes'."
            )

        identity_threshold, clustering_method = _ask_clustering_params(
            self.project_name, self.project_config, is_rerun=self.is_rerun
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
        output_csv_path  = clusters_output_dir / get_step_filename("CLUSTER", self.track_id)
        output_xlsx_path = clusters_output_dir / get_step_filename("CLUSTER", self.track_id, ext="xlsx")
        write_user_facing_csv(epitopes_dataframe, output_csv_path)
        _write_cluster_xlsx(epitopes_dataframe, output_xlsx_path)

        view_columns = [COLUMN_PEPTIDE, "cluster_id", "cluster_size", "cluster_method"]
        view_csv_path = clusters_output_dir / get_step_filename("CLUSTER_VIEW", self.track_id)
        write_user_facing_csv(epitopes_dataframe[view_columns], view_csv_path)

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
        console.print(
            f"\n[bold green]RESULT: {number_of_clusters} clusters from "
            f"{len(epitope_sequences)} epitopes "
            f"(method={clustering_method}, threshold={identity_threshold}).[/bold green]\n"
        )
        console.print(f"  CSV → {output_csv_path}")

        return {
            "cluster_csv":    str(output_csv_path),
            "n_clusters":     number_of_clusters,
            "n_singletons":   number_of_singletons,
        }
