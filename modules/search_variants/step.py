"""Orchestration for search_variants: searches UniProt for variants of a track's
reference protein, validates them and writes the VARIANTS FASTA + metadata."""

import csv
import datetime
import json
import time
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.fasta_utils import write_fasta
from utils.naming import get_step_filename
from utils.project_manager import save_project_config

from .core import _build_and_validate, _compute_identity, _search_uniprot_variants
from .io import _load_reference_accessions, _load_reference_sequence, _write_empty_outputs
from .prompts import _ask_redo, _ask_scope, _prompt_multi_selection
from .render import _display_variants_table

# ── Step class ────────────────────────────────────────────────────────────────

class SearchVariantsStep(BaseTrackStep):
    step_name   = "search_variants"
    description = (
        "Fetches variant sequences for the track's protein from UniProt "
        "(scope: intraspecific or interspecific) so the conservation step can "
        "compare each ★ epitope against real-world sequence variation."
    )
    long_description = (
        "Searches UniProt for additional sequences of the same protein, "
        "computes pairwise identity vs. the reference (the seed picked in "
        "fetch_sequences), filters near-identical (≥99%) and clearly "
        "unrelated (<30%) hits, and lets you pick which variants to keep. "
        "The resulting multi-FASTA is consumed by [bold]analyze_conservation[/bold] "
        "to flag epitopes that are conserved across real-world strains.\n\n"
        "Two scopes are available:\n"
        "  • [bold]intraspecific[/bold] — variants of the same species (isolates, strains).\n"
        "  • [bold]interspecific[/bold] — homologs in related species (family-level)."
    )
    methodology = (
        "1. UniProt search restricted to the target protein, scoped by tax ID "
        "(intraspecific) or virus family (interspecific).\n"
        "2. Pairwise global alignment (Biopython PairwiseAligner, match=1, gaps=0) "
        "computes percent identity against the reference seed.\n"
        "3. Cutoffs: <30% identity → likely unrelated (excluded); ≥99% identity AND "
        "≥80% length → near-duplicate of the reference (excluded). The rest passes.\n"
        "4. Reviewed (Swiss-Prot) entries are flagged so you can prioritise curated variants.\n"
        "5. After the interactive selection, sequences are revalidated (no ambiguous "
        "residues, minimum length) before being written to the variants FASTA."
    )
    references = [
        {
            'authors': 'The UniProt Consortium',
            'title':   'UniProt: the Universal Protein Knowledgebase in 2023',
            'journal': 'Nucleic Acids Research',
            'year':    2023,
            'doi':     '10.1093/nar/gkac1052',
        },
        {
            'authors': 'Cock PJ, Antao T, Chang JT, et al.',
            'title':   'Biopython: freely available Python tools for computational molecular biology and bioinformatics',
            'journal': 'Bioinformatics',
            'year':    2009,
            'doi':     '10.1093/bioinformatics/btp163',
        },
    ]
    data_format = (
        "Input is implicit — uses the reference FASTA produced by [bold]fetch_sequences[/bold] "
        "as the query and the same protein name you defined for this organism/protein pair. "
        "You will be asked once per pair:\n"
        "  • [bold]Scope[/bold]: intraspecific (default) or interspecific.\n"
        "  • [bold]Host filter[/bold] (optional): restrict to a host organism.\n"
        "  • [bold]Family taxid[/bold] (interspecific only): the higher-level taxon to search within."
    )
    outputs_overview = (
        "[bold]VARIANTS_{track_id}.fasta[/bold]       — multi-FASTA of selected variants (input for analyze_conservation).\n"
        "[bold]VARIANTS_VIEW_{track_id}.csv[/bold]    — slim per-step view (accession, organism, length, identity).\n"
        "[bold]VARIANTS_AUDIT_{track_id}.json[/bold]  — scope, filters, totals, rejected entries (full audit)."
    )
    tips = [
        "Intraspecific is usually what you want for vaccine work — captures isolate-level diversity.",
        "If UniProt returns few intraspecific variants (common for emerging viruses), consider switching "
        "to interspecific or supplying a curated FASTA directly in analyze_conservation.",
        "Identity is computed end-to-end — partial sequences will score lower than full-length variants.",
        "The variants FASTA is cached: rerunning the step asks before overwriting existing results.",
    ]

    def describe_outputs(self) -> dict:
        variants_dir = self.track_dir / "variants"
        return {
            variants_dir / get_step_filename("VARIANTS_VIEW", self.track_id):
                "Slim per-step view — accession + organism + protein + length + identity_to_seed for every selected variant.",
            variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta"):
                "Multi-FASTA of variant sequences — input for analyze_conservation.",
            variants_dir / get_step_filename("VARIANTS_AUDIT", self.track_id, ext="json"):
                "Run audit — scope (intraspecific/interspecific), filters, totals, rejected entries.",
        }

    def run(self, input_data=None):
        track_config = self.project_config["tracks"][self.track_id]
        protein_name = track_config.get("protein_name", "")
        tax_id       = track_config.get("tax_id")

        variants_dir = self.track_dir / "variants"
        variants_dir.mkdir(parents=True, exist_ok=True)

        fasta_path = variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta")
        audit_path = variants_dir / get_step_filename("VARIANTS_AUDIT", self.track_id, ext="json")

        # ── Cache check: existing FASTA → ask to keep or redo ─────────────────
        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            existing = list(SeqIO.parse(str(fasta_path), "fasta"))
            if not _ask_redo(len(existing)):
                console.print(f"[dim]→ Keeping existing FASTA ({len(existing)} sequences).[/dim]")
                return {
                    "variants_fasta": str(fasta_path),
                    "variants_audit": str(audit_path),
                    "total_variants": len(existing),
                    "cached": True,
                }
            # Clear cache and saved scope so user picks fresh
            fasta_path.unlink()
            if audit_path.exists():
                audit_path.unlink()
            track_config.pop("variants_scope", None)
            track_config.pop("variants_host_filter", None)
            track_config.pop("variants_family_taxid", None)
            save_project_config(self.project_name, self.project_config)
            console.print("[dim]→ Cache cleared. Restarting search.[/dim]")

        # ── Load reference ────────────────────────────────────────────────────
        reference      = _load_reference_sequence(self.input_dir, self.track_id)
        ref_accessions = _load_reference_accessions(self.input_dir, self.track_id)
        ref_seq        = str(reference.seq)

        console.print(f"\n[bold cyan]━━━ Track: {self.track_id} ━━━[/bold cyan]")
        console.print(f"[dim]Reference : {reference.id}  ({len(ref_seq)} aa)[/dim]")
        console.print(f"[dim]Protein   : {protein_name}[/dim]")
        console.print(f"[dim]Tax ID    : {tax_id}[/dim]")

        # ── Scope + host filter + family restriction ──────────────────────────
        scope, host_filter, family_taxid = _ask_scope(
            self.track_id, track_config,
            self.project_name, self.project_config,
        )

        # ── UniProt search ────────────────────────────────────────────────────
        candidates = _search_uniprot_variants(
            protein_name, tax_id, scope, host_filter, family_taxid
        )

        # ── Filter: exclude reference accessions ──────────────────────────────
        before_filter = len(candidates)
        candidates    = [c for c in candidates if c["accession"] not in ref_accessions]
        excluded_ref  = before_filter - len(candidates)
        if excluded_ref:
            console.print(f"[dim]→ Excluded {excluded_ref} reference accession(s).[/dim]")

        if not candidates:
            console.print("[yellow]⚠ No variants found after excluding references.[/yellow]")
            console.print("[dim]  UniProt has limited intraspecific coverage (few isolates per species).[/dim]")
            console.print("[dim]  For richer strain diversity, you can provide a pre-built FASTA[/dim]")
            console.print("[dim]  (e.g. from NCBI) in the analyze_conservation step.[/dim]")
            _write_empty_outputs(fasta_path, audit_path, self.track_id, scope, host_filter,
                                 family_taxid, ref_accessions,
                                 note="No variants found after excluding references.")
            return {"variants_fasta": str(fasta_path), "variants_audit": str(audit_path),
                    "total_variants": 0, "cached": False}

        # ── Compute identity & filter ─────────────────────────────────────────
        near_identical: list[str]  = []
        low_identity:   list[str]  = []
        passing:        list[dict] = []

        _MIN_IDENTITY_THRESHOLD = 30.0  # % — excludes unrelated proteins from other virus families

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
        ) as identity_progress_bar:
            identity_task_id = identity_progress_bar.add_task(
                "  Computing pairwise identity vs. reference",
                total=len(candidates),
            )
            for candidate in candidates:
                if not candidate["sequence"]:
                    identity_progress_bar.advance(identity_task_id)
                    continue
                identity_percent = _compute_identity(ref_seq, candidate["sequence"])
                candidate["identity"] = round(identity_percent, 2)
                if identity_percent < _MIN_IDENTITY_THRESHOLD:
                    low_identity.append(candidate["accession"])
                else:
                    is_full_length_candidate = len(candidate["sequence"]) >= 0.8 * len(ref_seq)
                    if identity_percent >= 99.0 and is_full_length_candidate:
                        near_identical.append(candidate["accession"])
                    else:
                        passing.append(candidate)
                identity_progress_bar.advance(identity_task_id)

        flush_stdin()

        if low_identity:
            sample = ", ".join(low_identity[:5])
            suffix = "..." if len(low_identity) > 5 else ""
            console.print(
                f"[dim]→ Excluded {len(low_identity)} candidate(s) with identity "
                f"< {_MIN_IDENTITY_THRESHOLD:.0f}% (unrelated proteins): {sample}{suffix}[/dim]"
            )

        if near_identical:
            sample = ", ".join(near_identical[:5])
            suffix = "..." if len(near_identical) > 5 else ""
            console.print(
                f"[dim]→ Excluded {len(near_identical)} near-identical (≥99%): {sample}{suffix}[/dim]"
            )

        candidates = sorted(passing, key=lambda c: c["identity"] or 0.0, reverse=True)

        if not candidates:
            console.print("[yellow]⚠ No variants remain after identity filtering.[/yellow]")
            console.print("[dim]  All candidates were near-identical to the reference or below the 30% identity threshold.[/dim]")
            _write_empty_outputs(fasta_path, audit_path, self.track_id, scope, host_filter,
                                 family_taxid, ref_accessions,
                                 note="No variants remain after identity filtering (min 30%).")
            return {"variants_fasta": str(fasta_path), "variants_audit": str(audit_path),
                    "total_variants": 0, "cached": False}

        # ── Display table & select ────────────────────────────────────────────
        _display_variants_table(candidates)
        selected = _prompt_multi_selection(candidates)

        # ── Validate + build SeqRecords ───────────────────────────────────────
        valid_records, rejected_log = _build_and_validate(selected)

        if rejected_log:
            console.print(f"[yellow]⚠ {len(rejected_log)} variant(s) failed validation:[/yellow]")
            for entry in rejected_log:
                console.print(f"  [yellow]{entry['id']}[/yellow] — {entry['reason']}")

        # ── Write FASTA (permanent) ───────────────────────────────────────────
        write_fasta(valid_records, fasta_path)

        # ── Write audit ───────────────────────────────────────────────────────
        audit = {
            "track_id":                self.track_id,
            "generated_at":            datetime.datetime.now().isoformat(),
            "scope":                   scope,
            "family_taxid":            family_taxid,
            "host_filter":             host_filter,
            "protein_name":            protein_name,
            "tax_id":                  tax_id,
            "reference_id":            reference.id,
            "reference_excluded":      list(ref_accessions),
            "near_identical_excluded": near_identical,
            "total_raw_results":       before_filter,
            "total_after_filter":      len(candidates),
            "total_selected":          len(selected),
            "total_valid":             len(valid_records),
            "total_rejected":          len(rejected_log),
            "rejected":                rejected_log,
            "selected_variants": [
                {
                    "accession": c["accession"],
                    "organism":  c["organism"],
                    "identity":  c["identity"],
                    "length":    c["length"],
                    "reviewed":  c["reviewed"],
                }
                for c in selected
            ],
        }
        with open(audit_path, "w", encoding="utf-8") as fh:
            json.dump(audit, fh, indent=2, ensure_ascii=False)

        view_path = variants_dir / get_step_filename("VARIANTS_VIEW", self.track_id)
        with open(view_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["accession", "organism", "protein", "length", "identity_to_seed"])
            for c in selected:
                writer.writerow([
                    c.get("accession", ""),
                    c.get("organism", ""),
                    protein_name,
                    c.get("length", ""),
                    c.get("identity", ""),
                ])

        console.print(f"\n[bold green]✓ {len(valid_records)} variant(s) saved.[/bold green]")
        console.print(f"  FASTA → {fasta_path}")
        console.print(f"  Audit → {audit_path}")

        return {
            "variants_fasta": str(fasta_path),
            "variants_audit": str(audit_path),
            "total_variants": len(valid_records),
            "cached": False,
        }
