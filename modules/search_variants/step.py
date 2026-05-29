"""Orchestration for search_variants: searches UniProt for variants of a track's
reference protein, validates them and writes the VARIANTS FASTA + metadata."""

import csv
import datetime
import json

from Bio import SeqIO
from rich.progress import BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn

from config import VARIANTS_NEAR_IDENTICAL_PERCENT, VARIANTS_UNRELATED_IDENTITY_PERCENT
from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.fasta_utils import write_fasta
from utils.naming import get_step_filename
from utils.project_manager import save_project_config

from .core import _build_and_validate, _compute_identity, _search_uniprot_variants
from .io import (
    _load_reference_accessions,
    _load_reference_sequence,
    _reference_has_chain_slice,
    _write_empty_outputs,
)
from .prompts import (
    _ask_redo,
    _ask_scope,
    _ask_variant_view_mode,
    _ensure_tax_id_for_variants,
    _prompt_multi_selection,
)
from .render import _display_genotype_grouped_table, _display_variants_table

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
        "  • [bold]intraspecific[/bold]: variants of the same species (isolates, strains).\n"
        "  • [bold]interspecific[/bold]: homologs in related species (family-level)."
    )
    methodology = (
        "1. UniProt search restricted to the target protein, scoped by tax ID "
        "(intraspecific) or virus family (interspecific).\n"
        "2. Pairwise global alignment (Biopython PairwiseAligner, match=1, gaps=0) "
        "computes percent identity against the reference seed.\n"
        "3. Identity is computed for every candidate and shown, but it never excludes: "
        "below 30% is flagged 'possibly unrelated' and at or above 99% is flagged "
        "'near-identical'. You decide what to keep during selection, and exact "
        "byte-identical duplicates are collapsed to one row.\n"
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
        "Input is implicit: it uses the reference FASTA produced by [bold]fetch_sequences[/bold] "
        "as the query and the same protein name you defined for this organism/protein pair. "
        "You will be asked once per pair:\n"
        "  • [bold]Scope[/bold]: intraspecific (default) or interspecific.\n"
        "  • [bold]Host filter[/bold] (optional): restrict to a host organism.\n"
        "  • [bold]Family taxid[/bold] (interspecific only): the higher-level taxon to search within."
    )
    outputs_overview = (
        "[bold]VARIANTS_{track_id}.fasta[/bold]       multi-FASTA of selected variants (input for analyze_conservation).\n"
        "[bold]VARIANTS_VIEW_{track_id}.csv[/bold]    slim per-step view (accession, organism, length, identity).\n"
        "[bold]VARIANTS_AUDIT_{track_id}.json[/bold]  scope, filters, totals, rejected entries (full audit)."
    )
    tips = [
        "Intraspecific is usually what you want for vaccine work; it captures isolate-level diversity.",
        "If UniProt returns few intraspecific variants (common for emerging viruses), consider switching "
        "to interspecific or supplying a curated FASTA directly in analyze_conservation.",
        "Identity is computed end-to-end, so partial sequences will score lower than full-length variants.",
        "The variants FASTA is cached: rerunning the step asks before overwriting existing results.",
    ]

    def clean_outputs(self) -> None:
        """On a rerun, delete the cached variant FASTA/audit/view so the search
        starts fresh — prevents a previous empty (0-variant) result from being
        kept and re-saved as zeros."""
        variants_dir = self.track_dir / "variants"
        for stale_path in (
            variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta"),
            variants_dir / get_step_filename("VARIANTS_AUDIT", self.track_id, ext="json"),
            variants_dir / get_step_filename("VARIANTS_VIEW", self.track_id),
        ):
            if stale_path.exists():
                stale_path.unlink()

    def describe_outputs(self) -> dict:
        variants_dir = self.track_dir / "variants"
        return {
            variants_dir / get_step_filename("VARIANTS_VIEW", self.track_id):
                "Slim per-step view: accession + organism + protein + length + identity_to_seed for every selected variant.",
            variants_dir / get_step_filename("VARIANTS", self.track_id, ext="fasta"):
                "Multi-FASTA of variant sequences, the input for analyze_conservation.",
            variants_dir / get_step_filename("VARIANTS_AUDIT", self.track_id, ext="json"):
                "Run audit: scope (intraspecific/interspecific), filters, totals, rejected entries.",
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

        # ── Guard local-FASTA references (no tax_id / possible polyprotein) ────
        tax_id = _ensure_tax_id_for_variants(
            self.track_id, track_config,
            self.project_name, self.project_config, len(ref_seq),
        )

        # ── Scope + host filter + family restriction ──────────────────────────
        scope, host_filter, family_taxid = _ask_scope(
            self.track_id, track_config,
            self.project_name, self.project_config,
            is_rerun=self.is_rerun,
        )

        # ── Polyprotein slicing gate ───────────────────────────────────────────
        # Slice variant candidates down to the mature chain only when the track
        # needs it: either the reference was itself sliced out of a polyprotein, or
        # the user supplied a local FASTA (a clean mature protein whose UniProt
        # variants may still be full polyproteins). Clean single-protein UniProt
        # tracks skip this — no ft_chain request, no per-candidate matching.
        input_source          = track_config.get("input_source", "uniprot")
        should_slice_variants = (
            _reference_has_chain_slice(self.input_dir, self.track_id)
            or input_source == "local"
        )

        # ── UniProt search ────────────────────────────────────────────────────
        candidates = _search_uniprot_variants(
            protein_name, tax_id, scope, host_filter, family_taxid,
            should_slice_variants=should_slice_variants,
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

        # ── Compute identity, collapse exact duplicates, flag (no exclusion) ───
        # Identity is computed for every candidate and shown in the table, but it
        # never excludes — <30% and ≥99% only set advisory flags so the user decides
        # during selection (intraspecific conservation *wants* the ~99% isolates).
        # The reference's own accessions are already gone; here we additionally
        # collapse byte-identical re-deposits (same sliced sequence = same crc64)
        # down to one row, preferring the reviewed (Swiss-Prot) copy.
        seen_sequence_keys:   dict[str, int] = {}
        unique_candidates:    list[dict]     = []
        duplicate_accessions: list[str]      = []

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
                identity_progress_bar.advance(identity_task_id)
                if not candidate["sequence"]:
                    continue

                identity_percent = _compute_identity(ref_seq, candidate["sequence"])
                candidate["identity"] = round(identity_percent, 2)

                candidate_flags: list[str] = []
                if identity_percent < VARIANTS_UNRELATED_IDENTITY_PERCENT:
                    candidate_flags.append("possibly_unrelated")
                if identity_percent >= VARIANTS_NEAR_IDENTICAL_PERCENT:
                    candidate_flags.append("near_identical")
                # Slicing was expected for this track but no Chain matched this candidate,
                # yet it is far longer than the reference — likely an uncut polyprotein whose
                # naming we could not match. Surface it instead of silently feeding the whole
                # polyprotein into conservation; the user can rename the protein or drop it.
                is_unsliced_oversized = (
                    should_slice_variants
                    and "chain_slice" not in candidate
                    and len(candidate["sequence"]) > 1.5 * len(ref_seq)
                )
                if is_unsliced_oversized:
                    candidate_flags.append("uncut_polyprotein")
                candidate["flags"] = candidate_flags

                sequence_key = candidate["sequence"]
                if sequence_key not in seen_sequence_keys:
                    seen_sequence_keys[sequence_key] = len(unique_candidates)
                    unique_candidates.append(candidate)
                else:
                    duplicate_accessions.append(candidate["accession"])
                    kept_index = seen_sequence_keys[sequence_key]
                    if candidate["reviewed"] and not unique_candidates[kept_index]["reviewed"]:
                        unique_candidates[kept_index] = candidate

        flush_stdin()

        if duplicate_accessions:
            sample = ", ".join(duplicate_accessions[:5])
            suffix = "..." if len(duplicate_accessions) > 5 else ""
            console.print(
                f"[dim]→ Collapsed {len(duplicate_accessions)} exact duplicate(s) "
                f"(identical sequence): {sample}{suffix}[/dim]"
            )

        candidates = sorted(unique_candidates, key=lambda c: c["identity"] or 0.0, reverse=True)

        flagged_count = sum(1 for c in candidates if c.get("flags"))
        if flagged_count:
            console.print(
                f"[dim]→ {flagged_count} candidate(s) flagged "
                f"(<{VARIANTS_UNRELATED_IDENTITY_PERCENT:.0f}% = possibly unrelated, "
                f"≥{VARIANTS_NEAR_IDENTICAL_PERCENT:.0f}% = near-identical), kept; you decide.[/dim]"
            )

        if not candidates:
            console.print("[yellow]⚠ No variants with a usable sequence were returned.[/yellow]")
            _write_empty_outputs(fasta_path, audit_path, self.track_id, scope, host_filter,
                                 family_taxid, ref_accessions,
                                 note="No variants with a usable sequence after dedup.")
            return {"variants_fasta": str(fasta_path), "variants_audit": str(audit_path),
                    "total_variants": 0, "cached": False}

        # ── Display table & select ────────────────────────────────────────────
        # Interspecific can span many genotypes; offer a per-genotype view so other
        # types (e.g. HPV18/31/33…) aren't buried under the reference's own isolates.
        if scope == "interspecific" and _ask_variant_view_mode() == "grouped":
            genotype_reps = _display_genotype_grouped_table(candidates)
            selected = _prompt_multi_selection(genotype_reps, unit_label="genotype")
        else:
            _display_variants_table(candidates)
            selected = _prompt_multi_selection(candidates)

        # ── Validate + build SeqRecords ───────────────────────────────────────
        valid_records, rejected_log = _build_and_validate(selected)

        if rejected_log:
            console.print(f"[yellow]⚠ {len(rejected_log)} variant(s) failed validation:[/yellow]")
            for entry in rejected_log:
                console.print(f"  [yellow]{entry['id']}[/yellow]: {entry['reason']}")

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
            "reference_id":             reference.id,
            "reference_excluded":       list(ref_accessions),
            "exact_duplicates_excluded": duplicate_accessions,
            "sliced_to_chain":          should_slice_variants,
            "total_raw_results":        before_filter,
            "total_after_dedup":        len(candidates),
            "total_selected":           len(selected),
            "total_valid":              len(valid_records),
            "total_rejected":           len(rejected_log),
            "rejected":                 rejected_log,
            "selected_variants": [
                {
                    "accession": selected_candidate["accession"],
                    "organism":  selected_candidate["organism"],
                    "identity":  selected_candidate["identity"],
                    "length":    selected_candidate["length"],
                    "reviewed":  selected_candidate["reviewed"],
                    "flags":     selected_candidate.get("flags", []),
                }
                for selected_candidate in selected
            ],
        }
        with open(audit_path, "w", encoding="utf-8") as fh:
            json.dump(audit, fh, indent=2, ensure_ascii=False)

        view_path = variants_dir / get_step_filename("VARIANTS_VIEW", self.track_id)
        with open(view_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["accession", "organism", "protein", "length", "identity_to_seed", "flags"])
            for selected_candidate in selected:
                writer.writerow([
                    selected_candidate.get("accession", ""),
                    selected_candidate.get("organism", ""),
                    protein_name,
                    selected_candidate.get("length", ""),
                    selected_candidate.get("identity", ""),
                    "; ".join(selected_candidate.get("flags", [])),
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
