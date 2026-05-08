"""
Step search_variants — retrieves protein sequence variants from UniProt.

Two search scopes:
  intraspecific  — UniProt restricted to the same taxonomy_id (isolates/strains of the
                   same species, e.g. different HPV16 submissions, SARS-CoV-2 variants)
  interspecific  — UniProt search by protein name across all organisms, with an optional
                   host filter (e.g. E5 from HPV16, HPV18, HPV31... infecting Homo sapiens)

Note: UniProt has limited intraspecific coverage. For richer strain diversity,
supplement the output with a pre-built FASTA in the analyze_conservation step.

The resulting multi-FASTA is permanent. When it already exists the user is prompted
to keep it or regenerate it. Non-interactive mode always keeps the existing file.
"""

import datetime
import json
import sys
import time
from pathlib import Path

import requests
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from utils.console import console
from rich.table import Table
from rich import box

from modules.base_step import BaseTrackStep
from utils.fasta_utils import write_fasta
from utils.naming import get_step_filename
from utils.project_manager import save_project_config

UNIPROT_SEARCH_URL   = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy/{tax_id}"
_MAX_VARIANTS        = 100
_REQUEST_TIMEOUT     = 20
_AMBIGUOUS_AAS       = set("XBZUO")
_INFORMATIVE_RANKS   = frozenset({
    "genus", "subgenus", "subfamily", "family", "superfamily",
    "suborder", "order",
})


# ── HTTP helper ───────────────────────────────────────────────────────────────

def _http_get(url: str, params: dict = None, max_attempts: int = 3) -> requests.Response:
    """GET with exponential-backoff retry for transient network errors."""
    last_err: Exception = RuntimeError("no attempts made")
    for attempt in range(max_attempts):
        try:
            response = requests.get(url, params=params, timeout=_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as err:
            last_err = err
            if attempt < max_attempts - 1:
                wait = 2 ** attempt
                console.print(
                    f"[yellow]⚠ UniProt request failed (attempt {attempt + 1}/{max_attempts}): "
                    f"{err} — retrying in {wait}s[/yellow]"
                )
                time.sleep(wait)
    raise last_err


# ── Taxonomy lineage ─────────────────────────────────────────────────────────

def _fetch_taxonomy_lineage(tax_id: int) -> list[dict]:
    """
    Returns genus/family/order ancestors for tax_id from UniProt, closest-first.
    Returns empty list on any failure (non-critical — skips the family prompt).
    """
    try:
        resp = _http_get(
            UNIPROT_TAXONOMY_URL.format(tax_id=tax_id), max_attempts=2
        )
        data    = resp.json()
        lineage = data.get("lineage", [])
        result  = []
        for entry in reversed(lineage):   # root→species → reverse to closest-first
            rank = (entry.get("rank") or "").lower().replace(" ", "")
            if rank in _INFORMATIVE_RANKS:
                result.append({
                    "name":  entry.get("scientificName", ""),
                    "rank":  rank,
                    "taxid": entry.get("taxonId", 0),
                })
        return result
    except Exception:
        return []


# ── Reference loading ─────────────────────────────────────────────────────────

def _load_reference_sequence(input_dir: Path, track_id: str) -> SeqRecord:
    fasta_path = input_dir / track_id / get_step_filename("SEQUENCES", track_id, ext="fasta")
    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in reference FASTA: {fasta_path}")
    return records[0]


def _load_reference_accessions(input_dir: Path, track_id: str) -> set[str]:
    registry_path = input_dir / track_id / get_step_filename("REGISTRY", track_id, ext="json")
    if not registry_path.exists():
        return set()
    with open(registry_path, encoding="utf-8") as fh:
        data = json.load(fh)
    return {s["accession_id"] for s in data.get("sequences", [])}


# ── Variant validation ────────────────────────────────────────────────────────

def _validate_variant(record: SeqRecord) -> tuple[bool, str]:
    """
    Lenient validation for variant sequences (conservation analysis, not binding prediction).
    Rejects only empty sequences or those composed entirely of ambiguous residues.
    Short sequences (< 50 aa) are accepted with a warning emitted by the caller.
    """
    seq = str(record.seq).upper().strip()
    if not seq:
        return False, "empty sequence"
    if all(aa in _AMBIGUOUS_AAS for aa in seq):
        return False, "sequence consists entirely of ambiguous residues (X/B/Z/U/O)"
    return True, ""


# ── Cache redo prompt ─────────────────────────────────────────────────────────

def _ask_redo(n_existing: int) -> bool:
    """
    When the FASTA cache already exists, asks whether to redo.
    Non-interactive mode always keeps the existing file.
    """
    if _is_non_interactive():
        return False

    console.print(f"\n[yellow]Variants FASTA already exists ({n_existing} sequences).[/yellow]")
    console.print("  [cyan]Enter[/cyan] / [cyan]n[/cyan] — keep existing and skip")
    console.print("  [cyan]y[/cyan]       — delete and redo the search\n")

    try:
        raw = input("Redo search? (y/N): ").strip().lower()
    except EOFError:
        raw = "n"

    return raw == "y"


# ── Scope + host filter selection ─────────────────────────────────────────────

def _ask_scope(
    track_id: str,
    track_config: dict,
    project_name: str,
    project_config: dict,
) -> tuple[str, str | None, int | None]:
    """
    Returns (scope, host_filter, family_taxid).

      scope        : "intraspecific" | "interspecific"
      host_filter  : e.g. "Homo sapiens" | None  (interspecific only)
      family_taxid : taxonomy_id of a virus family/genus | None  (interspecific only)
                     when set, query becomes (taxonomy_id:{family_taxid}) AND (protein_name:"…")
                     instead of an unrestricted protein-name search across all organisms.

    All three values are cached in project_config for reproducibility.
    When scope is intraspecific, family_taxid is forced to None (the species tax_id is used directly).
    """
    cached_scope = track_config.get("variants_scope")
    if cached_scope:
        cached_host   = track_config.get("variants_host_filter")
        cached_family = track_config.get("variants_family_taxid")
        parts = [f"[bold]{cached_scope}[/bold]"]
        if cached_family:
            parts.append(f"family tax_id: [bold]{cached_family}[/bold]")
        if cached_host:
            parts.append(f"host filter: [bold]{cached_host}[/bold]")
        console.print(f"[dim]→ Using saved scope: {'  '.join(parts)}[/dim]")
        return cached_scope, cached_host, cached_family

    console.print("\n[bold]Variant search scope[/bold]")
    console.print("  [cyan]1[/cyan] — Intraspecific  : variants within the same species/strain")
    console.print("                    (e.g. different HPV16 isolates, SARS-CoV-2 strains)")
    console.print("  [cyan]2[/cyan] — Interspecific  : same protein across different species")
    console.print("                    (e.g. E5 from HPV16/18/31; spike from SARS/MERS/OC43)")
    console.print("                    → will ask for optional family restriction + host filter\n")
    console.print("  [dim]Note: UniProt intraspecific coverage can be sparse.[/dim]")
    console.print("  [dim]Supplement with a local FASTA in analyze_conservation if needed.[/dim]\n")

    try:
        raw = input("Select scope (1/2, default=1): ").strip()
    except EOFError:
        raw = "1"

    scope = "interspecific" if raw == "2" else "intraspecific"
    console.print(f"[dim]→ Scope: [bold]{scope}[/bold][/dim]")

    family_taxid: int | None = None
    host_filter:  str | None = None

    if scope == "interspecific":
        # ── Family / genus restriction ────────────────────────────────────────
        console.print(
            "\n[bold]Taxonomic restriction[/bold] [dim](optional — recommended to avoid cross-family contamination)[/dim]"
        )
        console.print(
            "[dim]  Restricts search to a specific virus family or genus.[/dim]"
        )
        console.print(
            "[dim]  Without this, UniProt returns proteins from all organisms with the same name.[/dim]\n"
        )

        tax_id = track_config.get("tax_id")
        lineage: list[dict] = []
        if tax_id:
            console.print(f"[dim]  Fetching lineage for tax_id {tax_id}...[/dim]")
            lineage = _fetch_taxonomy_lineage(tax_id)

        if lineage:
            console.print("  Lineage options (closest ancestor first):\n")
            for i, entry in enumerate(lineage, start=1):
                console.print(
                    f"    [cyan]{i}[/cyan]  {entry['name']:<32} [{entry['rank']}]"
                    f"  [dim]tax_id: {entry['taxid']}[/dim]"
                )
            no_restrict_idx = len(lineage) + 1
            console.print(f"    [cyan]{no_restrict_idx}[/cyan]  No restriction (all organisms)")
            console.print()
            try:
                raw_f = input(
                    f"Choice (1-{no_restrict_idx} or custom tax_id, Enter={no_restrict_idx}): "
                ).strip()
            except EOFError:
                raw_f = ""

            if not raw_f or raw_f == str(no_restrict_idx):
                family_taxid = None
                console.print("[dim]→ No taxonomic restriction.[/dim]")
            else:
                try:
                    idx = int(raw_f) - 1
                    if 0 <= idx < len(lineage):
                        family_taxid = lineage[idx]["taxid"]
                        console.print(
                            f"[dim]→ Restricted to [bold]{lineage[idx]['name']}[/bold]"
                            f" ({family_taxid})[/dim]"
                        )
                    else:
                        family_taxid = int(raw_f)
                        console.print(f"[dim]→ Restricted to custom tax_id: [bold]{family_taxid}[/bold][/dim]")
                except ValueError:
                    console.print("[dim]→ Invalid input — no restriction applied.[/dim]")
        else:
            console.print("[dim]  (Lineage unavailable — enter a tax_id manually or press Enter to skip)[/dim]")
            try:
                raw_f = input("  Family/genus tax_id (Enter = no restriction): ").strip()
            except EOFError:
                raw_f = ""
            if raw_f:
                try:
                    family_taxid = int(raw_f)
                    console.print(f"[dim]→ Restricted to tax_id: [bold]{family_taxid}[/bold][/dim]")
                except ValueError:
                    console.print("[dim]→ Invalid input — no restriction applied.[/dim]")

        # ── Host filter ───────────────────────────────────────────────────────
        console.print("\n[bold]Host filter[/bold] [dim](optional)[/dim]")
        console.print("  Example: [cyan]Homo sapiens[/cyan]  — keeps only viruses that infect humans")
        console.print("  Press [cyan]Enter[/cyan] to skip\n")
        try:
            host_raw = input("Host filter (default = no filter): ").strip()
        except EOFError:
            host_raw = ""
        host_filter = host_raw if host_raw else None
        if host_filter:
            console.print(f"[dim]→ Host filter: [bold]{host_filter}[/bold][/dim]")
        else:
            console.print("[dim]→ No host filter applied.[/dim]")

    project_config["tracks"][track_id]["variants_scope"]        = scope
    project_config["tracks"][track_id]["variants_host_filter"]  = host_filter
    project_config["tracks"][track_id]["variants_family_taxid"] = family_taxid
    save_project_config(project_name, project_config)
    return scope, host_filter, family_taxid


# ── UniProt search ────────────────────────────────────────────────────────────

def _extract_protein_name(result: dict) -> str:
    desc = result.get("proteinDescription", {})
    recommended = desc.get("recommendedName", {})
    if recommended:
        return recommended.get("fullName", {}).get("value", "N/A")
    submitted = desc.get("submittedName", [])
    if submitted:
        return submitted[0].get("fullName", {}).get("value", "N/A")
    return "N/A"


def _build_candidates(results: list) -> list[dict]:
    candidates = []
    for r in results:
        reviewed = "Swiss-Prot" in r.get("entryType", "")
        organism = r.get("organism", {})
        seq_data = r.get("sequence", {})
        candidates.append({
            "accession":    r.get("primaryAccession", ""),
            "reviewed":     reviewed,
            "protein_name": _extract_protein_name(r),
            "length":       seq_data.get("length", 0),
            "organism":     organism.get("scientificName", ""),
            "tax_id":       organism.get("taxonId", 0),
            "sequence":     seq_data.get("value", ""),
            "identity":     None,
        })
    return candidates


def _search_uniprot_variants(
    protein_name: str,
    tax_id: int | None,
    scope: str,
    host_filter: str | None,
    family_taxid: int | None = None,
) -> list[dict]:
    """Searches UniProt for protein variants, including sequence data."""
    fields = "accession,reviewed,protein_name,organism_name,organism_id,length,sequence"

    if scope == "intraspecific" and tax_id:
        query = f'(taxonomy_id:{tax_id}) AND (protein_name:"{protein_name}")'
    elif scope == "interspecific" and family_taxid:
        query = f'(taxonomy_id:{family_taxid}) AND (protein_name:"{protein_name}")'
        if host_filter:
            query += f' AND (virus_host_name:"{host_filter}")'
    else:
        query = f'(protein_name:"{protein_name}")'
        if host_filter:
            query += f' AND (virus_host_name:"{host_filter}")'

    console.print(f"[yellow]Searching UniProt...[/yellow]")
    console.print(f"[dim]Query: {query}[/dim]")

    response = _http_get(UNIPROT_SEARCH_URL, params={
        "query":  query,
        "fields": fields,
        "format": "json",
        "size":   str(_MAX_VARIANTS),
    })
    results = response.json().get("results", [])
    console.print(f"[dim]→ {len(results)} raw results returned.[/dim]")
    return _build_candidates(results)


# ── Identity computation ──────────────────────────────────────────────────────

def _compute_identity(seq_a: str, seq_b: str) -> float:
    """
    Returns percent identity (0–100) via global alignment (match=1, gaps=0).

    Denominator is min(len_a, len_b) so that a short fragment that perfectly
    matches its corresponding region in a longer sequence scores 100%, rather
    than being penalised for the length difference (which max() would cause).
    This is the appropriate metric when comparing full proteins to fragments
    or to proteins of different lengths.
    """
    aligner = PairwiseAligner()
    aligner.mode             = "global"
    aligner.match_score      = 1
    aligner.mismatch_score   = 0
    aligner.open_gap_score   = 0
    aligner.extend_gap_score = 0
    score   = aligner.score(seq_a, seq_b)
    min_len = min(len(seq_a), len(seq_b))
    if min_len == 0:
        return 0.0
    return float(score) / min_len * 100.0


# ── Rich table display ────────────────────────────────────────────────────────

def _display_variants_table(candidates: list[dict]):
    """Renders a Rich table of variant candidates colour-coded by identity."""
    table = Table(
        box=box.ROUNDED, show_header=True, header_style="bold white",
        title=f"UniProt Variants — {len(candidates)} candidates",
        title_style="bold cyan",
    )
    table.add_column("#",           no_wrap=True, justify="right", min_width=3)
    table.add_column("Accession",   no_wrap=True, style="cyan",    min_width=12)
    table.add_column("Organism",    no_wrap=False,                  min_width=30, max_width=45)
    table.add_column("% Identity",  no_wrap=True, justify="right", min_width=10)
    table.add_column("Length",      no_wrap=True, justify="right", min_width=8)
    table.add_column("Status",      no_wrap=True,                  min_width=12)

    for i, c in enumerate(candidates, start=1):
        identity = c["identity"]

        if identity is None or identity < 50:
            row_style = "dim red"
        elif identity < 80:
            row_style = "dim"
        else:
            row_style = ""

        identity_str = f"{identity:.1f}%" if identity is not None else "—"
        status_str   = "[bold yellow]★ Swiss-Prot[/bold yellow]" if c["reviewed"] else "[dim]TrEMBL[/dim]"

        table.add_row(
            str(i),
            c["accession"],
            c["organism"][:45],
            identity_str,
            f"{c['length']} aa",
            status_str,
            style=row_style,
        )

    console.print(table)


# ── Selection prompt ──────────────────────────────────────────────────────────

def _parse_selection(raw: str, max_index: int) -> list[int]:
    """
    Parses "1,3,5-8" style selection into sorted 0-based index list.
    "all" → all indices, "none"/"0"/"" → empty list.
    """
    raw = raw.strip().lower()
    if raw in ("none", "0", ""):
        return []
    if raw == "all":
        return list(range(max_index))

    indices: set[int] = set()
    for part in raw.split(","):
        part = part.strip()
        if "-" in part:
            try:
                start_s, end_s = part.split("-", 1)
                for i in range(int(start_s) - 1, int(end_s)):
                    if 0 <= i < max_index:
                        indices.add(i)
            except ValueError:
                pass
        else:
            try:
                idx = int(part) - 1
                if 0 <= idx < max_index:
                    indices.add(idx)
            except ValueError:
                pass

    return sorted(indices)


def _is_non_interactive() -> bool:
    return not sys.stdin.isatty()


def _prompt_multi_selection(candidates: list[dict]) -> list[dict]:
    """Prompts user to select variants by index range. Non-interactive selects all."""
    if _is_non_interactive():
        console.print("[dim]→ Non-interactive mode: selecting all variants.[/dim]")
        return candidates

    console.print(
        f"\n[bold]Select variants to include[/bold] "
        f"[dim](e.g. 1,3,5-8 / all / none — Enter = all)[/dim]"
    )
    try:
        raw = input("> ").strip()
    except EOFError:
        raw = "all"

    if not raw:
        raw = "all"

    indices  = _parse_selection(raw, len(candidates))
    selected = [candidates[i] for i in indices]
    console.print(f"[dim]→ {len(selected)} variant(s) selected.[/dim]")
    return selected


# ── Empty output helper ───────────────────────────────────────────────────────

def _write_empty_outputs(
    fasta_path: Path,
    audit_path: Path,
    track_id: str,
    scope: str,
    host_filter: str | None,
    family_taxid: int | None,
    ref_accessions: set[str],
    note: str,
):
    fasta_path.touch()
    audit = {
        "track_id":           track_id,
        "generated_at":       datetime.datetime.now().isoformat(),
        "scope":              scope,
        "family_taxid":       family_taxid,
        "host_filter":        host_filter,
        "reference_excluded": list(ref_accessions),
        "total_valid":        0,
        "note":               note,
    }
    with open(audit_path, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, ensure_ascii=False)


# ── Validate + build SeqRecords (shared) ─────────────────────────────────────

def _build_and_validate(candidates: list[dict]) -> tuple[list[SeqRecord], list[dict]]:
    """
    Builds SeqRecords from candidate dicts and applies lenient variant validation.
    Short sequences (< 50 aa) are included with a console warning.
    Returns (valid_records, rejected_log).
    """
    valid_records:  list[SeqRecord] = []
    rejected_log:   list[dict]      = []
    short_warnings: list[str]       = []

    for c in candidates:
        seq_str = c.get("sequence", "")
        record  = SeqRecord(
            seq=Seq(seq_str),
            id=c.get("accession", "unknown"),
            description=(
                f"{c.get('organism', '')} | {c.get('protein_name', '')[:60]}"
                + (f" | identity={c['identity']:.1f}%" if c.get("identity") is not None else "")
            ).strip(" |"),
        )
        ok, reason = _validate_variant(record)
        if not ok:
            rejected_log.append({"id": record.id, "reason": reason})
            continue

        if len(seq_str) < 50:
            short_warnings.append(f"{record.id} ({len(seq_str)} aa)")

        valid_records.append(record)

    if short_warnings:
        console.print(
            f"[yellow]⚠ {len(short_warnings)} short sequence(s) included "
            f"(< 50 aa — valid for conservation, not for binding prediction):[/yellow]"
        )
        for w in short_warnings:
            console.print(f"  [yellow]{w}[/yellow]")

    return valid_records, rejected_log


# ── Step class ────────────────────────────────────────────────────────────────

class SearchVariantsStep(BaseTrackStep):
    step_name = "search_variants"

    def describe_outputs(self) -> dict:
        variants_dir = self.track_dir / "variants"
        return {
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
        console.print("[yellow]Computing pairwise identity vs. reference...[/yellow]")
        near_identical: list[str] = []
        low_identity:   list[str] = []
        passing: list[dict]       = []

        _MIN_IDENTITY_THRESHOLD = 30.0  # % — excludes unrelated proteins from other virus families

        for c in candidates:
            if not c["sequence"]:
                continue
            pct = _compute_identity(ref_seq, c["sequence"])
            c["identity"] = round(pct, 2)
            if pct < _MIN_IDENTITY_THRESHOLD:
                low_identity.append(c["accession"])
                continue
            # Near-identical filter: only full-length candidates (≥ 80% of ref length)
            is_full_length = len(c["sequence"]) >= 0.8 * len(ref_seq)
            if pct >= 99.0 and is_full_length:
                near_identical.append(c["accession"])
            else:
                passing.append(c)

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

        console.print(f"\n[bold green]✓ {len(valid_records)} variant(s) saved.[/bold green]")
        console.print(f"  FASTA → {fasta_path}")
        console.print(f"  Audit → {audit_path}")

        return {
            "variants_fasta": str(fasta_path),
            "variants_audit": str(audit_path),
            "total_variants": len(valid_records),
            "cached": False,
        }
