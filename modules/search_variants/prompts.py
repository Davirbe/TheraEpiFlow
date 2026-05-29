"""Interactive prompts for search_variants: cache-redo, scope/host-filter
selection and multi-selection of candidate variants."""

from utils.console import console, is_interactive_session
from utils.project_manager import save_project_config

from .core import _fetch_taxonomy_lineage

# ── Cache redo prompt ─────────────────────────────────────────────────────────

def _ask_redo(n_existing: int) -> bool:
    """When the FASTA cache exists, asks whether to redo; non-interactive mode keeps it."""
    if not is_interactive_session():
        return False

    console.print(f"\n[yellow]Variants FASTA already exists ({n_existing} sequences).[/yellow]")
    console.print("  [cyan]Enter[/cyan] / [cyan]n[/cyan]  keep existing and skip")
    console.print("  [cyan]y[/cyan]       delete and redo the search\n")

    try:
        raw = input("Redo search? (y/N): ").strip().lower()
    except EOFError:
        raw = "n"

    return raw == "y"


# Above this length a single viral "protein" is more likely an uncut polyprotein.
_POLYPROTEIN_LENGTH_HINT = 1500


def _ensure_tax_id_for_variants(
    track_id: str,
    track_config: dict,
    project_name: str,
    project_config: dict,
    ref_length: int,
) -> int | None:
    """Guards the variant search when the reference came from a local FASTA.

    Warns if the reference looks like an uncut polyprotein, and — when no tax_id is
    known (typical for a local FASTA) — explains that intraspecific search needs one and
    (interactively) offers to enter it, saving it to the track config. Returns the tax_id."""
    tax_id = track_config.get("tax_id")

    if ref_length and ref_length > _POLYPROTEIN_LENGTH_HINT:
        console.print(
            f"[yellow]⚠ Reference is {ref_length} aa; this may be an uncut polyprotein. "
            f"Variant identity can be skewed; consider providing the mature protein.[/yellow]"
        )

    if tax_id:
        return tax_id

    console.print(
        "[yellow]⚠ No tax_id for this track (common for a local FASTA). "
        "Intraspecific search needs a tax_id to restrict to the species; without it the "
        "search falls back to a protein-name query across all organisms.[/yellow]"
    )
    if not is_interactive_session():
        return None

    try:
        raw = input("  NCBI tax_id for species-restricted search (Enter to skip): ").strip()
    except EOFError:
        raw = ""
    if raw.isdigit():
        tax_id = int(raw)
        track_config["tax_id"] = tax_id
        save_project_config(project_name, project_config)
        console.print(f"[dim]→ tax_id set to {tax_id} and saved.[/dim]")
    else:
        console.print("[dim]→ No tax_id provided; proceeding without species restriction.[/dim]")
    return tax_id


# ── Scope + host filter selection ─────────────────────────────────────────────

def _ask_scope(
    track_id: str,
    track_config: dict,
    project_name: str,
    project_config: dict,
    is_rerun: bool = False,
) -> tuple[str, str | None, int | None]:
    """Returns (scope, host_filter, family_taxid); all three cached in project_config.
      scope        : "intraspecific" | "interspecific"
      host_filter  : "Homo sapiens" | None  (interspecific only)
      family_taxid : virus family/genus tax_id | None  (interspecific only) — when set,
                     query becomes (taxonomy_id:{family_taxid}) AND (protein_name:"…")
                     instead of an unrestricted protein-name search.
    On a rerun (is_rerun) the saved scope is offered for editing instead of reused silently.
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
        if not (is_rerun and is_interactive_session()):
            return cached_scope, cached_host, cached_family
        console.print(
            "  [cyan][1][/cyan] Keep saved scope   [cyan][2][/cyan] Change scope/filters"
        )
        try:
            edit_choice = input("> ").strip()
        except EOFError:
            edit_choice = "1"
        if edit_choice != "2":
            return cached_scope, cached_host, cached_family
        console.print("[dim]→ Re-selecting variant search scope…[/dim]")

    console.print("\n[bold]Variant search scope[/bold]")
    console.print("  [cyan]1[/cyan]  Intraspecific : variants within the same species/strain")
    console.print("                    (e.g. different HPV16 isolates, SARS-CoV-2 strains)")
    console.print("  [cyan]2[/cyan]  Interspecific : same protein across different species")
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
            "\n[bold]Taxonomic restriction[/bold] [dim](optional, recommended to avoid cross-family contamination)[/dim]"
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
                    console.print("[dim]→ Invalid input; no restriction applied.[/dim]")
        else:
            console.print("[dim]  (Lineage unavailable; enter a tax_id manually or press Enter to skip)[/dim]")
            try:
                raw_f = input("  Family/genus tax_id (Enter = no restriction): ").strip()
            except EOFError:
                raw_f = ""
            if raw_f:
                try:
                    family_taxid = int(raw_f)
                    console.print(f"[dim]→ Restricted to tax_id: [bold]{family_taxid}[/bold][/dim]")
                except ValueError:
                    console.print("[dim]→ Invalid input; no restriction applied.[/dim]")

        # ── Host filter ───────────────────────────────────────────────────────
        console.print("\n[bold]Host filter[/bold] [dim](optional)[/dim]")
        console.print("  Example: [cyan]Homo sapiens[/cyan]  keeps only viruses that infect humans")
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


# ── Selection prompt ──────────────────────────────────────────────────────────

def _parse_selection(raw: str, max_index: int) -> list[int]:
    """Parses "1,3,5-8" style selection into sorted 0-based indices.
    "all" → every index; "none"/"0"/"" → []."""
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


def _ask_variant_view_mode() -> str:
    """Interspecific only: asks how to view/select candidates.
    Returns 'grouped' (one ★ best per genotype) or 'flat' (full list by identity).
    Non-interactive defaults to 'flat' (preserves prior behaviour)."""
    if not is_interactive_session():
        return "flat"

    console.print("\n[bold]How do you want to pick variants?[/bold]")
    console.print("  [cyan]1[/cyan]  By genotype: one best representative (★) per genotype "
                  "[dim](HPV18, 31, 33, 45, 52…)[/dim]")
    console.print("  [cyan]2[/cyan]  Full list, sorted by identity [dim](current behaviour)[/dim]")
    try:
        raw = input("Select view (1/2, default=1): ").strip()
    except EOFError:
        raw = "1"
    return "flat" if raw == "2" else "grouped"


def _prompt_multi_selection(candidates: list[dict], unit_label: str = "variant") -> list[dict]:
    """Prompts user to select items by index range. Non-interactive selects all.
    unit_label tailors the hint (e.g. 'genotype' in grouped mode)."""
    if not is_interactive_session():
        console.print(f"[dim]→ Non-interactive mode: selecting all {unit_label}s.[/dim]")
        return candidates

    console.print(
        f"\n[bold]Select {unit_label}s to include[/bold] "
        f"[dim](e.g. 1,3,5-8 / all / none; Enter = all)[/dim]"
    )
    try:
        raw = input("> ").strip()
    except EOFError:
        raw = "all"

    if not raw:
        raw = "all"

    indices  = _parse_selection(raw, len(candidates))
    selected = [candidates[i] for i in indices]
    console.print(f"[dim]→ {len(selected)} {unit_label}(s) selected.[/dim]")
    return selected


