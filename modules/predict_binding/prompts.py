"""Interactive binding-parameter prompt for predict_binding (alleles + lengths)."""

from rich.prompt import Prompt

from utils.console import console, is_interactive_session
from utils.naming import parse_hla_allele

# ── Parameter setup ───────────────────────────────────────────────────────────

def _parse_peptide_lengths(raw_input: str, default_lengths: list) -> list:
    """Parses peptide length input — comma list or range.
    Examples: "9" → [9]; "9,10,11" → [9,10,11]; "9-11" → [9,10,11]."""
    raw_input = raw_input.strip()
    if not raw_input:
        return list(default_lengths)

    if '-' in raw_input and ',' not in raw_input:
        range_parts = raw_input.split('-')
        if len(range_parts) == 2:
            try:
                start_length = int(range_parts[0].strip())
                end_length   = int(range_parts[1].strip())
                if 1 <= start_length <= end_length <= 20:
                    return list(range(start_length, end_length + 1))
            except ValueError:
                pass

    try:
        parsed = [int(token.strip()) for token in raw_input.split(',') if token.strip()]
        return parsed if parsed else list(default_lengths)
    except ValueError:
        return list(default_lengths)


def _ask_binding_params(
    project_name: str, project_config: dict, is_rerun: bool = False,
) -> tuple[list[str], list[int]]:
    """Returns (hla_alleles, peptide_lengths) — cached on project_config or prompted once.
    Defaults: 27-allele MHC-I panel + peptide length 9.
    On a rerun (is_rerun) the saved values are offered for editing instead of being
    reused silently."""
    from config import DEFAULT_HLA_ALLELES, DEFAULT_PEPTIDE_LENGTHS

    if 'hla_alleles' in project_config and project_config['hla_alleles'] \
            and 'peptide_lengths' in project_config and project_config['peptide_lengths']:
        hla_alleles     = project_config['hla_alleles']
        peptide_lengths = project_config['peptide_lengths']
        console.print(f"[dim]  Alleles (saved): {len(hla_alleles)} alleles[/dim]")
        console.print(f"[dim]  Peptide lengths (saved): {peptide_lengths}[/dim]")
        if not (is_rerun and is_interactive_session()):
            return hla_alleles, peptide_lengths
        console.print(
            "  [cyan][1][/cyan] Keep saved values   [cyan][2][/cyan] Edit alleles/lengths"
        )
        try:
            edit_choice = input("> ").strip()
        except EOFError:
            edit_choice = "1"
        if edit_choice != "2":
            return hla_alleles, peptide_lengths
        console.print("[dim]  Re-entering binding setup…[/dim]")

    console.print("\n[bold]Binding prediction setup[/bold]")

    # ── HLA alleles ───────────────────────────────────────────────────────────
    console.print(
        f"\n[dim]Default: {len(DEFAULT_HLA_ALLELES)} standard MHC-I alleles "
        f"(A*01:01, A*02:01 ... B*57:01, B*58:01)[/dim]"
    )
    try:
        use_defaults = Prompt.ask("Use default 27 alleles?", choices=['y', 'n'], default='y')
    except Exception:
        use_defaults = 'y'

    if use_defaults.lower() == 'y':
        hla_alleles = list(DEFAULT_HLA_ALLELES)
        console.print(f"[dim]  → {len(hla_alleles)} default alleles selected.[/dim]")
    else:
        hla_alleles = _prompt_and_validate_alleles(default_alleles=DEFAULT_HLA_ALLELES)

    # ── Peptide lengths ───────────────────────────────────────────────────────
    console.print(
        "\n[dim]Peptide lengths to predict (MHC-I binding groove: 8–12 aa).[/dim]\n"
        "[dim]Examples:  9  |  9,10,11  |  9-11[/dim]"
    )
    try:
        length_input = Prompt.ask("Peptide lengths", default="9").strip()
    except EOFError:
        length_input = '9'

    peptide_lengths = _parse_peptide_lengths(length_input, DEFAULT_PEPTIDE_LENGTHS)

    project_config['hla_alleles']     = hla_alleles
    project_config['peptide_lengths'] = peptide_lengths

    from utils.project_manager import save_project_config
    save_project_config(project_name, project_config)

    return hla_alleles, peptide_lengths


def _prompt_and_validate_alleles(default_alleles: list[str]) -> list[str]:
    """Prompts the user for a comma-separated HLA list and validates each entry.
    Loops until every token is a parseable IMGT allele; auto-corrections (missing
    asterisk, lowercase) are shown to the user and require explicit confirmation.
    Returns the list of normalized IMGT alleles. Empty input falls back to defaults."""
    while True:
        console.print("[dim]Enter HLA alleles in IMGT format, comma-separated.[/dim]")
        console.print("[dim]Example: HLA-A*02:01,HLA-B*07:02[/dim]")
        try:
            allele_input = Prompt.ask("HLA alleles").strip()
        except EOFError:
            allele_input = ''

        raw_tokens = [token.strip() for token in allele_input.split(',') if token.strip()]
        if not raw_tokens:
            console.print("[yellow]  No alleles entered, using defaults.[/yellow]")
            return list(default_alleles)

        normalized_alleles: list[str] = []
        corrections: list[tuple[str, str, str]] = []  # (original, normalized, note)
        invalid_tokens: list[str] = []
        for token in raw_tokens:
            normalized, note = parse_hla_allele(token)
            if normalized is None:
                invalid_tokens.append(token)
                continue
            normalized_alleles.append(normalized)
            if note is not None:
                corrections.append((token, normalized, note))

        if invalid_tokens:
            console.print(
                f"[red]  Could not parse {len(invalid_tokens)} entr{'y' if len(invalid_tokens) == 1 else 'ies'}: "
                f"{', '.join(invalid_tokens)}[/red]"
            )
            console.print(
                "[dim]  HLA alleles must follow IMGT format like 'HLA-A*02:01' "
                "(letter A/B/C/E/G, then *, then NN:NN).[/dim]"
            )
            continue

        if corrections:
            console.print("[yellow]  Auto-corrections detected:[/yellow]")
            for original, fixed, note in corrections:
                console.print(f"    [dim]{original}[/dim] → [cyan]{fixed}[/cyan]  [dim]({note})[/dim]")
            try:
                accept = Prompt.ask(
                    "Apply these corrections?",
                    choices=['y', 'n'],
                    default='y',
                ).lower()
            except EOFError:
                accept = 'y'
            if accept != 'y':
                console.print("[dim]  Re-enter the allele list:[/dim]")
                continue

        console.print(f"[dim]  → {len(normalized_alleles)} alleles validated.[/dim]")
        return normalized_alleles

