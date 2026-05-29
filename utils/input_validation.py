"""Per-field-type input validation for interactive prompts.

Free-text identifier fields (organism, protein, description) block shell-,
glob- and quoting-breaking characters outright, because those values flow into
file paths and command lines downstream. Accented characters are not blocked —
instead an accent-free NFKD form is suggested, so a user typing "Vírus" is
offered "Virus" rather than having the accent silently mangled.

Path and peptide fields have their own policies (paths keep separators and
drive colons; peptides accept only the 20 amino-acid letters).

All `validate_*` functions are pure (no I/O) so they can be unit-tested in
isolation; `prompt_validated()` wires a validator to `input()` with re-ask and
suggestion handling.
"""

import unicodedata
from dataclasses import dataclass
from typing import Callable, Optional

# Characters that break shells, globs, redirection, or downstream tooling when
# they appear in a free-text identifier. Blocked hard — the user must retype.
_FORBIDDEN_NAME_CHARS = set('/\\()<>|&;"\'*?')

# Path fields legitimately contain separators and Windows drive colons, so only
# redirection, pipe, glob and quoting characters are blocked.
_FORBIDDEN_PATH_CHARS = set('<>|&;"\'*?')

_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


@dataclass
class ValidationResult:
    """Outcome of validating one raw input string.

    ok          — value is acceptable (possibly with an accent suggestion).
    value       — the trimmed input as the user typed it.
    suggestion  — an accent-stripped alternative to offer (None when identical).
    error       — why the value was hard-blocked (None when ok).
    """

    ok: bool
    value: str = ''
    suggestion: Optional[str] = None
    error: Optional[str] = None


def _has_control_chars(text: str) -> bool:
    return any(ord(character) < 32 for character in text)


def strip_accents(text: str) -> str:
    """Returns `text` with combining diacritics removed (NFKD)."""
    decomposed = unicodedata.normalize('NFKD', text)
    return ''.join(character for character in decomposed if not unicodedata.combining(character))


def _format_chars(chars: set) -> str:
    return ' '.join(sorted(chars))


def _validate_free_text(raw: str, field: str, allow_empty: bool = False) -> ValidationResult:
    value = raw.strip()
    if not value:
        if allow_empty:
            return ValidationResult(ok=True, value='')
        return ValidationResult(ok=False, error=f'{field} cannot be empty.')
    if _has_control_chars(value):
        return ValidationResult(ok=False, value=value, error=f'{field} contains control characters.')
    forbidden = {character for character in value if character in _FORBIDDEN_NAME_CHARS}
    if forbidden:
        return ValidationResult(
            ok=False, value=value,
            error=f'{field} cannot contain: {_format_chars(forbidden)}',
        )
    suggestion = strip_accents(value)
    if suggestion != value:
        return ValidationResult(ok=True, value=value, suggestion=suggestion)
    return ValidationResult(ok=True, value=value)


def validate_organism_name(raw: str) -> ValidationResult:
    return _validate_free_text(raw, field='Organism name')


def validate_protein_name(raw: str) -> ValidationResult:
    return _validate_free_text(raw, field='Protein name')


def validate_description(raw: str) -> ValidationResult:
    """Validates a free-prose description. Unlike identifier fields, a
    description never becomes a path or command argument, so only control
    characters are blocked — parentheses, slashes etc. are legitimate prose."""
    value = raw.strip()
    if not value:
        return ValidationResult(ok=True, value='')
    if _has_control_chars(value):
        return ValidationResult(ok=False, value=value, error='Description contains control characters.')
    return ValidationResult(ok=True, value=value)


def validate_local_path(raw: str) -> ValidationResult:
    """Validates a filesystem path string (existence is checked by the caller).

    Surrounding quotes — common when a path is dragged into the terminal — are
    stripped before validation.
    """
    value = raw.strip().strip('"').strip("'")
    if not value:
        return ValidationResult(ok=False, error='Path cannot be empty.')
    if _has_control_chars(value):
        return ValidationResult(ok=False, value=value, error='Path contains control characters.')
    forbidden = {character for character in value if character in _FORBIDDEN_PATH_CHARS}
    if forbidden:
        return ValidationResult(
            ok=False, value=value,
            error=f'Path cannot contain: {_format_chars(forbidden)}',
        )
    return ValidationResult(ok=True, value=value)


def validate_peptide(raw: str) -> ValidationResult:
    """Validates a peptide sequence; accepts only the 20 standard amino acids."""
    value = raw.strip().upper()
    if not value:
        return ValidationResult(ok=False, error='Peptide cannot be empty.')
    offending = set(value) - _AMINO_ACIDS
    if offending:
        return ValidationResult(
            ok=False, value=value,
            error=f'Not a valid amino-acid sequence (offending: {_format_chars(offending)}).',
        )
    return ValidationResult(ok=True, value=value)


def prompt_validated(
    validator: Callable[[str], ValidationResult],
    indent: str = '',
) -> str:
    """`input()` loop that enforces `validator`, re-asking on a hard block.

    Callers render their own header/example panel first; this only draws the
    `> ` input line. When the validator returns an accent suggestion, the user
    is offered the accent-free form
    ([y] accept / [n] keep original). Returns '' in non-interactive sessions so
    callers can gate on an empty result, matching the other prompt helpers in
    `utils.console`.
    """
    from utils.console import console, is_interactive_session

    if not is_interactive_session():
        return ''

    while True:
        try:
            raw = input(f'{indent}> ')
        except EOFError:
            return ''
        result = validator(raw)
        if not result.ok:
            console.print(f'{indent}[red]{result.error}[/red]')
            continue
        if result.suggestion and result.suggestion != result.value:
            console.print(
                f'{indent}[yellow]Did you mean[/yellow] '
                f'[cyan]{result.suggestion}[/cyan]? [dim](accents removed)[/dim]'
            )
            try:
                answer = input(f'{indent}  [y] use suggestion  [n] keep "{result.value}": ').strip().lower()
            except EOFError:
                answer = 'y'
            if answer in ('', 'y', 'yes'):
                return result.suggestion
        return result.value
