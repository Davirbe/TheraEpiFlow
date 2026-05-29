"""Project lifecycle management.

A project contains one or more sequence tracks (e.g. HPV16_E6, ZIKV_ENVELOPE).
Tracks share pipeline settings but progress independently through the per-track
steps; the global steps run once after all tracks complete.

Configuration is built incrementally — only what is needed is asked at each stage:
  create_project_interactive()        → project name + description
  setup_project_tracks_interactive()  → organisms, proteins, input source
  Each step adds its own config when it runs (alleles at predict_binding, etc.)

File structure per project:
  projects/{project_name}/
    project_config.json    → grows as each step adds its configuration
    pipeline.json          → per-track progress + global step states
    data/
      input/{track_id}/         — fetch_sequences outputs
      intermediate/{track_id}/  — per-step folders (predictions, consensus, clusters, …)
      output/                   — MASTER_TABLE_FULL/VIEW/AUDIT_{project}.*, REPORT_{project}.html
      downloads/                — optional destination for the [z] download menu archives
"""

import json
import re
import shutil
import datetime
from pathlib import Path

from rich.columns import Columns
from rich.panel import Panel
from rich.table import Table
from rich import box

from utils.console import console, confirm_value, show_recap_and_confirm
from utils.input_validation import (
    prompt_validated,
    validate_description,
    validate_local_path,
    validate_organism_name,
    validate_protein_name,
)
from utils.naming import build_track_id

PROJECTS_DIR = Path('projects')
REGISTRY_FILE = PROJECTS_DIR / 'projects_registry.json'

TRACK_INTERMEDIATE_FOLDERS = [
    'predictions',
    'consensus',
    'clusters',
    'toxicity',
    'variants',
    'conservation',
    'coverage',
    'murine',
]

# ── Registry helpers ───────────────────────────────────────────────────────────

def _ensure_projects_dir():
    """Creates projects/ and projects_registry.json if they do not exist."""
    PROJECTS_DIR.mkdir(exist_ok=True)
    if not REGISTRY_FILE.exists():
        with open(REGISTRY_FILE, 'w', encoding='utf-8') as registry_file:
            json.dump({'projects': {}}, registry_file, indent=2)


def _load_registry() -> dict:
    _ensure_projects_dir()
    with open(REGISTRY_FILE, encoding='utf-8') as registry_file:
        return json.load(registry_file)


def _save_registry(registry_data: dict):
    with open(REGISTRY_FILE, 'w', encoding='utf-8') as registry_file:
        json.dump(registry_data, registry_file, indent=2, ensure_ascii=False)


# ── Label suggestion ───────────────────────────────────────────────────────────

def _suggest_organism_label(organism_full_name: str) -> str:
    """Suggests the standard virus abbreviation for an organism name (track ID prefix).
    Falls back to acronym (first letters + numbers) when no pattern matches.

    Examples:
      'Human papillomavirus 16' → 'HPV16'  (numbered family)
      'SARS-CoV-2'              → 'SARS2'  (hyphenated abbreviation)
    """
    name_lowercase = organism_full_name.strip().lower()

    # Human papillomavirus N / HPV type N → HPVN
    hpv_number_match = re.search(r'papillomavirus\s+(?:type\s+)?(\d+)', name_lowercase)
    if hpv_number_match:
        return f'HPV{hpv_number_match.group(1)}'
    if 'papillomavirus' in name_lowercase:
        return 'HPV'

    # Zika virus → ZIKV
    if 'zika' in name_lowercase:
        return 'ZIKV'

    # Dengue virus N → DENVN
    dengue_number_match = re.search(r'dengue.*?(\d+)', name_lowercase)
    if dengue_number_match:
        return f'DENV{dengue_number_match.group(1)}'
    if 'dengue' in name_lowercase:
        return 'DENV'

    # SARS-CoV-2 / SARS-CoV → SARS2 / SARS1
    if 'sars-cov-2' in name_lowercase or 'sars cov 2' in name_lowercase:
        return 'SARS2'
    if 'sars-cov' in name_lowercase or 'sars cov' in name_lowercase:
        return 'SARS1'

    # HIV-1 / HIV-2 / Human immunodeficiency virus 1 → HIV1 / HIV2
    hiv_number_match = re.search(r'(?:hiv[-\s]?|immunodeficiency virus\s+)(\d+)', name_lowercase)
    if hiv_number_match:
        return f'HIV{hiv_number_match.group(1)}'
    if 'immunodeficiency virus' in name_lowercase or 'hiv' in name_lowercase:
        return 'HIV'

    # Influenza A / B / C → INFA / INFB / INFC
    for influenza_type in ['a', 'b', 'c']:
        if f'influenza {influenza_type}' in name_lowercase:
            return f'INF{influenza_type.upper()}'
    if 'influenza' in name_lowercase:
        return 'INF'

    # Ebola → EBOV
    if 'ebola' in name_lowercase:
        return 'EBOV'

    # Hepatitis B / C / E → HBVN / HCV / HEV
    hepatitis_match = re.search(r'hepatitis\s+([bce])', name_lowercase)
    if hepatitis_match:
        return f'H{hepatitis_match.group(1).upper()}V'

    # Fallback: first letter of each word + all numbers in the name
    words = re.split(r'[\s\-]+', organism_full_name)
    first_letters = ''.join(word[0] for word in words if word and not word[0].isdigit())
    numbers_found = ''.join(re.findall(r'\d+', organism_full_name))
    return (first_letters + numbers_found).upper()


def _suggest_protein_label(protein_full_name: str) -> str:
    """Suggests the standard abbreviation for a protein name (second part of a track ID).
    Falls back to uppercase of the input when no pattern matches.

    Examples:
      'Early 6'       → 'E6'  (descriptive form → abbreviation)
      'spike protein' → 'S'   (descriptive long name → single letter)
    """
    protein_name_stripped  = protein_full_name.strip()
    protein_name_lowercase = protein_name_stripped.lower()

    # Already a short label (2-5 alphanumeric chars, possibly with digit) → return as-is
    if re.match(r'^[A-Za-z]{1,4}\d{0,2}$', protein_name_stripped):
        return protein_name_stripped.upper()

    # HPV Early proteins: "early 6", "early protein 6", "early protein E6", "E6 protein" → E6
    early_protein_match = re.search(r'early\s+(?:protein\s+)?[eE]?(\d)', protein_name_lowercase)
    if early_protein_match:
        return f'E{early_protein_match.group(1)}'

    # HPV Late proteins: "late 1", "late protein L1" → L1
    late_protein_match = re.search(r'late\s+(?:protein\s+)?[lL]?(\d)', protein_name_lowercase)
    if late_protein_match:
        return f'L{late_protein_match.group(1)}'

    # Envelope: "envelope protein E", "envelope protein", "envelope" → E
    if 'envelope' in protein_name_lowercase:
        return 'E'

    # Non-structural proteins: "non-structural 1", "nonstructural protein 3", "NS5" → NS1/NS3/NS5
    ns_match = re.search(r'(?:non[-\s]?structural|ns)\s*(?:protein\s*)?(\d+)', protein_name_lowercase)
    if ns_match:
        return f'NS{ns_match.group(1)}'

    # Spike / surface → S
    if 'spike' in protein_name_lowercase or 'surface protein' in protein_name_lowercase:
        return 'S'

    # Nucleocapsid → N
    if 'nucleocapsid' in protein_name_lowercase:
        return 'N'

    # Membrane → M
    if protein_name_lowercase in ('membrane', 'membrane protein'):
        return 'M'

    # Capsid → C
    if 'capsid' in protein_name_lowercase and 'nucleocapsid' not in protein_name_lowercase:
        return 'C'

    # Polymerase (generic) → POL
    if 'polymerase' in protein_name_lowercase and 'rna' not in protein_name_lowercase:
        return 'POL'
    if 'rna-dependent rna polymerase' in protein_name_lowercase:
        return 'RDRP'

    # HIV/retrovirus proteins
    if 'gag' in protein_name_lowercase:
        return 'GAG'
    if 'pol' in protein_name_lowercase and len(protein_name_stripped) <= 6:
        return 'POL'
    if protein_name_lowercase in ('env', 'env glycoprotein', 'env protein'):
        return 'ENV'

    # Hemagglutinin / Neuraminidase (influenza)
    if 'hemagglutinin' in protein_name_lowercase:
        return 'HA'
    if 'neuraminidase' in protein_name_lowercase:
        return 'NA'

    # Fallback: uppercase, remove spaces and special chars, keep alphanumeric
    label_cleaned = re.sub(r'[^A-Za-z0-9]', '', protein_name_stripped).upper()
    return label_cleaned[:8] if label_cleaned else protein_name_stripped.upper()


# ── Project creation ───────────────────────────────────────────────────────────

def _render_prompt_with_example(title: str, hint: str, example_lines: list[str]):
    """Renders a prompt header in two columns: explanation on the left, a concrete
    filled-in example on the right. Caller draws the `> ` input line afterwards."""
    explanation_panel = Panel(
        f'[bold]{title}[/bold]\n[dim]{hint}[/dim]',
        box=box.ROUNDED,
        border_style='cyan',
        padding=(0, 1),
    )
    example_panel = Panel(
        '\n'.join(example_lines),
        title='[dim]example[/dim]',
        title_align='left',
        box=box.ROUNDED,
        border_style='dim',
        padding=(0, 1),
    )
    console.print()
    console.print(Columns([explanation_panel, example_panel], expand=False, equal=False, padding=(0, 1)))


def create_project(
    project_name: str,
    description: str = '',
) -> Path:
    """Creates a new project with the minimal required structure; returns the project directory.
    Tracks are added later by setup_project_tracks_interactive() before fetch_sequences runs."""
    _ensure_projects_dir()

    project_dir = PROJECTS_DIR / project_name
    if project_dir.exists():
        raise FileExistsError(f"Project '{project_name}' already exists.")

    # Shared folders (track folders created later when tracks are defined)
    (project_dir / 'data' / 'input').mkdir(parents=True, exist_ok=True)
    (project_dir / 'data' / 'output').mkdir(parents=True, exist_ok=True)
    (project_dir / 'data' / 'intermediate').mkdir(parents=True, exist_ok=True)

    creation_timestamp = datetime.datetime.now().isoformat()

    # Minimal project config — grows as steps add their configurations.
    # target_host is fixed at config.TARGET_HOST (Homo sapiens) — HLA-I is human-only.
    project_config = {
        'project_name': project_name,
        'description':  description,
        'created_at':   creation_timestamp,
        'tracks':       {},   # populated by setup_project_tracks_interactive()
    }

    with open(project_dir / 'project_config.json', 'w', encoding='utf-8') as config_file:
        json.dump(project_config, config_file, indent=2, ensure_ascii=False)

    # Pipeline state — tracks added when setup_project_tracks_interactive() runs.
    # Step keys are bare step names (no numeric prefix). See
    # utils/pipeline_state.py for the migration helper that handles legacy keys.
    # Global steps are pulled live from step_registry so a new global added later
    # is initialised automatically without touching this file.
    from step_registry import GLOBAL_STEPS
    pipeline_state = {
        'tracks': {},
        'global_steps': {
            step_name: {'status': 'pending'} for step_name in GLOBAL_STEPS
        },
    }

    with open(project_dir / 'pipeline.json', 'w', encoding='utf-8') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)

    # Register in global registry
    registry = _load_registry()
    registry['projects'][project_name] = {
        'created_at':  creation_timestamp,
        'last_used':   creation_timestamp,
        'description': description,
        'track_count': 0,
        'status':      'in_progress',
    }
    _save_registry(registry)

    return project_dir


def create_project_interactive() -> str:
    """Minimal project creation wizard — asks project name + optional description.
    Everything else is asked contextually later; returns the project name."""
    console.print(Panel.fit(
        '[bold cyan]New Project[/bold cyan]\n'
        '[dim]Minimal setup — organisms, proteins and parameters are asked later.[/dim]',
        box=box.ROUNDED,
    ))

    # ── 1. Project name ───────────────────────────────────────────────────────
    _render_prompt_with_example(
        title='Project name',
        hint='Lowercase letters, numbers, hyphens and underscores only. No spaces.',
        example_lines=[
            '[bold]Example values[/bold]',
            '  hpv_analysis',
            '  zikv_v1',
            '  dengue-2026',
        ],
    )
    while True:
        try:
            project_name_input = input('> ').strip().lower()
        except EOFError:
            project_name_input = ''
        if not project_name_input:
            console.print('[red]Empty input not allowed. Please type a project name.[/red]')
            continue
        if not re.match(r'^[a-z0-9][a-z0-9_-]*$', project_name_input):
            console.print('[red]Invalid name. Use: a-z, 0-9, _ and - only.[/red]')
            continue
        if (PROJECTS_DIR / project_name_input).exists():
            console.print(f'[red]Project "{project_name_input}" already exists.[/red]')
            continue
        break
    project_name = project_name_input

    # ── 2. Description ────────────────────────────────────────────────────────
    console.print('\n[bold]Description[/bold] [dim](optional — press Enter to skip)[/dim]')
    description = prompt_validated(validate_description)

    # ── Create project ────────────────────────────────────────────────────────
    project_dir = create_project(
        project_name=project_name,
        description=description,
    )

    console.print(f'\n[bold green]✓ Project "{project_name}" created at {project_dir}[/bold green]')
    console.print('[dim]Next: define organisms and proteins when fetch_sequences runs.[/dim]')

    return project_name


# ── Track setup ────────────────────────────────────────────────────────────────

def _collect_tracks_interactive(default_organism_count: int = 1) -> dict:
    """Walks the user through defining organisms and proteins; returns tracks_to_create.

    Critical free-text fields are followed by a y/n confirm so a typo is fixable in one
    keystroke; the per-pair recap lets the user redo a single pair. Nothing is persisted
    here — caller owns the final recap and saves only on confirmation, so the outer restart
    loop can throw away the whole result on 'n'."""
    # ── How many organisms? ───────────────────────────────────────────────────
    console.print('\n[bold]How many organisms (or genotypes) to analyze?[/bold] '
                  f'[dim](default: {default_organism_count})[/dim]')
    try:
        organism_count_raw = input('> ').strip()
    except EOFError:
        organism_count_raw = ''
    number_of_organisms = int(organism_count_raw) \
        if organism_count_raw.isdigit() and int(organism_count_raw) > 0 \
        else default_organism_count

    tracks_to_create: dict = {}

    try:
        from modules.fetch_sequences import ORGANISM_ALIASES
        known_aliases = ', '.join(sorted(ORGANISM_ALIASES.keys()))
        organism_hint = (
            f'Known aliases: {known_aliases}.\n'
            'You may type an alias (e.g. HPV16) or a full scientific name. '
            'Unknown names are searched as free text in UniProt.'
        )
    except ImportError:
        organism_hint = 'Type the organism name as it should be searched in UniProt.'

    for organism_index in range(1, number_of_organisms + 1):
        console.print(
            f'\n[bold cyan]── Organism {organism_index} of '
            f'{number_of_organisms} ──[/bold cyan]'
        )

        # Per-field confirm: organism full name
        while True:
            _render_prompt_with_example(
                title='Full organism name',
                hint=organism_hint,
                example_lines=[
                    '[bold]Example values[/bold]',
                    '  HPV16',
                    '  Human papillomavirus 16',
                    '  CHIKV',
                    '  Chikungunya virus',
                ],
            )
            organism_full_name = prompt_validated(validate_organism_name)
            if confirm_value('Organism', organism_full_name):
                break

        # Short label (default-allowed → no per-field confirm needed)
        suggested_label = _suggest_organism_label(organism_full_name)
        _render_prompt_with_example(
            title='Short label',
            hint=f'Used in file names and result identifiers. '
                 f'Press Enter to accept the default ({suggested_label}).',
            example_lines=[
                '[bold]Example values[/bold]',
                '  HPV16',
                '  CHIKV',
                '  ZIKV',
                f'  [dim]default →[/dim] {suggested_label}',
            ],
        )
        try:
            organism_label_input = input('> ').strip()
        except EOFError:
            organism_label_input = ''
        organism_label = organism_label_input.upper() \
            if organism_label_input else suggested_label

        # How many proteins for this organism?
        console.print(
            f'\n[bold]How many proteins for {organism_label}?[/bold] '
            f'[dim](default: 1)[/dim]'
        )
        try:
            protein_count_raw = input('> ').strip()
        except EOFError:
            protein_count_raw = ''
        number_of_proteins = int(protein_count_raw) \
            if protein_count_raw.isdigit() and int(protein_count_raw) > 0 else 1

        for protein_index in range(1, number_of_proteins + 1):
            # Per-protein retry loop — the user can redo a single pair without
            # restarting the whole organism.
            while True:
                console.print(
                    f'\n  [bold]Protein {protein_index} of {number_of_proteins} '
                    f'— {organism_label}[/bold]'
                )

                # Per-field confirm: protein full name
                while True:
                    _render_prompt_with_example(
                        title='Protein name',
                        hint='Name as it should be searched in UniProt — full or abbreviated.',
                        example_lines=[
                            '[bold]Example values[/bold]',
                            '  E6',
                            '  envelope protein',
                            '  nsP1',
                            '  spike glycoprotein',
                        ],
                    )
                    protein_full_name = prompt_validated(validate_protein_name, indent='  ')
                    if confirm_value('Protein name', protein_full_name, indent='  '):
                        break

                # Protein short label (default-allowed)
                suggested_protein_label = _suggest_protein_label(protein_full_name)
                console.print(
                    f'  [bold]Protein label[/bold] '
                    f'[dim](used in file names and identifiers — default: {suggested_protein_label})[/dim]'
                )
                try:
                    protein_label_input = input('  > ').strip()
                except EOFError:
                    protein_label_input = ''
                protein_label = protein_label_input.upper() \
                    if protein_label_input else suggested_protein_label

                # Per-field confirm: sequence source
                while True:
                    _render_prompt_with_example(
                        title='Sequence source',
                        hint='Type 1 to search UniProt (default), or 2 to provide a local FASTA file.',
                        example_lines=[
                            '[bold]Options[/bold]',
                            '  [cyan]1[/cyan]  Search UniProt   [dim](default)[/dim]',
                            '  [cyan]2[/cyan]  Local FASTA file',
                        ],
                    )
                    try:
                        source_input = input('  > ').strip()
                    except EOFError:
                        source_input = ''
                    input_source = 'local' if source_input == '2' else 'uniprot'
                    source_label_display = 'Local FASTA file' if input_source == 'local' else 'Search UniProt'
                    if confirm_value('Source', source_label_display, indent='  '):
                        break

                # Per-field confirm: local FASTA path (only if local)
                local_file_path = None
                if input_source == 'local':
                    while True:
                        console.print('  [bold]Path to FASTA file:[/bold]')
                        local_file_path = prompt_validated(validate_local_path, indent='  ')
                        if confirm_value('FASTA path', local_file_path, indent='  '):
                            break

                # Per-pair recap — chance to redo this protein before committing it
                recap_source = (
                    f'local: {local_file_path}' if input_source == 'local'
                    else 'Search UniProt'
                )
                if show_recap_and_confirm(
                    f'Recap: {organism_label} / {protein_label}',
                    [
                        ('Organism',       organism_full_name),
                        ('Organism label', organism_label),
                        ('Protein',        protein_full_name),
                        ('Protein label', protein_label),
                        ('Source',         recap_source),
                    ],
                    proceed_label='Keep this pair',
                ):
                    break
                console.print('\n  [yellow]Re-entering this protein…[/yellow]')

            # Build the internal identifier and store the pair
            track_id = build_track_id(organism_label, protein_label)
            if track_id in tracks_to_create:
                console.print(
                    f'  [yellow]Warning: "{organism_label} / {protein_label}" '
                    f'already defined. Overwriting.[/yellow]'
                )
            tracks_to_create[track_id] = {
                'organism_name':    organism_full_name,
                'organism_label':   organism_label,
                'protein_name':     protein_full_name,
                'protein_label':    protein_label,
                'input_source':     input_source,
                'local_file_path':  local_file_path,
            }
            console.print(f'  [green]✓ Added: {organism_label} / {protein_label}[/green]')

    return tracks_to_create


def setup_project_tracks_interactive(project_name: str) -> dict:
    """Asks which organisms, proteins and input sources to use, called once at the
    beginning of fetch_sequences when no tracks are defined yet. Creates one track entry
    per (organism, protein) pair in project_config.json + the corresponding folder
    structure. Returns the tracks dict (also saved to project_config.json)."""
    project_config = load_project_config(project_name)

    console.print(Panel.fit(
        '[bold cyan]Define organisms and proteins[/bold cyan]\n'
        '[dim]Each organism+protein pair becomes a [bold]track[/bold] — a unit of\n'
        'work the pipeline processes end-to-end and reports separately. All\n'
        'tracks in this project share the same HLA alleles and parameters.\n'
        'Sequences are fetched from UniProt (Swiss-Prot preferred over TrEMBL).[/dim]',
        box=box.ROUNDED,
    ))

    # Wizard runs inside an outer "restart loop". Nothing is persisted until
    # the final recap returns True; if the user types 'n' at the final recap,
    # the entire collection restarts and re-asks every value (including the
    # organism count).
    while True:
        tracks_to_create = _collect_tracks_interactive(default_organism_count=1)

        final_recap_fields = []
        for track_id, td in tracks_to_create.items():
            source_label = td['input_source']
            if source_label == 'local' and td.get('local_file_path'):
                source_label = f'local: {td["local_file_path"]}'
            final_recap_fields.append((
                track_id,
                f'{td["organism_name"]} / {td["protein_name"]}  ({source_label})',
            ))

        if show_recap_and_confirm(
            f'{len(tracks_to_create)} organism/protein pair(s) ready to save',
            final_recap_fields,
            proceed_label='Save and continue',
        ):
            break
        console.print('\n[yellow]Restarting the wizard — no values were saved.[/yellow]\n')

    # ── Save to project_config (only after final recap is confirmed) ──────────
    project_config['tracks'] = tracks_to_create
    save_project_config(project_name, project_config)

    # ── Create folder structure for each track ────────────────────────────────
    project_dir = PROJECTS_DIR / project_name
    for track_id in tracks_to_create:
        for intermediate_folder in TRACK_INTERMEDIATE_FOLDERS:
            (project_dir / 'data' / 'intermediate' / track_id / intermediate_folder).mkdir(
                parents=True, exist_ok=True
            )
        (project_dir / 'data' / 'input' / track_id).mkdir(parents=True, exist_ok=True)

    # ── Update pipeline.json ──────────────────────────────────────────────────
    pipeline_path = project_dir / 'pipeline.json'
    with open(pipeline_path, encoding='utf-8') as pipeline_file:
        pipeline_state = json.load(pipeline_file)

    for track_id in tracks_to_create:
        if track_id not in pipeline_state['tracks']:
            pipeline_state['tracks'][track_id] = {
                'current_step': 0,
                'steps':        {},
            }

    with open(pipeline_path, 'w', encoding='utf-8') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)

    # ── Update registry track count ───────────────────────────────────────────
    registry = _load_registry()
    if project_name in registry['projects']:
        registry['projects'][project_name]['track_count'] = len(tracks_to_create)
    _save_registry(registry)

    console.print('[green]Tracks saved. Ready to run step01.[/green]')
    return tracks_to_create


# ── Track editing ─────────────────────────────────────────────────────────────

def _clean_track_data(project_name: str, track_id: str):
    """Wipes input + intermediate folders of a track and resets every step status.
    Used after a track is edited — upstream FASTA may now point at a different sequence."""
    project_dir = PROJECTS_DIR / project_name

    track_input_dir = project_dir / 'data' / 'input' / track_id
    if track_input_dir.exists():
        shutil.rmtree(track_input_dir)
    track_input_dir.mkdir(parents=True, exist_ok=True)

    track_intermediate_dir = project_dir / 'data' / 'intermediate' / track_id
    if track_intermediate_dir.exists():
        shutil.rmtree(track_intermediate_dir)
    for intermediate_folder in TRACK_INTERMEDIATE_FOLDERS:
        (track_intermediate_dir / intermediate_folder).mkdir(parents=True, exist_ok=True)

    pipeline_path = project_dir / 'pipeline.json'
    with open(pipeline_path, encoding='utf-8') as pipeline_file:
        pipeline_state = json.load(pipeline_file)
    if track_id in pipeline_state.get('tracks', {}):
        pipeline_state['tracks'][track_id] = {'current_step': 0, 'steps': {}}
    # The global steps (integrate_data, generate_report) aggregate every track,
    # so editing one track makes their outputs stale — drop their 'done' status
    # so they rebuild on the next run.
    pipeline_state['global_steps'] = {}
    with open(pipeline_path, 'w', encoding='utf-8') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)


def _rename_track_on_disk(project_name: str, old_track_id: str, new_track_id: str):
    """Renames the folders and pipeline.json keys for a track when its identifier changes.
    Caller must guarantee new_track_id is free."""
    project_dir = PROJECTS_DIR / project_name

    for parent in ('input', 'intermediate'):
        old_folder = project_dir / 'data' / parent / old_track_id
        new_folder = project_dir / 'data' / parent / new_track_id
        if old_folder.exists():
            old_folder.rename(new_folder)

    pipeline_path = project_dir / 'pipeline.json'
    with open(pipeline_path, encoding='utf-8') as pipeline_file:
        pipeline_state = json.load(pipeline_file)
    if old_track_id in pipeline_state.get('tracks', {}):
        pipeline_state['tracks'][new_track_id] = pipeline_state['tracks'].pop(old_track_id)
    with open(pipeline_path, 'w', encoding='utf-8') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)


def edit_track_interactive(project_name: str, track_id: str) -> str:
    """Edits a track's configuration interactively (Enter to keep a field, type a value to change).
    Any change clears generated data + resets step statuses; label changes rename folders + pipeline keys.
    Returns the resulting track_id (new if labels changed, otherwise unchanged) or '' on cancel/no-op."""
    project_config = load_project_config(project_name)
    tracks_dict = project_config.get('tracks', {})

    if track_id not in tracks_dict:
        console.print(f'[red]Track "{track_id}" not found in project "{project_name}".[/red]')
        return ''

    current_track = tracks_dict[track_id]
    new_track = dict(current_track)

    console.print(Panel.fit(
        f'[bold cyan]Editing track[/bold cyan]  [white]{track_id}[/white]\n'
        '[dim]Press Enter to keep the current value. Type a new value to change it.[/dim]',
        box=box.ROUNDED,
    ))

    fields_to_edit = [
        ('organism_name',  'Full organism name', None),
        ('organism_label', 'Short organism label (used in track ID)', 'upper'),
        ('protein_name',   'Full protein name', None),
        ('protein_label',  'Short protein label (used in track ID)', 'upper'),
    ]
    for field_key, prompt_label, transform in fields_to_edit:
        current_value = current_track.get(field_key, '')
        console.print(f'[bold]{prompt_label}[/bold] [dim](current: {current_value or "(empty)"})[/dim]')
        try:
            user_input = input('> ').strip()
        except EOFError:
            user_input = ''
        if user_input:
            new_track[field_key] = user_input.upper() if transform == 'upper' else user_input

    console.print(
        '[bold]Input source[/bold]  '
        '[cyan][1][/cyan] uniprot  [cyan][2][/cyan] local  '
        f'[dim](current: {current_track.get("input_source", "uniprot")}; Enter to keep)[/dim]'
    )
    try:
        source_input = input('> ').strip()
    except EOFError:
        source_input = ''
    if source_input == '1':
        new_track['input_source'] = 'uniprot'
        new_track['local_file_path'] = None
    elif source_input == '2':
        new_track['input_source'] = 'local'

    if new_track.get('input_source') == 'local':
        current_path = current_track.get('local_file_path') or ''
        console.print(f'[bold]Local FASTA path[/bold] [dim](current: {current_path or "(empty)"})[/dim]')
        try:
            path_input = input('> ').strip()
        except EOFError:
            path_input = ''
        if path_input:
            new_track['local_file_path'] = path_input
    else:
        new_track['local_file_path'] = None

    if new_track == current_track:
        console.print('[dim]No changes made.[/dim]')
        return ''

    new_track_id = build_track_id(new_track['organism_label'], new_track['protein_label'])

    diff_table = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    diff_table.add_column('Field', style='cyan')
    diff_table.add_column('Before')
    diff_table.add_column('After')
    for key in current_track:
        if current_track.get(key) != new_track.get(key):
            diff_table.add_row(key, str(current_track.get(key)), str(new_track.get(key)))
    if track_id != new_track_id:
        diff_table.add_row('track_id', track_id, new_track_id)

    console.print('\n[bold yellow]Pending changes:[/bold yellow]')
    console.print(diff_table)
    console.print(
        '\n[bold]Apply changes?[/bold] '
        '[dim]All generated data for this track will be cleared and every step '
        'status will be reset. (y/n)[/dim]'
    )
    try:
        confirmation = input('> ').strip().lower()
    except EOFError:
        confirmation = 'n'
    if confirmation not in ('y', 'yes'):
        console.print('[dim]Edit cancelled.[/dim]')
        return ''

    if track_id != new_track_id:
        if new_track_id in tracks_dict:
            console.print(
                f'[red]Cannot rename to "{new_track_id}": that track ID already exists. '
                f'Aborting edit.[/red]'
            )
            return ''
        _rename_track_on_disk(project_name, track_id, new_track_id)
        del tracks_dict[track_id]

    tracks_dict[new_track_id] = new_track
    project_config['tracks'] = tracks_dict
    save_project_config(project_name, project_config)

    _clean_track_data(project_name, new_track_id)

    console.print(f'[green]✓ Track updated. Active track ID: {new_track_id}[/green]')
    return new_track_id


# ── Project loading ────────────────────────────────────────────────────────────

def load_project_config(project_name: str) -> dict:
    """Returns the project_config for an existing project."""
    config_path = PROJECTS_DIR / project_name / 'project_config.json'
    if not config_path.exists():
        raise FileNotFoundError(f"Project '{project_name}' not found.")
    with open(config_path, encoding='utf-8') as config_file:
        return json.load(config_file)


def save_project_config(project_name: str, updated_config: dict):
    """Saves updated project_config — used by each step to add its configuration."""
    config_path = PROJECTS_DIR / project_name / 'project_config.json'
    with open(config_path, 'w', encoding='utf-8') as config_file:
        json.dump(updated_config, config_file, indent=2, ensure_ascii=False)


def update_last_used(project_name: str):
    """Updates the last_used timestamp in the registry."""
    registry = _load_registry()
    if project_name in registry['projects']:
        registry['projects'][project_name]['last_used'] = (
            datetime.datetime.now().isoformat()
        )
        _save_registry(registry)


# ── Project listing ────────────────────────────────────────────────────────────

def list_projects(expected_track_step_names: list = None) -> list:
    """Returns all projects with metadata for the TUI list.
    A track is 'completed' when every step in expected_track_step_names is 'done' for it;
    caller (main.py) owns the canonical step list. Without it completed_tracks falls back to 0."""
    registry = _load_registry()
    projects_list = []
    expected_set = set(expected_track_step_names or [])

    for project_name, project_meta in registry['projects'].items():
        pipeline_path = PROJECTS_DIR / project_name / 'pipeline.json'
        completed_track_count = 0
        total_track_count     = project_meta.get('track_count', 0)

        if pipeline_path.exists() and expected_set:
            with open(pipeline_path, encoding='utf-8') as pipeline_file:
                pipeline_state = json.load(pipeline_file)
            for track_state in pipeline_state.get('tracks', {}).values():
                done_step_names = {
                    name for name, value in track_state.get('steps', {}).items()
                    if value.get('status') == 'done'
                }
                if expected_set.issubset(done_step_names):
                    completed_track_count += 1

        projects_list.append({
            'name':             project_name,
            'description':      project_meta.get('description', ''),
            'status':           project_meta.get('status', 'unknown'),
            'created_at':       project_meta.get('created_at', ''),
            'last_used':        project_meta.get('last_used', ''),
            'track_count':      total_track_count,
            'completed_tracks': completed_track_count,
        })

    return projects_list


# ── Track utilities ────────────────────────────────────────────────────────────

def tracks_are_defined(project_name: str) -> bool:
    """Returns True if the project has at least one track defined."""
    project_config = load_project_config(project_name)
    return bool(project_config.get('tracks'))


# ── Project deletion ───────────────────────────────────────────────────────────

def delete_project(project_name: str):
    """Permanently deletes a project folder and removes it from the registry."""
    project_dir = PROJECTS_DIR / project_name
    if not project_dir.exists():
        raise FileNotFoundError(f"Project '{project_name}' not found.")

    shutil.rmtree(project_dir)

    registry = _load_registry()
    registry['projects'].pop(project_name, None)
    _save_registry(registry)
