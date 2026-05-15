"""
Project lifecycle management.

A project contains one or more sequence tracks (e.g. HPV16_E6, ZIKV_ENVELOPE).
Tracks share the same pipeline settings but progress independently through steps 01-12.
Steps 13-14 run once after all tracks complete.

Configuration is built incrementally — only what is needed is asked at each stage:
  create_project_interactive()        → name, email, host (5 questions)
  setup_project_tracks_interactive()  → organisms, proteins, input source (at step01)
  Each step adds its own config when it runs (alleles at step03, etc.)

File structure per project:
  projects/{project_name}/
    project_config.json    → grows as each step adds its configuration
    pipeline.json          → per-track progress + global step states
    data/
      input/
        {track_id}/        → sequences_{track_id}.fasta + sequence_registry.json
      intermediate/
        {track_id}/
          predictions/
          consensus/
          clusters/
          toxicity/
          variants/
          conservation/
          coverage/
          murine/
      output/              → master_table.xlsx, report.html
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

from utils.console import console
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

# Track-level "completed" check counts how many per-track steps are marked
# 'done' for each track. The caller (main.py) owns the authoritative list of
# step names, so we receive a count; no step-name knowledge lives here.
TOTAL_GLOBAL_STEPS = 2


# ── Registry helpers ───────────────────────────────────────────────────────────

def _ensure_projects_dir():
    """Creates projects/ and projects_registry.json if they do not exist."""
    PROJECTS_DIR.mkdir(exist_ok=True)
    if not REGISTRY_FILE.exists():
        with open(REGISTRY_FILE, 'w') as registry_file:
            json.dump({'projects': {}}, registry_file, indent=2)


def _load_registry() -> dict:
    _ensure_projects_dir()
    with open(REGISTRY_FILE) as registry_file:
        return json.load(registry_file)


def _save_registry(registry_data: dict):
    with open(REGISTRY_FILE, 'w') as registry_file:
        json.dump(registry_data, registry_file, indent=2, ensure_ascii=False)


# ── Label suggestion ───────────────────────────────────────────────────────────

def _suggest_organism_label(organism_full_name: str) -> str:
    """
    Suggests the standard virus abbreviation label for a given organism name.
    Used as track ID prefix — e.g. HPV16_E6, ZIKV_ENVELOPE, DENV2_NS5.

    Handles the most common pathogens used in vaccine research.
    Falls back to acronym (first letters + numbers) for unknown organisms.

    Examples:
      'Human papillomavirus 16'        → 'HPV16'
      'Human papillomavirus type 18'   → 'HPV18'
      'Zika virus'                     → 'ZIKV'
      'Dengue virus 2'                 → 'DENV2'
      'SARS-CoV-2'                     → 'SARS2'
      'Human immunodeficiency virus 1' → 'HIV1'
      'Influenza A virus'              → 'INFA'
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
    """
    Suggests the standard abbreviation for a protein name.
    Used as the second part of a track ID — e.g. HPV16_E6, ZIKV_E, HIV1_GAG.

    Handles common naming patterns found in GenBank protein records.
    Falls back to uppercase of the input if no pattern matches.

    Examples:
      'E6'                    → 'E6'    (already a label — returned as-is)
      'Early 6'               → 'E6'
      'early protein E5'      → 'E5'
      'Late 1'                → 'L1'
      'envelope protein E'    → 'E'
      'envelope'              → 'E'
      'Non-structural 1'      → 'NS1'
      'NS5'                   → 'NS5'
      'spike protein'         → 'S'
      'nucleocapsid'          → 'N'
      'membrane protein'      → 'M'
      'capsid protein C'      → 'C'
      'polymerase'            → 'POL'
      'gag polyprotein'       → 'GAG'
      'env glycoprotein'      → 'ENV'
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
    """Renders a prompt header in a two-column layout: explanation on the left,
    a concrete filled-in example on the right.

    The `> ` input line is drawn afterwards by the caller. Designed for the
    project / track wizards where pesquisadores com pouco uso de CLI benefit
    from seeing an example side-by-side instead of a description on top.
    """
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


def _prompt_required_nonempty(prompt_label: str, indent: str = '') -> str:
    """input() with non-empty validation. Re-prompts until the user types a value.

    Handles EOFError as empty (so non-interactive runs do not crash; the
    surrounding code is responsible for not gating required fields in non-TTY).
    """
    while True:
        try:
            value = input(f'{indent}> ').strip()
        except EOFError:
            value = ''
        if value:
            return value
        console.print(f'{indent}[red]Empty input not allowed. Please type a value.[/red]')


def create_project(
    project_name: str,
    description: str = '',
    entrez_email: str = '',
) -> Path:
    """
    Creates a new project with the minimal required structure.
    Tracks are NOT defined here — they are added later by setup_project_tracks_interactive()
    when step01 is about to run.

    Returns the project directory path.
    """
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
        'entrez_email': entrez_email,
        'tracks':       {},   # populated by setup_project_tracks_interactive()
    }

    with open(project_dir / 'project_config.json', 'w') as config_file:
        json.dump(project_config, config_file, indent=2, ensure_ascii=False)

    # Pipeline state — tracks added when setup_project_tracks_interactive() runs.
    # Step keys are bare step names (no numeric prefix). See
    # utils/pipeline_state.py for the migration helper that handles legacy keys.
    pipeline_state = {
        'tracks': {},
        'global_steps': {
            'integrate_data':  {'status': 'pending'},
            'generate_report': {'status': 'pending'},
        },
    }

    with open(project_dir / 'pipeline.json', 'w') as pipeline_file:
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
    """
    Minimal project creation wizard — asks only what is needed to create
    the project structure. Everything else is asked contextually later.

    Questions asked:
      1. Project name
      2. Description (optional)
      3. Entrez email (for NCBI GenBank access)

    Target host is fixed at config.TARGET_HOST (Homo sapiens) — HLA-I is human-only.

    Returns the project name.
    """
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
    try:
        description = input('> ').strip()
    except EOFError:
        description = ''

    # ── 3. Entrez email ───────────────────────────────────────────────────────
    try:
        from config import ENTREZ_EMAIL
        default_email_from_config = ENTREZ_EMAIL or ''
    except ImportError:
        default_email_from_config = ''

    email_hint = f'[dim](default: {default_email_from_config})[/dim]' \
        if default_email_from_config else ''
    console.print(f'\n[bold]Entrez email[/bold] {email_hint}')
    console.print('[dim]Required for NCBI GenBank searches (Biopython Entrez).[/dim]')
    try:
        entrez_email_input = input('> ').strip()
    except EOFError:
        entrez_email_input = ''
    entrez_email = entrez_email_input if entrez_email_input else default_email_from_config

    # ── Create project ────────────────────────────────────────────────────────
    project_dir = create_project(
        project_name=project_name,
        description=description,
        entrez_email=entrez_email,
    )

    console.print(f'\n[bold green]✓ Project "{project_name}" created at {project_dir}[/bold green]')
    console.print('[dim]Next: define organisms and proteins when step01 runs.[/dim]')

    return project_name


# ── Track setup ────────────────────────────────────────────────────────────────

def setup_project_tracks_interactive(project_name: str) -> dict:
    """
    Asks which organisms, proteins and input sources to use for this project.
    Called once at the beginning of step01, when no tracks are defined yet.

    For each combination of organism + protein, creates a track entry in
    project_config.json and the corresponding folder structure.

    Questions asked per track:
      - Full organism name (for GenBank search)
      - Short label (for track ID and file names)
      - Protein name
      - Input source: GenBank or local FASTA file
      - If local file: path to the file

    Polyprotein handling is automatic — step01 runs direct + polyprotein
    searches in parallel, no user decision needed here.

    Returns the tracks dict (also saved to project_config.json).
    """
    project_config = load_project_config(project_name)

    console.print(Panel.fit(
        '[bold cyan]Define organisms and proteins[/bold cyan]\n'
        '[dim]Each organism+protein pair is analysed independently and gets its\n'
        'own set of results. All pairs share the same HLA alleles and parameters.\n'
        'Sequences are fetched from UniProt (Swiss-Prot preferred over TrEMBL).[/dim]',
        box=box.ROUNDED,
    ))

    # ── How many organisms? ───────────────────────────────────────────────────
    console.print('\n[bold]How many organisms (or genotypes) to analyze?[/bold] '
                  '[dim](default: 1)[/dim]')
    organism_count_raw = input('> ').strip()
    number_of_organisms = int(organism_count_raw) \
        if organism_count_raw.isdigit() and int(organism_count_raw) > 0 else 1

    tracks_to_create = {}

    for organism_index in range(1, number_of_organisms + 1):
        console.print(f'\n[bold cyan]── Organism {organism_index} of {number_of_organisms} ──[/bold cyan]')

        # Full organism name
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
        organism_full_name = _prompt_required_nonempty('Full organism name')

        # Short label used internally in file names and identifiers
        suggested_label = _suggest_organism_label(organism_full_name)
        _render_prompt_with_example(
            title='Short label',
            hint=f'Used in file names and result identifiers. Press Enter to accept the default ({suggested_label}).',
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
            console.print(
                f'\n  [bold]Protein {protein_index} of {number_of_proteins} '
                f'— {organism_label}[/bold]'
            )

            # Protein full name (used in UniProt search)
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
            protein_full_name = _prompt_required_nonempty('Protein name', indent='  ')

            # Protein label used internally in file names and identifiers
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

            # Use protein_full_name for GenBank search, protein_label for track ID
            protein_name = protein_full_name

            # Input source
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

            local_file_path = None

            if input_source == 'local':
                console.print('  [bold]Path to FASTA file:[/bold]')
                local_file_path = _prompt_required_nonempty('Path to FASTA file', indent='  ')

            # Build track_id from labels (abbreviated): e.g. HPV16_E6, ZIKV_E
            track_id = build_track_id(organism_label, protein_label)

            if track_id in tracks_to_create:
                console.print(
                    f'  [yellow]Warning: "{organism_label} / {protein_label}" already defined. '
                    f'Overwriting.[/yellow]'
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

    # ── Summary ───────────────────────────────────────────────────────────────
    console.print(f'\n[bold green]{len(tracks_to_create)} organism/protein pair(s) defined:[/bold green]')
    summary_table = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    summary_table.add_column('Identifier',   style='cyan', no_wrap=True)
    summary_table.add_column('Organism',   no_wrap=True)
    summary_table.add_column('Protein',    no_wrap=True)
    summary_table.add_column('Source',     no_wrap=True)
    for track_id, track_data in tracks_to_create.items():
        source_label = track_data['input_source']
        if source_label == 'local' and track_data.get('local_file_path'):
            source_label = f'local: {track_data["local_file_path"]}'
        summary_table.add_row(
            track_id,
            track_data['organism_name'],
            track_data['protein_name'],
            source_label,
        )
    console.print(summary_table)

    # ── Save to project_config ────────────────────────────────────────────────
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
    with open(pipeline_path) as pipeline_file:
        pipeline_state = json.load(pipeline_file)

    for track_id in tracks_to_create:
        if track_id not in pipeline_state['tracks']:
            pipeline_state['tracks'][track_id] = {
                'current_step': 0,
                'steps':        {},
            }

    with open(pipeline_path, 'w') as pipeline_file:
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
    """
    Wipes the input and intermediate folders of a track and resets every
    step status in pipeline.json. Used after a track is edited, since the
    upstream FASTA may now point at a different sequence.
    """
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
    with open(pipeline_path) as pipeline_file:
        pipeline_state = json.load(pipeline_file)
    if track_id in pipeline_state.get('tracks', {}):
        pipeline_state['tracks'][track_id] = {'current_step': 0, 'steps': {}}
    with open(pipeline_path, 'w') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)


def _rename_track_on_disk(project_name: str, old_track_id: str, new_track_id: str):
    """
    Renames the folders and pipeline.json keys associated with a track when
    its identifier changes. Caller must guarantee that new_track_id is free.
    """
    project_dir = PROJECTS_DIR / project_name

    for parent in ('input', 'intermediate'):
        old_folder = project_dir / 'data' / parent / old_track_id
        new_folder = project_dir / 'data' / parent / new_track_id
        if old_folder.exists():
            old_folder.rename(new_folder)

    pipeline_path = project_dir / 'pipeline.json'
    with open(pipeline_path) as pipeline_file:
        pipeline_state = json.load(pipeline_file)
    if old_track_id in pipeline_state.get('tracks', {}):
        pipeline_state['tracks'][new_track_id] = pipeline_state['tracks'].pop(old_track_id)
    with open(pipeline_path, 'w') as pipeline_file:
        json.dump(pipeline_state, pipeline_file, indent=2, ensure_ascii=False)


def edit_track_interactive(project_name: str, track_id: str) -> str:
    """
    Walks the user through editing a track's configuration (organism, protein,
    labels, input source). Press Enter to keep a field, type a new value to
    change it.

    Any actual change clears the track's generated data and resets every step
    status, since the input has shifted. If the labels change, the track folders
    and pipeline.json keys are renamed in place.

    Returns the resulting track_id (the new one if labels changed, the original
    otherwise) on success, or an empty string when the user cancels or nothing
    needed changing.
    """
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
        '[bold]Input source[/bold] '
        f'[dim](current: {current_track.get("input_source", "uniprot")}; '
        'options: 1=uniprot, 2=local; Enter to keep)[/dim]'
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
    with open(config_path) as config_file:
        return json.load(config_file)


def save_project_config(project_name: str, updated_config: dict):
    """Saves updated project_config — used by each step to add its configuration."""
    config_path = PROJECTS_DIR / project_name / 'project_config.json'
    with open(config_path, 'w') as config_file:
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
    """
    Returns all projects with metadata for display in the TUI.

    A track is "completed" when every step name in `expected_track_step_names`
    has status='done' for that track. The caller (main.py) owns the canonical
    step list and passes it in — project_manager itself is registry-agnostic.
    If no list is passed, completed_tracks defaults to 0 (display fallback).
    """
    registry = _load_registry()
    projects_list = []
    expected_set = set(expected_track_step_names or [])

    for project_name, project_meta in registry['projects'].items():
        pipeline_path = PROJECTS_DIR / project_name / 'pipeline.json'
        completed_track_count = 0
        total_track_count     = project_meta.get('track_count', 0)

        if pipeline_path.exists() and expected_set:
            with open(pipeline_path) as pipeline_file:
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
