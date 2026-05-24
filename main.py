"""
TheraEPIflow MHC-I epitope pipeline. Entry point for every pipeline operation.

Two modes of operation:

  1. Interactive session, recommended for normal use:
       python main.py --project NAME
       python main.py --project NAME --run

  2. One-shot commands, for automation and scripting:
       python main.py --list                              list all projects
       python main.py --new-project                       create a project
       python main.py --project NAME --step STEP_NAME     run one step, then exit
       python main.py --project NAME --status             show progress and exit
       python main.py --project NAME --delete             delete a project

Step identity is the step name itself (for example fetch_sequences or
consensus_filter). There are no numeric step IDs anywhere: folder names,
pipeline.json keys, and CLI flags all use the bare name. Adding, removing
or reordering steps is a STEP_REGISTRY-only edit; the order of entries in
that dict is the order of the pipeline.

Inside the interactive session the user can run the next step, run all
pending steps, retry the last one, jump back, edit a track configuration,
or quit, without the shell closing between commands.
"""

import argparse
import sys
from typing import Optional

from rich.align import Align
from rich.columns import Columns
from rich.console import Group
from rich.table import Table
from rich.panel import Panel
from rich.rule import Rule
from rich import box

from utils.console import console, press_enter_to_continue
from utils.project_manager import (
    create_project_interactive,
    setup_project_tracks_interactive,
    edit_track_interactive,
    load_project_config,
    list_projects,
    tracks_are_defined,
    delete_project,
    update_last_used,
)
from utils.pipeline_state import (
    load_pipeline_state,
    get_track_step_status,
    get_global_step_status,
    reset_track_step,
)
from step_registry import (
    GLOBAL_STEPS,
    STEP_REGISTRY,
    TRACK_STEPS,
    _import_step_class,
    _resolve_step_name_from_user_input,
    _step_is_implemented,
)

PIPELINE_VERSION = '1.0.0-dev'

WELCOME_BANNER = """[bold cyan]
  ████████╗██╗  ██╗███████╗██████╗  █████╗ ███████╗██████╗ ██╗███████╗██╗      ██████╗ ██╗    ██╗
     ██╔══╝██║  ██║██╔════╝██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔════╝██║     ██╔═══██╗██║    ██║
     ██║   ███████║█████╗  ██████╔╝███████║█████╗  ██████╔╝██║█████╗  ██║     ██║   ██║██║ █╗ ██║
     ██║   ██╔══██║██╔══╝  ██╔══██╗██╔══██║██╔══╝  ██╔═══╝ ██║██╔══╝  ██║     ██║   ██║██║███╗██║
     ██║   ██║  ██║███████╗██║  ██║██║  ██║███████╗██║     ██║██║     ███████╗╚██████╔╝╚███╔███╔╝
     ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝     ╚══════╝ ╚═════╝  ╚══╝╚══╝[/bold cyan]"""


# Short, human-friendly column headers for the per-track status table.
# Falls back to the step_name itself if absent.
_STEP_DISPLAY_LABELS: dict = {
    'fetch_sequences':         'fetch',
    'predict_binding':         'predict',
    'consensus_filter':        'consensus',
    'cluster_epitopes':        'cluster',
    'select_representatives':  'select',
    'screen_toxicity':         'toxicity',
    'search_variants':         'variants',
    'analyze_conservation':    'conserv',
    'population_coverage':     'coverage',
    'predict_murine':          'm-pred',
    'curate_murine':           'm-curate',
    'integrate_data':          'integrate',
    'generate_report':         'report',
}


# ── Welcome / header ──────────────────────────────────────────────────────────

def _print_header():
    """Prints the welcome banner and version, centred in the terminal."""
    console.print(Align.center(WELCOME_BANNER))
    console.print(Align.center(
        f'[dim]MHC-I Epitope Pipeline for Vaccine Design  '
        f'·  v{PIPELINE_VERSION}[/dim]\n'
    ))


def _render_pre_step_page(step_class) -> None:
    """Renders the rich pre-step page consumed by `_run_step_interactively`.

    Reads optional ClassVar attrs (`long_description`, `methodology`,
    `references`, `data_format`, `outputs_overview`, `tips`) from the step
    class and renders a two-column layout: explanation on the left,
    references + tips on the right. Falls back to the short `description`
    blurb when no rich attrs are set, preserving behaviour for steps that
    have not been ported yet.
    """
    long_description = (getattr(step_class, 'long_description', '') or '').strip()
    methodology      = (getattr(step_class, 'methodology', '') or '').strip()
    references       = getattr(step_class, 'references', []) or []
    data_format      = (getattr(step_class, 'data_format', '') or '').strip()
    outputs_overview = (getattr(step_class, 'outputs_overview', '') or '').strip()
    tips             = getattr(step_class, 'tips', []) or []
    short_description = (getattr(step_class, 'description', '') or '').strip()

    has_rich_content = any([long_description, methodology, references, data_format, outputs_overview, tips])
    if not has_rich_content:
        if short_description:
            console.print(Panel(
                short_description,
                box=box.SIMPLE,
                border_style='dim',
                padding=(0, 1),
            ))
        return

    # ── Left column: what it does / methodology / data / outputs ──────────────
    left_renderables = []
    if long_description:
        left_renderables.append(Panel(
            long_description,
            title='[bold cyan]What this step does[/bold cyan]',
            title_align='left',
            box=box.ROUNDED,
            border_style='cyan',
            padding=(0, 1),
        ))
    if methodology:
        left_renderables.append(Panel(
            methodology,
            title='[bold]Methodology[/bold]',
            title_align='left',
            box=box.ROUNDED,
            border_style='dim',
            padding=(0, 1),
        ))
    if data_format:
        left_renderables.append(Panel(
            data_format,
            title='[bold]Expected data / input[/bold]',
            title_align='left',
            box=box.ROUNDED,
            border_style='dim',
            padding=(0, 1),
        ))
    if outputs_overview:
        left_renderables.append(Panel(
            outputs_overview,
            title='[bold]What will be generated[/bold]',
            title_align='left',
            box=box.ROUNDED,
            border_style='dim',
            padding=(0, 1),
        ))

    # ── Right column: references + tips ───────────────────────────────────────
    right_renderables = []
    if references:
        reference_lines = []
        for index, reference in enumerate(references, start=1):
            authors = reference.get('authors', '').strip()
            title   = reference.get('title', '').strip()
            journal = reference.get('journal', '').strip()
            year    = reference.get('year', '')
            doi     = reference.get('doi', '').strip()
            url     = reference.get('url', '').strip()
            citation = f'[bold cyan]\\[{index}][/bold cyan] '
            if authors:
                citation += f'{authors}. '
            if title:
                citation += f'[italic]{title}[/italic]. '
            if journal:
                citation += f'{journal}. '
            if year:
                citation += f'{year}.'
            if doi:
                citation += f'\n    [dim]doi:[/dim] [link=https://doi.org/{doi}]{doi}[/link]'
            elif url:
                citation += f'\n    [link]{url}[/link]'
            reference_lines.append(citation)
        right_renderables.append(Panel(
            '\n\n'.join(reference_lines),
            title='[bold]References[/bold]',
            title_align='left',
            box=box.ROUNDED,
            border_style='dim',
            padding=(0, 1),
        ))
    if tips:
        tip_lines = '\n'.join(f'• {tip}' for tip in tips)
        right_renderables.append(Panel(
            tip_lines,
            title='[bold yellow]Tips[/bold yellow]',
            title_align='left',
            box=box.ROUNDED,
            border_style='yellow',
            padding=(0, 1),
        ))

    if left_renderables and right_renderables:
        console.print(Columns(
            [Group(*left_renderables), Group(*right_renderables)],
            expand=True,
            equal=False,
            padding=(0, 1),
        ))
    else:
        for renderable in (left_renderables + right_renderables):
            console.print(renderable)

    console.print(Rule(style='dim'))
    press_enter_to_continue('Press Enter to start the step')


def _print_welcome_page():
    """Renders the landing screen shown when `python main.py` is run with no args.

    Shows the ASCII banner, an About panel (author / lab / links / short
    description), and a quick-start cheatsheet. Pauses for the user to press
    Enter before transitioning to the project menu.
    """
    console.clear()
    console.print(Align.center(WELCOME_BANNER))
    console.print(Align.center(
        f'[dim]MHC-I Epitope Pipeline for Vaccine Design  '
        f'·  v{PIPELINE_VERSION}[/dim]\n'
    ))

    about_body = (
        '[bold]TheraEPIflow[/bold] is an automated CLI pipeline that identifies and selects\n'
        'MHC-I (CTL/CD8+) epitopes for therapeutic vaccine design. A project\n'
        'bundles one or more [bold]tracks[/bold] (each track = one organism × protein\n'
        'pair), processed end-to-end through: sequence fetching, binding\n'
        'prediction (NetMHCpan + MHCFlurry), consensus filtering with Calis\n'
        'immunogenicity, toxicity screening (ToxinPred3), clustering,\n'
        'conservation analysis, population coverage and murine cross-validation.\n\n'
        '[bold]Author[/bold]    Davi Ribeiro  ·  [dim]davi.ribeiro@ufpe.br[/dim]\n'
        '[bold]Lab[/bold]       [bold magenta]◆ LEMTE ◆[/bold magenta]  '
        '[dim]Laboratório de Estudos Moleculares e\n'
        '            Terapia Experimental — UFPE[/dim]\n'
        '[bold]GitHub[/bold]    [link]https://github.com/Davirbe[/link]\n'
        '[bold]LinkedIn[/bold]  [link]https://www.linkedin.com/in/davi-ribeiro-861588186/[/link]'
    )
    console.print(Panel(
        about_body,
        title='[bold cyan]About[/bold cyan]',
        title_align='left',
        box=box.ROUNDED,
        border_style='cyan',
        padding=(1, 2),
    ))

    commands_body = (
        '[cyan]On the next screen[/cyan] you can:\n'
        '  • [bold]Open[/bold] an existing project by typing its number\n'
        r'  • [bold]Create[/bold] a new project           [dim]key:[/dim] [cyan]\[n][/cyan]' '\n'
        r'  • [bold]Edit[/bold] a project track config    [dim]key:[/dim] [cyan]\[e][/cyan]' '\n'
        r'  • [bold]Delete[/bold] a project               [dim]key:[/dim] [cyan]\[d][/cyan]' '\n'
        r'  • [bold]Quit[/bold]                           [dim]key:[/dim] [cyan]\[q][/cyan]' '\n\n'
        '[dim]Inside a project, every step shows its own description and a Commands\n'
        'panel; you do not need to memorize anything — just follow the prompts.[/dim]'
    )
    console.print(Panel(
        commands_body,
        title='[bold]How to use[/bold]',
        title_align='left',
        box=box.ROUNDED,
        border_style='dim',
        padding=(1, 2),
    ))

    press_enter_to_continue('Press Enter to continue to the project menu')


# ── Project menu (interactive TUI) ────────────────────────────────────────────

def _run_project_menu_loop():
    """
    Interactive project menu shown after the welcome page when the user runs
    `python main.py` with no arguments. Loops until the user chooses [q]uit,
    so destructive actions (delete) return to the menu instead of exiting the
    program — which was the previous CLI-only behaviour.
    """
    while True:
        console.clear()
        _print_header()
        all_projects = list_projects(expected_track_step_names=TRACK_STEPS)
        _render_project_menu_table(all_projects)
        _render_project_menu_actions(any_projects=bool(all_projects))

        try:
            raw_input_value = input('  > ').strip()
        except EOFError:
            return
        raw_input_lower = raw_input_value.lower()

        if raw_input_lower in ('q', 'quit', 'exit'):
            console.print('\n[dim]Goodbye.[/dim]')
            return

        if raw_input_lower in ('n', 'new'):
            new_project_name = create_project_interactive()
            console.print(
                f'\n[bold green]Project "{new_project_name}" created. '
                f'Opening interactive session…[/bold green]'
            )
            command_interactive_session(project_name=new_project_name)
            continue

        if raw_input_lower in ('e', 'edit'):
            picked_name = _pick_project_from_menu(all_projects, action_label='edit')
            if picked_name is not None:
                _edit_project_tracks_from_menu(project_name=picked_name)
            continue

        if raw_input_lower in ('d', 'delete'):
            picked_name = _pick_project_from_menu(all_projects, action_label='delete')
            if picked_name is not None:
                command_delete_project(picked_name)
                press_enter_to_continue('Press Enter to return to the project menu')
            continue

        if raw_input_value.isdigit() and all_projects:
            index_chosen = int(raw_input_value)
            if 1 <= index_chosen <= len(all_projects):
                selected_name = all_projects[index_chosen - 1]['name']
                command_interactive_session(project_name=selected_name)
                continue
            console.print(
                f'[red]Out of range. Pick 1..{len(all_projects)}, or n / e / d / q.[/red]'
            )
            press_enter_to_continue()
            continue

        console.print(
            f'[red]Unrecognized: "{raw_input_value}". Pick a project number, '
            f'or use n / e / d / q.[/red]'
        )
        press_enter_to_continue()


def _render_project_menu_table(all_projects: list):
    """Renders the numbered project list shown in the menu loop."""
    if not all_projects:
        console.print(Panel(
            '[dim]No projects yet. Press [cyan]n[/cyan] to create your first one.[/dim]',
            box=box.ROUNDED,
            title='[bold]Projects[/bold]',
            title_align='left',
            border_style='dim',
            padding=(1, 2),
        ))
        return

    table = Table(
        box=box.ROUNDED, show_header=True, header_style='bold white',
        title='Projects', title_style='bold',
        row_styles=['', 'dim'],
    )
    table.add_column('#',           no_wrap=True,  justify='right', min_width=3)
    table.add_column('Project',     no_wrap=True,  style='cyan',    min_width=18)
    table.add_column('Description', no_wrap=False, max_width=32)
    table.add_column('Tracks',      no_wrap=True,  justify='center')
    table.add_column('Status',      no_wrap=True)
    table.add_column('Last used',   no_wrap=True,  style='dim')

    for index, project_data in enumerate(all_projects, start=1):
        track_progress = (
            f"{project_data['completed_tracks']}/{project_data['track_count']}"
            if project_data['track_count'] > 0
            else '—'
        )
        last_used_date = project_data['last_used'][:10] \
            if project_data['last_used'] else '—'

        table.add_row(
            str(index),
            project_data['name'],
            project_data['description'] or '—',
            track_progress,
            project_data['status'],
            last_used_date,
        )

    console.print(table)


def _render_project_menu_actions(any_projects: bool):
    """Renders the action panel shown below the project list."""
    action_lines = []
    if any_projects:
        action_lines.append(
            r'[cyan]\[1..N][/cyan]  open project   '
            r'[cyan]\[n][/cyan]  new   '
            r'[cyan]\[e][/cyan]  edit tracks   '
            r'[cyan]\[d][/cyan]  delete   '
            r'[cyan]\[q][/cyan]  quit'
        )
    else:
        action_lines.append(
            r'[cyan]\[n][/cyan]  new project   '
            r'[cyan]\[q][/cyan]  quit'
        )
    console.print(Panel(
        '\n'.join(action_lines),
        box=box.HEAVY_EDGE,
        border_style='cyan',
        title='[bold cyan]Commands[/bold cyan]',
        title_align='left',
        padding=(0, 2),
    ))
    console.print('  [dim]Type a key from above:[/dim]')


def _pick_project_from_menu(all_projects: list, action_label: str):
    """Asks the user to pick a project by number for an edit/delete action.

    Returns the project name or None if cancelled / invalid.
    """
    if not all_projects:
        console.print(f'[yellow]No projects to {action_label}.[/yellow]')
        press_enter_to_continue()
        return None

    console.print(
        f'\n[bold]Which project do you want to {action_label}?[/bold] '
        f'[dim](number 1..{len(all_projects)}, Enter to cancel)[/dim]'
    )
    try:
        selection_input = input('  > ').strip()
    except EOFError:
        return None
    if not selection_input:
        return None
    if not selection_input.isdigit():
        console.print(f'[red]Pick a number 1..{len(all_projects)}.[/red]')
        press_enter_to_continue()
        return None
    index_chosen = int(selection_input)
    if not (1 <= index_chosen <= len(all_projects)):
        console.print(f'[red]Out of range. Pick 1..{len(all_projects)}.[/red]')
        press_enter_to_continue()
        return None
    return all_projects[index_chosen - 1]['name']


def _edit_project_tracks_from_menu(project_name: str):
    """Edits one of the project's tracks (delegates to _edit_track_from_menu)."""
    if not tracks_are_defined(project_name):
        console.print(
            f'[yellow]Project "{project_name}" has no tracks defined yet — '
            f'open it ([cyan]select its number[/cyan]) and the wizard will run.[/yellow]'
        )
        press_enter_to_continue()
        return
    _edit_track_from_menu(project_name=project_name)
    press_enter_to_continue('Press Enter to return to the project menu')


# ── Project listing ───────────────────────────────────────────────────────────

def command_list_projects(show_header: bool = False):
    """Displays all projects with their current status."""
    if show_header:
        _print_header()

    all_projects = list_projects(expected_track_step_names=TRACK_STEPS)

    if not all_projects:
        console.print(Panel(
            '[dim]No projects found.\n\n'
            'Run [bold cyan]python main.py --new-project[/bold cyan] to create your first project.[/dim]',
            box=box.ROUNDED,
            title='[bold]Projects[/bold]',
            title_align='left',
        ))
        return

    table = Table(
        box=box.ROUNDED, show_header=True, header_style='bold white',
        title='Projects', title_style='bold',
        row_styles=['', 'dim'],
    )
    table.add_column('#',           no_wrap=True,  justify='right', min_width=3)
    table.add_column('Project',     no_wrap=True,  style='cyan',    min_width=18)
    table.add_column('Description', no_wrap=False, max_width=32)
    table.add_column('Tracks',      no_wrap=True,  justify='center')
    table.add_column('Status',      no_wrap=True)
    table.add_column('Last used',   no_wrap=True,  style='dim')

    for index, project_data in enumerate(all_projects, start=1):
        track_progress = (
            f"{project_data['completed_tracks']}/{project_data['track_count']}"
            if project_data['track_count'] > 0
            else '—'
        )
        last_used_date = project_data['last_used'][:10] \
            if project_data['last_used'] else '—'

        table.add_row(
            str(index),
            project_data['name'],
            project_data['description'] or '—',
            track_progress,
            project_data['status'],
            last_used_date,
        )

    console.print(table)


# ── Project status ────────────────────────────────────────────────────────────

def _build_status_table(project_name: str) -> Optional[Table]:
    """
    Builds the per-track per-step status table grouped by organism.
    Returns None if the project has no tracks yet.
    """
    project_config  = load_project_config(project_name)
    pipeline_state  = load_pipeline_state(project_name)
    defined_tracks  = project_config.get('tracks', {})

    if not defined_tracks:
        return None

    table = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    table.add_column('Track', style='cyan', no_wrap=True, min_width=18)

    for step_name in TRACK_STEPS:
        column_label = _STEP_DISPLAY_LABELS.get(step_name, step_name)
        table.add_column(column_label, no_wrap=True, justify='center')

    last_organism_label = None
    for track_id, track_config in defined_tracks.items():
        organism_label = track_config.get('organism_label', track_id)
        protein_label  = track_config.get('protein_label', '')
        display_name   = f'{organism_label} / {protein_label}' if protein_label else track_id

        if last_organism_label is not None and organism_label != last_organism_label:
            table.add_section()
        last_organism_label = organism_label

        track_pipeline_state = pipeline_state.get('tracks', {}).get(track_id, {})
        track_steps_state    = track_pipeline_state.get('steps', {})

        row_values = [display_name]
        for step_name in TRACK_STEPS:
            step_status_record = track_steps_state.get(step_name, {})
            recorded_status    = step_status_record.get('status', 'pending')

            if not _step_is_implemented(step_name):
                status_display = '[dim]—[/dim]'
            elif recorded_status == 'done':
                status_display = '[green]✓[/green]'
            elif recorded_status == 'error':
                status_display = '[red]✗[/red]'
            else:
                status_display = '[dim]○[/dim]'

            row_values.append(status_display)

        table.add_row(*row_values)

    return table


def command_show_status(project_name: str):
    """Shows per-track, per-step status for a project."""
    project_config = load_project_config(project_name)
    status_table   = _build_status_table(project_name)

    if status_table is None:
        console.print(
            f'[yellow]Project "{project_name}" has no tracks defined yet.[/yellow]'
        )
        console.print('[dim]Run with --run to define tracks and start the pipeline.[/dim]')
        return

    console.print(Panel.fit(
        f'[bold cyan]{project_name}[/bold cyan]  '
        f'[dim]{project_config.get("description", "")}[/dim]',
        box=box.ROUNDED,
    ))
    console.print(status_table)
    console.print('[dim]✓ done  ○ pending  ✗ error  — not implemented[/dim]')


# ── Step runners ──────────────────────────────────────────────────────────────

def _run_track_step_for_track(
    step_name: str,
    project_name: str,
    project_config: dict,
    track_id: str,
    force_rerun: bool = False,
    reconfigure: bool = False,
    preflight_config: Optional[dict] = None,
) -> dict:
    """
    Instantiates and executes a per-track step for one track.
    Returns the status dict from execute(): {'status': ..., ...}.
    If the step is not implemented, returns {'status': 'not_implemented'}.

    `preflight_config` (when provided) is forwarded to the step instance so
    it can read overrides decided during the step's pre-iteration hook.
    `reconfigure` runs the step regardless of cached status and re-offers config edits.
    """
    step_class = _import_step_class(step_name)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
        track_id=track_id,
    )
    if preflight_config is not None:
        step_instance.preflight_config = preflight_config
    return step_instance.execute(force_rerun=force_rerun, reconfigure=reconfigure)


def _run_global_step(
    step_name: str,
    project_name: str,
    project_config: dict,
    force_rerun: bool = False,
    reconfigure: bool = False,
) -> dict:
    """Same as above, but for global steps (no track_id)."""
    step_class = _import_step_class(step_name)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
    )
    return step_instance.execute(force_rerun=force_rerun, reconfigure=reconfigure)


# ── Interactive session state inspection ──────────────────────────────────────

def _find_next_pending_step(project_name: str) -> Optional[str]:
    """
    Returns the next step_name with at least one track not yet 'done'.
    Returns None if everything implemented is done.
    """
    project_config    = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    if not defined_track_ids:
        return None

    for step_name in TRACK_STEPS:
        if not _step_is_implemented(step_name):
            continue
        for track_id in defined_track_ids:
            if get_track_step_status(project_name, track_id, step_name) != 'done':
                return step_name

    for step_name in GLOBAL_STEPS:
        if not _step_is_implemented(step_name):
            continue
        if get_global_step_status(project_name, step_name) != 'done':
            return step_name

    return None


def _find_last_completed_step(project_name: str) -> Optional[str]:
    """Returns the last step_name (in registry order) with at least one track 'done'."""
    project_config    = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    last_completed_step_name = None
    for step_name in TRACK_STEPS:
        if not _step_is_implemented(step_name):
            continue
        for track_id in defined_track_ids:
            if get_track_step_status(project_name, track_id, step_name) == 'done':
                last_completed_step_name = step_name
                break
    return last_completed_step_name


# ── Interactive session — core loop ───────────────────────────────────────────

def _run_step_interactively(
    step_name: str,
    project_name: str,
    force_rerun: bool = False,
    reconfigure: bool = False,
    only_tracks: Optional[list] = None,
) -> str:
    """
    Runs one step (for all tracks if track-type, once if global).
    Returns an outcome keyword used by the REPL:
      'completed'        — ran successfully for everything
      'not_implemented'  — step has no class yet
      'had_errors'       — at least one track/global failed
      'aborted'          — user aborted after an error

    `reconfigure` runs regardless of cached status and re-offers config edits.
    `only_tracks` (when given) restricts a track step to those track ids — used by
    the 'retry' command to target one or more specific tracks.
    """
    project_config = load_project_config(project_name)
    _, _, step_type = STEP_REGISTRY[step_name]

    if not _step_is_implemented(step_name):
        console.print(
            f'[yellow]Step "{step_name}" is not implemented yet.[/yellow]'
        )
        return 'not_implemented'

    console.print()
    console.print(Panel.fit(
        f"[bold cyan]{step_name}[/bold cyan]",
        box=box.HEAVY_EDGE,
        border_style="cyan",
        padding=(0, 2),
    ))

    step_class_for_blurb = _import_step_class(step_name)
    if step_class_for_blurb is not None:
        _render_pre_step_page(step_class_for_blurb)

    if step_type == 'global':
        outcome = _run_global_step(
            step_name=step_name,
            project_name=project_name,
            project_config=project_config,
            force_rerun=force_rerun,
            reconfigure=reconfigure,
        )
        if outcome['status'] == 'error':
            return _handle_step_failure(
                step_name=step_name,
                project_name=project_name,
                failed_entity='global',
                error_message=outcome.get('error_message', 'unknown error'),
            )
        return 'completed'

    # Track step — iterate all tracks (or only the ones requested by 'retry')
    defined_track_ids = list(project_config.get('tracks', {}).keys())
    if only_tracks:
        defined_track_ids = [t for t in defined_track_ids if t in only_tracks]
    any_failed = False
    track_outcomes: dict[str, dict] = {}

    # Optional pre-iteration hook: lets the step ask the user a single global
    # question (e.g. "got a local FASTA?") before the loop instead of
    # prompting per-track inside run().
    step_class_for_hooks = _import_step_class(step_name)
    preflight_config: Optional[dict] = None
    if step_class_for_hooks is not None:
        try:
            preflight_config = step_class_for_hooks.preflight(project_name, project_config, defined_track_ids)
        except Exception as preflight_exception:
            console.print(f'[yellow]preflight() raised: {preflight_exception} — continuing without it.[/yellow]')
            preflight_config = None

    for track_id in defined_track_ids:
        outcome = _run_track_step_for_track(
            step_name=step_name,
            project_name=project_name,
            project_config=project_config,
            track_id=track_id,
            force_rerun=force_rerun,
            reconfigure=reconfigure,
            preflight_config=preflight_config,
        )
        track_outcomes[track_id] = outcome
        if outcome['status'] == 'error':
            any_failed = True
            recovery = _handle_step_failure(
                step_name=step_name,
                project_name=project_name,
                failed_entity=track_id,
                error_message=outcome.get('error_message', 'unknown error'),
            )
            if recovery == 'aborted':
                return 'aborted'
            # 'retried_ok' / 'skipped' → continue to next track

        # Reload config in case the step updated it
        project_config = load_project_config(project_name)

    # Optional post-iteration hook: lets the step offer recovery actions
    # (e.g. "want to re-run track X with a local FASTA?") after seeing the
    # full set of outcomes.
    if step_class_for_hooks is not None:
        try:
            step_class_for_hooks.postflight(project_name, project_config, track_outcomes)
        except Exception as postflight_exception:
            console.print(f'[yellow]postflight() raised: {postflight_exception} — ignoring.[/yellow]')

    return 'had_errors' if any_failed else 'completed'


def _retry_failed_step(
    step_name: str,
    project_name: str,
    target_entity: str,
    is_track: bool,
) -> dict:
    """Reloads project_config and reruns the step for the given entity."""
    project_config = load_project_config(project_name)
    if is_track:
        return _run_track_step_for_track(
            step_name=step_name,
            project_name=project_name,
            project_config=project_config,
            track_id=target_entity,
            force_rerun=False,
        )
    return _run_global_step(
        step_name=step_name,
        project_name=project_name,
        project_config=project_config,
        force_rerun=False,
    )


def _handle_step_failure(
    step_name: str,
    project_name: str,
    failed_entity: str,
    error_message: str,
) -> str:
    """
    Interactive recovery prompt after a step fails. Returns one of
    'retried_ok', 'skipped', or 'aborted'. The [e] option opens the track
    configuration editor and retries against the resulting (possibly
    renamed) track. This is the typical fix when fetch_sequences fails
    because the protein name was wrong or a local FASTA path was bad.
    """
    console.print(Panel(
        f'[bold red]Step "{step_name}" failed[/bold red]\n'
        f'[dim]Entity:[/dim] {failed_entity}\n'
        f'[dim]Error: [/dim] {error_message}',
        box=box.ROUNDED, border_style='red',
    ))

    is_track_failure = failed_entity != 'global'
    target_entity = failed_entity

    while True:
        edit_option = '[cyan][e][/cyan] edit track config and retry  ' if is_track_failure else ''
        console.print(
            '\n[bold]What to do?[/bold]  '
            '[cyan][r][/cyan] retry  '
            f'{edit_option}'
            '[cyan][s][/cyan] skip and continue  '
            '[cyan][a][/cyan] abort (back to menu)'
        )
        try:
            user_choice = input('> ').strip().lower()
        except EOFError:
            user_choice = 's'

        if user_choice in ('e', 'edit'):
            if not is_track_failure:
                console.print('[yellow]Cannot edit configuration of a global step.[/yellow]')
                continue
            new_track_id = edit_track_interactive(project_name, target_entity)
            if not new_track_id:
                continue
            target_entity = new_track_id
            user_choice = 'r'

        if user_choice in ('r', 'retry'):
            retry_outcome = _retry_failed_step(
                step_name=step_name,
                project_name=project_name,
                target_entity=target_entity,
                is_track=is_track_failure,
            )
            if retry_outcome['status'] == 'done':
                return 'retried_ok'
            if retry_outcome['status'] == 'error':
                console.print('[red]Retry also failed.[/red]')
                continue
            return 'retried_ok'

        if user_choice in ('s', 'skip'):
            return 'skipped'

        if user_choice in ('a', 'abort', 'q'):
            return 'aborted'

        valid_keys = 'r / e / s / a' if is_track_failure else 'r / s / a'
        console.print(f'[dim]Unrecognized option. Try {valid_keys}.[/dim]')


# ── Interactive session — main REPL ───────────────────────────────────────────

def command_interactive_session(project_name: str):
    """
    Main interactive REPL: stays open until the user quits.
    Auto-advances through pending steps, with retry/skip/jump controls
    between each one.
    """
    update_last_used(project_name)
    _print_header()

    # Ensure tracks exist before entering the loop
    if not tracks_are_defined(project_name):
        console.print(
            f'\n[bold yellow]Project "{project_name}" has no tracks defined.[/bold yellow]'
        )
        console.print('[dim]Defining tracks now...[/dim]\n')
        setup_project_tracks_interactive(project_name)

    while True:
        _print_interactive_status(project_name)
        user_choice = _prompt_interactive_menu(project_name)

        if user_choice == 'quit':
            console.print('\n[dim]Session closed.[/dim]')
            return

        if user_choice == 'status':
            command_show_status(project_name)
            continue

        if user_choice == 'run_next':
            next_step_name = _find_next_pending_step(project_name)
            if next_step_name is None:
                console.print(
                    '\n[bold green]Nothing pending — all implemented steps are done.[/bold green]'
                )
                continue
            _run_step_interactively(step_name=next_step_name, project_name=project_name)
            continue

        if user_choice == 'run_all':
            _run_all_pending(project_name)
            continue

        if user_choice == 'rerun_last':
            last_completed_name = _find_last_completed_step(project_name)
            if last_completed_name is None:
                console.print('[yellow]No completed step to rerun.[/yellow]')
                continue
            console.print(
                f'[dim]Rerunning step "{last_completed_name}" (force).[/dim]'
            )
            _run_step_interactively(
                step_name=last_completed_name,
                project_name=project_name,
                force_rerun=True,
            )
            continue

        if user_choice == 'retry':
            _retry_step_interactively(project_name=project_name)
            continue

        if user_choice == 'edit_track':
            _edit_track_from_menu(project_name=project_name)
            continue

        if user_choice == 'browse':
            from utils.file_browser import run_project_browser
            run_project_browser(project_name=project_name)
            continue

        if isinstance(user_choice, tuple) and user_choice[0] == 'jump':
            target_step_name = user_choice[1]
            _jump_to_step(project_name=project_name, target_step_name=target_step_name)
            continue


def _print_interactive_status(project_name: str):
    """Compact header shown every iteration of the REPL."""
    project_config  = load_project_config(project_name)
    description      = project_config.get('description', '')
    next_pending_name = _find_next_pending_step(project_name)

    if next_pending_name is None:
        next_summary = '[green]All implemented steps done[/green]'
    else:
        next_summary = f'[cyan]{next_pending_name}[/cyan]'

    console.print(Panel.fit(
        f'[bold cyan]{project_name}[/bold cyan]  '
        f'[dim]{description}[/dim]\n'
        f'[dim]Next pending:[/dim] {next_summary}',
        box=box.ROUNDED, border_style='cyan',
    ))

    status_table = _build_status_table(project_name)
    if status_table is not None:
        console.print(status_table)


def _prompt_interactive_menu(project_name: str):
    """
    Shows the REPL menu and returns one of:
      'run_next' | 'run_all' | 'rerun_last' | 'retry' | ('jump', step_name) |
      'edit_track' | 'status' | 'quit'
    """
    menu_table = Table(box=box.SIMPLE_HEAD, show_header=False, padding=(0, 1), expand=False)
    menu_table.add_column(style='bold cyan', no_wrap=True, justify='right')
    menu_table.add_column(style='white')
    menu_table.add_row(r'\[Enter]',      'run next step')
    menu_table.add_row(r'\[a]',          'run all pending steps')
    menu_table.add_row(r'\[r]',          'repeat last step (force re-run)')
    menu_table.add_row(r'\[x]',          'retry a step on chosen track(s), editing config')
    menu_table.add_row(r'\[j <step>]',   'jump to a step (prefix accepted, e.g. "j consensus")')
    menu_table.add_row(r'\[b]',          'browse intermediate files')
    menu_table.add_row(r'\[t]',          'edit track configuration')
    menu_table.add_row(r'\[s]',          'show full status')
    menu_table.add_row(r'\[q]',          'quit')
    console.print()
    console.print(Panel(
        menu_table,
        title='[bold cyan]Commands[/bold cyan]',
        title_align='left',
        box=box.HEAVY_EDGE,
        border_style='cyan',
        padding=(1, 2),
    ))
    console.print('  [dim]Type a key from above (Enter runs the next pending step):[/dim]')
    raw_input_value = input('  > ').strip()
    raw_input_lower = raw_input_value.lower()

    if raw_input_value == '':
        return 'run_next'
    if raw_input_lower in ('q', 'quit', 'exit'):
        return 'quit'
    if raw_input_lower in ('a', 'all'):
        return 'run_all'
    if raw_input_lower in ('r', 'rerun'):
        return 'rerun_last'
    if raw_input_lower in ('x', 'retry'):
        return 'retry'
    if raw_input_lower in ('s', 'status'):
        return 'status'
    if raw_input_lower in ('t', 'track', 'edit'):
        return 'edit_track'
    if raw_input_lower.startswith('j'):
        # 'j fetch' or 'jconsensus_filter'
        remaining_text = raw_input_value[1:].strip().lower()
        if not remaining_text:
            console.print('[red]Provide a step name after "j" (prefix is enough). E.g. "j consensus".[/red]')
            return 'status'
        resolved_step_name = _resolve_step_name_from_user_input(remaining_text)
        if resolved_step_name is None:
            ambiguous = [s for s in STEP_REGISTRY if s.startswith(remaining_text)]
            if len(ambiguous) > 1:
                console.print(f'[red]Ambiguous prefix "{remaining_text}". Matches: {ambiguous}[/red]')
            else:
                console.print(
                    f'[red]No step matches "{remaining_text}". '
                    f'Available: {list(STEP_REGISTRY.keys())}[/red]'
                )
            return 'status'
        return ('jump', resolved_step_name)

    if raw_input_lower in ('b', 'browse'):
        return 'browse'

    console.print(f'[dim]Unrecognized: "{raw_input_value}". Try Enter / a / r / j <step> / b / s / q.[/dim]')
    return 'status'


def _run_all_pending(project_name: str):
    """
    Runs every pending step back-to-back until all are done or a step is aborted.
    User still gets the error-recovery prompt on failures (retry/skip/abort).
    """
    while True:
        next_step_name = _find_next_pending_step(project_name)
        if next_step_name is None:
            console.print('\n[bold green]All done.[/bold green]')
            return
        outcome = _run_step_interactively(step_name=next_step_name, project_name=project_name)
        if outcome == 'aborted':
            console.print('[yellow]Auto-run aborted — back to menu.[/yellow]')
            return
        if outcome == 'not_implemented':
            console.print(
                f'[yellow]Reached step "{next_step_name}" which is not implemented. '
                f'Stopping auto-run.[/yellow]'
            )
            return


def _edit_track_from_menu(project_name: str):
    """
    Lists the project's tracks and lets the user pick one to edit.
    Delegates the actual editing to edit_track_interactive.
    """
    project_config = load_project_config(project_name)
    defined_tracks = project_config.get('tracks', {})
    if not defined_tracks:
        console.print('[yellow]No tracks defined yet.[/yellow]')
        return

    track_id_list = list(defined_tracks.keys())

    console.print('\n[bold]Tracks in this project[/bold]')
    listing = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    listing.add_column('#', style='dim', no_wrap=True, justify='right')
    listing.add_column('Track ID', style='cyan', no_wrap=True)
    listing.add_column('Organism', no_wrap=True)
    listing.add_column('Protein', no_wrap=True)
    listing.add_column('Source', no_wrap=True)
    for index, track_id in enumerate(track_id_list, start=1):
        track_data = defined_tracks[track_id]
        listing.add_row(
            str(index),
            track_id,
            track_data.get('organism_name', ''),
            track_data.get('protein_name', ''),
            track_data.get('input_source', ''),
        )
    console.print(listing)

    console.print(
        '\n[bold]Pick a track to edit[/bold] [dim](number or track ID, '
        'Enter to cancel)[/dim]'
    )
    try:
        selection_input = input('> ').strip()
    except EOFError:
        selection_input = ''
    if not selection_input:
        console.print('[dim]Edit cancelled.[/dim]')
        return

    if selection_input.isdigit():
        index_chosen = int(selection_input)
        if not (1 <= index_chosen <= len(track_id_list)):
            console.print(f'[red]Out of range. Pick 1..{len(track_id_list)}.[/red]')
            return
        track_to_edit = track_id_list[index_chosen - 1]
    elif selection_input in defined_tracks:
        track_to_edit = selection_input
    else:
        console.print(f'[red]No track matches "{selection_input}".[/red]')
        return

    edit_track_interactive(project_name=project_name, track_id=track_to_edit)


def _jump_to_step(project_name: str, target_step_name: str):
    """
    Resets the target step (and all later track steps in registry order) for all
    tracks, so the pipeline effectively rewinds. Output files are not deleted —
    only the pipeline state is cleared.
    """
    project_config    = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    if target_step_name not in TRACK_STEPS:
        console.print(f'[red]"{target_step_name}" is not a per-track step.[/red]')
        return

    target_index = TRACK_STEPS.index(target_step_name)
    steps_to_reset = [
        name for name in TRACK_STEPS[target_index:]
        if _step_is_implemented(name)
    ]
    if not steps_to_reset:
        console.print(
            f'[yellow]No implemented steps from "{target_step_name}" onwards to reset.[/yellow]'
        )
        return

    console.print(
        f'[bold yellow]Jump to step "{target_step_name}"[/bold yellow]\n'
        f'[dim]This will clear the "done" state for: '
        f'{", ".join(steps_to_reset)} across all tracks. '
        f'Output files are NOT deleted — only the pipeline state is reset.[/dim]'
    )
    console.print('[bold]Confirm? (y/n)[/bold]')
    confirmation_input = input('> ').strip().lower()
    if confirmation_input not in ('y', 'yes'):
        console.print('[dim]Jump cancelled.[/dim]')
        return

    for step_name in steps_to_reset:
        for track_id in defined_track_ids:
            reset_track_step(
                project_name=project_name,
                track_id=track_id,
                step_key=step_name,
            )

    console.print(
        f'[green]Reset complete. Next run will start at "{target_step_name}".[/green]'
    )


def _retry_step_interactively(project_name: str):
    """Retry one step on chosen track(s), re-running regardless of cached status and
    re-offering config edits (alleles, strain, scope, length filter, …).

    Unlike 'r' (last completed step only) and 'j' (resets state, then a normal run keeps
    saved config), this targets any step + any track(s) and turns reconfigure on."""
    project_config    = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())
    if not defined_track_ids:
        console.print('[yellow]No tracks defined.[/yellow]')
        return

    console.print('\n[bold]Retry which step?[/bold] [dim](name or prefix, e.g. "predict_murine" / "cons")[/dim]')
    try:
        step_text = input('  step > ').strip().lower()
    except EOFError:
        step_text = ''
    target_step_name = _resolve_step_name_from_user_input(step_text)
    if target_step_name is None or not _step_is_implemented(target_step_name):
        console.print(f'[red]No implemented step matches "{step_text}".[/red]')
        return
    if target_step_name not in TRACK_STEPS:
        console.print('[yellow]Retry currently supports per-track steps only.[/yellow]')
        return

    console.print(
        f'\n[bold]Which track(s)?[/bold] [dim](comma-separated ids, or "all")[/dim]\n'
        f'  [dim]{", ".join(defined_track_ids)}[/dim]'
    )
    try:
        track_text = input('  tracks > ').strip()
    except EOFError:
        track_text = ''
    if not track_text or track_text.lower() == 'all':
        chosen_tracks = list(defined_track_ids)
    else:
        chosen_tracks = [t.strip() for t in track_text.split(',') if t.strip() in defined_track_ids]
    if not chosen_tracks:
        console.print('[yellow]No matching track ids — cancelled.[/yellow]')
        return

    console.print(
        f'[dim]Retrying "{target_step_name}" on: {", ".join(chosen_tracks)} '
        f'(config editable).[/dim]'
    )
    _run_step_interactively(
        step_name=target_step_name,
        project_name=project_name,
        reconfigure=True,
        only_tracks=chosen_tracks,
    )


# ── One-shot commands (for scripting / automation) ────────────────────────────

def command_run_single_step(project_name: str, step_name: str):
    """Runs ONE step for all tracks (or once globally) and exits."""
    update_last_used(project_name)

    if not tracks_are_defined(project_name):
        console.print(
            f'[yellow]Project "{project_name}" has no tracks defined. '
            f'Run without --step first to define them.[/yellow]'
        )
        return

    _run_step_interactively(step_name=step_name, project_name=project_name)


def command_delete_project(project_name: str):
    """Asks for confirmation then permanently deletes the project."""
    console.print(
        f'[bold red]Delete project "{project_name}"?[/bold red] '
        f'[dim]This cannot be undone. Type the project name to confirm:[/dim]'
    )
    confirmation_input = input('> ').strip()
    if confirmation_input == project_name:
        delete_project(project_name)
        console.print(f'[green]Project "{project_name}" deleted.[/green]')
    else:
        console.print('[dim]Deletion cancelled.[/dim]')


# ── CLI entry point ───────────────────────────────────────────────────────────

def main():
    available_step_names = ', '.join(STEP_REGISTRY.keys())

    argument_parser = argparse.ArgumentParser(
        prog='TheraEPIflow',
        description='MHC-I epitope pipeline for vaccine design.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  python main.py --new-project\n'
            '  python main.py --project hpv_analysis                     # interactive session\n'
            '  python main.py --project hpv_analysis --step consensus_filter\n'
            '  python main.py --project hpv_analysis --status\n'
            '  python main.py --list\n\n'
            f'Available steps: {available_step_names}'
        ),
    )

    argument_parser.add_argument('--new-project', action='store_true',
                                 help='Create a new project (interactive wizard).')
    argument_parser.add_argument('--project', metavar='NAME',
                                 help='Name of the project to work with.')
    argument_parser.add_argument('--run', action='store_true',
                                 help='Start the interactive session (same as passing only --project).')
    argument_parser.add_argument('--step', metavar='STEP_NAME', type=str,
                                 help='Run a specific step by name (e.g. consensus_filter), then exit. '
                                      'Unique prefixes accepted.')
    argument_parser.add_argument('--status', action='store_true',
                                 help='Show step-by-step progress, then exit.')
    argument_parser.add_argument('--list', action='store_true',
                                 help='List all projects.')
    argument_parser.add_argument('--delete', action='store_true',
                                 help='Delete the project (requires --project).')

    parsed_args = argument_parser.parse_args()

    # No arguments → show welcome page, then the interactive project menu loop
    if len(sys.argv) == 1:
        _print_welcome_page()
        _run_project_menu_loop()
        return

    if parsed_args.list:
        command_list_projects(show_header=True)
        return

    if parsed_args.new_project:
        _print_header()
        new_project_name = create_project_interactive()
        console.print(f'\n[bold green]Project "{new_project_name}" created. Starting session...[/bold green]')
        command_interactive_session(project_name=new_project_name)
        return

    # All commands below require --project
    if not parsed_args.project:
        console.print('[red]--project NAME is required for this command.[/red]')
        argument_parser.print_help()
        sys.exit(1)

    selected_project_name = parsed_args.project

    if parsed_args.delete:
        command_delete_project(selected_project_name)
        return

    if parsed_args.status:
        command_show_status(selected_project_name)
        return

    if parsed_args.step:
        resolved_step_name = _resolve_step_name_from_user_input(parsed_args.step.lower())
        if resolved_step_name is None:
            console.print(
                f'[red]Step "{parsed_args.step}" is not a known step name. '
                f'Available: {available_step_names}[/red]'
            )
            sys.exit(1)
        command_run_single_step(
            project_name=selected_project_name,
            step_name=resolved_step_name,
        )
        return

    # --project NAME (with or without --run) → interactive session
    command_interactive_session(project_name=selected_project_name)


if __name__ == '__main__':
    main()
