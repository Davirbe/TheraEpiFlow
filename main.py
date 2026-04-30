"""
TheraEPIflow — MHC-I Epitope Pipeline
Entry point for all pipeline operations.

Two modes of operation:

  1. Interactive session (recommended — stays open between steps):
       python main.py --project NAME
       python main.py --project NAME --run        (starts in interactive session)

  2. One-shot commands (for automation / scripting):
       python main.py --list                            list all projects
       python main.py --new-project                     create a project
       python main.py --project NAME --step STEP_NAME   run ONE step, then exit
       python main.py --project NAME --status           show progress and exit
       python main.py --project NAME --delete           delete a project

Step identity is the step NAME (e.g. fetch_sequences, consensus_filter) — there
are no numeric step IDs anywhere. Adding/removing/reordering steps is a
registry-only edit; folder names, pipeline.json keys and CLI flags all use the
bare step name. The order of steps is the order of entries in STEP_REGISTRY.

The interactive session auto-advances through steps and lets you retry, skip,
jump back, or run everything unattended — without the shell closing between
commands.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box

from utils.project_manager import (
    create_project_interactive,
    setup_project_tracks_interactive,
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

console = Console(width=120)

PIPELINE_VERSION = '1.0.0-dev'

WELCOME_BANNER = """[bold cyan]
  ████████╗██╗  ██╗███████╗██████╗  █████╗ ███████╗██████╗ ██╗███████╗██╗      ██████╗ ██╗    ██╗
     ██╔══╝██║  ██║██╔════╝██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔════╝██║     ██╔═══██╗██║    ██║
     ██║   ███████║█████╗  ██████╔╝███████║█████╗  ██████╔╝██║█████╗  ██║     ██║   ██║██║ █╗ ██║
     ██║   ██╔══██║██╔══╝  ██╔══██╗██╔══██║██╔══╝  ██╔═══╝ ██║██╔══╝  ██║     ██║   ██║██║███╗██║
     ██║   ██║  ██║███████╗██║  ██║██║  ██║███████╗██║     ██║██║     ███████╗╚██████╔╝╚███╔███╔╝
     ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝     ╚══════╝ ╚═════╝  ╚══╝╚══╝[/bold cyan]"""


# ── Step registry ─────────────────────────────────────────────────────────────
# Ordered dict (insertion order = pipeline order). Each entry maps:
#   step_name → (module_path, class_name, step_type)
# step_type: 'track'  = runs once per track
#            'global' = runs once after all tracks have completed every track step
#
# To add/remove/reorder steps, edit ONLY this dict.
# step_name is also the folder name under modules/ and the key in pipeline.json.

STEP_REGISTRY: dict = {
    # Per-track steps (run once per sequence track)
    'fetch_sequences':         ('modules.fetch_sequences',         'FetchSequencesStep',         'track'),
    'predict_binding':         ('modules.predict_binding',         'PredictBindingStep',         'track'),
    'consensus_filter':        ('modules.consensus_filter',        'ConsensusFilterStep',        'track'),
    'screen_toxicity':         ('modules.screen_toxicity',         'ScreenToxicityStep',         'track'),
    'cluster_epitopes':        ('modules.cluster_epitopes',        'ClusterEpitopesStep',        'track'),
    'select_representatives':  ('modules.select_representatives',  'SelectRepresentativesStep', 'track'),
    'search_variants':         ('modules.search_variants',         'SearchVariantsStep',         'track'),
    'analyze_conservation':    ('modules.analyze_conservation',    'AnalyzeConservationStep',    'track'),
    'population_coverage':     ('modules.population_coverage',     'PopulationCoverageStep',     'track'),
    'predict_murine':          ('modules.predict_murine',          'PredictMurineStep',          'track'),
    'curate_murine':           ('modules.curate_murine',           'CurateMurineStep',           'track'),
    # Global steps (run once after all per-track steps complete on every track)
    'integrate_data':          ('modules.integrate_data',          'IntegrateDataStep',          'global'),
    'generate_report':         ('modules.generate_report',         'GenerateReportStep',         'global'),
}

TRACK_STEPS:  list = [name for name, entry in STEP_REGISTRY.items() if entry[2] == 'track']
GLOBAL_STEPS: list = [name for name, entry in STEP_REGISTRY.items() if entry[2] == 'global']

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


# ── Step loading ──────────────────────────────────────────────────────────────

def _import_step_class(step_name: str):
    """
    Dynamically imports and returns the step class for a given step name.
    Returns None if the module exists but the class is not yet implemented
    (empty __init__.py) or the step name is not in the registry.
    """
    if step_name not in STEP_REGISTRY:
        return None
    module_path, class_name, _ = STEP_REGISTRY[step_name]
    try:
        import importlib
        module = importlib.import_module(module_path)
        return getattr(module, class_name, None)
    except ImportError:
        return None


def _step_is_implemented(step_name: str) -> bool:
    return _import_step_class(step_name) is not None


def _resolve_step_name_from_user_input(user_text: str) -> Optional[str]:
    """
    Allows users to type a unique prefix instead of the full step name.
    Returns the matched step_name, or None if no match / ambiguous.
    """
    candidates = [s for s in STEP_REGISTRY if s.startswith(user_text)]
    if len(candidates) == 1:
        return candidates[0]
    if user_text in STEP_REGISTRY:
        return user_text
    return None


# ── Welcome / header ──────────────────────────────────────────────────────────

def _print_header():
    """Prints the welcome banner and version."""
    console.print(WELCOME_BANNER)
    console.print(
        f'  [dim]MHC-I Epitope Pipeline for Vaccine Design  '
        f'·  v{PIPELINE_VERSION}[/dim]\n'
    )


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
    )
    table.add_column('#',           no_wrap=True,  justify='right', min_width=3)
    table.add_column('Project',     no_wrap=True,  style='cyan',    min_width=18)
    table.add_column('Description', no_wrap=False, max_width=32)
    table.add_column('Tracks',      no_wrap=True,  justify='center')
    table.add_column('Status',      no_wrap=True)
    table.add_column('Last used',   no_wrap=True)

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
                status_display = '[dim]·[/dim]'

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
    console.print('[dim]✓ done  · pending  ✗ error  — not implemented[/dim]')


# ── Step runners ──────────────────────────────────────────────────────────────

def _run_track_step_for_track(
    step_name: str,
    project_name: str,
    project_config: dict,
    track_id: str,
    force_rerun: bool = False,
) -> dict:
    """
    Instantiates and executes a per-track step for one track.
    Returns the status dict from execute(): {'status': ..., ...}.
    If the step is not implemented, returns {'status': 'not_implemented'}.
    """
    step_class = _import_step_class(step_name)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
        track_id=track_id,
    )
    return step_instance.execute(force_rerun=force_rerun)


def _run_global_step(
    step_name: str,
    project_name: str,
    project_config: dict,
    force_rerun: bool = False,
) -> dict:
    """Same as above, but for global steps (no track_id)."""
    step_class = _import_step_class(step_name)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
    )
    return step_instance.execute(force_rerun=force_rerun)


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
) -> str:
    """
    Runs one step (for all tracks if track-type, once if global).
    Returns an outcome keyword used by the REPL:
      'completed'        — ran successfully for everything
      'not_implemented'  — step has no class yet
      'had_errors'       — at least one track/global failed
      'aborted'          — user aborted after an error
    """
    project_config = load_project_config(project_name)
    _, _, step_type = STEP_REGISTRY[step_name]

    if not _step_is_implemented(step_name):
        console.print(
            f'[yellow]Step "{step_name}" is not implemented yet.[/yellow]'
        )
        return 'not_implemented'

    console.print(f'\n[bold cyan]══ {step_name} ══[/bold cyan]')

    if step_type == 'global':
        outcome = _run_global_step(
            step_name=step_name,
            project_name=project_name,
            project_config=project_config,
            force_rerun=force_rerun,
        )
        if outcome['status'] == 'error':
            return _handle_step_failure(
                step_name=step_name,
                project_name=project_name,
                failed_entity='global',
                error_message=outcome.get('error_message', 'unknown error'),
            )
        return 'completed'

    # Track step — iterate all tracks
    defined_track_ids = list(project_config.get('tracks', {}).keys())
    any_failed = False

    for track_id in defined_track_ids:
        outcome = _run_track_step_for_track(
            step_name=step_name,
            project_name=project_name,
            project_config=project_config,
            track_id=track_id,
            force_rerun=force_rerun,
        )
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

    return 'had_errors' if any_failed else 'completed'


def _handle_step_failure(
    step_name: str,
    project_name: str,
    failed_entity: str,
    error_message: str,
) -> str:
    """
    Interactive recovery prompt after a step fails.
    Returns one of: 'retried_ok', 'retried_failed', 'skipped', 'aborted'.
    """
    console.print(Panel(
        f'[bold red]Step "{step_name}" failed[/bold red]\n'
        f'[dim]Entity:[/dim] {failed_entity}\n'
        f'[dim]Error: [/dim] {error_message}',
        box=box.ROUNDED, border_style='red',
    ))

    while True:
        console.print(
            '\n[bold]What to do?[/bold]  '
            '[cyan][r][/cyan] retry  '
            '[cyan][s][/cyan] skip (mark pending, continue)  '
            '[cyan][a][/cyan] abort (back to menu)'
        )
        try:
            user_choice = input('> ').strip().lower()
        except EOFError:
            user_choice = 's'

        if user_choice in ('r', 'retry'):
            project_config = load_project_config(project_name)
            if failed_entity == 'global':
                retry_outcome = _run_global_step(
                    step_name=step_name,
                    project_name=project_name,
                    project_config=project_config,
                    force_rerun=False,
                )
            else:
                retry_outcome = _run_track_step_for_track(
                    step_name=step_name,
                    project_name=project_name,
                    project_config=project_config,
                    track_id=failed_entity,
                    force_rerun=False,
                )
            if retry_outcome['status'] == 'done':
                return 'retried_ok'
            if retry_outcome['status'] == 'error':
                console.print('[red]Retry also failed.[/red]')
                continue  # ask again
            return 'retried_ok'

        if user_choice in ('s', 'skip'):
            return 'skipped'

        if user_choice in ('a', 'abort', 'q'):
            return 'aborted'

        console.print('[dim]Unrecognized option. Try r / s / a.[/dim]')


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
      'run_next' | 'run_all' | 'rerun_last' | ('jump', step_name) | 'status' | 'quit'
    """
    console.print(
        '\n[bold]Actions:[/bold]  '
        '[cyan][Enter][/cyan] run next step  '
        '[cyan][a][/cyan] run all pending  '
        '[cyan][r][/cyan] rerun last completed  '
        '[cyan][j NAME][/cyan] jump to step (prefix OK)  '
        '[cyan][s][/cyan] full status  '
        '[cyan][q][/cyan] quit'
    )
    raw_input_value = input('> ').strip()
    raw_input_lower = raw_input_value.lower()

    if raw_input_value == '':
        return 'run_next'
    if raw_input_lower in ('q', 'quit', 'exit'):
        return 'quit'
    if raw_input_lower in ('a', 'all'):
        return 'run_all'
    if raw_input_lower in ('r', 'rerun'):
        return 'rerun_last'
    if raw_input_lower in ('s', 'status'):
        return 'status'
    if raw_input_lower.startswith('j'):
        # 'j fetch' or 'jconsensus_filter'
        remaining_text = raw_input_value[1:].strip().lower()
        if not remaining_text:
            console.print('[red]Provide a step name (or unique prefix). E.g. "j consensus".[/red]')
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

    console.print(f'[dim]Unrecognized: "{raw_input_value}". Try Enter / a / r / j NAME / s / q.[/dim]')
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
        prog='theraEPIflow',
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

    # No arguments → show banner + project list + quick-start hint
    if len(sys.argv) == 1:
        command_list_projects(show_header=True)
        console.print(Panel(
            '[bold]Quick start:[/bold]\n'
            '  [cyan]python main.py --new-project[/cyan]                       Create a new project\n'
            '  [cyan]python main.py --project NAME[/cyan]                      Start interactive session\n'
            '  [cyan]python main.py --project NAME --step STEP_NAME[/cyan]     Run one specific step\n'
            '  [cyan]python main.py --project NAME --status[/cyan]             Show progress\n'
            '  [cyan]python main.py --help[/cyan]                              Full command reference',
            box=box.ROUNDED, border_style='dim',
        ))
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
