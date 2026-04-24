"""
TheraEPIflow — MHC-I Epitope Pipeline
Entry point for all pipeline operations.

Two modes of operation:

  1. Interactive session (recommended — stays open between steps):
       python main.py --project NAME
       python main.py --project NAME --run       (starts in interactive session)

  2. One-shot commands (for automation / scripting):
       python main.py --list                     list all projects
       python main.py --new-project              create a project
       python main.py --project NAME --step N    run ONE step, then exit
       python main.py --project NAME --status    show progress and exit
       python main.py --project NAME --delete    delete a project

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

# ── Step registry ──────────────────────────────────────────────────────────────
# Maps step number to (module_path, class_name, step_type)
# step_type: 'track'  = runs per track (steps 01-12)
#            'global' = runs once after all tracks (steps 13-14)

STEP_REGISTRY = {
    1:  ('modules.step01_fetch_sequences',        'FetchSequencesStep',         'track'),
    2:  ('modules.step02_validate_input',         'ValidateInputStep',          'track'),
    3:  ('modules.step03_predict_binding',        'PredictBindingStep',         'track'),
    4:  ('modules.step04_consensus_filter',       'ConsensusFilterStep',        'track'),
    5:  ('modules.step05_cluster_epitopes',       'ClusterEpitopesStep',        'track'),
    6:  ('modules.step06_select_representatives', 'SelectRepresentativesStep', 'track'),
    7:  ('modules.step07_screen_toxicity',        'ScreenToxicityStep',         'track'),
    8:  ('modules.step08_search_variants',        'SearchVariantsStep',         'track'),
    9:  ('modules.step09_analyze_conservation',   'AnalyzeConservationStep',    'track'),
    10: ('modules.step10_population_coverage',    'PopulationCoverageStep',     'track'),
    11: ('modules.step11_predict_murine',         'PredictMurineStep',          'track'),
    12: ('modules.step12_curate_murine',          'CurateMurineStep',           'track'),
    13: ('modules.step13_integrate_data',         'IntegrateDataStep',          'global'),
    14: ('modules.step14_generate_report',        'GenerateReportStep',         'global'),
}

TRACK_STEPS  = [number for number, entry in STEP_REGISTRY.items() if entry[2] == 'track']
GLOBAL_STEPS = [number for number, entry in STEP_REGISTRY.items() if entry[2] == 'global']


# ── Step loading ──────────────────────────────────────────────────────────────

def _import_step_class(step_number: int):
    """
    Dynamically imports and returns the step class for a given step number.
    Returns None if the module exists but the class is not yet implemented
    (empty __init__.py).
    """
    module_path, class_name, _ = STEP_REGISTRY[step_number]
    try:
        import importlib
        module = importlib.import_module(module_path)
        return getattr(module, class_name, None)
    except ImportError:
        return None


def _step_is_implemented(step_number: int) -> bool:
    return _import_step_class(step_number) is not None


# ── Step key helpers (pipeline.json uses step_key like "step01_fetch_sequences") ──

def _build_step_key(step_number: int) -> str:
    """Builds the pipeline.json step key for a given step number."""
    step_class = _import_step_class(step_number)
    if step_class is None:
        return f'step{step_number:02d}_unknown'
    return f'step{step_number:02d}_{step_class.step_name}'


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

    all_projects = list_projects()

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
    Builds the per-track per-step status table. Returns None if the project
    has no tracks yet.
    """
    project_config  = load_project_config(project_name)
    pipeline_state  = load_pipeline_state(project_name)
    defined_tracks  = project_config.get('tracks', {})

    if not defined_tracks:
        return None

    table = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    table.add_column('Track', style='cyan', no_wrap=True, min_width=14)

    for step_number in TRACK_STEPS:
        table.add_column(f'S{step_number:02d}', no_wrap=True, justify='center')

    for track_id in defined_tracks:
        track_pipeline_state = pipeline_state.get('tracks', {}).get(track_id, {})
        track_steps_state    = track_pipeline_state.get('steps', {})

        row_values = [track_id]
        for step_number in TRACK_STEPS:
            matching_step_status = 'pending'
            for saved_key, saved_state in track_steps_state.items():
                if saved_key.startswith(f'step{step_number:02d}'):
                    matching_step_status = saved_state.get('status', 'pending')
                    break

            if not _step_is_implemented(step_number):
                status_display = '[dim]—[/dim]'
            elif matching_step_status == 'done':
                status_display = '[green]✓[/green]'
            elif matching_step_status == 'error':
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
    step_number: int,
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
    step_class = _import_step_class(step_number)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
        track_id=track_id,
    )
    return step_instance.execute(force_rerun=force_rerun)


def _run_global_step(
    step_number: int,
    project_name: str,
    project_config: dict,
    force_rerun: bool = False,
) -> dict:
    """Same as above, but for global steps (no track_id)."""
    step_class = _import_step_class(step_number)
    if step_class is None:
        return {'status': 'not_implemented'}

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
    )
    return step_instance.execute(force_rerun=force_rerun)


# ── Interactive session state inspection ──────────────────────────────────────

def _find_next_pending_step(project_name: str) -> Optional[int]:
    """
    Returns the next step number that has at least one track with status
    other than 'done'. Returns None if everything is done (or not implemented).
    Only looks at implemented steps.
    """
    project_config   = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    if not defined_track_ids:
        return None

    for step_number in TRACK_STEPS:
        if not _step_is_implemented(step_number):
            continue
        step_key = _build_step_key(step_number)
        for track_id in defined_track_ids:
            track_status = get_track_step_status(project_name, track_id, step_key)
            if track_status != 'done':
                return step_number

    # All track steps done → check global steps
    for step_number in GLOBAL_STEPS:
        if not _step_is_implemented(step_number):
            continue
        step_key = _build_step_key(step_number)
        if get_global_step_status(project_name, step_key) != 'done':
            return step_number

    return None


def _find_last_completed_step(project_name: str) -> Optional[int]:
    """Returns the highest step number with at least one track marked 'done'."""
    project_config   = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    last_completed_step = None
    for step_number in TRACK_STEPS:
        if not _step_is_implemented(step_number):
            continue
        step_key = _build_step_key(step_number)
        for track_id in defined_track_ids:
            if get_track_step_status(project_name, track_id, step_key) == 'done':
                last_completed_step = step_number
                break
    return last_completed_step


# ── Interactive session — core loop ───────────────────────────────────────────

def _run_step_interactively(
    step_number: int,
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
    _, _, step_type = STEP_REGISTRY[step_number]

    if not _step_is_implemented(step_number):
        console.print(
            f'[yellow]Step {step_number:02d} is not implemented yet.[/yellow]'
        )
        return 'not_implemented'

    console.print(f'\n[bold cyan]══ Step {step_number:02d} ══[/bold cyan]')

    if step_type == 'global':
        outcome = _run_global_step(
            step_number=step_number,
            project_name=project_name,
            project_config=project_config,
            force_rerun=force_rerun,
        )
        if outcome['status'] == 'error':
            return _handle_step_failure(
                step_number=step_number,
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
            step_number=step_number,
            project_name=project_name,
            project_config=project_config,
            track_id=track_id,
            force_rerun=force_rerun,
        )
        if outcome['status'] == 'error':
            any_failed = True
            recovery = _handle_step_failure(
                step_number=step_number,
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
    step_number: int,
    project_name: str,
    failed_entity: str,
    error_message: str,
) -> str:
    """
    Interactive recovery prompt after a step fails.
    Returns one of: 'retried_ok', 'retried_failed', 'skipped', 'aborted'.
    """
    console.print(Panel(
        f'[bold red]Step {step_number:02d} failed[/bold red]\n'
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
        user_choice = input('> ').strip().lower()

        if user_choice in ('r', 'retry'):
            project_config = load_project_config(project_name)
            if failed_entity == 'global':
                retry_outcome = _run_global_step(
                    step_number=step_number,
                    project_name=project_name,
                    project_config=project_config,
                    force_rerun=False,
                )
            else:
                retry_outcome = _run_track_step_for_track(
                    step_number=step_number,
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
            next_step = _find_next_pending_step(project_name)
            if next_step is None:
                console.print(
                    '\n[bold green]Nothing pending — all implemented steps are done.[/bold green]'
                )
                continue
            _run_step_interactively(step_number=next_step, project_name=project_name)
            continue

        if user_choice == 'run_all':
            _run_all_pending(project_name)
            continue

        if user_choice == 'rerun_last':
            last_completed = _find_last_completed_step(project_name)
            if last_completed is None:
                console.print('[yellow]No completed step to rerun.[/yellow]')
                continue
            console.print(
                f'[dim]Rerunning step {last_completed:02d} (force).[/dim]'
            )
            _run_step_interactively(
                step_number=last_completed,
                project_name=project_name,
                force_rerun=True,
            )
            continue

        if isinstance(user_choice, tuple) and user_choice[0] == 'jump':
            target_step = user_choice[1]
            _jump_to_step(project_name=project_name, target_step=target_step)
            continue


def _print_interactive_status(project_name: str):
    """Compact header shown every iteration of the REPL."""
    project_config  = load_project_config(project_name)
    description      = project_config.get('description', '')
    next_pending    = _find_next_pending_step(project_name)

    if next_pending is None:
        next_summary = '[green]All implemented steps done[/green]'
    else:
        step_class = _import_step_class(next_pending)
        step_label = step_class.step_name if step_class else 'unknown'
        next_summary = f'Step {next_pending:02d} ([cyan]{step_label}[/cyan])'

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
      'run_next' | 'run_all' | 'rerun_last' | ('jump', N) | 'status' | 'quit'
    """
    console.print(
        '\n[bold]Actions:[/bold]  '
        '[cyan][Enter][/cyan] run next step  '
        '[cyan][a][/cyan] run all pending  '
        '[cyan][r][/cyan] rerun last completed  '
        '[cyan][j N][/cyan] jump to step N  '
        '[cyan][s][/cyan] full status  '
        '[cyan][q][/cyan] quit'
    )
    raw_input_value = input('> ').strip().lower()

    if raw_input_value == '':
        return 'run_next'
    if raw_input_value in ('q', 'quit', 'exit'):
        return 'quit'
    if raw_input_value in ('a', 'all'):
        return 'run_all'
    if raw_input_value in ('r', 'rerun'):
        return 'rerun_last'
    if raw_input_value in ('s', 'status'):
        return 'status'
    if raw_input_value.startswith('j'):
        # 'j 3' or 'j3'
        remaining_text = raw_input_value[1:].strip()
        if remaining_text.isdigit():
            target_step = int(remaining_text)
            if target_step in STEP_REGISTRY:
                return ('jump', target_step)
        console.print('[red]Invalid jump target. Use "j N" with N between 1 and 14.[/red]')
        return 'status'

    console.print(f'[dim]Unrecognized: "{raw_input_value}". Try Enter / a / r / j N / s / q.[/dim]')
    return 'status'


def _run_all_pending(project_name: str):
    """
    Runs every pending step back-to-back until all are done or a step is aborted.
    User still gets the error-recovery prompt on failures (retry/skip/abort).
    """
    while True:
        next_step = _find_next_pending_step(project_name)
        if next_step is None:
            console.print('\n[bold green]All done.[/bold green]')
            return
        outcome = _run_step_interactively(step_number=next_step, project_name=project_name)
        if outcome == 'aborted':
            console.print('[yellow]Auto-run aborted — back to menu.[/yellow]')
            return
        if outcome == 'not_implemented':
            console.print(
                f'[yellow]Reached step {next_step:02d} which is not implemented. '
                f'Stopping auto-run.[/yellow]'
            )
            return


def _jump_to_step(project_name: str, target_step: int):
    """
    Resets the target step (and all later track steps) for all tracks, so the
    pipeline effectively rewinds. Asks for confirmation — this discards state.
    """
    project_config   = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    steps_to_reset = [n for n in TRACK_STEPS if n >= target_step and _step_is_implemented(n)]
    if not steps_to_reset:
        console.print(f'[yellow]No implemented steps from {target_step:02d} onwards to reset.[/yellow]')
        return

    console.print(
        f'[bold yellow]Jump to step {target_step:02d}[/bold yellow]\n'
        f'[dim]This will clear the "done" state for steps '
        f'{", ".join(f"{n:02d}" for n in steps_to_reset)} across all tracks. '
        f'Output files are NOT deleted — only the pipeline state is reset.[/dim]'
    )
    console.print('[bold]Confirm? (y/n)[/bold]')
    confirmation_input = input('> ').strip().lower()
    if confirmation_input not in ('y', 'yes'):
        console.print('[dim]Jump cancelled.[/dim]')
        return

    for step_number in steps_to_reset:
        step_key = _build_step_key(step_number)
        for track_id in defined_track_ids:
            reset_track_step(
                project_name=project_name,
                track_id=track_id,
                step_key=step_key,
            )

    console.print(
        f'[green]Reset complete. Next run will start at step {target_step:02d}.[/green]'
    )


# ── One-shot commands (for scripting / automation) ────────────────────────────

def command_run_single_step(project_name: str, step_number: int):
    """Runs ONE step for all tracks (or once globally) and exits."""
    update_last_used(project_name)

    if not tracks_are_defined(project_name):
        console.print(
            f'[yellow]Project "{project_name}" has no tracks defined. '
            f'Run without --step first to define them.[/yellow]'
        )
        return

    _run_step_interactively(step_number=step_number, project_name=project_name)


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
    argument_parser = argparse.ArgumentParser(
        prog='theraEPIflow',
        description='MHC-I epitope pipeline for vaccine design.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  python main.py --new-project\n'
            '  python main.py --project hpv_analysis            # interactive session\n'
            '  python main.py --project hpv_analysis --step 1   # run one step, then exit\n'
            '  python main.py --project hpv_analysis --status\n'
            '  python main.py --list\n'
        ),
    )

    argument_parser.add_argument('--new-project', action='store_true',
                                 help='Create a new project (interactive wizard).')
    argument_parser.add_argument('--project', metavar='NAME',
                                 help='Name of the project to work with.')
    argument_parser.add_argument('--run', action='store_true',
                                 help='Start the interactive session (same as passing only --project).')
    argument_parser.add_argument('--step', metavar='N', type=int,
                                 help='Run a specific step number (1-14), then exit.')
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
            '  [cyan]python main.py --new-project[/cyan]               Create a new project\n'
            '  [cyan]python main.py --project NAME[/cyan]              Start interactive session\n'
            '  [cyan]python main.py --project NAME --step N[/cyan]     Run one specific step\n'
            '  [cyan]python main.py --project NAME --status[/cyan]     Show progress\n'
            '  [cyan]python main.py --help[/cyan]                      Full command reference',
            box=box.ROUNDED, border_style='dim',
        ))
        return

    if parsed_args.list:
        command_list_projects(show_header=True)
        return

    if parsed_args.new_project:
        _print_header()
        new_project_name = create_project_interactive()
        console.print(
            f'\n[dim]Run "python main.py --project {new_project_name}" '
            f'to start the interactive session.[/dim]'
        )
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
        step_number = parsed_args.step
        if step_number not in STEP_REGISTRY:
            console.print(f'[red]Step {step_number} does not exist. Valid range: 1-14.[/red]')
            sys.exit(1)
        command_run_single_step(
            project_name=selected_project_name,
            step_number=step_number,
        )
        return

    # --project NAME (with or without --run) → interactive session
    command_interactive_session(project_name=selected_project_name)


if __name__ == '__main__':
    main()
