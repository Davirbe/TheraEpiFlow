"""
TheraEPIflow — MHC-I Epitope Pipeline
Entry point for all pipeline operations.

Usage:
  python main.py                              # list projects
  python main.py --new-project               # create project + define tracks
  python main.py --project NAME --run        # run all pending steps
  python main.py --project NAME --step N     # run one specific step (all tracks)
  python main.py --project NAME --status     # show step-by-step progress
  python main.py --list                      # list all projects
  python main.py --project NAME --delete     # delete a project (asks confirmation)
"""

import argparse
import sys
from pathlib import Path

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
from utils.pipeline_state import load_pipeline_state

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
# step_type: 'track' = runs per track (steps 01-12)
#            'global' = runs once after all tracks (steps 13-14)

STEP_REGISTRY = {
    1:  ('modules.step01_fetch_sequences',      'FetchSequencesStep',         'track'),
    2:  ('modules.step02_validate_input',       'ValidateInputStep',          'track'),
    3:  ('modules.step03_predict_binding',      'PredictBindingStep',         'track'),
    4:  ('modules.step04_consensus_filter',     'ConsensusFilterStep',        'track'),
    5:  ('modules.step05_cluster_epitopes',     'ClusterEpitopesStep',        'track'),
    6:  ('modules.step06_select_representatives', 'SelectRepresentativesStep','track'),
    7:  ('modules.step07_screen_toxicity',      'ScreenToxicityStep',         'track'),
    8:  ('modules.step08_search_variants',      'SearchVariantsStep',         'track'),
    9:  ('modules.step09_analyze_conservation', 'AnalyzeConservationStep',    'track'),
    10: ('modules.step10_population_coverage',  'PopulationCoverageStep',     'track'),
    11: ('modules.step11_predict_murine',       'PredictMurineStep',          'track'),
    12: ('modules.step12_curate_murine',        'CurateMurineStep',           'track'),
    13: ('modules.step13_integrate_data',       'IntegrateDataStep',          'global'),
    14: ('modules.step14_generate_report',      'GenerateReportStep',         'global'),
}

TRACK_STEPS  = [number for number, entry in STEP_REGISTRY.items() if entry[2] == 'track']
GLOBAL_STEPS = [number for number, entry in STEP_REGISTRY.items() if entry[2] == 'global']


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
        step_class = getattr(module, class_name, None)
        return step_class
    except ImportError:
        return None


# ── Project listing ────────────────────────────────────────────────────────────

def _print_header():
    """Prints the welcome banner and version."""
    console.print(WELCOME_BANNER)
    console.print(
        f'  [dim]MHC-I Epitope Pipeline for Vaccine Design  '
        f'·  v{PIPELINE_VERSION}[/dim]\n'
    )


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
        box=box.ROUNDED,
        show_header=True,
        header_style='bold white',
        title='Projects',
        title_style='bold',
    )
    table.add_column('#',           no_wrap=True,  justify='right', min_width=3)
    table.add_column('Project',     no_wrap=True,  style='cyan',   min_width=18)
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


# ── Project status ─────────────────────────────────────────────────────────────

def command_show_status(project_name: str):
    """Shows per-track, per-step status for a project."""
    project_config  = load_project_config(project_name)
    pipeline_state  = load_pipeline_state(project_name)
    defined_tracks  = project_config.get('tracks', {})

    if not defined_tracks:
        console.print(f'[yellow]Project "{project_name}" has no tracks defined yet.[/yellow]')
        console.print('[dim]Run with --run to define tracks and start the pipeline.[/dim]')
        return

    console.print(Panel.fit(
        f'[bold cyan]{project_name}[/bold cyan]  '
        f'[dim]{project_config.get("description", "")}[/dim]',
        box=box.ROUNDED,
    ))

    table = Table(box=box.SIMPLE, show_header=True, header_style='bold')
    table.add_column('Track', style='cyan', no_wrap=True, min_width=14)

    # Add one column per step
    for step_number in TRACK_STEPS:
        table.add_column(f'S{step_number:02d}', no_wrap=True, justify='center')

    for track_id in defined_tracks:
        track_pipeline_state = pipeline_state.get('tracks', {}).get(track_id, {})
        track_steps_state    = track_pipeline_state.get('steps', {})

        row_values = [track_id]
        for step_number in TRACK_STEPS:
            module_path, class_name, _ = STEP_REGISTRY[step_number]
            step_key = f'step{step_number:02d}_{class_name.replace("Step", "").lower()}'

            # Find the actual step key in state (key format may vary)
            matching_step_status = 'pending'
            for saved_key, saved_state in track_steps_state.items():
                if saved_key.startswith(f'step{step_number:02d}'):
                    matching_step_status = saved_state.get('status', 'pending')
                    break

            step_class_exists = _import_step_class(step_number) is not None
            if not step_class_exists:
                status_display = '[dim]—[/dim]'
            elif matching_step_status == 'done':
                status_display = '[green]✓[/green]'
            elif matching_step_status == 'error':
                status_display = '[red]✗[/red]'
            else:
                status_display = '[dim]·[/dim]'

            row_values.append(status_display)

        table.add_row(*row_values)

    console.print(table)
    console.print('[dim]✓ done  · pending  ✗ error  — not implemented[/dim]')


# ── Step runner ────────────────────────────────────────────────────────────────

def _run_track_step(
    step_number: int,
    project_name: str,
    project_config: dict,
    track_id: str,
) -> bool:
    """
    Instantiates and executes a single track step.
    Returns True if the step completed successfully, False otherwise.
    """
    step_class = _import_step_class(step_number)
    if step_class is None:
        console.print(
            f'[dim]Step {step_number:02d} not yet implemented — skipping.[/dim]'
        )
        return False

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
        track_id=track_id,
    )
    try:
        step_instance.execute()
        return True
    except Exception as step_error:
        console.print(f'[bold red]Step {step_number:02d} failed for {track_id}: {step_error}[/bold red]')
        return False


def _run_global_step(
    step_number: int,
    project_name: str,
    project_config: dict,
) -> bool:
    """
    Instantiates and executes a single global step (runs once for the whole project).
    Returns True if the step completed successfully, False otherwise.
    """
    step_class = _import_step_class(step_number)
    if step_class is None:
        console.print(
            f'[dim]Step {step_number:02d} not yet implemented — skipping.[/dim]'
        )
        return False

    step_instance = step_class(
        project_name=project_name,
        project_config=project_config,
    )
    try:
        step_instance.execute()
        return True
    except Exception as step_error:
        console.print(f'[bold red]Step {step_number:02d} failed: {step_error}[/bold red]')
        return False


def command_run_project(project_name: str, only_step: int = None):
    """
    Runs all pending steps for a project, or one specific step if only_step is set.

    If no tracks are defined, calls setup_project_tracks_interactive() first.
    Track steps (01-12) run per track; global steps (13-14) run once.
    Stops on error and reports which step failed.
    """
    update_last_used(project_name)

    # If tracks not yet defined, ask now (beginning of step01 flow)
    if not tracks_are_defined(project_name):
        console.print(
            f'\n[bold yellow]Project "{project_name}" has no tracks defined.[/bold yellow]'
        )
        console.print('[dim]Defining tracks now before running step01...[/dim]\n')
        setup_project_tracks_interactive(project_name)

    # Reload config after possible track setup
    project_config = load_project_config(project_name)
    defined_track_ids = list(project_config.get('tracks', {}).keys())

    if not defined_track_ids:
        console.print('[red]No tracks defined. Cannot run.[/red]')
        return

    steps_to_run_track  = [only_step] if only_step and only_step in TRACK_STEPS \
                           else TRACK_STEPS
    steps_to_run_global = [only_step] if only_step and only_step in GLOBAL_STEPS \
                          else GLOBAL_STEPS

    # ── Per-track steps (01-12) ───────────────────────────────────────────────
    for step_number in steps_to_run_track:
        console.print(
            f'\n[bold cyan]══ Step {step_number:02d} ══[/bold cyan]'
        )
        for track_id in defined_track_ids:
            step_succeeded = _run_track_step(
                step_number=step_number,
                project_name=project_name,
                project_config=project_config,
                track_id=track_id,
            )
            if not step_succeeded and _import_step_class(step_number) is not None:
                console.print(
                    f'[red]Stopping at step {step_number:02d} — '
                    f'fix the error before continuing.[/red]'
                )
                return

        # Reload config after each step in case a step updated it
        project_config = load_project_config(project_name)

    # ── Global steps (13-14) ──────────────────────────────────────────────────
    if not only_step or only_step in GLOBAL_STEPS:
        for step_number in steps_to_run_global:
            console.print(f'\n[bold cyan]══ Step {step_number:02d} (global) ══[/bold cyan]')
            step_succeeded = _run_global_step(
                step_number=step_number,
                project_name=project_name,
                project_config=project_config,
            )
            if not step_succeeded and _import_step_class(step_number) is not None:
                console.print(f'[red]Stopping at step {step_number:02d}.[/red]')
                return
            project_config = load_project_config(project_name)

    console.print('\n[bold green]Run complete.[/bold green]')


# ── Delete ─────────────────────────────────────────────────────────────────────

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


# ── CLI entry point ────────────────────────────────────────────────────────────

def main():
    argument_parser = argparse.ArgumentParser(
        prog='theraEPIflow',
        description='MHC-I epitope pipeline for vaccine design.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  python main.py --new-project\n'
            '  python main.py --project hpv_analysis --run\n'
            '  python main.py --project hpv_analysis --step 1\n'
            '  python main.py --project hpv_analysis --status\n'
            '  python main.py --list\n'
        ),
    )

    argument_parser.add_argument(
        '--new-project',
        action='store_true',
        help='Create a new project (interactive wizard).',
    )
    argument_parser.add_argument(
        '--project',
        metavar='NAME',
        help='Name of the project to work with.',
    )
    argument_parser.add_argument(
        '--run',
        action='store_true',
        help='Run all pending steps for the project.',
    )
    argument_parser.add_argument(
        '--step',
        metavar='N',
        type=int,
        help='Run a specific step number (1-14) for the project.',
    )
    argument_parser.add_argument(
        '--status',
        action='store_true',
        help='Show step-by-step progress for the project.',
    )
    argument_parser.add_argument(
        '--list',
        action='store_true',
        help='List all projects.',
    )
    argument_parser.add_argument(
        '--delete',
        action='store_true',
        help='Delete the project (requires --project).',
    )

    parsed_args = argument_parser.parse_args()

    # No arguments → show banner + project list + quick-start hint
    if len(sys.argv) == 1:
        command_list_projects(show_header=True)
        console.print(Panel(
            '[bold]Quick start:[/bold]\n'
            '  [cyan]python main.py --new-project[/cyan]               Create a new project\n'
            '  [cyan]python main.py --project NAME --run[/cyan]        Run all pending steps\n'
            '  [cyan]python main.py --project NAME --step N[/cyan]     Run a specific step\n'
            '  [cyan]python main.py --project NAME --status[/cyan]     Show progress\n'
            '  [cyan]python main.py --help[/cyan]                      Full command reference',
            box=box.ROUNDED,
            border_style='dim',
        ))
        return

    if parsed_args.list:
        command_list_projects(show_header=True)
        return

    if parsed_args.new_project:
        _print_header()
        new_project_name = create_project_interactive()
        console.print(
            f'\n[dim]Run "python main.py --project {new_project_name} --run" '
            f'to start the pipeline.[/dim]'
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

    if parsed_args.run:
        command_run_project(project_name=selected_project_name)
        return

    if parsed_args.step:
        step_number = parsed_args.step
        if step_number not in STEP_REGISTRY:
            console.print(f'[red]Step {step_number} does not exist. Valid range: 1-14.[/red]')
            sys.exit(1)
        command_run_project(
            project_name=selected_project_name,
            only_step=step_number,
        )
        return

    # --project NAME with no action → show status
    command_show_status(selected_project_name)


if __name__ == '__main__':
    main()
