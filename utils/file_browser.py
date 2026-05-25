"""Interactive file browser (Rich-only, no external viewer).

Two entry points share the same per-file renderers:

  - `browse_step_outputs(output_descriptions)` — the **per-step popup**
    invoked from `BaseTrackStep.execute()` and `BaseGlobalStep.execute()`
    right after a step finishes. Shows a numbered table of artifacts the
    step declared via `describe_outputs()`.

  - `run_project_browser(project_name)` — the **standalone browser** bound
    to the `[b]` REPL key. Opens a top-level chooser (Project outputs /
    Downloads / Tracks) so master tables, REPORT html and tar.gz archives
    are reachable without scripting paths by hand.

Per-file renderers dispatch on extension:
  - CSV / TSV → Rich Table (first N rows, truncated cells)
  - XLSX → Rich Table (first sheet, first N rows)
  - JSON → Rich Syntax highlighting
  - FASTA / FA / FAA / FNA → record snippets
  - HTML / HTM → summary panel + `webbrowser.open()` prompt
  - tar.gz / tgz → archive contents preview
  - everything else → plain text head
"""

from __future__ import annotations

import json
from pathlib import Path

from rich import box
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from utils.console import console, format_file_size_human

MAX_INLINE_CSV_ROWS = 15
MAX_INLINE_XLSX_ROWS = 10
MAX_INLINE_JSON_LINES = 30
MAX_INLINE_FASTA_RECORDS = 5
MAX_INLINE_FASTA_SEQ_CHARS = 80
MAX_INLINE_TEXT_LINES = 20
MAX_DISPLAYED_CELL_CHARS = 60
MAX_DISPLAYED_COLUMNS = 8


def _truncate_cell_value(raw_value: object) -> str:
    text_value = "" if raw_value is None else str(raw_value)
    if len(text_value) > MAX_DISPLAYED_CELL_CHARS:
        return text_value[: MAX_DISPLAYED_CELL_CHARS - 1] + "…"
    return text_value


def _render_csv_or_tsv_header(file_path: Path) -> None:
    import pandas as pd

    separator = "\t" if file_path.suffix.lower() == ".tsv" else ","
    try:
        dataframe_head = pd.read_csv(file_path, sep=separator, nrows=MAX_INLINE_CSV_ROWS)
    except Exception as csv_read_exception:
        console.print(f"[red]Could not read {file_path.name}: {csv_read_exception}[/red]")
        return

    total_columns = len(dataframe_head.columns)
    columns_to_show = list(dataframe_head.columns[:MAX_DISPLAYED_COLUMNS])
    has_hidden_columns = total_columns > MAX_DISPLAYED_COLUMNS
    visible_dataframe = dataframe_head[columns_to_show]

    inline_table = Table(box=box.SIMPLE, header_style="bold cyan", show_lines=False)
    for column_name in columns_to_show:
        inline_table.add_column(str(column_name), overflow="fold")
    if has_hidden_columns:
        inline_table.add_column(f"+{total_columns - MAX_DISPLAYED_COLUMNS} more", style="dim")

    for _, row_series in visible_dataframe.iterrows():
        row_cells = [_truncate_cell_value(cell_value) for cell_value in row_series.tolist()]
        if has_hidden_columns:
            row_cells.append("…")
        inline_table.add_row(*row_cells)

    if has_hidden_columns:
        hidden_column_names = list(dataframe_head.columns[MAX_DISPLAYED_COLUMNS:])
        column_summary_text = (
            f"[dim]Hidden columns ({len(hidden_column_names)}): "
            + ", ".join(hidden_column_names[:12])
            + ("…" if len(hidden_column_names) > 12 else "")
            + "[/dim]"
        )
    else:
        column_summary_text = ""

    title_text = (
        f"[cyan]{file_path.name}[/cyan]  "
        f"[dim](first {len(dataframe_head)} rows × {total_columns} cols)[/dim]"
    )

    console.print(Panel(inline_table, title=title_text, border_style="dim", box=box.ROUNDED))
    if column_summary_text:
        console.print(f"  {column_summary_text}")


def _render_xlsx_header(file_path: Path) -> None:
    import pandas as pd

    try:
        excel_file = pd.ExcelFile(file_path)
    except Exception as excel_read_exception:
        console.print(f"[red]Could not open {file_path.name}: {excel_read_exception}[/red]")
        return

    sheets_overview_table = Table(box=box.SIMPLE, header_style="bold cyan")
    sheets_overview_table.add_column("Sheet", style="cyan")
    sheets_overview_table.add_column("Rows", justify="right")
    sheets_overview_table.add_column("Cols", justify="right")
    for sheet_name in excel_file.sheet_names:
        sheet_dataframe = excel_file.parse(sheet_name)
        sheets_overview_table.add_row(sheet_name, str(len(sheet_dataframe)), str(len(sheet_dataframe.columns)))

    first_sheet_dataframe = excel_file.parse(excel_file.sheet_names[0]).head(MAX_INLINE_XLSX_ROWS)
    first_sheet_table = Table(box=box.SIMPLE, header_style="bold cyan")
    for column_name in first_sheet_dataframe.columns:
        first_sheet_table.add_column(str(column_name), overflow="fold")
    for _, row_series in first_sheet_dataframe.iterrows():
        first_sheet_table.add_row(*[_truncate_cell_value(cell_value) for cell_value in row_series.tolist()])

    console.print(
        Panel(
            sheets_overview_table,
            title=f"[cyan]{file_path.name}[/cyan]  [dim](sheets)[/dim]",
            border_style="dim",
            box=box.ROUNDED,
        )
    )
    console.print(
        Panel(
            first_sheet_table,
            title=f"[cyan]{excel_file.sheet_names[0]}[/cyan]  [dim](first {len(first_sheet_dataframe)} rows)[/dim]",
            border_style="dim",
            box=box.ROUNDED,
        )
    )


def _render_json_header(file_path: Path) -> None:
    try:
        loaded_object = json.loads(file_path.read_text(encoding="utf-8"))
    except Exception as json_read_exception:
        console.print(f"[red]Could not parse {file_path.name}: {json_read_exception}[/red]")
        return

    pretty_text = json.dumps(loaded_object, indent=2, ensure_ascii=False)
    line_list = pretty_text.splitlines()
    if len(line_list) > MAX_INLINE_JSON_LINES:
        truncated_text = "\n".join(line_list[:MAX_INLINE_JSON_LINES]) + f"\n… ({len(line_list) - MAX_INLINE_JSON_LINES} more lines)"
    else:
        truncated_text = pretty_text

    console.print(
        Panel(
            Syntax(truncated_text, "json", theme="ansi_dark", background_color="default"),
            title=f"[cyan]{file_path.name}[/cyan]",
            border_style="dim",
            box=box.ROUNDED,
        )
    )


def _render_fasta_header(file_path: Path) -> None:
    record_lines: list[Text] = []
    record_count = 0
    current_header: str | None = None
    current_sequence_chars: list[str] = []

    def flush_current_record() -> None:
        nonlocal record_count, current_header, current_sequence_chars
        if current_header is None:
            return
        sequence_string = "".join(current_sequence_chars)
        sequence_preview = sequence_string[:MAX_INLINE_FASTA_SEQ_CHARS]
        sequence_suffix = " …" if len(sequence_string) > MAX_INLINE_FASTA_SEQ_CHARS else ""
        record_lines.append(Text.from_markup(f"[bold green]>{current_header}[/bold green]"))
        record_lines.append(Text(f"  {sequence_preview}{sequence_suffix}  [{len(sequence_string)} aa]"))
        record_count += 1
        current_header = None
        current_sequence_chars = []

    with file_path.open("r", encoding="utf-8") as fasta_file_handle:
        for raw_line in fasta_file_handle:
            stripped_line = raw_line.rstrip()
            if stripped_line.startswith(">"):
                flush_current_record()
                if record_count >= MAX_INLINE_FASTA_RECORDS:
                    break
                current_header = stripped_line[1:]
                current_sequence_chars = []
            elif current_header is not None:
                current_sequence_chars.append(stripped_line)
        else:
            flush_current_record()

    body_renderable = Text("\n").join(record_lines) if record_lines else Text("(empty FASTA)")
    console.print(
        Panel(
            body_renderable,
            title=f"[cyan]{file_path.name}[/cyan]  [dim](first {record_count} records)[/dim]",
            border_style="dim",
            box=box.ROUNDED,
        )
    )


def _render_plain_text_header(file_path: Path) -> None:
    try:
        all_lines = file_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception as text_read_exception:
        console.print(f"[red]Could not read {file_path.name}: {text_read_exception}[/red]")
        return

    head_lines = all_lines[:MAX_INLINE_TEXT_LINES]
    body_text = "\n".join(head_lines)
    if len(all_lines) > MAX_INLINE_TEXT_LINES:
        body_text += f"\n… ({len(all_lines) - MAX_INLINE_TEXT_LINES} more lines)"

    syntax_language = "json" if file_path.suffix.lower() == ".json" else "text"
    console.print(
        Panel(
            Syntax(body_text, syntax_language, theme="ansi_dark", background_color="default"),
            title=f"[cyan]{file_path.name}[/cyan]",
            border_style="dim",
            box=box.ROUNDED,
        )
    )


def _render_html_header(file_path: Path) -> None:
    """For HTML files (typically the final REPORT_*.html), show a summary panel
    and offer to open the file in the user's default browser."""
    import webbrowser

    file_size_kb = file_path.stat().st_size / 1024
    summary_text = (
        f"[bold]Interactive HTML report[/bold] — {file_size_kb:.1f} KB, fully offline.\n"
        f"[dim]Path:[/dim] {file_path}\n\n"
        "Open in your default browser? [bold]\\[Y][/bold]es / [bold]\\[n][/bold]o"
    )
    console.print(Panel(
        summary_text,
        title=f"[cyan]{file_path.name}[/cyan]",
        border_style="cyan",
        box=box.ROUNDED,
        padding=(1, 2),
    ))
    try:
        user_response = console.input("  > ").strip().lower()
    except EOFError:
        return
    if user_response in ("", "y", "yes"):
        try:
            opened = webbrowser.open(file_path.resolve().as_uri())
        except Exception as browser_exception:
            console.print(f"  [yellow]Could not launch browser: {browser_exception}[/yellow]")
            console.print(f"  [dim]Open this path manually:[/dim] [cyan]{file_path.resolve()}[/cyan]")
            return
        if not opened:
            console.print("  [yellow]No browser available on this system.[/yellow]")
            console.print(f"  [dim]Open this path manually:[/dim] [cyan]{file_path.resolve()}[/cyan]")


def _render_archive_header(file_path: Path) -> None:
    """For tar.gz archives, list the top-level contents and the size."""
    import tarfile

    try:
        with tarfile.open(file_path, 'r:gz') as archive:
            archive_members = archive.getnames()
    except Exception as archive_exception:
        console.print(f"[red]Could not open archive {file_path.name}: {archive_exception}[/red]")
        return

    total_members = len(archive_members)
    members_preview = archive_members[:MAX_INLINE_TEXT_LINES]
    body_lines = [
        f"[bold]{total_members:,}[/bold] entries  ·  "
        f"{format_file_size_human(file_path.stat().st_size)}",
        "",
        *members_preview,
    ]
    if total_members > MAX_INLINE_TEXT_LINES:
        body_lines.append(f"… ({total_members - MAX_INLINE_TEXT_LINES} more entries)")

    console.print(Panel(
        "\n".join(body_lines),
        title=f"[cyan]{file_path.name}[/cyan]",
        border_style="dim",
        box=box.ROUNDED,
    ))


def _render_header_for_path(file_path: Path) -> None:
    """Dispatch to the right inline-header renderer based on file extension."""
    suffix_lower = file_path.suffix.lower()
    if suffix_lower in {".csv", ".tsv"}:
        _render_csv_or_tsv_header(file_path)
    elif suffix_lower in {".xlsx", ".xls"}:
        _render_xlsx_header(file_path)
    elif suffix_lower == ".json":
        _render_json_header(file_path)
    elif suffix_lower in {".fasta", ".fa", ".faa", ".fna"}:
        _render_fasta_header(file_path)
    elif suffix_lower in {".html", ".htm"}:
        _render_html_header(file_path)
    elif file_path.name.endswith(".tar.gz") or suffix_lower in {".tgz"}:
        _render_archive_header(file_path)
    else:
        _render_plain_text_header(file_path)


def browse_step_outputs(output_descriptions: dict[Path, str]) -> None:
    """Show a numbered listing of step artifacts and let the user peek at headers.

    Parameters
    ----------
    output_descriptions:
        Mapping of artifact path to short description (single short sentence).
        Paths that do not exist on disk are silently filtered out. If nothing
        is left to show, the function returns immediately.
    """
    existing_path_description_pairs: list[tuple[Path, str]] = [
        (artifact_path, description_text)
        for artifact_path, description_text in output_descriptions.items()
        if artifact_path.exists()
    ]
    if not existing_path_description_pairs:
        return

    listing_table = Table(box=box.SIMPLE, header_style="bold dim", pad_edge=False)
    listing_table.add_column("#", justify="right", style="bold cyan", width=3)
    listing_table.add_column("File", style="cyan", overflow="fold")
    listing_table.add_column("Description", overflow="fold")
    listing_table.add_column("Size", justify="right", style="dim")
    for one_based_index, (artifact_path, description_text) in enumerate(existing_path_description_pairs, start=1):
        listing_table.add_row(
            str(one_based_index),
            artifact_path.name,
            description_text,
            format_file_size_human(artifact_path.stat().st_size),
        )

    console.print(
        Panel(
            listing_table,
            title="[bold]Files written by this step[/bold]",
            border_style="cyan",
            box=box.ROUNDED,
            padding=(0, 2),
        )
    )

    while True:
        try:
            user_input = console.input(
                "  [dim]Type a number to preview the file header, or press Enter to continue:[/dim] "
            ).strip()
        except EOFError:
            return
        if not user_input:
            return

        try:
            chosen_index = int(user_input)
        except ValueError:
            console.print("  [yellow]Please type a number from the list, or just press Enter to continue.[/yellow]")
            continue

        if not (1 <= chosen_index <= len(existing_path_description_pairs)):
            console.print(f"  [yellow]Number must be between 1 and {len(existing_path_description_pairs)}.[/yellow]")
            continue

        chosen_path, _ = existing_path_description_pairs[chosen_index - 1]
        _render_header_for_path(chosen_path)


# ── On-demand project-wide file browser ──────────────────────────────────────


# Subfolders that BaseTrackStep subclasses write into inside
# projects/{project}/data/intermediate/{track_id}/. Ordered to match pipeline
# execution.
_INTERMEDIATE_PHASES: list[str] = [
    "predictions",
    "consensus",
    "toxicity",
    "clusters",
    "variants",
    "conservation",
    "coverage",
    "murine",
]


def _list_project_tracks(project_root: Path) -> list[str]:
    """Returns track folder names found under intermediate/, sorted."""
    intermediate_root = project_root / "data" / "intermediate"
    if not intermediate_root.exists():
        return []
    return sorted(
        candidate.name for candidate in intermediate_root.iterdir()
        if candidate.is_dir() and not candidate.name.startswith("_")
    )


def _list_phase_folders(track_dir: Path) -> list[str]:
    """Returns existing phase folder names under a track dir, in pipeline order."""
    return [phase for phase in _INTERMEDIATE_PHASES if (track_dir / phase).is_dir()]


def _list_phase_files(phase_dir: Path) -> list[Path]:
    """Returns files inside a phase folder (no subdirectories), sorted."""
    return sorted(
        candidate for candidate in phase_dir.iterdir()
        if candidate.is_file() and not candidate.name.startswith(".")
    )


def _render_numbered_menu(title: str, breadcrumb: str, options: list[str], extras: list[str]) -> None:
    """Prints a Rich Panel listing numbered options plus 'back' / 'quit'."""
    body_lines: list[str] = [f"[dim]{breadcrumb}[/dim]", ""]
    for one_based_index, option_label in enumerate(options, start=1):
        body_lines.append(f"  [bold cyan][{one_based_index}][/bold cyan] {option_label}")
    if extras:
        body_lines.append("")
        for extra_label in extras:
            # Escape any literal [..] tokens so Rich doesn't eat them as markup.
            safe_label = extra_label.replace("[", r"\[")
            body_lines.append(f"  [bold yellow]{safe_label}[/bold yellow]")
    console.print(
        Panel(
            "\n".join(body_lines),
            title=f"[bold]{title}[/bold]",
            border_style="cyan",
            box=box.HEAVY_EDGE,
            padding=(1, 2),
        )
    )


def _prompt_menu_choice(allowed_numbers: range, extra_keys: set[str]) -> str | int | None:
    """
    Returns:
      - int N when the user picks a numbered option in `allowed_numbers`
      - one of the strings in `extra_keys` (e.g. 'b', 'q') when picked
      - None when the user hits Enter / EOF / quits
    """
    hint_parts: list[str] = ["a number"]
    if "b" in extra_keys:
        hint_parts.append("'b' for back")
    if "q" in extra_keys:
        hint_parts.append("'q' to quit")
    console.print(f"  [dim]Type {', '.join(hint_parts)}:[/dim]")
    try:
        user_text = console.input("  > ").strip().lower()
    except EOFError:
        return None
    if not user_text or user_text == "q":
        return "q"
    if user_text == "b":
        if "b" in extra_keys:
            return "b"
    if user_text.isdigit():
        as_int = int(user_text)
        if as_int in allowed_numbers:
            return as_int
    console.print("  [yellow]Invalid choice. Type a number, 'b' for back, 'q' to quit.[/yellow]")
    return None


def _list_files_in(folder: Path) -> list[Path]:
    """Lists regular files in a folder, sorted, hidden files excluded."""
    if not folder.is_dir():
        return []
    return sorted(
        candidate for candidate in folder.iterdir()
        if candidate.is_file() and not candidate.name.startswith(".")
    )


def _browse_flat_folder(
    project_name: str,
    folder_label: str,
    folder_path:  Path,
) -> str | None:
    """Numbered list of files in `folder_path`, preview on selection.

    Returns 'q' if the user quit the whole browser, None to step back.
    """
    while True:
        files_in_folder = _list_files_in(folder_path)
        if not files_in_folder:
            console.print(f"  [yellow]No files in {folder_label} yet.[/yellow]")
            return None

        file_labels = [
            f"{file_path.name}  [dim]({format_file_size_human(file_path.stat().st_size)})[/dim]"
            for file_path in files_in_folder
        ]
        _render_numbered_menu(
            title=f"File browser — {folder_label}",
            breadcrumb=f"{project_name} / {folder_label} /",
            options=file_labels,
            extras=["[b] back", "[q] quit browser"],
        )
        choice = _prompt_menu_choice(range(1, len(files_in_folder) + 1), {"b", "q"})
        if choice == "q":
            return "q"
        if choice in (None, "b"):
            return None

        console.print()
        _render_header_for_path(files_in_folder[int(choice) - 1])
        console.print()


def run_project_browser(project_name: str) -> None:
    """
    Interactive on-demand browser for everything under `projects/{project_name}/`.

    Top-level chooser:
        [1] Project outputs  (data/output/ — master tables + REPORT html)
        [2] Downloads        (downloads/   — tar.gz archives, if any)
        [3] Tracks           (data/intermediate/ — per-track per-phase files)

    Inside Tracks: Track → Phase folder → Files → Header preview.
    Press 'b' to step back, 'q' (or Enter at the top) to leave.
    """
    project_root = Path("projects") / project_name
    if not project_root.exists():
        console.print(f"[red]Project not found: {project_name}[/red]")
        return

    output_dir    = project_root / "data" / "output"
    downloads_dir = project_root / "downloads"

    while True:
        # ── Level 0: top-level chooser ───────────────────────────────────────
        track_names      = _list_project_tracks(project_root)
        top_level_items: list[tuple[str, str, Path | None]] = []
        if output_dir.is_dir() and any(output_dir.iterdir()):
            top_level_items.append((
                f"Project outputs/  [dim]({len(_list_files_in(output_dir))} files — "
                "master tables, REPORT html)[/dim]",
                "outputs", output_dir,
            ))
        if downloads_dir.is_dir() and any(downloads_dir.iterdir()):
            top_level_items.append((
                f"Downloads/  [dim]({len(_list_files_in(downloads_dir))} archive(s))[/dim]",
                "downloads", downloads_dir,
            ))
        if track_names:
            top_level_items.append((
                f"Tracks (intermediate)/  [dim]({len(track_names)} tracks)[/dim]",
                "tracks", None,
            ))

        if not top_level_items:
            console.print(
                f"[yellow]No data yet for '{project_name}'. "
                "Run at least one step first.[/yellow]"
            )
            return

        _render_numbered_menu(
            title="File browser",
            breadcrumb=f"{project_name} /",
            options=[label for label, _, _ in top_level_items],
            extras=["[q] quit browser"],
        )
        section_choice = _prompt_menu_choice(range(1, len(top_level_items) + 1), {"q"})
        if section_choice in (None, "q"):
            return

        _, section_kind, section_path = top_level_items[int(section_choice) - 1]

        if section_kind == "outputs":
            if _browse_flat_folder(project_name, "outputs", section_path) == "q":
                return
            continue
        if section_kind == "downloads":
            if _browse_flat_folder(project_name, "downloads", section_path) == "q":
                return
            continue

        # ── section_kind == "tracks" — pick a track ──────────────────────────
        _render_numbered_menu(
            title="File browser — tracks",
            breadcrumb=f"{project_name} / intermediate /",
            options=[
                f"{tid}  [dim]({len(_list_phase_folders(project_root / 'data' / 'intermediate' / tid))} phase folders)[/dim]"
                for tid in track_names
            ],
            extras=["[b] back", "[q] quit browser"],
        )
        track_choice = _prompt_menu_choice(range(1, len(track_names) + 1), {"b", "q"})
        if track_choice == "q":
            return
        if track_choice in (None, "b"):
            continue

        track_id  = track_names[int(track_choice) - 1]
        track_dir = project_root / "data" / "intermediate" / track_id

        # ── Level 2: pick a phase ────────────────────────────────────────────
        while True:
            phase_folders = _list_phase_folders(track_dir)
            if not phase_folders:
                console.print(f"  [yellow]No phase data yet for {track_id}.[/yellow]")
                break

            _render_numbered_menu(
                title=f"File browser — {track_id}",
                breadcrumb=f"{project_name} / {track_id} /",
                options=[
                    f"{phase}/  [dim]({len(_list_phase_files(track_dir / phase))} files)[/dim]"
                    for phase in phase_folders
                ],
                extras=["[b] back to tracks", "[q] quit browser"],
            )
            phase_choice = _prompt_menu_choice(range(1, len(phase_folders) + 1), {"b", "q"})
            if phase_choice == "q":
                return
            if phase_choice in (None, "b"):
                break

            phase_name = phase_folders[int(phase_choice) - 1]
            phase_dir  = track_dir / phase_name

            # ── Level 3: pick a file ────────────────────────────────────────
            while True:
                phase_files = _list_phase_files(phase_dir)
                if not phase_files:
                    console.print(f"  [yellow]Empty folder: {phase_name}.[/yellow]")
                    break

                file_labels: list[str] = []
                for file_path in phase_files:
                    size_label = format_file_size_human(file_path.stat().st_size)
                    file_labels.append(f"{file_path.name}  [dim]({size_label})[/dim]")

                _render_numbered_menu(
                    title=f"File browser — {track_id} / {phase_name}",
                    breadcrumb=f"{project_name} / {track_id} / {phase_name} /",
                    options=file_labels,
                    extras=["[b] back to phases", "[q] quit browser"],
                )
                file_choice = _prompt_menu_choice(range(1, len(phase_files) + 1), {"b", "q"})
                if file_choice == "q":
                    return
                if file_choice in (None, "b"):
                    break

                chosen_file = phase_files[int(file_choice) - 1]
                console.print()
                _render_header_for_path(chosen_file)
                console.print()
