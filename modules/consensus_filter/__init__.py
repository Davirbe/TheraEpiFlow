"""
Step "consensus_filter".

Receives prediction CSVs from the "predict_binding" step and applies:
  Stage 1 — presentation filter (percentile threshold + intersection between tools)
  Stage 2 — immunogenicity filter (Calis 2013, score > 0)

Inputs (track_dir/predictions/):
  PRED_NET_{track_id}.csv
  PRED_FLURRY_{track_id}.csv

Outputs in track_dir/consensus/:
  0a_PRED_NET_no_nan.csv            rows after NaN removal (audit)
  0a_PRED_FLURRY_no_nan.csv
  0b_PRED_NET_thresholded.csv       after percentile cutoff (audit)
  0b_PRED_FLURRY_thresholded.csv
  1_PRED_NET_consolidated.csv       1 row per peptide (audit)
  1_PRED_FLURRY_consolidated.csv
  2_intersection.csv                peptides present in both tools (audit)
  CONSENSUS_{track_id}.csv          final table (prefix-only column names)
  CONSENSUS_IMMUNOGENIC_{track_id}.csv  Calis survivors (score > 0)
  consensus_audit_summary.json      per-stage counts
"""

import contextlib
import datetime
import io
import json
from pathlib import Path

import pandas as pd
from rich import box
from utils.console import console
from rich.panel import Panel
from rich.table import Table
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from modules.base_step import BaseTrackStep
from utils.project_manager import save_project_config
from utils.naming import get_prediction_filename, get_step_filename


# Substrings used to locate the percentile column in each tool's output
_NET_PERCENTILE_SUBSTRINGS    = ['netmhcpan_el', 'percentile']
_FLURRY_PERCENTILE_SUBSTRINGS = ['mhcflurry', 'presentation', 'percentile']


# ── Threshold ─────────────────────────────────────────────────────────────────

def _ask_threshold(project_name: str, project_config: dict) -> float:
    """
    Returns the consensus threshold from project_config if already saved,
    otherwise prompts the user once and saves it.
    """
    existing_threshold = project_config.get('consensus_threshold')
    if existing_threshold is not None:
        return float(existing_threshold)

    console.print(Panel.fit(
        '[bold cyan]Consensus threshold[/bold cyan]\n'
        '[dim]Keeps peptides with percentile rank ≤ threshold in BOTH tools[/dim]\n\n'
        '  [cyan][1][/cyan] Strong + weak binders  (≤ 2.0)  [dim]default[/dim]\n'
        '  [cyan][2][/cyan] Strong binders only    (≤ 0.5)\n'
        '  [cyan][3][/cyan] Custom value',
        box=box.ROUNDED, border_style='cyan',
    ))

    while True:
        try:
            user_choice = input('> ').strip()
        except EOFError:
            user_choice = '1'
        if user_choice in ('', '1'):
            chosen_threshold = 2.0
            break
        if user_choice == '2':
            chosen_threshold = 0.5
            break
        if user_choice == '3':
            try:
                raw_threshold = input('Enter threshold (e.g. 1.5): ').strip().replace(',', '.')
            except EOFError:
                raw_threshold = '2.0'
            try:
                chosen_threshold = float(raw_threshold)
                if chosen_threshold <= 0:
                    raise ValueError
                break
            except ValueError:
                console.print('[red]Invalid value — enter a number greater than 0.[/red]')
                continue
        console.print('[red]Choose 1, 2 or 3.[/red]')

    project_config['consensus_threshold'] = chosen_threshold
    save_project_config(project_name, project_config)
    console.print(f'[dim]Threshold set to ≤ {chosen_threshold} — saved to project_config.[/dim]')
    return chosen_threshold


# ── DataFrame helpers ──────────────────────────────────────────────────────────

def _find_col(dataframe: pd.DataFrame, substrings: list) -> str | None:
    """Returns the first column name that contains all given substrings."""
    for column_name in dataframe.columns:
        if all(substring in column_name for substring in substrings):
            return column_name
    return None


def _normalize_col_name(raw_name: str) -> str:
    """Lowercase + underscores. 'seq #' → 'sequence_number'."""
    cleaned = raw_name.strip()
    if cleaned.lower() == 'seq #':
        return 'sequence_number'
    return cleaned.replace(' ', '_').replace('#', 'num').lower()


# ── Stage 1 — Presentation filter ─────────────────────────────────────────────

def _load_and_filter(csv_path: Path, percentile_substrings: list, threshold: float) -> dict:
    """
    Loads a prediction CSV and applies two sub-phases:
      0a — drops rows with NaN in peptide, allele, or percentile
      0b — keeps only rows with percentile <= threshold

    Returns a dict with both intermediate DataFrames and audit counts.
    """
    dataframe = pd.read_csv(csv_path, sep=',', dtype=str, keep_default_na=False,
                            na_values=['', 'NA', 'N/A'], engine='python')
    if dataframe.shape[1] == 1:
        dataframe = pd.read_csv(csv_path, sep=';', dtype=str, keep_default_na=False,
                                na_values=['', 'NA', 'N/A'], engine='python')

    dataframe.columns = [_normalize_col_name(col) for col in dataframe.columns]
    if 'peptide' in dataframe.columns:
        dataframe['peptide'] = dataframe['peptide'].str.strip()

    numeric_columns = [col for col in dataframe.columns
                       if any(keyword in col for keyword in ('percentile', 'score', 'ic50', 'rank', 'affinity'))]
    for numeric_col in numeric_columns:
        dataframe[numeric_col] = pd.to_numeric(
            dataframe[numeric_col].str.replace(',', '.', regex=False), errors='coerce'
        )

    percentile_column = _find_col(dataframe, percentile_substrings)
    if percentile_column is None:
        raise ValueError(
            f'Percentile column not found in {csv_path}. '
            f'Expected substrings: {percentile_substrings}. '
            f'Available columns: {list(dataframe.columns)}'
        )

    # 0a — drop NaN
    n_raw       = len(dataframe)
    df_0a       = dataframe.dropna(subset=['peptide', 'allele', percentile_column]).copy()
    n_after_nan = len(df_0a)

    # 0b — apply percentile threshold
    df_0b               = df_0a[df_0a[percentile_column] <= threshold].copy()
    n_after_threshold   = len(df_0b)

    # Reference counts for comparison display (strong and weak binder breakpoints)
    n_strong_binders = int((df_0a[percentile_column] <= 0.5).sum())
    n_weak_binders   = int((df_0a[percentile_column] <= 2.0).sum())

    return {
        'df_0a':       df_0a,
        'df_0b':       df_0b,
        'pct_col':     percentile_column,
        'n_raw':       n_raw,
        'n_0a':        n_after_nan,
        'dropped_nan': n_raw - n_after_nan,
        'n_0b':        n_after_threshold,
        'dropped_thr': n_after_nan - n_after_threshold,
        'n_strong':    n_strong_binders,
        'n_weak':      n_weak_binders,
    }


def _consolidate(
    df: pd.DataFrame,
    pct_col: str,
    tool_prefix: str,
    metric_name: str,
) -> pd.DataFrame:
    """
    Phase 1 — collapses N rows per allele into 1 row per peptide.

    Produces a DataFrame with only the columns needed for consensus,
    using the prefix-only naming convention (no redundant tool suffix):

      peptide
      {tool_prefix}_best_allele
      {tool_prefix}_{metric_name}_percentile       <- best (minimum) across alleles
      {tool_prefix}_alleles                        <- all alleles alphabetically, semicolon-separated
      {tool_prefix}_{metric_name}_percentiles_all  <- percentile for each allele, same order
      {tool_prefix}_num_alleles

    Args:
        tool_prefix:  "netmhcpan" or "mhcflurry"
        metric_name:  "el" for NetMHCpan, "presentation" for MHCFlurry
    """
    if df.empty:
        return pd.DataFrame()

    aggregated_rows = []

    for peptide_sequence, allele_group in df.groupby('peptide'):
        best_percentile_per_allele = (
            allele_group.groupby('allele')[pct_col]
            .min()
            .sort_index()
        )
        alleles_joined     = ';'.join(best_percentile_per_allele.index.tolist())
        percentiles_joined = ';'.join(f"{value:.4g}" for value in best_percentile_per_allele.values)
        number_of_alleles  = len(best_percentile_per_allele)

        best_row_index  = allele_group[pct_col].idxmin()
        best_percentile = float(allele_group.loc[best_row_index, pct_col])
        best_allele     = str(allele_group.loc[best_row_index, 'allele'])

        aggregated_rows.append({
            'peptide':                                          peptide_sequence,
            f'{tool_prefix}_best_allele':                      best_allele,
            f'{tool_prefix}_{metric_name}_percentile':         best_percentile,
            f'{tool_prefix}_alleles':                          alleles_joined,
            f'{tool_prefix}_{metric_name}_percentiles_all':    percentiles_joined,
            f'{tool_prefix}_num_alleles':                      number_of_alleles,
        })

    return pd.DataFrame(aggregated_rows)


def _intersect(net_consolidated: pd.DataFrame, flurry_consolidated: pd.DataFrame) -> dict:
    """
    Phase 2 — keeps only peptides present in BOTH tools.
    Returns a dict with the common-peptide DataFrame and audit counts.
    """
    if net_consolidated.empty or flurry_consolidated.empty:
        return {
            'common':       pd.DataFrame(columns=['peptide']),
            'net_only':     len(net_consolidated),
            'flurry_only':  len(flurry_consolidated),
            'common_count': 0,
        }

    peptides_in_net    = set(net_consolidated['peptide'])
    peptides_in_flurry = set(flurry_consolidated['peptide'])
    peptides_in_both   = peptides_in_net & peptides_in_flurry

    return {
        'common':       pd.DataFrame({'peptide': sorted(peptides_in_both)}),
        'net_only':     len(peptides_in_net - peptides_in_flurry),
        'flurry_only':  len(peptides_in_flurry - peptides_in_net),
        'common_count': len(peptides_in_both),
    }


def _finalize(
    common_peptides: pd.DataFrame,
    net_consolidated: pd.DataFrame,
    flurry_consolidated: pd.DataFrame,
) -> pd.DataFrame:
    """
    Phase 3 — assembles the final table.

    Columns are already named correctly by _consolidate() (prefix-only convention,
    no redundant suffixes). Filters each tool to common peptides and merges.
    """
    if common_peptides.empty:
        return pd.DataFrame(columns=['peptide'])

    net_filtered    = net_consolidated[net_consolidated['peptide'].isin(common_peptides['peptide'])].copy()
    flurry_filtered = flurry_consolidated[flurry_consolidated['peptide'].isin(common_peptides['peptide'])].copy()

    merged = pd.merge(common_peptides, net_filtered,    on='peptide', how='left')
    merged = pd.merge(merged,          flurry_filtered, on='peptide', how='left')
    return merged


# ── Stage 2 — Calis 2013 ─────────────────────────────────────────────────────

# Allele → anchor positions (1-based), reproduced from predict_immunogenicity.py
# to resolve the mask without calling validate() (which uses sys.exit and reads a file).
_CALIS_ALLELES = {
    "H-2-Db":"2,5,9","H-2-Dd":"2,3,5","H-2-Kb":"2,3,9","H-2-Kd":"2,5,9",
    "H-2-Kk":"2,8,9","H-2-Ld":"2,5,9","HLA-A0101":"2,3,9","HLA-A0201":"1,2,9",
    "HLA-A0202":"1,2,9","HLA-A0203":"1,2,9","HLA-A0206":"1,2,9","HLA-A0211":"1,2,9",
    "HLA-A0301":"1,2,9","HLA-A1101":"1,2,9","HLA-A2301":"2,7,9","HLA-A2402":"2,7,9",
    "HLA-A2601":"1,2,9","HLA-A2902":"2,7,9","HLA-A3001":"1,3,9","HLA-A3002":"2,7,9",
    "HLA-A3101":"1,2,9","HLA-A3201":"1,2,9","HLA-A3301":"1,2,9","HLA-A6801":"1,2,9",
    "HLA-A6802":"1,2,9","HLA-A6901":"1,2,9","HLA-B0702":"1,2,9","HLA-B0801":"2,5,9",
    "HLA-B1501":"1,2,9","HLA-B1502":"1,2,9","HLA-B1801":"1,2,9","HLA-B2705":"2,3,9",
    "HLA-B3501":"1,2,9","HLA-B3901":"1,2,9","HLA-B4001":"1,2,9","HLA-B4002":"1,2,9",
    "HLA-B4402":"2,3,9","HLA-B4403":"2,3,9","HLA-B4501":"1,2,9","HLA-B4601":"1,2,9",
    "HLA-B5101":"1,2,9","HLA-B5301":"1,2,9","HLA-B5401":"1,2,9","HLA-B5701":"1,2,9",
    "HLA-B5801":"1,2,9",
}


def _score_calis(peptide_list: list, allele_imgt: str | None) -> dict:
    """
    Scores a list of peptides via Calis 2013 (predict_immunogenicity.Prediction).
    Captures stdout because Prediction.predict() prints instead of returning.

    Args:
        peptide_list: list of uppercase amino acid sequences
        allele_imgt:  IMGT format e.g. 'HLA-A*02:01', or None → default mask

    Returns dict {peptide: score}.
    """
    from .predict_immunogenicity import Prediction

    allele_calis = None
    anchor_mask  = None
    if allele_imgt:
        allele_normalized = allele_imgt.replace('*', '').replace(':', '')
        if allele_normalized in _CALIS_ALLELES:
            allele_calis = allele_normalized
            anchor_mask  = _CALIS_ALLELES[allele_normalized]

    captured_output = io.StringIO()
    with contextlib.redirect_stdout(captured_output):
        Prediction().predict((list(peptide_list), anchor_mask, allele_calis))

    # Parse output lines with format "PEPTIDE,length,score"
    peptide_scores = {}
    for output_line in captured_output.getvalue().splitlines():
        output_line = output_line.strip()
        if not output_line or ',' not in output_line:
            continue
        if output_line.startswith(('peptide,', 'masking:', 'masked', 'allele:')):
            continue
        line_parts = output_line.split(',')
        if len(line_parts) == 3:
            try:
                peptide_scores[line_parts[0]] = float(line_parts[2])
            except ValueError:
                continue
    return peptide_scores


def _apply_calis(consensus_df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Applies Calis 2013 grouped by netmhcpan_best_allele.
    Returns (DataFrame filtered to score > 0, audit dict).
    """
    if consensus_df.empty:
        return consensus_df.assign(calis_score=pd.Series(dtype=float)), {
            'input': 0, 'scored': 0, 'survived': 0, 'unsupported_allele': 0,
        }

    peptide_score_map  = {}
    peptide_allele_map = {}
    unsupported_count  = 0

    calis_groups_by_allele = list(
        consensus_df.groupby('netmhcpan_best_allele', dropna=False)
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as calis_progress_bar:
        calis_task_id = calis_progress_bar.add_task(
            "  Scoring Calis 2013 immunogenicity (per allele)",
            total=len(calis_groups_by_allele),
        )
        for allele_name, allele_group_df in calis_groups_by_allele:
            peptides_in_allele_group = allele_group_df['peptide'].tolist()
            if pd.isna(allele_name) or not allele_name:
                calis_score_result      = _score_calis(peptides_in_allele_group, None)
                allele_label_for_audit  = 'default_mask'
                unsupported_count      += len(peptides_in_allele_group)
            else:
                allele_normalized = allele_name.replace('*', '').replace(':', '')
                if allele_normalized not in _CALIS_ALLELES:
                    calis_score_result      = _score_calis(peptides_in_allele_group, None)
                    allele_label_for_audit  = f'{allele_name} (default_mask)'
                    unsupported_count      += len(peptides_in_allele_group)
                else:
                    calis_score_result      = _score_calis(peptides_in_allele_group, allele_name)
                    allele_label_for_audit  = allele_name

            for peptide_sequence in peptides_in_allele_group:
                peptide_score_map[peptide_sequence]  = calis_score_result.get(peptide_sequence)
                peptide_allele_map[peptide_sequence] = allele_label_for_audit

            calis_progress_bar.advance(calis_task_id)

    scored_df = consensus_df.copy()
    scored_df['calis_score']       = scored_df['peptide'].map(peptide_score_map)
    scored_df['calis_allele_used'] = scored_df['peptide'].map(peptide_allele_map)

    immunogenic_df = scored_df[
        scored_df['calis_score'].notna() & (scored_df['calis_score'] > 0)
    ].copy()

    return immunogenic_df, {
        'input':               len(consensus_df),
        'scored':              int(scored_df['calis_score'].notna().sum()),
        'survived':            len(immunogenic_df),
        'unsupported_allele':  unsupported_count,
    }


# ── Progressive display ───────────────────────────────────────────────────────

def _print_stage1_filtering(net_data: dict, flu_data: dict, threshold: float):
    """Stage 1 — loading, NaN removal, and threshold with a cutoff comparison."""
    display_table = Table(
        box=box.SIMPLE, show_header=True, header_style='bold',
        title='[bold cyan]Stage 1 — Presentation threshold filter[/bold cyan]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('NetMHCpan', justify='right')
    display_table.add_column('MHCFlurry', justify='right')

    display_table.add_row('Prediction rows (input)',  str(net_data['n_raw']),        str(flu_data['n_raw']))
    display_table.add_row('  removed (NaN)',          f'-{net_data["dropped_nan"]}', f'-{flu_data["dropped_nan"]}')
    display_table.add_row('  remaining',              str(net_data['n_0a']),         str(flu_data['n_0a']))
    display_table.add_section()

    def _mark_active(cutoff):
        return ' [bold yellow]★[/bold yellow]' if cutoff == threshold else ''

    display_table.add_row(
        f'≤ 0.5  (strong binders only){_mark_active(0.5)}',
        str(net_data['n_strong']), str(flu_data['n_strong'])
    )
    display_table.add_row(
        f'≤ 2.0  (strong + weak binders){_mark_active(2.0)}',
        str(net_data['n_weak']), str(flu_data['n_weak'])
    )
    if threshold not in (0.5, 2.0):
        display_table.add_row(
            f'≤ {threshold}  (custom) [bold yellow]★[/bold yellow]',
            str(net_data['n_0b']), str(flu_data['n_0b'])
        )

    console.print(display_table)


def _print_stage2_consolidation(
    net_data: dict,
    flu_data: dict,
    net_consolidated: pd.DataFrame,
    flurry_consolidated: pd.DataFrame,
):
    """Stage 2 — consolidation of allele rows into unique peptides."""
    n_net_peptides    = len(net_consolidated)
    n_flurry_peptides = len(flurry_consolidated)

    display_table = Table(
        box=box.SIMPLE, show_header=True, header_style='bold',
        title='[bold green]Stage 2 — Consolidation (allele×peptide rows → unique peptides)[/bold green]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('NetMHCpan', justify='right')
    display_table.add_column('MHCFlurry', justify='right')

    display_table.add_row('Rows after threshold',     str(net_data['n_0b']),                     str(flu_data['n_0b']))
    display_table.add_row('Unique peptides',          str(n_net_peptides),                       str(n_flurry_peptides))
    display_table.add_row('Collapsed rows',           str(net_data['n_0b'] - n_net_peptides),    str(flu_data['n_0b'] - n_flurry_peptides))

    console.print(display_table)


def _print_stage3_intersection(intersection_data: dict):
    """Stage 3 — intersection between NetMHCpan and MHCFlurry."""
    display_table = Table(
        box=box.SIMPLE, show_header=False,
        title='[bold yellow]Stage 3 — Tool intersection[/bold yellow]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('', justify='right')

    display_table.add_row('NetMHCpan only',  str(intersection_data['net_only']))
    display_table.add_row('MHCFlurry only',  str(intersection_data['flurry_only']))
    display_table.add_row(
        'In consensus (both tools) [bold yellow]★[/bold yellow]',
        f'[bold green]{intersection_data["common_count"]}[/bold green]'
    )

    console.print(display_table)


def _print_stage4_immunogenicity(n_input: int, n_survivors: int):
    """Stage 4 — result of the Calis 2013 immunogenicity filter."""
    n_discarded  = n_input - n_survivors
    pct_survived = (n_survivors / n_input * 100) if n_input > 0 else 0.0
    pct_discarded = 100.0 - pct_survived

    display_table = Table(
        box=box.SIMPLE, show_header=False,
        title='[bold magenta]Stage 4 — Calis 2013 immunogenicity (score > 0)[/bold magenta]',
    )
    display_table.add_column('', style='bold cyan', no_wrap=True)
    display_table.add_column('', justify='right')

    display_table.add_row('Input (consensus)',  str(n_input))
    display_table.add_row(
        'Survivors',
        f'[bold green]{n_survivors}[/bold green]  ({pct_survived:.0f}%)'
    )
    display_table.add_row(
        'Discarded',
        f'[dim]{n_discarded}[/dim]  ({pct_discarded:.0f}%)'
    )

    console.print(display_table)
    console.print(f'\n[bold green]FINAL RESULT: {n_survivors} immunogenic epitopes[/bold green]\n')


# ── Step ──────────────────────────────────────────────────────────────────────

class ConsensusFilterStep(BaseTrackStep):
    step_name   = 'consensus_filter'
    description = (
        "Keeps only peptides predicted as binders by BOTH NetMHCpan and "
        "MHCFlurry (percentile threshold), consolidates one row per peptide "
        "with the joined HLA list, then drops anything Calis 2013 "
        "immunogenicity scores at zero or below."
    )

    def describe_outputs(self) -> dict:
        consensus_dir = self.track_dir / 'consensus'
        return {
            consensus_dir / get_step_filename("CONSENSUS_VIEW", self.track_id):
                "Slim per-step view — peptide + best allele/percentile per tool + Calis score.",
            consensus_dir / get_step_filename("CONSENSUS_IMMUNOGENIC", self.track_id):
                "Final consensus output — peptides that passed BOTH NetMHCpan and MHCFlurry "
                "AND have a positive Calis 2013 immunogenicity score. Feeds screen_toxicity / cluster_epitopes.",
            consensus_dir / get_step_filename("CONSENSUS", self.track_id):
                "Same peptide set BEFORE the Calis immunogenicity filter (Stage 1 only).",
            consensus_dir / 'consensus_audit_summary.json':
                "Funnel audit — counts at each stage of the consensus pipeline.",
        }

    def run(self, input_data=None):
        # 1. Threshold
        threshold = _ask_threshold(self.project_name, self.project_config)

        # 2. Locate prediction CSVs from the previous step
        pred_dir = self.track_dir / 'predictions'
        net_csv  = pred_dir / get_prediction_filename("NET_PRED",   self.track_id)
        flu_csv  = pred_dir / get_prediction_filename("FLURRY_PRED", self.track_id)
        for csv_path in (net_csv, flu_csv):
            if not csv_path.exists():
                raise FileNotFoundError(
                    f'{csv_path} not found. Run "predict_binding" first.'
                )

        output_dir = self.track_dir / 'consensus'
        output_dir.mkdir(parents=True, exist_ok=True)

        console.print(Panel.fit(
            f'[bold cyan]Consensus filter — {self.track_id}[/bold cyan]\n'
            f'[dim]Threshold:[/dim] [bold]≤ {threshold}[/bold]',
            box=box.ROUNDED, border_style='cyan',
        ))

        # 3. Phase 0 — load + drop NaN + apply threshold
        net_data = _load_and_filter(net_csv, _NET_PERCENTILE_SUBSTRINGS,    threshold)
        flu_data = _load_and_filter(flu_csv, _FLURRY_PERCENTILE_SUBSTRINGS, threshold)

        net_data['df_0a'].to_csv(output_dir / '0a_PRED_NET_no_nan.csv',         index=False, sep=';', decimal=',')
        net_data['df_0b'].to_csv(output_dir / '0b_PRED_NET_thresholded.csv',    index=False, sep=';', decimal=',')
        flu_data['df_0a'].to_csv(output_dir / '0a_PRED_FLURRY_no_nan.csv',      index=False, sep=';', decimal=',')
        flu_data['df_0b'].to_csv(output_dir / '0b_PRED_FLURRY_thresholded.csv', index=False, sep=';', decimal=',')

        _print_stage1_filtering(net_data, flu_data, threshold)

        # 4. Phase 1 — consolidate per peptide
        net_consolidated    = _consolidate(net_data['df_0b'], net_data['pct_col'], 'netmhcpan', 'el')
        flurry_consolidated = _consolidate(flu_data['df_0b'], flu_data['pct_col'], 'mhcflurry', 'presentation')

        net_consolidated.to_csv(output_dir / '1_PRED_NET_consolidated.csv',    index=False, sep=';', decimal=',')
        flurry_consolidated.to_csv(output_dir / '1_PRED_FLURRY_consolidated.csv', index=False, sep=';', decimal=',')

        _print_stage2_consolidation(net_data, flu_data, net_consolidated, flurry_consolidated)

        # 5. Phase 2 — intersect between tools
        intersection = _intersect(net_consolidated, flurry_consolidated)
        intersection['common'].to_csv(output_dir / '2_intersection.csv', index=False, sep=';', decimal=',')

        _print_stage3_intersection(intersection)

        # 6. Phase 3 — build final table
        consensus_df  = _finalize(intersection['common'], net_consolidated, flurry_consolidated)
        consensus_csv = output_dir / get_step_filename("CONSENSUS", self.track_id)
        consensus_df.to_csv(consensus_csv, index=False)

        # 7. Stage 2 — Calis 2013 immunogenicity filter
        immunogenic_df, calis_audit = _apply_calis(consensus_df)
        immunogenic_csv = output_dir / get_step_filename("CONSENSUS_IMMUNOGENIC", self.track_id)
        immunogenic_df.to_csv(immunogenic_csv, index=False)

        view_columns = [
            'peptide',
            'netmhcpan_best_allele',
            'netmhcpan_el_percentile',
            'mhcflurry_best_allele',
            'mhcflurry_presentation_percentile',
            'calis_score',
        ]
        view_present_columns = [c for c in view_columns if c in immunogenic_df.columns]
        view_csv = output_dir / get_step_filename("CONSENSUS_VIEW", self.track_id)
        immunogenic_df[view_present_columns].to_csv(view_csv, index=False)

        _print_stage4_immunogenicity(len(consensus_df), len(immunogenic_df))

        # 8. Audit JSON
        audit_data = {
            'track_id':     self.track_id,
            'generated_at': datetime.datetime.now().isoformat(),
            'threshold':    threshold,
            'netmhcpan':    {k: net_data[k] for k in ('n_raw', 'dropped_nan', 'n_0a', 'dropped_thr', 'n_0b')},
            'mhcflurry':    {k: flu_data[k] for k in ('n_raw', 'dropped_nan', 'n_0a', 'dropped_thr', 'n_0b')},
            'intersection': {k: intersection[k] for k in ('net_only', 'flurry_only', 'common_count')},
            'phase3_count': len(consensus_df),
            'calis':        calis_audit,
            'outputs': {
                'audit_0a_net':       str(output_dir / '0a_PRED_NET_no_nan.csv'),
                'audit_0a_flurry':    str(output_dir / '0a_PRED_FLURRY_no_nan.csv'),
                'audit_0b_net':       str(output_dir / '0b_PRED_NET_thresholded.csv'),
                'audit_0b_flurry':    str(output_dir / '0b_PRED_FLURRY_thresholded.csv'),
                'audit_1_net':        str(output_dir / '1_PRED_NET_consolidated.csv'),
                'audit_1_flurry':     str(output_dir / '1_PRED_FLURRY_consolidated.csv'),
                'audit_2_intersect':  str(output_dir / '2_intersection.csv'),
                'consensus_csv':      str(consensus_csv),
                'immunogenic_csv':    str(immunogenic_csv),
            },
        }
        audit_path = output_dir / 'consensus_audit_summary.json'
        with open(audit_path, 'w', encoding='utf-8') as audit_file:
            json.dump(audit_data, audit_file, indent=2, ensure_ascii=False)
        console.print(f'[dim]Audit → {audit_path}[/dim]')

        return {
            'consensus_csv':        str(consensus_csv),
            'immunogenic_csv':      str(immunogenic_csv),
            'audit_summary_json':   str(audit_path),
            'phase3_peptide_count': len(consensus_df),
            'stage2_peptide_count': len(immunogenic_df),
        }
