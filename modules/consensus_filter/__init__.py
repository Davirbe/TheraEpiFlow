"""
Step "consensus_filter" — orchestrator.

Adapted from existing_scripts/analise_consenso_predicao.py with two adjustments:
  1) Reduced from 4 sources (PRED_NET, PRED_FLURRY, PROC_NET, PROC_FLURRY) to
     2 sources (PRED_NET = NetMHCpan EL %Rank, PRED_FLURRY = MHCflurry
     presentation_percentile). The PROC_NET/PROC_FLURRY processing data adds
     no signal beyond what EL / presentation already capture.
  2) Added a Phase 0b threshold filter (the original consumed pre-filtered
     tables; here we run on raw step "predict_binding" output).

The 4-phase structure of the original is preserved and visible in the console
log so the user can audit every drop:

  FASE 0a  load + drop NaN in essential columns               → 0a_*_no_nan.csv
  FASE 0b  apply percentile ≤ threshold                       → 0b_*_thresholded.csv
  FASE 1   consolidate per peptide (HLAs_agregados, Num_HLAs) → 1_*_consolidated.csv
  FASE 2   intersection of peptides between the two tools     → 2_Intermediario_PRED.csv
  FASE 3   final table with _net_pred / _flurry_pred suffixes → 3_FINAL_consensus.csv
  CONTROL  Rich panel + JSON audit (BEFORE Stage 2)
  STAGE 2  Calis 2013 immunogenicity > 0                      → 4_consensus_immunogenic.csv

Threshold is asked the first time step runs in a project, then persisted in
project_config.json["consensus_filter"]. Two presets (≤ 2.0 / ≤ 0.5) plus a
custom-value option.
"""

import datetime
import json
from pathlib import Path

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import box

from modules.base_step import BaseTrackStep
from utils.project_manager import save_project_config

from .consensus_phases import (
    load_prediction_csv_with_decimal_normalization,
    phase0_clean_and_threshold,
    phase1_consolidate_per_peptide,
    phase2_intersect_common_peptides,
    phase3_finalize_with_source_prefixes,
)
from .calis_wrapper import (
    score_peptides_calis_immunogenicity,
    is_allele_supported_by_calis,
)

console = Console(width=120)


# ── Column-name substring matchers ────────────────────────────────────────────
# Step "predict_binding" hasn't been implemented yet, so we accept several
# plausible column names (substring-style match — see find_column_by_substring_match).

NETMHCPAN_EL_PERCENTILE_SUBSTRINGS = ['netmhcpan_el', 'percentile']
MHCFLURRY_PRESENTATION_PERCENTILE_SUBSTRINGS = ['mhcflurry', 'presentation', 'percentile']


# ── Threshold profile presets ────────────────────────────────────────────────

THRESHOLD_PROFILE_BINDERS_STRONG_AND_WEAK = {
    'profile_id':     'binders_strong_and_weak',
    'percentile_max': 2.0,
    'description_pt': 'Ligantes fortes + fracos (≤ 2.0)',
}
THRESHOLD_PROFILE_BINDERS_STRONG_ONLY = {
    'profile_id':     'binders_strong_only',
    'percentile_max': 0.5,
    'description_pt': 'Ligantes fortes apenas (≤ 0.5)',
}


# ── Step class ────────────────────────────────────────────────────────────────

class ConsensusFilterStep(BaseTrackStep):
    """Consensus filter — runs Stage 1 (presentation) + Stage 2 (immunogenicity)."""
    step_name = 'consensus_filter'

    def run(self, input_data=None):
        # 1. Resolve threshold (project_config or interactive prompt)
        threshold_config = self._resolve_or_ask_threshold()
        percentile_threshold_max = threshold_config['percentile_max']

        # 2. Locate input CSVs from the previous step
        netmhcpan_input_csv, mhcflurry_input_csv = self._locate_predict_binding_outputs()

        consensus_output_dir = self.track_dir / 'consensus'
        consensus_output_dir.mkdir(parents=True, exist_ok=True)

        console.print(Panel.fit(
            f'[bold cyan]Consensus filter — {self.track_id}[/bold cyan]\n'
            f'[dim]Threshold profile:[/dim] [bold]{threshold_config["profile_id"]}[/bold]  '
            f'(percentile ≤ {percentile_threshold_max})',
            box=box.ROUNDED, border_style='cyan',
        ))

        # 3. Load raw CSVs
        console.print('\n[bold]FASE 0 — Carregamento + limpeza + threshold[/bold]')
        netmhcpan_raw_dataframe = load_prediction_csv_with_decimal_normalization(netmhcpan_input_csv)
        mhcflurry_raw_dataframe = load_prediction_csv_with_decimal_normalization(mhcflurry_input_csv)

        # 4. Phase 0 — clean + threshold (per source)
        netmhcpan_phase0 = phase0_clean_and_threshold(
            raw_dataframe=netmhcpan_raw_dataframe,
            source_label='NetMHCpan',
            percentile_column_substrings=NETMHCPAN_EL_PERCENTILE_SUBSTRINGS,
            percentile_threshold_max=percentile_threshold_max,
        )
        mhcflurry_phase0 = phase0_clean_and_threshold(
            raw_dataframe=mhcflurry_raw_dataframe,
            source_label='MHCflurry',
            percentile_column_substrings=MHCFLURRY_PRESENTATION_PERCENTILE_SUBSTRINGS,
            percentile_threshold_max=percentile_threshold_max,
        )
        netmhcpan_phase0['phase0a_dataframe'].to_csv(
            consensus_output_dir / '0a_PRED_NET_no_nan.csv',
            index=False, sep=';', decimal=',',
        )
        netmhcpan_phase0['phase0b_dataframe'].to_csv(
            consensus_output_dir / '0b_PRED_NET_thresholded.csv',
            index=False, sep=';', decimal=',',
        )
        mhcflurry_phase0['phase0a_dataframe'].to_csv(
            consensus_output_dir / '0a_PRED_FLURRY_no_nan.csv',
            index=False, sep=';', decimal=',',
        )
        mhcflurry_phase0['phase0b_dataframe'].to_csv(
            consensus_output_dir / '0b_PRED_FLURRY_thresholded.csv',
            index=False, sep=';', decimal=',',
        )

        # 5. Phase 1 — consolidate per peptide
        console.print('[bold]FASE 1 — Consolidação por peptídeo[/bold]')
        netmhcpan_consolidated = phase1_consolidate_per_peptide(
            thresholded_dataframe=netmhcpan_phase0['phase0b_dataframe'],
            percentile_column_name=netmhcpan_phase0['percentile_column_name'],
        )
        mhcflurry_consolidated = phase1_consolidate_per_peptide(
            thresholded_dataframe=mhcflurry_phase0['phase0b_dataframe'],
            percentile_column_name=mhcflurry_phase0['percentile_column_name'],
        )
        netmhcpan_consolidated.to_csv(
            consensus_output_dir / '1_PRED_NET_consolidated.csv',
            index=False, sep=';', decimal=',',
        )
        mhcflurry_consolidated.to_csv(
            consensus_output_dir / '1_PRED_FLURRY_consolidated.csv',
            index=False, sep=';', decimal=',',
        )

        # 6. Phase 2 — intersect peptides between tools
        console.print('[bold]FASE 2 — Intersecção entre ferramentas[/bold]')
        phase2_result = phase2_intersect_common_peptides(
            netmhcpan_consolidated=netmhcpan_consolidated,
            mhcflurry_consolidated=mhcflurry_consolidated,
        )
        phase2_result['common_peptides_dataframe'].to_csv(
            consensus_output_dir / '2_Intermediario_PRED.csv',
            index=False, sep=';', decimal=',',
        )

        # 7. Phase 3 — finalize with source prefixes
        console.print('[bold]FASE 3 — Tabela final com prefixos[/bold]')
        consensus_dataframe = phase3_finalize_with_source_prefixes(
            common_peptides_dataframe=phase2_result['common_peptides_dataframe'],
            netmhcpan_consolidated=netmhcpan_consolidated,
            mhcflurry_consolidated=mhcflurry_consolidated,
        )
        consensus_csv_path = consensus_output_dir / '3_FINAL_consensus.csv'
        consensus_dataframe.to_csv(consensus_csv_path, index=False, sep=';', decimal=',')

        # 8. CONTROL panel — surface losses BEFORE running Stage 2
        audit_summary = self._build_audit_summary(
            threshold_config=threshold_config,
            netmhcpan_phase0=netmhcpan_phase0,
            mhcflurry_phase0=mhcflurry_phase0,
            netmhcpan_consolidated=netmhcpan_consolidated,
            mhcflurry_consolidated=mhcflurry_consolidated,
            phase2_result=phase2_result,
            consensus_dataframe=consensus_dataframe,
        )
        self._print_control_panel(audit_summary)

        # 9. Stage 2 — Calis 2013 immunogenicity
        console.print('\n[bold]STAGE 2 — Imunogenicidade Calis 2013 (score > 0)[/bold]')
        immunogenic_dataframe, calis_audit = self._apply_stage2_calis(consensus_dataframe)
        immunogenic_csv_path = consensus_output_dir / '4_consensus_immunogenic.csv'
        immunogenic_dataframe.to_csv(
            immunogenic_csv_path, index=False, sep=';', decimal=',',
        )
        audit_summary['stage2_calis'] = calis_audit
        audit_summary['outputs'] = {
            'phase0a_no_nan_netmhcpan':       str(consensus_output_dir / '0a_PRED_NET_no_nan.csv'),
            'phase0a_no_nan_mhcflurry':       str(consensus_output_dir / '0a_PRED_FLURRY_no_nan.csv'),
            'phase0b_thresholded_netmhcpan':  str(consensus_output_dir / '0b_PRED_NET_thresholded.csv'),
            'phase0b_thresholded_mhcflurry':  str(consensus_output_dir / '0b_PRED_FLURRY_thresholded.csv'),
            'phase1_consolidated_netmhcpan':  str(consensus_output_dir / '1_PRED_NET_consolidated.csv'),
            'phase1_consolidated_mhcflurry':  str(consensus_output_dir / '1_PRED_FLURRY_consolidated.csv'),
            'phase2_intersection':            str(consensus_output_dir / '2_Intermediario_PRED.csv'),
            'phase3_consensus':               str(consensus_csv_path),
            'stage2_immunogenic':             str(immunogenic_csv_path),
        }
        audit_summary_path = consensus_output_dir / 'consensus_audit_summary.json'
        with open(audit_summary_path, 'w', encoding='utf-8') as audit_file:
            json.dump(audit_summary, audit_file, indent=2, ensure_ascii=False)

        console.print(
            f'\n[dim]Audit summary → {audit_summary_path}[/dim]'
        )
        console.print(
            f'[bold green]Stage 2 survivors (Calis > 0): {len(immunogenic_dataframe)}[/bold green]'
        )

        return {
            'consensus_csv':         str(consensus_csv_path),
            'immunogenic_csv':       str(immunogenic_csv_path),
            'audit_summary_json':    str(audit_summary_path),
            'phase3_peptide_count':  len(consensus_dataframe),
            'stage2_peptide_count':  len(immunogenic_dataframe),
        }

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _locate_predict_binding_outputs(self):
        """Returns (netmhcpan_csv_path, mhcflurry_csv_path); raises if missing."""
        track_predictions_dir = self.track_dir / 'predictions'
        netmhcpan_csv = (
            track_predictions_dir / f'netmhcpan_predictions_{self.track_id}.csv'
        )
        mhcflurry_csv = (
            track_predictions_dir / f'mhcflurry_predictions_{self.track_id}.csv'
        )
        if not netmhcpan_csv.exists():
            raise FileNotFoundError(
                f'NetMHCpan predictions not found at {netmhcpan_csv}. '
                f'Run step "predict_binding" first.'
            )
        if not mhcflurry_csv.exists():
            raise FileNotFoundError(
                f'MHCflurry predictions not found at {mhcflurry_csv}. '
                f'Run step "predict_binding" first.'
            )
        return netmhcpan_csv, mhcflurry_csv

    def _resolve_or_ask_threshold(self) -> dict:
        """Reads consensus_filter from project_config or asks the user once."""
        existing_threshold_config = self.project_config.get('consensus_filter')
        if existing_threshold_config:
            return existing_threshold_config

        chosen_threshold_config = _ask_threshold_interactively()
        chosen_threshold_config['asked_at'] = datetime.datetime.now().isoformat()
        self.project_config['consensus_filter'] = chosen_threshold_config
        save_project_config(self.project_name, self.project_config)
        return chosen_threshold_config

    def _apply_stage2_calis(self, consensus_dataframe: pd.DataFrame) -> tuple:
        """
        Score each consensus peptide via Calis using the best NetMHCpan allele
        (melhor_allele_net_pred). Group by allele to apply the right
        anchor-position mask in one batched call per allele.

        Returns (filtered_dataframe, audit_dict).
        """
        if consensus_dataframe.empty:
            empty_with_score_columns = consensus_dataframe.assign(
                calis_immunogenicity_score=pd.Series(dtype=float),
                calis_allele_used=pd.Series(dtype=str),
            )
            return empty_with_score_columns, {
                'input_count':                  0,
                'scored_count':                 0,
                'rejected_below_or_equal_zero': 0,
                'survived_above_zero':          0,
                'unsupported_allele_count':     0,
            }

        peptide_to_score       = {}
        peptide_to_allele_used = {}
        unsupported_allele_count = 0

        for best_allele, peptide_group in consensus_dataframe.groupby(
            'melhor_allele_net_pred', dropna=False,
        ):
            group_peptide_list = peptide_group['peptide'].tolist()

            if pd.isna(best_allele) or not best_allele:
                # No allele for this group → run with default mask
                group_calis_scores = score_peptides_calis_immunogenicity(
                    group_peptide_list, None,
                )
                allele_used_label = 'default_mask'
                unsupported_allele_count += len(group_peptide_list)
            elif not is_allele_supported_by_calis(best_allele):
                # Allele exists but is not in the Calis mask table → default mask
                group_calis_scores = score_peptides_calis_immunogenicity(
                    group_peptide_list, None,
                )
                allele_used_label = f'{best_allele} (default_mask — unsupported)'
                unsupported_allele_count += len(group_peptide_list)
            else:
                group_calis_scores = score_peptides_calis_immunogenicity(
                    group_peptide_list, best_allele,
                )
                allele_used_label = best_allele

            for peptide_sequence in group_peptide_list:
                peptide_to_score[peptide_sequence]       = group_calis_scores.get(peptide_sequence)
                peptide_to_allele_used[peptide_sequence] = allele_used_label

        scored_consensus_dataframe = consensus_dataframe.copy()
        scored_consensus_dataframe['calis_immunogenicity_score'] = (
            scored_consensus_dataframe['peptide'].map(peptide_to_score)
        )
        scored_consensus_dataframe['calis_allele_used'] = (
            scored_consensus_dataframe['peptide'].map(peptide_to_allele_used)
        )

        survivors_dataframe = scored_consensus_dataframe[
            scored_consensus_dataframe['calis_immunogenicity_score'].notna()
            & (scored_consensus_dataframe['calis_immunogenicity_score'] > 0.0)
        ].copy()

        return survivors_dataframe, {
            'input_count':                  len(consensus_dataframe),
            'scored_count':                 int(scored_consensus_dataframe['calis_immunogenicity_score'].notna().sum()),
            'rejected_below_or_equal_zero': len(consensus_dataframe) - len(survivors_dataframe),
            'survived_above_zero':          len(survivors_dataframe),
            'unsupported_allele_count':     unsupported_allele_count,
        }

    def _build_audit_summary(self, **kwargs) -> dict:
        threshold_config       = kwargs['threshold_config']
        netmhcpan_phase0       = kwargs['netmhcpan_phase0']
        mhcflurry_phase0       = kwargs['mhcflurry_phase0']
        netmhcpan_consolidated = kwargs['netmhcpan_consolidated']
        mhcflurry_consolidated = kwargs['mhcflurry_consolidated']
        phase2_result          = kwargs['phase2_result']
        consensus_dataframe    = kwargs['consensus_dataframe']

        return {
            'track_id':           self.track_id,
            'generated_at':       datetime.datetime.now().isoformat(),
            'threshold_profile':  threshold_config,
            'netmhcpan': {
                'percentile_column':           netmhcpan_phase0['percentile_column_name'],
                'phase0a_input_count':         netmhcpan_phase0['phase0a_input_count'],
                'phase0a_dropped_count':       netmhcpan_phase0['phase0a_dropped_count'],
                'phase0a_output_count':        netmhcpan_phase0['phase0a_output_count'],
                'phase0b_dropped_count':       netmhcpan_phase0['phase0b_dropped_count'],
                'phase0b_output_count':        netmhcpan_phase0['phase0b_output_count'],
                'phase1_unique_peptide_count': len(netmhcpan_consolidated),
            },
            'mhcflurry': {
                'percentile_column':           mhcflurry_phase0['percentile_column_name'],
                'phase0a_input_count':         mhcflurry_phase0['phase0a_input_count'],
                'phase0a_dropped_count':       mhcflurry_phase0['phase0a_dropped_count'],
                'phase0a_output_count':        mhcflurry_phase0['phase0a_output_count'],
                'phase0b_dropped_count':       mhcflurry_phase0['phase0b_dropped_count'],
                'phase0b_output_count':        mhcflurry_phase0['phase0b_output_count'],
                'phase1_unique_peptide_count': len(mhcflurry_consolidated),
            },
            'phase2_intersection': {
                'netmhcpan_only_count': phase2_result['netmhcpan_only_count'],
                'mhcflurry_only_count': phase2_result['mhcflurry_only_count'],
                'common_count':         phase2_result['common_count'],
            },
            'phase3_consensus_count': len(consensus_dataframe),
        }

    def _print_control_panel(self, audit_summary: dict):
        """Renders the Rich panel + per-phase counts before Stage 2 starts."""
        net = audit_summary['netmhcpan']
        flu = audit_summary['mhcflurry']
        intersect = audit_summary['phase2_intersection']

        per_phase_table = Table(
            box=box.SIMPLE, show_header=True, header_style='bold',
            title='Controle de triagem (antes do Calis)',
        )
        per_phase_table.add_column('Fase', style='bold cyan', no_wrap=True)
        per_phase_table.add_column('NetMHCpan', justify='right')
        per_phase_table.add_column('MHCflurry', justify='right')
        per_phase_table.add_row(
            'FASE 0a — input',
            str(net['phase0a_input_count']),
            str(flu['phase0a_input_count']),
        )
        per_phase_table.add_row(
            'FASE 0a — removido (NaN)',
            f'-{net["phase0a_dropped_count"]}',
            f'-{flu["phase0a_dropped_count"]}',
        )
        per_phase_table.add_row(
            'FASE 0a — sobrevive',
            str(net['phase0a_output_count']),
            str(flu['phase0a_output_count']),
        )
        per_phase_table.add_row(
            'FASE 0b — removido (threshold)',
            f'-{net["phase0b_dropped_count"]}',
            f'-{flu["phase0b_dropped_count"]}',
        )
        per_phase_table.add_row(
            'FASE 0b — sobrevive',
            str(net['phase0b_output_count']),
            str(flu['phase0b_output_count']),
        )
        per_phase_table.add_row(
            'FASE 1 — peptídeos únicos',
            str(net['phase1_unique_peptide_count']),
            str(flu['phase1_unique_peptide_count']),
        )
        console.print(per_phase_table)

        intersection_table = Table(box=box.SIMPLE, show_header=False)
        intersection_table.add_column('label', style='bold cyan')
        intersection_table.add_column('value', justify='right')
        intersection_table.add_row(
            'Apenas NetMHCpan', str(intersect['netmhcpan_only_count']),
        )
        intersection_table.add_row(
            'Apenas MHCflurry', str(intersect['mhcflurry_only_count']),
        )
        intersection_table.add_row(
            'FASE 2 — comuns (consenso)',
            f'[bold green]{intersect["common_count"]}[/bold green]',
        )
        intersection_table.add_row(
            'FASE 3 — total no consenso final',
            f'[bold green]{audit_summary["phase3_consensus_count"]}[/bold green]',
        )
        console.print(Panel(
            intersection_table, title='Intersecção', border_style='cyan',
        ))


# ── Threshold prompt (asked once per project, then persisted) ─────────────────

def _ask_threshold_interactively() -> dict:
    """Asks the user which binding threshold to apply. Returns a profile dict."""
    console.print(Panel.fit(
        '[bold cyan]Threshold de consenso[/bold cyan]\n'
        '[dim]Mantém peptídeos com percentile ≤ threshold em AMBAS as ferramentas:\n'
        '  • NetMHCpan 4.1 EL %Rank\n'
        '  • MHCflurry 2.2 presentation_percentile[/dim]\n\n'
        '[bold]Escolha o perfil:[/bold]\n'
        '  [cyan][1][/cyan] Ligantes fortes + fracos    (≤ 2.0)   [dim]padrão[/dim]\n'
        '  [cyan][2][/cyan] Ligantes fortes apenas      (≤ 0.5)\n'
        '  [cyan][3][/cyan] Customizado',
        box=box.ROUNDED, border_style='cyan',
    ))

    while True:
        user_choice = input('> ').strip()

        if user_choice in ('', '1'):
            return dict(THRESHOLD_PROFILE_BINDERS_STRONG_AND_WEAK)
        if user_choice == '2':
            return dict(THRESHOLD_PROFILE_BINDERS_STRONG_ONLY)
        if user_choice == '3':
            console.print('[bold]Digite o threshold customizado (ex: 1.5):[/bold]')
            custom_value_text = input('> ').strip().replace(',', '.')
            try:
                custom_threshold_value = float(custom_value_text)
            except ValueError:
                console.print('[red]Valor inválido — digite um número.[/red]')
                continue
            if custom_threshold_value <= 0:
                console.print('[red]Threshold deve ser maior que 0.[/red]')
                continue
            return {
                'profile_id':     'custom',
                'percentile_max': custom_threshold_value,
                'description_pt': f'Customizado (≤ {custom_threshold_value})',
            }

        console.print('[red]Opção inválida — escolha 1, 2 ou 3.[/red]')
