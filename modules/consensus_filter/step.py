"""Orchestration for consensus_filter: runs the two-stage filter on the
NetMHCpan + MHCFlurry predictions and writes the consensus CSV."""

import datetime
import json

import pandas as pd
from rich import box
from rich.panel import Panel

from modules.base_step import BaseTrackStep
from utils.console import console, flush_stdin
from utils.naming import get_prediction_filename, get_step_filename

from .core import (
    _FLURRY_PERCENTILE_SUBSTRINGS,
    _NET_PERCENTILE_SUBSTRINGS,
    _apply_calis,
    _consolidate,
    _finalize,
    _intersect,
    _load_and_filter,
)
from .prompts import _ask_threshold
from .render import (
    _print_stage1_filtering,
    _print_stage2_consolidation,
    _print_stage3_intersection,
    _print_stage4_immunogenicity,
)

# ── Step ──────────────────────────────────────────────────────────────────────

class ConsensusFilterStep(BaseTrackStep):
    step_name   = 'consensus_filter'
    description = (
        "Keeps only peptides predicted as binders by BOTH NetMHCpan and "
        "MHCFlurry (percentile threshold), consolidates one row per peptide "
        "with the joined HLA list, then drops anything Calis 2013 "
        "immunogenicity scores at zero or below."
    )
    long_description = (
        "Two filtering stages back-to-back:\n\n"
        "  • [bold]Stage 1 — Cross-tool consensus[/bold]: a peptide only passes "
        "if its NetMHCpan EL %rank ≤ 2 AND its MHCFlurry presentation "
        "percentile ≤ 2. Eliminates ~95% of input volume while keeping the "
        "biologically plausible MHC-I binders.\n\n"
        "  • [bold]Stage 2 — Calis immunogenicity[/bold]: each surviving peptide "
        "is scored by Calis 2013 against its bound HLA allele(s). Peptides "
        "with score ≤ 0 are dropped — those are predicted not to elicit a "
        "T-cell response even though they bind."
    )
    methodology = (
        "1. Load NetMHCpan + MHCFlurry CSVs from predict_binding.\n"
        "2. Drop NaNs and rows where the percentile column is missing.\n"
        "3. Threshold both at ≤ 2% (configurable via project_config).\n"
        "4. Consolidate per tool: one row per peptide with the union of HLA "
        "alleles that bound it + the best percentile.\n"
        "5. Intersect the two peptide sets — only peptides flagged by both "
        "tools survive.\n"
        "6. Calis: anchor-mask scoring per peptide; peptides receive the "
        "score from the best of their bound alleles (or default mask if none "
        "of the bound alleles is in the Calis allele table)."
    )
    references = [
        {
            'authors': 'Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Keşmir C, Peters B.',
            'title':   'Properties of MHC class I presented peptides that enhance immunogenicity',
            'journal': 'PLoS Computational Biology',
            'year':    2013,
            'doi':     '10.1371/journal.pcbi.1003266',
        },
    ]
    data_format = (
        "Input is automatic — uses both prediction CSVs from predict_binding. "
        "You will be asked once for the percentile threshold (default 2.0%)."
    )
    outputs_overview = (
        "[bold]CONSENSUS_VIEW_{track_id}.csv[/bold]     — slim per-step view (peptide + best percentile per tool + Calis score).\n"
        "[bold]CONSENSUS_IMMUNOGENIC_{track_id}.csv[/bold] — final output: peptides passing Stage 1 + Stage 2.\n"
        "[bold]CONSENSUS_{track_id}.csv[/bold]          — Stage 1 only (before Calis filter).\n"
        "[bold]0a/0b/1/2_*.csv[/bold]                   — intermediate audit CSVs (one per processing stage).\n"
        "[bold]consensus_audit_summary.json[/bold]      — funnel: counts at each filtering stage."
    )
    tips = [
        "The 2% percentile cutoff is conservative — relaxing to 5% triples volume; tightening to 0.5% can miss true binders.",
        "Calis scoring uses anchor-position immunogenicity weights — peptides binding alleles outside the Calis table fall back to the default mask.",
        "If volume after Stage 1 looks low, check predict_binding's audit JSON — it may be an upstream prediction failure, not over-filtering here.",
    ]

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

        # 7. Calis 2013 immunogenicity
        immunogenic_df, calis_audit = _apply_calis(consensus_df)
        flush_stdin()

        if 'calis_score' in immunogenic_df.columns:
            immunogenic_df['calis_score'] = pd.to_numeric(
                immunogenic_df['calis_score'], errors='coerce'
            ).round(3)
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
