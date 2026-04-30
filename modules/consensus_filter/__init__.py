"""
Step "consensus_filter".

Recebe os CSVs de predição gerados pelo step "predict_binding" e aplica:
  Stage 1 — filtro de apresentação (percentile threshold + interseção entre tools)
  Stage 2 — filtro de imunogenicidade (Calis 2013, score > 0)

Entradas (track_dir/predictions/):
  PRED_NET_{track_id}.csv
  PRED_FLURRY_{track_id}.csv

Saídas em track_dir/consensus/:
  0a_PRED_NET_no_nan.csv            linhas sem NaN (auditoria)
  0a_PRED_FLURRY_no_nan.csv
  0b_PRED_NET_thresholded.csv       após corte de percentile (auditoria)
  0b_PRED_FLURRY_thresholded.csv
  1_PRED_NET_consolidated.csv       1 linha por peptídeo (auditoria)
  1_PRED_FLURRY_consolidated.csv
  2_Intermediario_PRED.csv          interseção entre os dois tools (auditoria)
  CONSENSUS_{track_id}.csv          tabela final com prefixos _net_pred/_flurry_pred
  CONSENSUS_IMMUNOGENIC_{track_id}.csv  sobreviventes do Calis (score > 0)
  consensus_audit_summary.json      contagens por fase
"""

import contextlib
import datetime
import io
import json
from pathlib import Path

import pandas as pd
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from modules.base_step import BaseTrackStep
from utils.project_manager import save_project_config
from utils.naming import get_prediction_filename, get_step_filename

console = Console(width=120)

# Substrings usadas para localizar as colunas de percentile em cada tool
_NET_PERCENTILE_SUBS    = ['netmhcpan_el', 'percentile']
_FLURRY_PERCENTILE_SUBS = ['mhcflurry', 'presentation', 'percentile']


# ── Threshold ─────────────────────────────────────────────────────────────────

def _ask_threshold(project_name: str, project_config: dict) -> float:
    """
    Lê o threshold do project_config ou pergunta uma vez ao usuário.
    Salva no project_config para não perguntar novamente.
    """
    existing = project_config.get('consensus_threshold')
    if existing is not None:
        return float(existing)

    console.print(Panel.fit(
        '[bold cyan]Threshold de consenso[/bold cyan]\n'
        '[dim]Mantém peptídeos com percentile ≤ threshold em AMBAS as ferramentas[/dim]\n\n'
        '  [cyan][1][/cyan] Ligantes fortes + fracos  (≤ 2.0)  [dim]padrão[/dim]\n'
        '  [cyan][2][/cyan] Ligantes fortes apenas    (≤ 0.5)\n'
        '  [cyan][3][/cyan] Valor customizado',
        box=box.ROUNDED, border_style='cyan',
    ))

    while True:
        try:
            choice = input('> ').strip()
        except EOFError:
            choice = '1'
        if choice in ('', '1'):
            threshold = 2.0
            break
        if choice == '2':
            threshold = 0.5
            break
        if choice == '3':
            try:
                raw = input('Digite o threshold (ex: 1.5): ').strip().replace(',', '.')
            except EOFError:
                raw = '2.0'
            try:
                threshold = float(raw)
                if threshold <= 0:
                    raise ValueError
                break
            except ValueError:
                console.print('[red]Valor inválido — digite um número maior que 0.[/red]')
                continue
        console.print('[red]Escolha 1, 2 ou 3.[/red]')

    project_config['consensus_threshold'] = threshold
    save_project_config(project_name, project_config)
    console.print(f'[dim]Threshold definido: ≤ {threshold} — salvo no project_config.[/dim]')
    return threshold


# ── Funções auxiliares de DataFrame ──────────────────────────────────────────

def _find_col(df: pd.DataFrame, substrings: list) -> str | None:
    """Retorna o primeiro nome de coluna que contém todas as substrings."""
    for col in df.columns:
        if all(s in col for s in substrings):
            return col
    return None


def _normalize_col_name(name: str) -> str:
    """Lowercase + underscores. 'seq #' → 'sequence_number'."""
    cleaned = name.strip()
    if cleaned.lower() == 'seq #':
        return 'sequence_number'
    return cleaned.replace(' ', '_').replace('#', 'num').lower()


# ── Stage 1 — Apresentação ────────────────────────────────────────────────────

def _load_and_filter(csv_path: Path, percentile_subs: list, threshold: float) -> dict:
    """
    Carrega um CSV de predição e aplica duas sub-fases:
      0a — remove linhas com NaN em peptide, allele ou percentile
      0b — mantém só linhas com percentile ≤ threshold

    Retorna dict com os dois DataFrames intermediários e metadados de auditoria.
    """
    # Carregamento com detecção de separador e normalização de decimais
    df = pd.read_csv(csv_path, sep=',', dtype=str, keep_default_na=False,
                     na_values=['', 'NA', 'N/A'], engine='python')
    if df.shape[1] == 1:
        df = pd.read_csv(csv_path, sep=';', dtype=str, keep_default_na=False,
                         na_values=['', 'NA', 'N/A'], engine='python')

    df.columns = [_normalize_col_name(c) for c in df.columns]
    if 'peptide' in df.columns:
        df['peptide'] = df['peptide'].str.strip()

    # Converte colunas numéricas
    numeric_cols = [c for c in df.columns
                    if any(k in c for k in ('percentile', 'score', 'ic50', 'rank', 'affinity'))]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col].str.replace(',', '.', regex=False), errors='coerce')

    # Localiza coluna de percentile
    pct_col = _find_col(df, percentile_subs)
    if pct_col is None:
        raise ValueError(f'Coluna de percentile não encontrada em {csv_path}. '
                         f'Substrings buscadas: {percentile_subs}. '
                         f'Colunas disponíveis: {list(df.columns)}')

    # 0a — remove NaN
    n_raw = len(df)
    df_0a = df.dropna(subset=['peptide', 'allele', pct_col]).copy()
    n_0a  = len(df_0a)

    # 0b — aplica threshold
    df_0b = df_0a[df_0a[pct_col] <= threshold].copy()
    n_0b  = len(df_0b)

    return {
        'df_0a':      df_0a,
        'df_0b':      df_0b,
        'pct_col':    pct_col,
        'n_raw':      n_raw,
        'n_0a':       n_0a,
        'dropped_nan': n_raw - n_0a,
        'n_0b':       n_0b,
        'dropped_thr': n_0a - n_0b,
    }


def _consolidate(df: pd.DataFrame, pct_col: str) -> pd.DataFrame:
    """
    Fase 1 — colapsa N linhas por allele em 1 linha por peptídeo.

    Colunas adicionadas:
      melhor_allele   — allele com menor percentile
      HLAs_agregados  — todos os alleles únicos separados por ';'
      Num_HLAs        — contagem de alleles únicos
    """
    if df.empty:
        return pd.DataFrame()

    grouped = df.groupby('peptide')

    # Agrega alleles
    hla_agg = (grouped['allele']
               .apply(lambda x: ';'.join(sorted(x.unique())))
               .reset_index()
               .rename(columns={'allele': 'HLAs_agregados'}))
    hla_agg['Num_HLAs'] = hla_agg['HLAs_agregados'].apply(lambda x: len(x.split(';')))

    # Linha com melhor (menor) percentile
    best = df.loc[grouped[pct_col].idxmin()].copy()
    best = best.rename(columns={'allele': 'melhor_allele'})

    return pd.merge(best, hla_agg, on='peptide', how='left')


def _intersect(net: pd.DataFrame, flurry: pd.DataFrame) -> dict:
    """
    Fase 2 — peptídeos presentes nos DOIS tools.
    Retorna dict com DataFrame de peptídeos comuns e contagens de auditoria.
    """
    if net.empty or flurry.empty:
        return {
            'common': pd.DataFrame(columns=['peptide']),
            'net_only':    len(net),
            'flurry_only': len(flurry),
            'common_count': 0,
        }

    set_net    = set(net['peptide'])
    set_flurry = set(flurry['peptide'])
    common     = set_net & set_flurry

    return {
        'common':       pd.DataFrame({'peptide': sorted(common)}),
        'net_only':     len(set_net - set_flurry),
        'flurry_only':  len(set_flurry - set_net),
        'common_count': len(common),
    }


def _finalize(common: pd.DataFrame, net: pd.DataFrame, flurry: pd.DataFrame) -> pd.DataFrame:
    """
    Fase 3 — monta tabela final.
    Renomeia colunas de cada tool com sufixos _net_pred / _flurry_pred
    e faz merge na coluna peptide.
    """
    if common.empty:
        return pd.DataFrame(columns=['peptide'])

    net_f = net[net['peptide'].isin(common['peptide'])].copy()
    flu_f = flurry[flurry['peptide'].isin(common['peptide'])].copy()

    net_f = net_f.rename(columns={c: f'{c}_net_pred'    for c in net_f.columns   if c != 'peptide'})
    flu_f = flu_f.rename(columns={c: f'{c}_flurry_pred' for c in flu_f.columns   if c != 'peptide'})

    df = pd.merge(common, net_f, on='peptide', how='left')
    df = pd.merge(df,     flu_f, on='peptide', how='left')
    return df


# ── Stage 2 — Calis 2013 ─────────────────────────────────────────────────────

# Tabela allele → posições âncora (1-based). Reproduzida do predict_immunogenicity.py
# para resolver a máscara sem chamar validate() (que usa sys.exit e lê arquivo).
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


def _score_calis(peptides: list, allele_imgt: str | None) -> dict:
    """
    Pontua uma lista de peptídeos via Calis 2013 (predict_immunogenicity.Prediction).
    Captura stdout porque Prediction.predict() imprime em vez de retornar.

    Args:
        peptides:    lista de sequências em maiúsculas
        allele_imgt: formato IMGT ex: 'HLA-A*02:01', ou None → máscara default

    Retorna dict {peptide: score}.
    """
    from .predict_immunogenicity import Prediction

    # Resolve máscara e allele no formato Calis (sem * e :)
    allele_calis = None
    mask         = None
    if allele_imgt:
        normalized = allele_imgt.replace('*', '').replace(':', '')
        if normalized in _CALIS_ALLELES:
            allele_calis = normalized
            mask         = _CALIS_ALLELES[normalized]

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        Prediction().predict((list(peptides), mask, allele_calis))

    # Parseia saída: linhas com formato "PEPTIDE,length,score"
    scores = {}
    for line in buf.getvalue().splitlines():
        line = line.strip()
        if not line or ',' not in line:
            continue
        if line.startswith(('peptide,', 'masking:', 'masked', 'allele:')):
            continue
        parts = line.split(',')
        if len(parts) == 3:
            try:
                scores[parts[0]] = float(parts[2])
            except ValueError:
                continue
    return scores


def _apply_calis(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Aplica Calis agrupando por melhor_allele_net_pred.
    Retorna (DataFrame filtrado score > 0, dict de auditoria).
    """
    if df.empty:
        return df.assign(calis_score=pd.Series(dtype=float)), {
            'input': 0, 'scored': 0, 'survived': 0, 'unsupported_allele': 0,
        }

    peptide_scores  = {}
    peptide_alleles = {}
    unsupported     = 0

    for allele, group in df.groupby('melhor_allele_net_pred', dropna=False):
        peps = group['peptide'].tolist()
        if pd.isna(allele) or not allele:
            result        = _score_calis(peps, None)
            allele_label  = 'default_mask'
            unsupported  += len(peps)
        else:
            normalized = allele.replace('*', '').replace(':', '')
            if normalized not in _CALIS_ALLELES:
                result        = _score_calis(peps, None)
                allele_label  = f'{allele} (default_mask)'
                unsupported  += len(peps)
            else:
                result       = _score_calis(peps, allele)
                allele_label = allele

        for p in peps:
            peptide_scores[p]  = result.get(p)
            peptide_alleles[p] = allele_label

    scored = df.copy()
    scored['calis_score']        = scored['peptide'].map(peptide_scores)
    scored['calis_allele_used']  = scored['peptide'].map(peptide_alleles)

    survivors = scored[scored['calis_score'].notna() & (scored['calis_score'] > 0)].copy()

    return survivors, {
        'input':               len(df),
        'scored':              int(scored['calis_score'].notna().sum()),
        'survived':            len(survivors),
        'unsupported_allele':  unsupported,
    }


# ── Painel Rich ───────────────────────────────────────────────────────────────

def _print_panel(net: dict, flu: dict, intersect: dict, n_phase3: int):
    """Imprime tabela de contagens por fase antes do Calis."""
    t = Table(box=box.SIMPLE, show_header=True, header_style='bold',
              title='Controle de triagem (antes do Calis)')
    t.add_column('Fase', style='bold cyan', no_wrap=True)
    t.add_column('NetMHCpan', justify='right')
    t.add_column('MHCflurry', justify='right')

    rows = [
        ('0a — input',              str(net['n_raw']),         str(flu['n_raw'])),
        ('0a — removido (NaN)',     f'-{net["dropped_nan"]}',  f'-{flu["dropped_nan"]}'),
        ('0a — sobrevive',          str(net['n_0a']),          str(flu['n_0a'])),
        ('0b — removido (thresh.)', f'-{net["dropped_thr"]}',  f'-{flu["dropped_thr"]}'),
        ('0b — sobrevive',          str(net['n_0b']),          str(flu['n_0b'])),
        ('1  — peptídeos únicos',   str(net['n_consolidated']),str(flu['n_consolidated'])),
    ]
    for row in rows:
        t.add_row(*row)
    console.print(t)

    i = Table(box=box.SIMPLE, show_header=False)
    i.add_column('label', style='bold cyan')
    i.add_column('value', justify='right')
    i.add_row('Apenas NetMHCpan',          str(intersect['net_only']))
    i.add_row('Apenas MHCflurry',          str(intersect['flurry_only']))
    i.add_row('2 — comuns (consenso)',      f'[bold green]{intersect["common_count"]}[/bold green]')
    i.add_row('3 — total no consenso final',f'[bold green]{n_phase3}[/bold green]')
    console.print(Panel(i, title='Intersecção', border_style='cyan'))


# ── Step ──────────────────────────────────────────────────────────────────────

class ConsensusFilterStep(BaseTrackStep):
    step_name = 'consensus_filter'

    def run(self, input_data=None):
        # 1. Threshold
        threshold = _ask_threshold(self.project_name, self.project_config)

        # 2. Localiza CSVs do step anterior
        pred_dir  = self.track_dir / 'predictions'
        net_csv   = pred_dir / get_prediction_filename("NET_PRED", self.track_id)
        flu_csv   = pred_dir / get_prediction_filename("FLURRY_PRED", self.track_id)
        for p in (net_csv, flu_csv):
            if not p.exists():
                raise FileNotFoundError(f'{p} não encontrado. Execute "predict_binding" primeiro.')

        out = self.track_dir / 'consensus'
        out.mkdir(parents=True, exist_ok=True)

        console.print(Panel.fit(
            f'[bold cyan]Consensus filter — {self.track_id}[/bold cyan]\n'
            f'[dim]Threshold:[/dim] [bold]≤ {threshold}[/bold]',
            box=box.ROUNDED, border_style='cyan',
        ))

        # 3. Fase 0 — carrega + NaN + threshold
        console.print('\n[bold]FASE 0 — Carregamento + limpeza + threshold[/bold]')
        net_data = _load_and_filter(net_csv,   _NET_PERCENTILE_SUBS,    threshold)
        flu_data = _load_and_filter(flu_csv,   _FLURRY_PERCENTILE_SUBS, threshold)

        net_data['df_0a'].to_csv(out / '0a_PRED_NET_no_nan.csv',          index=False, sep=';', decimal=',')
        net_data['df_0b'].to_csv(out / '0b_PRED_NET_thresholded.csv',     index=False, sep=';', decimal=',')
        flu_data['df_0a'].to_csv(out / '0a_PRED_FLURRY_no_nan.csv',       index=False, sep=';', decimal=',')
        flu_data['df_0b'].to_csv(out / '0b_PRED_FLURRY_thresholded.csv',  index=False, sep=';', decimal=',')

        # 4. Fase 1 — consolida por peptídeo
        console.print('[bold]FASE 1 — Consolidação por peptídeo[/bold]')
        net_cons = _consolidate(net_data['df_0b'], net_data['pct_col'])
        flu_cons = _consolidate(flu_data['df_0b'], flu_data['pct_col'])

        net_cons.to_csv(out / '1_PRED_NET_consolidated.csv',    index=False, sep=';', decimal=',')
        flu_cons.to_csv(out / '1_PRED_FLURRY_consolidated.csv', index=False, sep=';', decimal=',')

        # 5. Fase 2 — intersecta entre tools
        console.print('[bold]FASE 2 — Intersecção entre ferramentas[/bold]')
        intersect = _intersect(net_cons, flu_cons)
        intersect['common'].to_csv(out / '2_Intermediario_PRED.csv', index=False, sep=';', decimal=',')

        # 6. Fase 3 — tabela final com prefixos
        console.print('[bold]FASE 3 — Tabela final com prefixos[/bold]')
        consensus_df = _finalize(intersect['common'], net_cons, flu_cons)
        consensus_csv = out / get_step_filename("CONSENSUS", self.track_id)
        consensus_df.to_csv(consensus_csv, index=False)

        # 7. Painel Rich
        net_data['n_consolidated'] = len(net_cons)
        flu_data['n_consolidated'] = len(flu_cons)
        _print_panel(net_data, flu_data, intersect, len(consensus_df))

        # 8. Stage 2 — Calis
        console.print('\n[bold]STAGE 2 — Imunogenicidade Calis 2013 (score > 0)[/bold]')
        immuno_df, calis_audit = _apply_calis(consensus_df)
        immuno_csv = out / get_step_filename("CONSENSUS_IMMUNOGENIC", self.track_id)
        immuno_df.to_csv(immuno_csv, index=False)

        console.print(f'[bold green]Sobreviventes Calis > 0: {len(immuno_df)}[/bold green]')

        # 9. Audit JSON
        audit = {
            'track_id':      self.track_id,
            'generated_at':  datetime.datetime.now().isoformat(),
            'threshold':     threshold,
            'netmhcpan':     {k: net_data[k] for k in ('n_raw','dropped_nan','n_0a','dropped_thr','n_0b')},
            'mhcflurry':     {k: flu_data[k] for k in ('n_raw','dropped_nan','n_0a','dropped_thr','n_0b')},
            'intersection':  {k: intersect[k] for k in ('net_only','flurry_only','common_count')},
            'phase3_count':  len(consensus_df),
            'calis':         calis_audit,
            'outputs': {
                'audit_0a_net':       str(out / '0a_PRED_NET_no_nan.csv'),
                'audit_0a_flurry':    str(out / '0a_PRED_FLURRY_no_nan.csv'),
                'audit_0b_net':       str(out / '0b_PRED_NET_thresholded.csv'),
                'audit_0b_flurry':    str(out / '0b_PRED_FLURRY_thresholded.csv'),
                'audit_1_net':        str(out / '1_PRED_NET_consolidated.csv'),
                'audit_1_flurry':     str(out / '1_PRED_FLURRY_consolidated.csv'),
                'audit_2_intersect':  str(out / '2_Intermediario_PRED.csv'),
                'consensus_csv':      str(consensus_csv),
                'immunogenic_csv':    str(immuno_csv),
            },
        }
        audit_path = out / 'consensus_audit_summary.json'
        with open(audit_path, 'w', encoding='utf-8') as f:
            json.dump(audit, f, indent=2, ensure_ascii=False)
        console.print(f'[dim]Audit → {audit_path}[/dim]')

        return {
            'consensus_csv':        str(consensus_csv),
            'immunogenic_csv':      str(immuno_csv),
            'audit_summary_json':   str(audit_path),
            'phase3_peptide_count': len(consensus_df),
            'stage2_peptide_count': len(immuno_df),
        }
