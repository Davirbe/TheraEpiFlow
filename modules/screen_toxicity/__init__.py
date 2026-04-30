"""
Step "screen_toxicity".

Recebe os epítopos imunogênicos do consensus_filter e filtra os tóxicos
usando ToxinPred3 (Model 1: AAC + DPC, Extra Trees, threshold padrão 0.38).

Entrada (track_dir/consensus/):
  CONSENSUS_IMMUNOGENIC_{track_id}.csv

Saídas em track_dir/toxicity/:
  TOXICITY_ALL_{track_id}.csv      todos os peptídeos + scores
  TOXICITY_SAFE_{track_id}.csv     apenas não-tóxicos
  toxicity_audit_{track_id}.json   contagens e parâmetros
"""

import importlib.util
import json
import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rich import box
from rich.console import Console
from rich.table import Table

import config
from modules.base_step import BaseTrackStep
from utils.naming import get_step_filename, COLUMN_PEPTIDE

console = Console(width=120)

_STD_AAS = list("ACDEFGHIKLMNPQRSTVWY")


# ── ToxinPred3 core (in-memory, sem arquivos temporários) ─────────────────────

def _get_model_path() -> Path:
    spec = importlib.util.find_spec('toxinpred3')
    if spec is None:
        raise ImportError(
            "toxinpred3 não encontrado. Instale com: pip install toxinpred3>=1.4"
        )
    pkg_dir = Path(spec.submodule_search_locations[0])
    return pkg_dir / 'model' / 'toxinpred3.0_model.pkl'


def _aac_vector(seq: str) -> list:
    n = len(seq)
    return [(seq.count(aa) / n) * 100.0 for aa in _STD_AAS]


def _dpc_vector(seq: str) -> list:
    n = len(seq)
    denom = n - 1
    result = []
    for aa1 in _STD_AAS:
        for aa2 in _STD_AAS:
            count = sum(1 for i in range(n - 1) if seq[i] == aa1 and seq[i + 1] == aa2)
            result.append(count / denom * 100.0 if denom > 0 else 0.0)
    return result


def _predict(peptides: list, model_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Retorna (scores, ppv) — arrays 1-D, um valor por peptídeo."""
    clf = joblib.load(str(model_path))
    X = np.array([_aac_vector(p) + _dpc_vector(p) for p in peptides], dtype=np.float32)
    scores = clf.predict_proba(X)[:, -1]
    ppv = np.clip(scores * 1.2341 - 0.1182, 0.0, 1.0)
    return scores, ppv


# ── Display ───────────────────────────────────────────────────────────────────

def _print_result(n_input: int, n_toxic: int, n_safe: int, threshold: float):
    t = Table(box=box.SIMPLE, show_header=False,
              title=f'Triagem de toxicidade — ToxinPred3  (threshold={threshold:.2f})')
    t.add_column('', style='bold cyan', no_wrap=True)
    t.add_column('', justify='right')

    t.add_row('Entrada (epítopos imunogênicos)', str(n_input))
    t.add_row('Tóxicos removidos',
              f'[bold red]-{n_toxic}[/bold red]' if n_toxic else '[dim]0[/dim]')
    t.add_row('Seguros [bold green]✓[/bold green]',
              f'[bold green]{n_safe}[/bold green]')

    console.print(t)
    console.print(f'\n[bold green]RESULTADO: {n_safe} epítopos não-tóxicos passaram para o próximo step.[/bold green]\n')


# ── Step ──────────────────────────────────────────────────────────────────────

class ScreenToxicityStep(BaseTrackStep):
    step_name = 'screen_toxicity'

    def run(self, input_data=None):
        threshold = self.project_config.get(
            'toxicity_threshold', config.TOXICITY_SCORE_THRESHOLD
        )

        # Localiza input
        consensus_dir = self.track_dir / 'consensus'
        input_csv = consensus_dir / get_step_filename("CONSENSUS_IMMUNOGENIC", self.track_id)
        if not input_csv.exists():
            raise FileNotFoundError(
                f"Saída do consensus_filter não encontrada: {input_csv}\n"
                "Execute o step 'consensus_filter' antes de 'screen_toxicity'."
            )

        df = pd.read_csv(input_csv)
        if COLUMN_PEPTIDE not in df.columns:
            raise ValueError(
                f"Coluna '{COLUMN_PEPTIDE}' não encontrada em {input_csv.name}. "
                f"Colunas disponíveis: {list(df.columns)}"
            )

        peptides = df[COLUMN_PEPTIDE].tolist()
        n_input  = len(peptides)

        console.print(f'\n[bold]Carregando modelo ToxinPred3...[/bold]')
        model_path = _get_model_path()
        scores, ppv = _predict(peptides, model_path)

        df['toxinpred3_score'] = np.round(scores, 3)
        df['toxinpred3_ppv']   = np.round(ppv,    3)
        df['toxinpred3_label'] = [
            'Toxin' if s >= threshold else 'Non-Toxin' for s in scores
        ]

        out = self.track_dir / 'toxicity'
        out.mkdir(parents=True, exist_ok=True)

        all_csv  = out / get_step_filename("TOXICITY_ALL",  self.track_id)
        safe_csv = out / get_step_filename("TOXICITY_SAFE", self.track_id)

        df.to_csv(all_csv, index=False)

        df_safe = df[df['toxinpred3_label'] == 'Non-Toxin'].copy()
        df_safe.to_csv(safe_csv, index=False)

        n_toxic = int((df['toxinpred3_label'] == 'Toxin').sum())
        n_safe  = len(df_safe)

        _print_result(n_input, n_toxic, n_safe, threshold)

        audit = {
            'timestamp':      datetime.datetime.now().isoformat(),
            'track_id':       self.track_id,
            'tool':           'ToxinPred3',
            'model':          'AAC+DPC Extra Trees (Model 1)',
            'threshold':      threshold,
            'n_input':        n_input,
            'n_toxic':        n_toxic,
            'n_safe':         n_safe,
            'pct_removed':    round(n_toxic / n_input * 100, 1) if n_input else 0,
            'output_all':     str(all_csv),
            'output_safe':    str(safe_csv),
        }
        audit_path = out / f'toxicity_audit_{self.track_id}.json'
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False))

        return {
            'output_all':  str(all_csv),
            'output_safe': str(safe_csv),
            'n_input':     n_input,
            'n_toxic':     n_toxic,
            'n_safe':      n_safe,
        }
