"""
Validação da abordagem in-memory para ToxinPred3 (Model 1, AAC+DPC based ET).
Rodar antes de usar o step screen_toxicity para confirmar que o modelo funciona.
"""

import importlib.util
from pathlib import Path

import joblib
import numpy as np

_STD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

TEST_PEPTIDES = [
    # sequência              descrição esperada
    ("SLYNTVATLY",  "epítopo HIV gag — NÃO tóxico"),
    ("GILGFVFTL",   "epítopo influenza M1 — NÃO tóxico"),
    ("MKTIIALSYI",  "peptídeo sintético neutro — NÃO tóxico"),
    ("FLPIIAGAII",  "peptídeo antimicrobiano / membrana — POSSIVELMENTE tóxico"),
    ("KWKLFKKIEK",  "magainin-like (carga+, anfipático) — POSSIVELMENTE tóxico"),
]


def _get_model_path() -> Path:
    spec = importlib.util.find_spec('toxinpred3')
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


def predict(peptides: list, model_path: Path, threshold: float = 0.38):
    clf = joblib.load(str(model_path))
    X = np.array([_aac_vector(p) + _dpc_vector(p) for p in peptides], dtype=np.float32)
    scores = clf.predict_proba(X)[:, -1]
    ppv = np.clip(scores * 1.2341 - 0.1182, 0.0, 1.0)
    labels = ["Toxin" if s >= threshold else "Non-Toxin" for s in scores]
    return scores, ppv, labels


if __name__ == '__main__':
    model_path = _get_model_path()
    print(f"Modelo: {model_path}\n")

    threshold = 0.38
    seqs = [p[0] for p in TEST_PEPTIDES]
    scores, ppvs, labels = predict(seqs, model_path, threshold)

    print(f"{'Sequência':<15} {'Score':>6}  {'PPV':>6}  {'Label':<12}  Descrição")
    print("-" * 75)
    for (seq, desc), score, ppv, label in zip(TEST_PEPTIDES, scores, ppvs, labels):
        print(f"{seq:<15} {score:>6.3f}  {ppv:>6.3f}  {label:<12}  {desc}")
