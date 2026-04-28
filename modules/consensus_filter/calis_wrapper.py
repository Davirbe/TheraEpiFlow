"""
Stage 2 of consensus_filter — IEDB Calis 2013 immunogenicity scoring.

Wrapper around predict_immunogenicity.Prediction (the official IEDB tool
copied to this folder on 27/04/2026 from
https://nextgen-tools.iedb.org/download-all). We call Prediction.predict()
directly as a Python object, capturing stdout to parse the result, instead of
running the script via subprocess. This avoids file I/O and process overhead
while still using the EXACT official implementation (same coefficients, same
algorithm, same masks) — important for scientific reproducibility and
citability ("IEDB Class I Immunogenicity standalone tool, Calis et al. 2013").

The allele→anchor-position mask table is reproduced here as a module constant
so we can resolve the mask without going through Prediction.validate(), which
expects a file path and calls sys.exit() on errors.

Public API:
    score_peptides_calis_immunogenicity(peptide_list, allele_imgt_format=None)
        → dict mapping peptide_sequence → calis_score (float)

Score interpretation per Calis et al. 2013:
    score > 0  → predicted immunogenic
    score ≤ 0  → predicted non-immunogenic
"""

import contextlib
import io


# Allele → anchor positions (1-based, comma-separated). Reproduced verbatim
# from predict_immunogenicity.py, Prediction.validate.allele_dict.
# When the allele is in this table, those positions are masked out of the
# Calis score (they are MHC-binding contacts, not TCR contacts and therefore
# uninformative about T-cell recognition). For alleles NOT in this table, the
# tool's default mask is used: positions 1, 2, C-terminal.
CALIS_ALLELE_ANCHOR_POSITION_TABLE: dict = {
    "H-2-Db":   "2,5,9", "H-2-Dd":   "2,3,5", "H-2-Kb":   "2,3,9",
    "H-2-Kd":   "2,5,9", "H-2-Kk":   "2,8,9", "H-2-Ld":   "2,5,9",
    "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9", "HLA-A0202": "1,2,9",
    "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9",
    "HLA-A0301": "1,2,9", "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9",
    "HLA-A2402": "2,7,9", "HLA-A2601": "1,2,9", "HLA-A2902": "2,7,9",
    "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9", "HLA-A3101": "1,2,9",
    "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9",
    "HLA-A6802": "1,2,9", "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9",
    "HLA-B0801": "2,5,9", "HLA-B1501": "1,2,9", "HLA-B1502": "1,2,9",
    "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9", "HLA-B3501": "1,2,9",
    "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9",
    "HLA-B4402": "2,3,9", "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9",
    "HLA-B4601": "1,2,9", "HLA-B5101": "1,2,9", "HLA-B5301": "1,2,9",
    "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9", "HLA-B5801": "1,2,9",
}


def normalize_allele_to_calis_format(allele_imgt_format: str) -> str:
    """HLA-A*02:01 → HLA-A0201 (Calis tool format — no asterisks, no colons)."""
    return allele_imgt_format.replace("*", "").replace(":", "")


def is_allele_supported_by_calis(allele_imgt_format: str) -> bool:
    """
    Returns True if Calis has a specific anchor-position mask for this allele.
    For unsupported alleles the tool falls back to default mask (1, 2, C-term),
    which is a valid scoring path but slightly less precise.
    """
    return (
        normalize_allele_to_calis_format(allele_imgt_format)
        in CALIS_ALLELE_ANCHOR_POSITION_TABLE
    )


def score_peptides_calis_immunogenicity(
    peptide_list: list,
    allele_imgt_format: str = None,
) -> dict:
    """
    Score a list of peptides using the IEDB Calis 2013 standalone tool.

    Args:
        peptide_list:        sequences (uppercase, standard amino acids only)
        allele_imgt_format:  IMGT-format allele name (e.g. "HLA-A*02:01").
                             If None or not in the supported set, default mask
                             is used.

    Returns:
        dict mapping peptide_sequence → immunogenicity_score (float).
        Peptides that error out during scoring (e.g. invalid amino acid) are
        absent from the dict.

    Implementation notes:
        - We bypass Prediction.validate() entirely (it reads from a file path
          and calls sys.exit on error).
        - cleaned_data is the tuple shape that Prediction.predict() expects:
          (sequence_text_list, custom_mask_string_or_None, allele_or_None).
        - stdout is captured because Prediction.predict() prints its results
          rather than returning them.
    """
    from .predict_immunogenicity import Prediction

    allele_in_calis_format = None
    anchor_position_mask   = None
    if allele_imgt_format:
        normalized_allele = normalize_allele_to_calis_format(allele_imgt_format)
        if normalized_allele in CALIS_ALLELE_ANCHOR_POSITION_TABLE:
            allele_in_calis_format = normalized_allele
            anchor_position_mask   = CALIS_ALLELE_ANCHOR_POSITION_TABLE[normalized_allele]

    cleaned_data_for_predict = (
        list(peptide_list),
        anchor_position_mask,
        allele_in_calis_format,
    )

    prediction_engine = Prediction()
    captured_stdout_buffer = io.StringIO()

    with contextlib.redirect_stdout(captured_stdout_buffer):
        prediction_engine.predict(cleaned_data_for_predict)

    return _parse_calis_stdout_into_score_dict(captured_stdout_buffer.getvalue())


def _parse_calis_stdout_into_score_dict(raw_stdout_text: str) -> dict:
    """
    Prediction.predict() prints, in order:
      (optional)  allele: HLA-A0201
                  masking: custom (or default)
                  masked variables: [1, 2, 'cterm']
                  (blank line)
                  peptide,length,score
                  PEPTIDE,9,0.21243
                  PEPTIDE,9,0.15000

    Walk lines, skip metadata + header, parse data rows.
    """
    peptide_to_score = {}
    for raw_line in raw_stdout_text.splitlines():
        stripped_line = raw_line.strip()
        if not stripped_line:
            continue
        if stripped_line == 'peptide,length,score':
            continue
        if (stripped_line.startswith('masking:')
                or stripped_line.startswith('masked variables:')
                or stripped_line.startswith('allele:')):
            continue

        line_fields = stripped_line.split(',')
        if len(line_fields) != 3:
            continue

        peptide_sequence, _length_text, score_text = line_fields
        try:
            peptide_to_score[peptide_sequence] = float(score_text)
        except ValueError:
            continue

    return peptide_to_score
