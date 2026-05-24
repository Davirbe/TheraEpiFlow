"""Locale-robust CSV writer for every step's output.

`write_user_facing_csv` wraps pandas' `to_csv` with two settings that make the
resulting files open correctly on Excel / LibreOffice / Google Sheets in any
locale (US, UK, pt-BR, de-DE, fr-FR, ja-JP, …):

  - `encoding='utf-8-sig'` — emits the UTF-8 BOM so spreadsheet apps detect
    the encoding automatically on double-click instead of falling back to a
    locale-specific guess (which mangles accents and the ★ marker).
  - `quoting=csv.QUOTE_NONNUMERIC` — wraps every string cell in double quotes.
    The pipeline routinely emits cells whose contents include a semicolon
    (e.g. `HLA-A*02:01;HLA-B*07:02;…`); without quoting those get misread as
    column separators in locales where ';' is the default CSV delimiter,
    splitting one cell across many columns.

Numeric cells stay unquoted so spreadsheets parse them as numbers.

Usage:

    from utils.csv_write import write_user_facing_csv
    write_user_facing_csv(df, path)        # replaces df.to_csv(path, index=False)
    write_user_facing_csv(df, path, index=True)   # opt-in if you need the index

`utils.csv_write` is the single source of truth for CSV portability — every
new step that writes a CSV should import this helper instead of calling
`df.to_csv` directly.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


def write_user_facing_csv(
    dataframe: pd.DataFrame,
    output_path: str | Path,
    index: bool = False,
) -> None:
    """Writes `dataframe` to `output_path` with locale-robust settings.

    See module docstring for the rationale behind utf-8-sig + QUOTE_NONNUMERIC.
    """
    dataframe.to_csv(
        output_path,
        index    = index,
        quoting  = csv.QUOTE_NONNUMERIC,
        encoding = 'utf-8-sig',
    )
