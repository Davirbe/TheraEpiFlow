"""
Utility for flexible column name resolution in DataFrames.

Adapted from analise_consenso_predicao.py (find_column_name function).
Used across multiple steps when column names may vary between prediction tools.
"""

import pandas as pd


def find_column_name(dataframe: pd.DataFrame, possible_names: list[str]) -> str | None:
    """
    Finds the first matching column name from a list of candidates.

    Args:
        dataframe:      The DataFrame to search in.
        possible_names: Ordered list of column name candidates to try.

    Returns:
        The first matching column name found, or None if none match.

    Example:
        col = find_column_name(df, ["EL_Rank", "el_rank", "%Rank_EL"])
    """
    for name in possible_names:
        if name in dataframe.columns:
            return name
    return None


def require_column(dataframe: pd.DataFrame, possible_names: list[str], context: str = "") -> str:
    """
    Like find_column_name, but raises ValueError if no match is found.

    Args:
        dataframe:      The DataFrame to search in.
        possible_names: Ordered list of column name candidates to try.
        context:        Optional description for the error message (e.g. step name).

    Returns:
        The first matching column name found.

    Raises:
        ValueError if no matching column is found.
    """
    result = find_column_name(dataframe, possible_names)
    if result is None:
        prefix = f"[{context}] " if context else ""
        raise ValueError(
            f"{prefix}None of the expected columns found: {possible_names}. "
            f"Available columns: {list(dataframe.columns)}"
        )
    return result
