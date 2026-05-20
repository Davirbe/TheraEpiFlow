"""Shared HTTP helper: a GET with exponential-backoff retry for transient
network errors. Used by fetch_sequences and search_variants (both UniProt)."""

import time

import requests

from utils.console import console

DEFAULT_TIMEOUT = 15


def http_get(url: str, params: dict = None, max_attempts: int = 3,
             timeout: int = DEFAULT_TIMEOUT) -> requests.Response:
    """GET with exponential-backoff retry for transient network errors."""
    last_err: Exception = RuntimeError("no attempts made")
    for attempt in range(max_attempts):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as err:
            last_err = err
            if attempt < max_attempts - 1:
                wait = 2 ** attempt
                console.print(
                    f"[yellow]⚠ HTTP request failed (attempt {attempt + 1}/{max_attempts}): "
                    f"{err} — retrying in {wait}s[/yellow]"
                )
                time.sleep(wait)
    raise last_err
