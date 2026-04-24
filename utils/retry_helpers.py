"""
Retry helpers for flaky network calls.

Used by any step that hits external APIs (NCBI Entrez, IEDB, ToxinPred, etc.).
The idea is to make transient failures recoverable without crashing the pipeline:
  - Network hiccups
  - HTTP 429 (rate limit)
  - Server timeouts
  - Partial read errors from Bio.Entrez

Content errors (zero results, invalid input) are NOT retried — those are user-facing
problems the step must handle interactively.
"""

import time
import urllib.error
import http.client
import socket
from typing import Callable, TypeVar

from rich.console import Console

console = Console()

# ── Exceptions that are worth retrying ────────────────────────────────────────
# These are transient (network, rate limit, server hiccup) — retrying often works.
# Do NOT add ValueError, RuntimeError, etc. here — those are usually real bugs
# or content-level problems the step should surface to the user.

TRANSIENT_NETWORK_EXCEPTIONS = (
    urllib.error.URLError,
    urllib.error.HTTPError,
    http.client.HTTPException,
    socket.timeout,
    socket.gaierror,
    ConnectionError,
    TimeoutError,
)

ReturnType = TypeVar('ReturnType')


def retry_network_call(
    callable_to_execute: Callable[[], ReturnType],
    operation_description: str = 'network call',
    max_attempts: int = 3,
    initial_backoff_seconds: float = 2.0,
    backoff_multiplier: float = 2.0,
) -> ReturnType:
    """
    Executes a callable with retry on transient network failures.

    Retries are spaced with exponential backoff:
      attempt 1 fails → wait 2s → attempt 2 fails → wait 4s → attempt 3 fails → raise

    Args:
      callable_to_execute     — zero-argument callable (use lambda to bind arguments)
      operation_description   — short label shown in retry messages (e.g. "NCBI esearch")
      max_attempts            — total attempts including the first (default: 3)
      initial_backoff_seconds — first wait duration after the first failure
      backoff_multiplier      — each subsequent wait is multiplied by this factor

    Returns:
      Whatever the callable returns.

    Raises:
      The last caught transient exception if all attempts fail.
      Non-transient exceptions are raised immediately (no retry).
    """
    current_backoff_seconds = initial_backoff_seconds
    last_caught_transient_error: Exception = RuntimeError(
        f'{operation_description} failed with no attempts'
    )

    for current_attempt_number in range(1, max_attempts + 1):
        try:
            return callable_to_execute()

        except TRANSIENT_NETWORK_EXCEPTIONS as transient_network_error:
            last_caught_transient_error = transient_network_error
            is_final_attempt = current_attempt_number == max_attempts

            if is_final_attempt:
                console.print(
                    f'[red]✗ {operation_description} failed after {max_attempts} '
                    f'attempts: {transient_network_error}[/red]'
                )
                raise

            console.print(
                f'[yellow]⚠ {operation_description} failed '
                f'(attempt {current_attempt_number}/{max_attempts}): '
                f'{transient_network_error}[/yellow]'
            )
            console.print(
                f'[dim]   Retrying in {current_backoff_seconds:.0f}s...[/dim]'
            )
            time.sleep(current_backoff_seconds)
            current_backoff_seconds *= backoff_multiplier

    # Unreachable — loop above either returns or raises
    raise last_caught_transient_error
