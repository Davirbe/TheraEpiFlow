"""File-descriptor level stdout/stderr capture.

Python's `contextlib.redirect_stdout` / `redirect_stderr` only redirect the
Python-level `sys.stdout` / `sys.stderr` objects. They do NOT intercept
output written by C extensions (e.g., TensorFlow's INFO logs from the C++
runtime). This module provides a context manager that swaps the underlying
OS file descriptors 1 and 2 for pipes, capturing everything regardless of
the layer that emitted it.

Typical use: silence MHCFlurry / TensorFlow / Keras during prediction and
display the captured output only when something goes wrong.

Caveats:
- Not thread-safe across the whole process: `os.dup2` mutates the global
  process file descriptor table. Use it around a serialized region (or run
  the noisy work in a dedicated thread that joins inside the `with` block).
- The captured stream may include ANSI color escape codes from the inner
  library; callers can pass it through `rich.text.Text.from_ansi` if they
  want to render it in the terminal.
"""

from __future__ import annotations

import contextlib
import os
import sys
import threading
from typing import Iterator


class CapturedOutput:
    """Holds the captured stdout and stderr text after the context exits."""

    def __init__(self) -> None:
        self.stdout: str = ""
        self.stderr: str = ""

    def is_empty(self) -> bool:
        return not self.stdout and not self.stderr


def _drain_pipe_into_list(read_fd: int, byte_chunks: list[bytes]) -> None:
    """Read until EOF from `read_fd`, appending raw bytes into `byte_chunks`."""
    while True:
        try:
            chunk = os.read(read_fd, 4096)
        except OSError:
            break
        if not chunk:
            break
        byte_chunks.append(chunk)


@contextlib.contextmanager
def capture_fd_output() -> Iterator[CapturedOutput]:
    """Capture both Python-level and C-level stdout/stderr via dup2 + pipes.

    Yields a `CapturedOutput` whose `.stdout` and `.stderr` attributes are
    populated when the context block exits.
    """
    captured_output = CapturedOutput()

    sys.stdout.flush()
    sys.stderr.flush()

    saved_stdout_fd: int = os.dup(1)
    saved_stderr_fd: int = os.dup(2)

    stdout_read_fd, stdout_write_fd = os.pipe()
    stderr_read_fd, stderr_write_fd = os.pipe()

    os.dup2(stdout_write_fd, 1)
    os.dup2(stderr_write_fd, 2)
    os.close(stdout_write_fd)
    os.close(stderr_write_fd)

    stdout_byte_chunks: list[bytes] = []
    stderr_byte_chunks: list[bytes] = []

    stdout_drainer = threading.Thread(
        target=_drain_pipe_into_list, args=(stdout_read_fd, stdout_byte_chunks), daemon=True
    )
    stderr_drainer = threading.Thread(
        target=_drain_pipe_into_list, args=(stderr_read_fd, stderr_byte_chunks), daemon=True
    )
    stdout_drainer.start()
    stderr_drainer.start()

    try:
        yield captured_output
    finally:
        sys.stdout.flush()
        sys.stderr.flush()

        os.dup2(saved_stdout_fd, 1)
        os.dup2(saved_stderr_fd, 2)
        os.close(saved_stdout_fd)
        os.close(saved_stderr_fd)

        stdout_drainer.join(timeout=2.0)
        stderr_drainer.join(timeout=2.0)

        try:
            os.close(stdout_read_fd)
        except OSError:
            pass
        try:
            os.close(stderr_read_fd)
        except OSError:
            pass

        captured_output.stdout = b"".join(stdout_byte_chunks).decode("utf-8", errors="replace")
        captured_output.stderr = b"".join(stderr_byte_chunks).decode("utf-8", errors="replace")
