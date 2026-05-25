# export_bundle

The final pipeline step. Packages the project as a tar.gz that can be handed
off (collaborator, deposit, paper supplement). Wraps `utils/archive.py` so the
bundling logic stays in one place.

> **Code layout:** single-file module (`__init__.py`) — small enough that the
> per-role split (see top-level README) would just fragment it. Splits along
> the standard seams if it grows.

## Where it sits in the pipeline

```
... → curate_murine → integrate_data → generate_report → export_bundle
```

It is the very last step: it consumes nothing the pipeline computes — it
serialises everything already on disk under `projects/{project}/` into a single
archive file.

## What it does

On every invocation the step asks interactively:

1. **Scope** — bundle the whole project, bundle one earlier step's outputs
   across every track, or cancel.
2. **Predictions opt-in** (full scope only) — include the heavy
   `data/intermediate/*/predictions/` folders. Excluded by default because
   the raw NetMHCpan/MHCflurry CSVs are large and trivially regeneratable.
3. **Destination** — always offers the in-project `downloads/` folder; adds
   `~/Downloads` when that folder exists; adds the Windows-side Downloads
   folder when running under WSL (auto-detected via `/proc/version` +
   `/mnt/c/Users/{user}` — falls back to the only real user dir if the
   Linux/Windows usernames differ).

It then calls `archive_project` or `archive_step` from `utils/archive.py`,
which uses the stdlib `tarfile` module (no extra dependency).

Re-running always generates a new timestamped archive — previous archives stay
on disk so the user can keep the exact build they tested.

## Inputs

None. Reads `projects/{project}/` directly.

## Outputs

- `{project}_full_{stamp}.tar.gz` — full project bundle (default scope).
- `{project}_{step_name}_{stamp}.tar.gz` — single-step bundle (per-step scope).

`{stamp}` is `YYYYMMDD_HHMMSS` (no colons → safe on Windows). Output path is
the destination folder the user chose.

## REPL shortcut

The `[z]` key in the interactive REPL menu reroutes to running this step
via `force_rerun=True` — same flow without having to navigate to the step in
sequence.

## Verification

Smoke-tested on `hpv16`:
- Full project bundle ≈ 3.7 MB (predictions excluded).
- Per-step bundles in the 50–2000 KB range depending on the step.
- WSL Windows Downloads resolved to `/mnt/c/Users/davi/Downloads` on the
  development box.

```bash
python main.py --project hpv16 --step export_bundle
```
