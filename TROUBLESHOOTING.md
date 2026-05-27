# TheraEpiFlow — Troubleshooting

Known issues observed when colleagues installed and ran the tool, with
diagnosed cause and the workaround that actually unblocked them. Open an
issue on GitHub if you hit something not listed here.

---

## 1. `ModuleNotFoundError: No module named 'requests'` (or IEDB calls fail to connect)

**Symptom.** `predict_binding`, `consensus_filter` or the validation
suite's IEDB query module raises `ModuleNotFoundError: requests`, or the
HTTP request succeeds against `http://` but is rejected when the URL is
upgraded to `https://`.

**Cause.** The `requests` library was missing in the active environment, or
the user installed dependencies into the wrong conda env. A separate
historical bug had a few IEDB endpoints reachable only via `http://`; today
all endpoints used by TheraEpiFlow require `https://`.

**Workaround.**
```bash
conda activate TheraEpiFlow
pip install 'requests>=2.31'
```
Then re-run the failing step.

**Permanent fix.** `environment.yml` already pins `requests>=2.31` under
the `pip:` block. If the error returns, it means an unrelated step's
install dropped the package — re-create the env from scratch with
`bash setup.sh`.

---

## 2. `conda: command not found` after running `setup.sh`

**Symptom.** Immediately after `bash setup.sh` finishes, running
`conda activate TheraEpiFlow` returns `conda: command not found`.

**Cause.** `setup.sh` installs conda under `~/miniconda3/` but the shell
that called the script never sourced the new conda init script, so
`PATH` does not include miniconda's `bin/`. Common on shells other than
bash (zsh, fish) where `~/.bashrc` is never read.

**Workaround.**
```bash
# bash users
~/miniconda3/bin/conda init bash
exec bash
# zsh users
~/miniconda3/bin/conda init zsh
exec zsh
```

**Permanent fix on the roadmap.** `setup.sh` should detect `$SHELL` and
invoke the corresponding `conda init <shell>` automatically.

---

## 3. `setup.sh` installs a second conda when one already exists

**Symptom.** User already had Anaconda installed under `/opt/anaconda3`
or via Homebrew, but `setup.sh` still installed Miniconda into
`~/miniconda3/`. Now `which conda` resolves to the new install and the
old envs are invisible.

**Cause.** `setup.sh` looks for conda in four paths (`$PATH`,
`~/miniconda3`, `~/anaconda3`, `/opt/conda`) but not in the additional
locations Anaconda/Homebrew use.

**Workaround.** Tell setup explicitly where conda lives:
```bash
export CONDA_PATH=/opt/anaconda3   # or /usr/local/Caskroom/miniconda/base
bash setup.sh
```

**Permanent fix on the roadmap.** `setup.sh` should first try
`which conda` and fall back to the hard-coded paths only when that
returns nothing.

---

## 4. WSL2: env not detected until the WSL distro is restarted

**Symptom.** On the very first `setup.sh` run inside WSL Ubuntu,
`conda activate TheraEpiFlow` works in a fresh terminal but fails in
the one that ran setup. Re-sourcing `~/.bashrc` does not help; only
closing every WSL window and reopening one does.

**Cause.** WSL2 caches the per-session shell environment more
aggressively than native Linux. The conda init lines appended to
`~/.bashrc` are not re-read by an already-running session.

**Workaround.**
1. Close every WSL Ubuntu terminal.
2. From a PowerShell or CMD window: `wsl --shutdown`.
3. Wait ~5 seconds, reopen WSL.
4. `conda activate TheraEpiFlow` now works.

**Permanent fix.** This is a WSL behavior, not something we can change.
We document it as a mandatory post-install step in this file and in the
WSL section of the README.

---

## 5. WSL2: `~/Downloads` is invisible to Windows Explorer

**Symptom.** User exports a tar.gz to `~/Downloads/` using
`export_bundle`, then can't find it in Windows Explorer (only sees the
Windows `C:\Users\<name>\Downloads` folder).

**Cause.** In WSL, `~` resolves to `/home/<user>` on the Linux side of
the file system. Windows Explorer reads from `C:\` (mounted in WSL as
`/mnt/c/`). The two trees do not overlap unless you cross the bridge
explicitly.

**Workaround.** Use the built-in WSL detection in `export_bundle` —
it auto-detects WSL via `/proc/version` and offers
`/mnt/c/Users/<user>/Downloads/` as a destination so the archive lands
where Windows can see it. From the REPL, the `[z]` key shortcut runs
export_bundle directly and shows the destination chooser.

If you're scripting around this, build the path yourself:
```bash
WINDOWS_DOWNLOADS=/mnt/c/Users/$(cmd.exe /c "echo %USERNAME%" 2>/dev/null | tr -d '\r')/Downloads
```

**Permanent fix.** The README's WSL section explicitly flags this.
`export_bundle` already does the right thing automatically.

---

## 6. WSL2 only sees ~4 GiB RAM even though the host has 8 GiB

**Symptom.** Pipeline crashes with out-of-memory during
`predict_binding` or `consensus_filter` on a host advertised as having
8 GiB RAM. `free -h` inside WSL shows only `~3.8 Gi` total.

**Cause.** WSL2 by default allocates about half of host RAM to the
Linux VM. The 8 GiB host gives WSL ~4 GiB; large NetMHCpan+MHCFlurry
batches can outrun that.

**Workaround.** Create or edit `C:\Users\<you>\.wslconfig` from the
Windows side:
```ini
[wsl2]
memory=7GB
swap=4GB
```
Then run `wsl --shutdown` and reopen WSL. `free -h` should now show the
new total.

**Permanent fix.** Not something the tool can change — it is a per-host
Windows/WSL setting. Documented here so future users know to bump it
before running the heavier benchmark presets (`sars_nucleo_tr3`).

---

## Reporting a new issue

If you hit a crash not described here, please file a GitHub issue with:

- Output of `python -m tests.validation.hardware_probe --out /tmp/hw.json`
  and paste `/tmp/hw.json` into the issue.
- Last ~50 lines of the failing run.
- Exact command you ran.
- Whether you are on native Linux, macOS, or WSL.

This helps us extend this document instead of you tracking down the
problem alone the next time.
