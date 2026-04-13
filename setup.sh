#!/bin/bash
set -e

ENV_NAME="theraEPIflow"
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
MINICONDA_DEFAULT_PATH="$HOME/miniconda3"

echo ""
echo "=== TheraEPIflow — Setup ==="
echo ""

# ── Step 1: Find or install conda ────────────────────────────────────────────

# Locate conda binary — check common install locations
CONDA_BIN=""
for candidate in \
    "$(command -v conda 2>/dev/null)" \
    "$HOME/miniconda3/bin/conda" \
    "$HOME/anaconda3/bin/conda" \
    "/opt/conda/bin/conda"; do
    if [ -x "$candidate" ]; then
        CONDA_BIN="$candidate"
        break
    fi
done

if [ -z "$CONDA_BIN" ]; then
    echo "conda not found. Installing Miniconda..."
    echo ""

    curl -fsSL "$MINICONDA_URL" -o /tmp/miniconda_installer.sh
    bash /tmp/miniconda_installer.sh -b -p "$MINICONDA_DEFAULT_PATH"
    rm /tmp/miniconda_installer.sh

    CONDA_BIN="$MINICONDA_DEFAULT_PATH/bin/conda"

    # Initialize conda for future shell sessions
    "$CONDA_BIN" init bash
    echo ""
    echo "Miniconda installed at: $MINICONDA_DEFAULT_PATH"
else
    echo "conda found: $("$CONDA_BIN" --version)"
fi

# Source conda for this session
CONDA_BASE="$("$CONDA_BIN" info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"

# ── Step 2: Create or update conda environment ───────────────────────────────

echo ""
if "$CONDA_BIN" env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' already exists. Updating..."
    "$CONDA_BIN" env update -f environment.yml --prune
else
    echo "Creating environment '${ENV_NAME}'..."
    "$CONDA_BIN" env create -f environment.yml
fi

# ── Step 3: Done ─────────────────────────────────────────────────────────────

echo ""
echo "=== Setup complete! ==="
echo ""
echo "Activate the environment with:"
echo "    conda activate ${ENV_NAME}"
echo ""
echo "Then run the API tests with:"
echo "    python tests/test_apis.py"
echo ""
