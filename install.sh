#!/bin/bash
set -e

echo ""
echo "=== TheraEPIflow — Bootstrap Installer ==="
echo ""

# ── Step 1: Check OS ──────────────────────────────────────────────────────────
if [[ "$OSTYPE" != "linux-gnu"* ]]; then
    echo "ERROR: This installer requires Linux or WSL2 on Windows."
    echo "On Windows, please install WSL2 first: https://aka.ms/wsl"
    exit 1
fi

# ── Step 2: Install git if missing ───────────────────────────────────────────
if ! command -v git &>/dev/null; then
    echo "Git not found. Installing Git..."
    echo ""
    if command -v apt-get &>/dev/null; then
        sudo apt-get update -qq && sudo apt-get install -y git
    elif command -v dnf &>/dev/null; then
        sudo dnf install -y git
    elif command -v yum &>/dev/null; then
        sudo yum install -y git
    else
        echo "ERROR: Could not install Git automatically."
        echo "Please install Git manually and re-run this script."
        exit 1
    fi
    echo ""
    echo "Git installed successfully."
else
    echo "Git found: $(git --version)"
fi

# ── Step 3: Clone repository ─────────────────────────────────────────────────
REPO_URL="https://github.com/Davirbe/TheraEpiFlow.git"
REPO_DIR="TheraEpiFlow"

if [ -d "$REPO_DIR" ]; then
    echo ""
    echo "Folder '$REPO_DIR' already exists. Skipping clone."
    echo "If you want a fresh install, delete the folder and run again."
else
    echo ""
    echo "Cloning TheraEPIflow repository..."
    git clone "$REPO_URL"
    echo "Repository cloned successfully."
fi

# ── Step 4: Run setup ─────────────────────────────────────────────────────────
echo ""
cd "$REPO_DIR"
bash setup.sh

# ── Step 5: Done ─────────────────────────────────────────────────────────────
echo ""
echo "=== Installation complete! ==="
echo ""
echo "To start using TheraEPIflow, run:"
echo "    cd TheraEpiFlow"
echo "    conda activate TheraEPIflow"
echo "    python main.py"
echo ""
