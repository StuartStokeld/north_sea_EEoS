#!/usr/bin/env bash
# Clone the EEoS Python implementation (biomass.py) used via reticulate.
# Run from repo root: bash local_env_upload_28-7/download_scripts/02_clone_equation_of_state.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TARGET="$ROOT/equation_of_state"

if [ -f "$TARGET/biomass.py" ]; then
  echo "equation_of_state already present at $TARGET"
  exit 0
fi

if [ -d "$TARGET/.git" ]; then
  echo "Updating existing equation_of_state clone..."
  git -C "$TARGET" pull --ff-only
else
  echo "Cloning micbru/equation_of_state..."
  git clone --depth 1 https://github.com/micbru/equation_of_state.git "$TARGET"
fi

if [ ! -f "$TARGET/biomass.py" ]; then
  echo "ERROR: biomass.py not found in $TARGET"
  exit 1
fi
echo "OK: $TARGET/biomass.py"
