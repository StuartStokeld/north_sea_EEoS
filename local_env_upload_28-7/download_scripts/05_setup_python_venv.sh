#!/usr/bin/env bash
# Recreate the Python venv for reticulate / EEoS (not uploaded — 217 MB locally).
# Run from repo root: bash local_env_upload_28-7/download_scripts/05_setup_python_venv.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REQ="$ROOT/local_env_upload_28-7/python_requirements.txt"

python3 -m venv "$ROOT/.venv"
"$ROOT/.venv/bin/pip" install --upgrade pip

if [ -f "$REQ" ]; then
  "$ROOT/.venv/bin/pip" install -r "$REQ"
else
  "$ROOT/.venv/bin/pip" install numpy scipy pandas
fi

echo "OK: .venv ready. In R: reticulate::use_virtualenv('$ROOT/.venv', required = TRUE)"
