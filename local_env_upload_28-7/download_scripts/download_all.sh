#!/usr/bin/env bash
# Fetch all external dependencies that fit in git + document the rest.
# Run from repo root: bash local_env_upload_28-7/download_scripts/download_all.sh
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$DIR/../.." && pwd)"
cd "$ROOT"

echo "=== 1/4 FishGlob ==="
bash "$DIR/01_clone_fishglob.sh"

echo "=== 2/4 equation_of_state ==="
bash "$DIR/02_clone_equation_of_state.sh"

echo "=== 3/4 Couce fishing effort ==="
bash "$DIR/03_download_couce_fishing_effort.sh"

echo "=== 4/4 DATRAS HL (manual) ==="
echo "See $DIR/04_download_datras_hl.md — CSV not auto-downloaded (251 MB, ICES portal)."
if [ -f "$ROOT/outputs/datras_hl_raw.rds" ]; then
  echo "  outputs/datras_hl_raw.rds is present — H1 can run without the raw CSV."
else
  echo "  WARNING: datras_hl_raw.rds missing — download DATRAS or restore RDS from backup."
fi

echo ""
echo "Optional: bash $DIR/05_setup_python_venv.sh"
echo "Optional: R — open north_sea_eeos.Rproj and renv::restore()"
