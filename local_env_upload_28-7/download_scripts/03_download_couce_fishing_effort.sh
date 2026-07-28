#!/usr/bin/env bash
# Download Couce et al. (2020) North Sea trawling effort CSV (~2.5 MB).
# Source: Cefas Data Hub DOI 10.14466/CefasDataHub.61 (recordset 10953)
# Run from repo root: bash local_env_upload_28-7/download_scripts/03_download_couce_fishing_effort.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="$ROOT/data/external/couce_trawling_effort"
OUT_FILE="$OUT_DIR/NorthSea_trawling_effort_1985to2015_REVIEW_v2.csv"
URL="https://www.cefas.co.uk/data-and-publications/data-hub/api/export/10953?format=csv"

mkdir -p "$OUT_DIR"

if [ -f "$OUT_FILE" ]; then
  echo "Couce CSV already present: $OUT_FILE"
  exit 0
fi

echo "Downloading Couce fishing effort CSV..."
curl -fsSL "$URL" -o "$OUT_FILE"
echo "OK: $OUT_FILE ($(du -h "$OUT_FILE" | cut -f1))"
