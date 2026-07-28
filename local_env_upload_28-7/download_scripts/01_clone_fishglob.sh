#!/usr/bin/env bash
# Clone FishGlob and place NS-IBTS cleaned data where the pipeline expects it.
# Run from repo root: bash local_env_upload_28-7/download_scripts/01_clone_fishglob.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TARGET="$ROOT/FishGlob_data"

if [ -d "$TARGET/outputs/Cleaned_data/NS-IBTS_clean.RData" ] || \
   [ -f "$TARGET/outputs/Cleaned_data/NS-IBTS_clean.RData" ]; then
  echo "FishGlob NS-IBTS data already present at $TARGET"
  exit 0
fi

if [ -d "$TARGET/.git" ]; then
  echo "Updating existing FishGlob clone..."
  git -C "$TARGET" pull --ff-only
else
  echo "Cloning FishGlob_data..."
  git clone --depth 1 https://github.com/fishglob/FishGlob_data.git "$TARGET"
fi

RDATA="$TARGET/outputs/Cleaned_data/NS-IBTS_clean.RData"
if [ ! -f "$RDATA" ]; then
  echo "ERROR: Expected file missing after clone: $RDATA"
  exit 1
fi
echo "OK: $RDATA ($(du -h "$RDATA" | cut -f1))"
