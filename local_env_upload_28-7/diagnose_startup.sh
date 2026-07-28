#!/usr/bin/env bash
# Diagnostic script - does NOT modify .Rprofile, renv.lock, or move data.
# Includes timeout guards (no `timeout`/`gtimeout` binary available on this
# machine) so a hang in any single step doesn't hang the whole script.
set +e

run_timeout() {
  local secs="$1"; shift
  "$@" &
  local cmd_pid=$!
  (
    sleep "$secs"
    if kill -0 "$cmd_pid" 2>/dev/null; then
      echo ">>> TIMEOUT after ${secs}s - killing pid $cmd_pid <<<"
      kill -9 "$cmd_pid" 2>/dev/null
    fi
  ) &
  local watcher_pid=$!
  wait "$cmd_pid" 2>/dev/null
  local status=$?
  kill -9 "$watcher_pid" 2>/dev/null
  wait "$watcher_pid" 2>/dev/null
  return $status
}

echo "== 1. PROJECT SIZE BREAKDOWN =="
du -sh ./* 2>/dev/null | sort -rh | head -20

echo ""
echo "== 2. NESTED GIT REPOS INSIDE PROJECT ROOT =="
find . -maxdepth 3 -name ".git" -not -path "./.git" 2>/dev/null

echo ""
echo "== 3. .Rprofile CONTENTS =="
if [ -f .Rprofile ]; then cat .Rprofile; else echo "none"; fi

echo ""
echo "== 4. renv SETTINGS =="
if [ -f renv/settings.json ]; then cat renv/settings.json; fi
echo "renv.lock package count:"
grep -c '"Package"' renv.lock 2>/dev/null

echo ""
echo "== 5. RAW SHELL STARTUP TIME =="
time run_timeout 20 zsh -i -c 'echo shell ready'

echo ""
echo "== 6. BARE R STARTUP (no renv, no reticulate) =="
time run_timeout 20 Rscript --vanilla -e 'invisible(1)'

echo ""
echo "== 7. R STARTUP WITH renv ACTIVATED =="
time run_timeout 60 Rscript -e 'invisible(1)'

echo ""
echo "== 8. RETICULATE / PYTHON INIT TIME (isolated) =="
time run_timeout 60 Rscript -e '
  t0 <- Sys.time(); library(reticulate)
  cat("reticulate loaded:", Sys.time() - t0, "\n")
  t1 <- Sys.time()
  tryCatch({ py_config(); cat("py_config():", Sys.time() - t1, "\n") },
           error = function(e) cat("py_config() FAILED:", conditionMessage(e), "\n"))
'

echo ""
echo "== 9. LARGEST .RData / .rds FILES =="
find . -name "*.RData" -o -name "*.rds" 2>/dev/null | xargs du -h 2>/dev/null | sort -rh | head -10

echo ""
echo "== 10. renv CACHE HEALTH =="
run_timeout 20 Rscript -e 'cat("cache:", renv::paths$cache(), "\n"); cat("root:", renv::paths$root(), "\n")' 2>/dev/null

echo ""
echo "== DONE =="
