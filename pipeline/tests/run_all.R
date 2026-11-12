#!/usr/bin/env Rscript
# Run the whole test suite from the repo root:
#   Rscript pipeline/tests/run_all.R                    # fast tests
#   EEOS_SLOW_TESTS=1 Rscript pipeline/tests/run_all.R  # + slow tests (raw-CSV checks, all-haul solver check)
#
# Tests tagged [F-..] encode audit findings (audit/findings.md) and fail until those are fixed.
# Exit status: 0 if every test passes, 1 otherwise. Mutation testing: audit/scripts/mutation_test.R.

if (!requireNamespace("testthat", quietly = TRUE)) stop("Install testthat: renv::install('testthat')")
if (!file.exists("pipeline/tests/testthat")) stop("Run from the repo root: Rscript pipeline/tests/run_all.R")

res <- testthat::test_dir("pipeline/tests/testthat", reporter = testthat::SummaryReporter$new(max_reports = Inf),
                          stop_on_failure = FALSE)
df <- as.data.frame(res)
failed <- df[df$failed > 0 | df$error, c("file", "test")]
known <- grepl("^\\[F-", failed$test)
cat(sprintf("\nTests: %d | expectations failed in %d tests (%d tagged audit findings, %d untagged) | skipped: %d\n",
            nrow(df), nrow(failed), sum(known), sum(!known), sum(df$skipped)))
if (any(!known)) {
  cat("UNTAGGED failures (regressions or new problems):\n")
  cat(paste0("  ", failed$file[!known], ": ", failed$test[!known]), sep = "\n")
}
quit(status = if (nrow(failed) > 0) 1L else 0L)
