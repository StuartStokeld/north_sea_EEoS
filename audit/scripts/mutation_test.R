# Mutation testing (audit layer 4.8, scoped to 4 mutants): does the suite catch known bugs?
#
# Each mutant is applied to a temporary copy of pipeline/ (never to the repo). The suite runs against
# the copy via EEOS_PIPELINE_DIR. A mutant is KILLED when at least one test that passes on the
# reference copy fails on the mutant.
#
# Two of the mutants describe the current code (LngtCode "5" read as mm; signed residual). For those
# the reference copy is the CORRECTED code (the fix applied in the temporary copy only), so the
# question is whether the suite would catch the bug being (re)introduced.
#
# Run from the repo root: Rscript audit/scripts/mutation_test.R   (about 10-15 min)

suppressPackageStartupMessages(library(testthat))

edit_file <- function(dir, rel, old, new) {
  p <- file.path(dir, rel)
  x <- paste(readLines(p, warn = FALSE), collapse = "\n")
  if (!grepl(old, x, fixed = TRUE)) stop("mutation site not found in ", rel, ": ", old)
  writeLines(sub(old, new, x, fixed = TRUE), p)
}

fix_lngt <- list(rel = "R/datras_hl_helpers.R",
  old = '  is_mm <- code %in% c("0", ".", "5")\n  out[is_mm] <- (cls[is_mm] + 0.5) / 10',
  new = '  is_dot <- code == "."\n  out[is_dot] <- (cls[is_dot] + 0.5) / 10\n  is0 <- code == "0"\n  out[is0] <- (cls[is0] + 2.5) / 10\n  is5 <- code == "5"\n  out[is5] <- cls[is5] + 2.5')
fix_resid <- list(
  list(rel = "run_h2h3_phase_v2_refit.R", old = "formula_v2 <- residual ~", new = "formula_v2 <- abs_residual ~"),
  list(rel = "run_h2h3_phase_v2_reporting.R", old = "formula_car <- residual ~", new = "formula_car <- abs_residual ~")
)

mutants <- list(
  list(id = "M1", what = "remove the C x HaulDur/60 conversion", reference_fixes = list(),
       edits = list(list(rel = "R/datras_hh_helpers.R", old = 'if_else(DataType == "C", HaulDur / 60, 1)', new = "1"))),
  list(id = "M2", what = 'treat LngtCode "5" as mm (reference: ICES-correct midpoints)', reference_fixes = list(fix_lngt),
       edits = list(list(rel = "R/datras_hl_helpers.R", old = "  is5 <- code == \"5\"\n  out[is5] <- cls[is5] + 2.5",
                         new = "  is5 <- code == \"5\"\n  out[is5] <- (cls[is5] + 0.5) / 10"))),
  list(id = "M3", what = "signed instead of absolute residual (reference: absolute-residual formulas)",
       reference_fixes = fix_resid,
       edits = list(list(rel = "run_h2h3_phase_v2_reporting.R", old = "formula_car <- abs_residual ~", new = "formula_car <- residual ~"))),
  list(id = "M4", what = "FP_within not centred", reference_fixes = list(),
       edits = list(list(rel = "R/h2h3_within_between_helpers.R", old = "FP_within = log_hours_total - FP_between",
                         new = "FP_within = log_hours_total")))
)

run_suite <- function(pipeline_dir) {
  Sys.setenv(EEOS_PIPELINE_DIR = pipeline_dir)
  on.exit(Sys.unsetenv("EEOS_PIPELINE_DIR"))
  r <- as.data.frame(test_dir("pipeline/tests/testthat", reporter = SilentReporter$new(), stop_on_failure = FALSE))
  setNames(r$failed == 0 & !r$error, paste(r$file, r$test, sep = " :: "))
}

copy_pipeline <- function(fixes) {
  d <- file.path(tempfile("pipeline_mut_"), "pipeline")
  dir.create(d, recursive = TRUE)
  file.copy(list.files("pipeline", full.names = TRUE, all.files = FALSE), d, recursive = TRUE)
  unlink(file.path(d, "tests"), recursive = TRUE)
  for (f in fixes) edit_file(d, f$rel, f$old, f$new)
  d
}

ref_cache <- list()
rows <- list()
only <- strsplit(Sys.getenv("EEOS_MUTANTS", ""), ",")[[1]]
if (length(only)) mutants <- Filter(function(m) m$id %in% only, mutants)
for (m in mutants) {
  key <- paste(c("ref", vapply(m$reference_fixes, `[[`, "", "rel")), collapse = "+")
  if (is.null(ref_cache[[key]])) {
    cat("Reference run:", key, "\n")
    ref_cache[[key]] <- run_suite(copy_pipeline(m$reference_fixes))
  }
  ref <- ref_cache[[key]]
  d <- copy_pipeline(m$reference_fixes)
  for (e in m$edits) edit_file(d, e$rel, e$old, e$new)
  cat("Mutant", m$id, ":", m$what, "\n")
  mut <- run_suite(d)
  common <- intersect(names(ref), names(mut))
  killers <- common[ref[common] & !mut[common]]
  rows[[m$id]] <- data.frame(id = m$id, mutant = m$what, killed = length(killers) > 0,
                             n_killing_tests = length(killers),
                             killing_tests = paste(head(killers, 6), collapse = " | "))
}
res <- do.call(rbind, rows)
print(res[, c("id", "killed", "n_killing_tests")])
dir.create("audit/outputs/mutation", recursive = TRUE, showWarnings = FALSE)
out_csv <- "audit/outputs/mutation/mutation_results.csv"
if (length(only) && file.exists(out_csv)) {
  prev <- read.csv(out_csv)
  res <- rbind(prev[!prev$id %in% res$id, ], res)
  res <- res[order(res$id), ]
}
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("Kill rate: %d / %d\n", sum(res$killed), nrow(res)))
