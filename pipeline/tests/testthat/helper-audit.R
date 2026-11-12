# Shared helpers for the production-audit test suite (auto-loaded by testthat).
#
# Paths: tests run with the working directory at pipeline/tests/testthat. The pipeline code
# directory can be redirected with EEOS_PIPELINE_DIR (used by audit/scripts/mutation_test.R to
# run the suite against a mutated temporary copy); data and outputs always come from the repo.
#
# Test naming: tests that encode a confirmed audit finding carry its ID ("[F-A1]" etc., see
# audit/findings.md). They are written to the correct specification and are EXPECTED TO FAIL
# until the finding is fixed.

audit_repo_root <- function() {
  normalizePath(file.path(testthat::test_path(), "..", "..", ".."), mustWork = TRUE)
}

audit_pipeline_dir <- function() {
  d <- Sys.getenv("EEOS_PIPELINE_DIR", "")
  if (nzchar(d)) normalizePath(d, mustWork = TRUE) else file.path(audit_repo_root(), "pipeline")
}

#' Source pipeline/R/<file> (from EEOS_PIPELINE_DIR when set) into the calling environment.
source_pipeline <- function(file, env = parent.frame()) {
  sys.source(file.path(audit_pipeline_dir(), "R", file), envir = env)
}

#' Path to a live script (respects EEOS_PIPELINE_DIR).
pipeline_script <- function(name) file.path(audit_pipeline_dir(), name)

repo_path <- function(...) file.path(audit_repo_root(), ...)

skip_if_no_file <- function(...) {
  p <- repo_path(...)
  testthat::skip_if_not(file.exists(p), paste("missing input:", p))
  p
}

skip_if_not_slow <- function() {
  testthat::skip_if_not(identical(Sys.getenv("EEOS_SLOW_TESTS"), "1"),
                        "slow test: set EEOS_SLOW_TESTS=1")
}

#' Evaluate one top-level `name <- function(...)` (or `name <- value`) definition from a script
#' without running the script. Lets tests check functions defined inside live scripts.
extract_definition <- function(script, name, env = parent.frame()) {
  exprs <- as.list(parse(script, keep.source = FALSE))
  is_def <- function(e) {
    is.call(e) && length(e) == 3L && as.character(e[[1]])[1] %in% c("<-", "=") &&
      is.name(e[[2]]) && identical(as.character(e[[2]]), name)
  }
  hit <- Filter(is_def, exprs)
  if (length(hit) == 0L) {
    # search nested blocks (if / {), but not inside function bodies
    walk <- function(e) {
      if (is_def(e)) return(list(e))
      if (!is.call(e) || identical(e[[1]], as.name("function"))) return(list())
      unlist(lapply(as.list(e)[-1], function(x) if (is.call(x)) walk(x) else list()), recursive = FALSE)
    }
    hit <- unlist(lapply(exprs, walk), recursive = FALSE)
  }
  if (length(hit) == 0L) stop("No definition of ", name, " in ", script)
  eval(hit[[length(hit)]], envir = env)
  invisible(get(name, envir = env))
}

#' Cached readers for the large committed inputs (one read per test run).
.audit_cache <- new.env(parent = emptyenv())
cached_rds <- function(...) {
  p <- skip_if_no_file(...)
  if (is.null(.audit_cache[[p]])) .audit_cache[[p]] <- readRDS(p)
  .audit_cache[[p]]
}

#' Haul key (Survey_Year_Quarter_Country_Ship_HaulNo) from a FishGlob haul_id.
fishglob_key <- function(haul_id) {
  vapply(strsplit(haul_id, " "), function(p) paste(p[c(1:5, 8)], collapse = "_"), "")
}

#' Queen-contiguity neighbours of ICES rectangle codes (e.g. "38F1"): rows differ by <= 1 and
#' columns (letter + digit, with the A-M letter blocks of 10 columns; no "I") differ by <= 1.
ices_rect_rowcol <- function(code) {
  row <- as.integer(substr(code, 1, 2))
  letter <- substr(code, 3, 3)
  digit <- as.integer(substr(code, 4, 4))
  letters_ices <- c("A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L", "M")
  # A has columns 0-3 only (A0-A3); B..M have 0-9. Column index: A0 = 0 ... A3 = 3, B0 = 4 ...
  li <- match(letter, letters_ices)
  col <- ifelse(li == 1L, digit, 4L + (li - 2L) * 10L + digit)
  cbind(row = row, col = col)
}
