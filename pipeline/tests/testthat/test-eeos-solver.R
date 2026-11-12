# Layer 4.4: EEoS solver (equation_of_state/biomass.py) through the pipeline's reticulate path.
#
# Invariants relied on (METE units: each individual has metabolic rate eps in [1, E], mass eps^(4/3)):
#   lower: B >= N * (E/N)^(4/3)  (Jensen: E[eps^(4/3)] >= E[eps]^(4/3), with E[eps] = E/N)
#   upper: B <= E^(4/3)          (eps^(4/3) <= eps * E^(1/3), summed over individuals)
# Both hold only when the solver satisfies its constraints, so they also detect failed solves.

suppressPackageStartupMessages(library(dplyr))

eeos_ready <- function() {
  py <- repo_path(".venv", "bin", "python")
  bm <- repo_path("equation_of_state", "biomass.py")
  testthat::skip_if_not(file.exists(py) && file.exists(bm), "EEoS python env or biomass.py missing")
  if (is.null(.audit_cache$eeos)) {
    reticulate::use_python(py, required = TRUE)
    reticulate::source_python(bm, envir = globalenv())
    # pipeline wrapper, copied verbatim from pipeline/build_eeos_predictions.R
    reticulate::py_run_string("
def predict_biomass(S, N, E):
    import pandas as pd
    s = pd.Series({'S': float(S), 'N': float(N), 'E': float(E)})
    return biomass(s)
def constraint_dev(S, N, E):
    import pandas as pd, numpy as np
    s = pd.Series({'S': float(S), 'N': float(N), 'E': float(E)})
    return float(np.max(np.abs(constraints(mete_lambdas(s, thresh=False), s))))
")
    .audit_cache$eeos <- TRUE
  }
  reticulate::py
}

test_that("PROVENANCE (change detector): biomass.py is byte-identical to upstream a96f406", {
  p <- skip_if_no_file("equation_of_state", "biomass.py")
  expect_identical(digest::digest(file = p, algo = "sha256"),
                   "66069c2763b310d3454a5b3730303abfcc6432781d15ef2995f24d12cb7d9a05")
})

test_that("pipeline wrapper matches build_eeos_predictions.R", {
  src <- readLines(pipeline_script("build_eeos_predictions.R"))
  expect_true(any(grepl("def predict_biomass(S, N, E):", src, fixed = TRUE)))
  expect_true(any(grepl("return biomass(s)", src, fixed = TRUE)))
})

test_that("GOLDEN: 42 Harte et al. (2022) communities reproduce upstream pBnum (rel. tol 1e-8)", {
  p <- skip_if_no_file("equation_of_state", "data_statevariables.csv")
  expect_identical(digest::digest(file = p, algo = "sha256"),
                   "cf32ab0c5f65fd61985072ddab170df255a90adcd97d1b6c76cdbe52f910b1d5")
  py <- eeos_ready()
  g <- read.csv(p, check.names = FALSE)
  expect_equal(nrow(g), 42L)
  b <- vapply(seq_len(nrow(g)), function(i) py$predict_biomass(g$S[i], g$N[i], g$E[i]), 0)
  expect_lt(max(abs(b / g$pBnum - 1)), 1e-8)
})

test_that("R (reticulate) and standalone Python give identical results", {
  py <- eeos_ready()
  cases <- data.frame(S = c(5, 20, 40), N = c(50, 800, 5000), E = c(400, 9e4, 2e6))
  r_side <- vapply(seq_len(nrow(cases)), function(i) py$predict_biomass(cases$S[i], cases$N[i], cases$E[i]), 0)
  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys, pandas as pd",
    sprintf("sys.path.insert(0, %s)", shQuote(repo_path("equation_of_state"))),
    "from biomass import biomass",
    sprintf("cases = [(%s)]", paste(sprintf("(%g, %g, %g)", cases$S, cases$N, cases$E), collapse = ", ")),
    "for S, N, E in cases[0]:",
    "    print(repr(float(biomass(pd.Series({'S': float(S), 'N': float(N), 'E': float(E)})))))"
  ), script)
  out <- system2(repo_path(".venv", "bin", "python"), c("-I", script), stdout = TRUE, stderr = FALSE)
  expect_equal(as.numeric(tail(out, 3)), r_side, tolerance = 1e-12)
})

test_that("boundary inputs (S = 2; N = S; E just above N) give finite output within the invariants", {
  py <- eeos_ready()
  cases <- data.frame(S = c(2, 10, 10), N = c(40, 10, 100), E = c(400, 1000, 100 * (1 + 1e-3)))
  for (i in seq_len(nrow(cases))) {
    b <- py$predict_biomass(cases$S[i], cases$N[i], cases$E[i])
    dev <- py$constraint_dev(cases$S[i], cases$N[i], cases$E[i])
    lo <- cases$N[i] * (cases$E[i] / cases$N[i])^(4 / 3)
    info <- sprintf("S=%g N=%g E=%g B=%g dev=%g", cases$S[i], cases$N[i], cases$E[i], b, dev)
    expect_true(is.finite(b) && b > 0, info = info)
    if (dev <= 0.05) {
      expect_gte(b, lo * (1 - 1e-6))
      expect_lte(b, cases$E[i]^(4 / 3) * (1 + 1e-6))
    }
  }
})

test_that("[F-A16] solver failure is detectable: the five known failed solves violate the constraints", {
  py <- eeos_ready()
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  bad <- h %>% filter(grepl(
    "^NS-IBTS (1989 1 NL 64SS .* 13|2008 1 DK 26D4 GOV 135 33|2006 1 NL 64T2 .* 34|1998 1 DE 06NI .* 38|1997 1 NL 64T2 .* 29)$",
    haul_id))
  expect_equal(nrow(bad), 5L)
  dev <- vapply(seq_len(nrow(bad)), function(i) py$constraint_dev(bad$S[i], bad$N[i], bad$E[i]), 0)
  expect_true(all(dev > 0.05)) # the solver really did fail on these
})

test_that("[F-A16] the pipeline never keeps a prediction whose solve failed (|constraint dev| <= 5%)", {
  py <- eeos_ready()
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  bad <- h %>% filter(grepl("^NS-IBTS (1989 1 NL 64SS .* 13|1997 1 NL 64T2 .* 29)$", haul_id))
  dev <- vapply(seq_len(nrow(bad)), function(i) py$constraint_dev(bad$S[i], bad$N[i], bad$E[i]), 0)
  expect_true(all(dev <= 0.05))
})

test_that("[F-A16] biomass.py's own check flags negative constraint deviations (|dev|, not dev)", {
  p <- skip_if_no_file("equation_of_state", "biomass.py")
  src <- readLines(p)
  expect_false(any(grepl("np.any(con > thresh)", src, fixed = TRUE)))
})

test_that("[F-A16] every pipeline prediction respects the solver invariants", {
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  lo <- h$N * (h$E / h$N)^(4 / 3)
  hi <- h$E^(4 / 3)
  bad <- h$haul_id[h$B_pred_norm < lo * (1 - 1e-9) | h$B_pred_norm > hi * (1 + 1e-9)]
  expect_length(bad, 0L)
})

test_that("SLOW: every H1 haul's solve satisfies the constraints within 5%", {
  skip_if_not_slow()
  py <- eeos_ready()
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  dev <- vapply(seq_len(nrow(h)), function(i) py$constraint_dev(h$S[i], h$N[i], h$E[i]), 0)
  expect_equal(sum(dev > 0.05), 0L)
})
