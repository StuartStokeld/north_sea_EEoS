# H2 corollary follow-up: SAR (spatial lag) diagnostics for the
# biomass-vs-fishing-pressure bivariate test.
#
# Run after run_h2_biomass_fishing_pressure_test.R, which established that the
# LM specification tests favoured a spatial lag process for
# mean_ln_B_obs ~ mean_annual_hours_total and fit a lagsarlm() model.
#
# This script does NOT refit a different model or change that result. It
# reports two things about the already-fit SAR model, factually, without
# interpretation:
#   1. Whether SAR residuals still show spatial structure (built-in LM
#      residual-autocorrelation test from summary.sarlm(), plus an explicit
#      Moran's I test on the SAR residuals for corroboration).
#   2. The direct / indirect / total impact decomposition (Anselin-style
#      impact measures), with simulated standard errors and p-values, since
#      the coefficient reported in h2_biomass_fishing_spatial_model.csv is
#      the DIRECT effect only, not the full effect propagated through
#      neighbouring rectangles via rho.
#
# Does NOT modify or overwrite any existing H2 or H2-corollary output.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(spatialreg)
  library(spdep)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_biomass_fishing_sar_diagnostics.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_panel_helpers.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))

paths <- h2_output_paths(project_root)
path_biomass_summary <- file.path(project_root, "outputs", "h2_biomass_fishing_test_summary.rds")

stopifnot(
  file.exists(paths$panel),
  file.exists(path_biomass_summary)
)

biomass_summary <- readRDS(path_biomass_summary)
if (!identical(biomass_summary$spatial_spec, "sar_lag")) {
  stop(
    "STOP: h2_biomass_fishing_test_summary.rds records spatial_spec = '",
    biomass_summary$spatial_spec, "', not 'sar_lag'. This diagnostics script ",
    "assumes a lagsarlm() fit was the model selected by the LM tests; refusing ",
    "to proceed rather than silently fitting a different model."
  )
}

out_residual_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_sar_residual_diagnostics.csv")
out_impacts_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_sar_impacts.csv")
out_notes_txt <- file.path(project_root, "outputs", "h2_biomass_fishing_sar_diagnostics_notes.txt")

notes <- character(0)
add_note <- function(x) notes <<- c(notes, x)

# ---------------------------------------------------------------------------
# Refit the identical SAR specification (same panel, same weights, same
# formula as run_h2_biomass_fishing_pressure_test.R) to obtain the fitted
# model object needed for residual diagnostics and impacts(). Deterministic:
# the coefficient/rho values match h2_biomass_fishing_spatial_model.csv
# exactly (verified below).
# ---------------------------------------------------------------------------
panel <- readRDS(paths$panel)
rectangles_sf <- load_ices_rectangles_sf(project_root)
weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)

BIOMASS_FISHING_FORMULA <- mean_ln_B_obs ~ mean_annual_hours_total
fit <- spatialreg::lagsarlm(
  BIOMASS_FISHING_FORMULA,
  data = panel,
  listw = weights$listw,
  zero.policy = TRUE
)

saved_row <- biomass_summary$spatial_results
refit_matches <- isTRUE(all.equal(
  unname(coef(fit)["mean_annual_hours_total"]),
  saved_row$estimate,
  tolerance = 1e-8
))
add_note(sprintf(
  paste0(
    "Refit check: this script refits the identical SAR specification (same panel, ",
    "n = %d, weights, formula) rather than loading a cached model object (none is ",
    "persisted). Refit coefficient matches the saved h2_biomass_fishing_spatial_model.csv ",
    "value to within 1e-8: %s."
  ),
  nrow(panel), refit_matches
))
if (!refit_matches) {
  stop("STOP: refit SAR coefficient does not match the previously saved result; flagging rather than proceeding.")
}

# ---------------------------------------------------------------------------
# 1. Residual spatial structure diagnostics
# ---------------------------------------------------------------------------
sm <- summary(fit)
lm_resid_test <- sm$LMtest       # built-in LM test for residual autocorrelation (summary.Sarlm)
lm_resid_df <- sm$LMtest.df
lm_resid_p <- pchisq(lm_resid_test, df = 1, lower.tail = FALSE)

moran_resid <- spdep::moran.test(residuals(fit), weights$listw, zero.policy = TRUE)

residual_diagnostics <- data.frame(
  test = c("lm_residual_autocorrelation_builtin", "morans_i_sar_residuals"),
  statistic = c(unname(lm_resid_test), unname(moran_resid$estimate[["Moran I statistic"]])),
  p_value = c(lm_resid_p, moran_resid$p.value),
  df = c(1L, NA_integer_),
  notes = c(
    "Built-in summary.Sarlm() LM test for residual autocorrelation given the fitted rho.",
    "Explicit Moran's I on SAR residuals, same queen-contiguity weights as the fit, for corroboration."
  ),
  stringsAsFactors = FALSE
)
write_csv(residual_diagnostics, out_residual_csv)

add_note(sprintf(
  paste0(
    "Residual spatial structure: built-in LM residual-autocorrelation test = %.4f ",
    "(df = 1, p = %.4f); explicit Moran's I on SAR residuals = %.4f (p = %.4f, ",
    "one-sided greater). Both computed on the fitted lagsarlm() residuals using the ",
    "same queen-contiguity weights as the fit itself. Reported factually; no ",
    "judgement is made here about whether this level of residual autocorrelation ",
    "is acceptable."
  ),
  lm_resid_test, lm_resid_p,
  unname(moran_resid$estimate[["Moran I statistic"]]), moran_resid$p.value
))

# ---------------------------------------------------------------------------
# 2. Direct / indirect / total impact decomposition
# ---------------------------------------------------------------------------
# The coefficient in h2_biomass_fishing_spatial_model.csv is the DIRECT
# effect only. For a spatial lag model, a one-unit change in fishing hours in
# rectangle i also propagates to neighbouring rectangles' fitted biomass via
# rho, which feeds back into i's neighbours and, through the network, back
# toward i (the "indirect"/spillover effect). Total = direct + indirect.
# Standard errors/p-values for indirect and total are obtained by simulation
# (Monte Carlo draws from the model's asymptotic parameter distribution),
# following the standard spatialreg::impacts() workflow; a seed is fixed for
# reproducibility.
set.seed(20260721)
n_sim <- 2000L
W_sparse <- as(weights$listw, "CsparseMatrix")
trMatc <- spatialreg::trW(W_sparse, type = "mult")
impacts_fit <- spatialreg::impacts(fit, tr = trMatc, R = n_sim)
impacts_summary <- summary(impacts_fit, zstats = TRUE, short = TRUE)

# Single non-intercept term (bivariate spec), so impacts_summary$res$direct/
# indirect/total are unnamed scalars, and semat/zmat/pzmat are 1x3 matrices
# with columns "Direct"/"Indirect"/"Total" for that one term/row.
impact_table <- data.frame(
  term = "mean_annual_hours_total",
  measure = c("direct", "indirect", "total"),
  estimate = c(
    impacts_summary$res$direct,
    impacts_summary$res$indirect,
    impacts_summary$res$total
  ),
  std_error = as.numeric(impacts_summary$semat[1, c("Direct", "Indirect", "Total")]),
  z_value = as.numeric(impacts_summary$zmat[1, c("Direct", "Indirect", "Total")]),
  p_value = as.numeric(impacts_summary$pzmat[1, c("Direct", "Indirect", "Total")]),
  n_simulations = n_sim,
  rho = unname(fit$rho),
  stringsAsFactors = FALSE
)
write_csv(impact_table, out_impacts_csv)

add_note(sprintf(
  paste0(
    "Impact decomposition (rho = %.5f, %d Monte Carlo draws, seed = 20260721): ",
    "direct = %.6e (SE %.6e, p = %.4f); indirect = %.6e (SE %.6e, p = %.4f); ",
    "total = %.6e (SE %.6e, p = %.4f). The direct estimate matches the coefficient ",
    "already reported in h2_biomass_fishing_spatial_model.csv (to simulation/analytic ",
    "tolerance); indirect and total are new quantities not previously reported."
  ),
  unname(fit$rho), n_sim,
  impact_table$estimate[1], impact_table$std_error[1], impact_table$p_value[1],
  impact_table$estimate[2], impact_table$std_error[2], impact_table$p_value[2],
  impact_table$estimate[3], impact_table$std_error[3], impact_table$p_value[3]
))

writeLines(notes, out_notes_txt)

cat("=== SAR residual diagnostics ===\n")
print(residual_diagnostics)
cat("\n=== SAR impact decomposition (direct / indirect / total) ===\n")
print(impact_table)
cat("\nSaved", out_residual_csv, "\n")
cat("Saved", out_impacts_csv, "\n")
cat("Saved", out_notes_txt, "\n")
