# H2 spatial lag (SAR) models — compare with existing SEM fits
# Run after build_h2_rectangle_panel.R (uses min10 panel + same weights as primary SEM)
#
# Does NOT overwrite h2_sem_results.csv, h2_sem_biomass_covariate.*, or OLS outputs

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_sar_lag_models.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_robustness_helpers.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h2_provenance_helpers.R"))

paths <- h2_output_paths(project_root)
path_panel <- paths$panel

stopifnot(file.exists(path_panel))

provenance <- h2_collect_provenance(project_root)
panel10 <- readRDS(path_panel)
rectangles_sf <- load_ices_rectangles_sf(project_root)

# Single weights build on min10 panel (identical spec to primary SEM in run_h2_models.R)
weights10 <- build_h2_spatial_weights(panel10, rectangles_sf, project_root)
listw10 <- weights10$listw

cat("Panel min10: N =", nrow(panel10), " n_isolated =", weights10$n_isolated, "\n")
cat("Weights: queen contiguity, style W, zero.policy = TRUE\n\n")

# ---------------------------------------------------------------------------
# SAR lag models
# ---------------------------------------------------------------------------
sar_primary <- fit_h2_sar_lag(
  panel10,
  listw10,
  formula = SEM_PRIMARY_FORMULA,
  model_id = "sar_lag_primary_abs"
)

sar_biomass <- fit_h2_sar_lag(
  panel10,
  listw10,
  formula = SEM_BIOMASS_FORMULA,
  model_id = "sar_lag_biomass_covariate"
)

sar_coefficients <- bind_rows(sar_primary, sar_biomass)

# ---------------------------------------------------------------------------
# SEM fit stats for AIC / log-likelihood comparison (refit, same panel/weights)
# ---------------------------------------------------------------------------
sem_primary_stats <- fit_h2_sem_fit_stats(
  panel10,
  listw10,
  formula = SEM_PRIMARY_FORMULA,
  model_id = "sem_primary_abs"
)

sem_biomass_stats <- fit_h2_sem_fit_stats(
  panel10,
  listw10,
  formula = SEM_BIOMASS_FORMULA,
  model_id = "sem_biomass_covariate"
)

sar_primary_stats <- sar_primary %>%
  distinct(model_id, log_likelihood, aic, rho, rho_p_value, n) %>%
  mutate(spatial_spec = "sar_lag", spatial_param = rho, spatial_param_name = "rho")

sar_biomass_stats <- sar_biomass %>%
  distinct(model_id, log_likelihood, aic, rho, rho_p_value, n) %>%
  mutate(spatial_spec = "sar_lag", spatial_param = rho, spatial_param_name = "rho")

model_comparison <- bind_rows(
  sem_primary_stats %>%
    transmute(
      comparison_id = "primary_hours_only",
      spatial_spec,
      model_id,
      log_likelihood,
      aic,
      spatial_param,
      spatial_param_name,
      n
    ),
  sar_primary_stats %>%
    transmute(
      comparison_id = "primary_hours_only",
      spatial_spec,
      model_id,
      log_likelihood,
      aic,
      spatial_param,
      spatial_param_name,
      n
    ),
  sem_biomass_stats %>%
    transmute(
      comparison_id = "biomass_covariate",
      spatial_spec,
      model_id,
      log_likelihood,
      aic,
      spatial_param,
      spatial_param_name,
      n
    ),
  sar_biomass_stats %>%
    transmute(
      comparison_id = "biomass_covariate",
      spatial_spec,
      model_id,
      log_likelihood,
      aic,
      spatial_param,
      spatial_param_name,
      n
    )
) %>%
  group_by(comparison_id) %>%
  mutate(
    sem_aic = aic[spatial_spec == "sem"][1L],
    sem_loglik = log_likelihood[spatial_spec == "sem"][1L],
    aic_delta_sar_minus_sem = ifelse(spatial_spec == "sar_lag", aic - sem_aic, NA_real_),
    loglik_delta_sar_minus_sem = ifelse(
      spatial_spec == "sar_lag",
      log_likelihood - sem_loglik,
      NA_real_
    )
  ) %>%
  ungroup()

weights_info <- data.frame(
  panel = "min10",
  n_rectangles = nrow(panel10),
  n_isolated = weights10$n_isolated,
  zero_policy = TRUE,
  nb_style = "W",
  nb_type = "queen",
  stringsAsFactors = FALSE
)

out <- h2_stamp_result(
  list(
    sar_coefficients = sar_coefficients,
    model_comparison = model_comparison,
    weights_info = weights_info
  ),
  provenance
)

dir.create(dirname(paths$sar_lag_results), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, paths$sar_lag_results)
write_csv(sar_coefficients, sub("\\.rds$", ".csv", paths$sar_lag_results))
write_csv(model_comparison, sub("\\.rds$", "_comparison.csv", paths$sar_lag_results))

cat("=== SAR primary (hours only) ===\n")
print(sar_primary %>% select(term, estimate, p_value, rho, rho_p_value, aic))

cat("\n=== SAR biomass covariate ===\n")
print(sar_biomass %>% select(term, estimate, p_value, rho, rho_p_value, aic))

cat("\n=== Model comparison (AIC) ===\n")
print(model_comparison)

cat("\nSaved", paths$sar_lag_results, "\n")
cat("Saved", sub("\\.rds$", ".csv", paths$sar_lag_results), "\n")
