# H2 robustness checks and provenance (diagnostic layer — not core pipeline)
# Run after build_h2_rectangle_panel.R and run_h2_models.R
#
# Does NOT overwrite h2_ols_results.csv, h2_sem_results.csv, h2_model_summary.rds

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2_robustness_checks.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_panel_helpers.R"))
source(file.path(script_dir, "R", "h2_model_helpers.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h2_provenance_helpers.R"))
source(file.path(script_dir, "R", "h2_robustness_helpers.R"))

paths <- h2_output_paths(project_root)
path_panel <- paths$panel
path_couce <- paths$couce_rectangle
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")

stopifnot(
  file.exists(path_panel),
  file.exists(path_couce),
  file.exists(path_haul)
)

provenance <- h2_collect_provenance(project_root)
dir.create(dirname(paths$provenance_stamp), recursive = TRUE, showWarnings = FALSE)

panel10 <- readRDS(path_panel)
couce_rectangle <- readRDS(path_couce)
haul <- readRDS(path_haul)
residual_panel <- build_h2_rectangle_residuals(haul)

panel5 <- build_h2_analysis_panel(
  residual_panel, couce_rectangle, min_hauls = 5L, require_fishing = TRUE
)
panel20 <- build_h2_analysis_panel(
  residual_panel, couce_rectangle, min_hauls = 20L, require_fishing = TRUE
)

rectangles_sf <- load_ices_rectangles_sf(project_root)
weights10 <- build_h2_spatial_weights(panel10, rectangles_sf, project_root)

# ---------------------------------------------------------------------------
# Check 1 — SEM with biomass covariate (min10 weights)
# ---------------------------------------------------------------------------
cat("=== Check 1: SEM with biomass covariate ===\n")

sem_biomass <- fit_h2_sem(
  panel10,
  weights10$listw,
  formula = SEM_BIOMASS_FORMULA,
  model_id = "sem_biomass_covariate"
) %>%
  mutate(panel = "min10")

sem_biomass_out <- h2_stamp_result(
  list(
    sem_coefficients = sem_biomass,
    primary_sem_lambda_reference = 0.8584177740031406
  ),
  provenance
)
saveRDS(sem_biomass_out, paths$sem_biomass_covariate)
write_csv(sem_biomass, sub("\\.rds$", ".csv", paths$sem_biomass_covariate))

# ---------------------------------------------------------------------------
# Check 2 — Spatial models on min5 / min10 / min20 panels
# ---------------------------------------------------------------------------
cat("=== Check 2: Spatial sensitivity across haul thresholds ===\n")

run_min10 <- h2_spatial_panel_summary(
  panel10, rectangles_sf, project_root, "min10"
)
run_min5 <- h2_spatial_panel_summary(
  panel5, rectangles_sf, project_root, "min5"
)
run_min20 <- h2_spatial_panel_summary(
  panel20, rectangles_sf, project_root, "min20"
)

spatial_sensitivity_table <- h2_spatial_sensitivity_table(list(run_min5, run_min10, run_min20))

spatial_sensitivity_out <- h2_stamp_result(
  list(
    summary_table = spatial_sensitivity_table,
    sem_by_panel = bind_rows(
      run_min5$sem_results,
      run_min10$sem_results,
      run_min20$sem_results
    ),
    spatial_diagnostics_by_panel = bind_rows(
      run_min5$spatial_diagnostics,
      run_min10$spatial_diagnostics,
      run_min20$spatial_diagnostics
    )
  ),
  provenance
)
saveRDS(spatial_sensitivity_out, paths$spatial_sensitivity)
write_csv(spatial_sensitivity_table, sub("\\.rds$", ".csv", paths$spatial_sensitivity))

# ---------------------------------------------------------------------------
# Check 3 — Leverage / Cook's distance on primary OLS (min10)
# ---------------------------------------------------------------------------
cat("=== Check 3: Leverage diagnostics ===\n")

leverage <- h2_leverage_diagnostics(panel10)
leverage_out <- h2_stamp_result(
  list(
    cook_threshold = leverage$cook_threshold,
    n_flagged = leverage$n_flagged,
    flagged_rectangles = leverage$flagged_rectangles,
    all_rectangles = leverage$diagnostics,
    full_panel_ols = leverage$full_panel,
    trimmed_panel_ols = leverage$trimmed_panel
  ),
  provenance
)
saveRDS(leverage_out, paths$leverage_diagnostics)
write_csv(leverage$diagnostics, sub("\\.rds$", ".csv", paths$leverage_diagnostics))

# ---------------------------------------------------------------------------
# Check 4 — Lagrange Multiplier tests (min10 OLS + weights)
# ---------------------------------------------------------------------------
cat("=== Check 4: LM spatial specification tests ===\n")

primary_ols_min10 <- lm(PRIMARY_OLS_FORMULA, data = panel10)
lm_tests <- tryCatch(
  h2_lm_spatial_tests(primary_ols_min10, weights10$listw),
  error = function(e) {
    data.frame(
      test = "error",
      statistic = NA_real_,
      p_value = NA_real_,
      df = NA_real_,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  }
)

lm_out <- h2_stamp_result(
  list(lm_tests = lm_tests),
  provenance
)
saveRDS(lm_out, paths$lm_tests)
write_csv(lm_tests, sub("\\.rds$", ".csv", paths$lm_tests))

# ---------------------------------------------------------------------------
# Provenance stamp (standalone)
# ---------------------------------------------------------------------------
saveRDS(provenance, paths$provenance_stamp)
prov_df <- as.data.frame(provenance, stringsAsFactors = FALSE)
write_csv(prov_df, sub("\\.rds$", ".csv", paths$provenance_stamp))

cat("\n=== Outputs saved ===\n")
cat(paths$sem_biomass_covariate, "\n")
cat(paths$spatial_sensitivity, "\n")
cat(paths$leverage_diagnostics, "\n")
cat(paths$lm_tests, "\n")
cat(paths$provenance_stamp, "\n")
