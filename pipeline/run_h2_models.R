# H2 statistical models: OLS, spatial autocorrelation, spatial error model
# Run after build_h2_rectangle_panel.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2_models.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_panel_helpers.R"))
source(file.path(script_dir, "R", "h2_model_helpers.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))

paths <- h2_output_paths(project_root)
path_panel <- paths$panel
path_couce <- paths$couce_rectangle
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")

stopifnot(
  file.exists(path_panel),
  file.exists(path_couce),
  file.exists(path_haul)
)

dir.create(paths$fig_dir, recursive = TRUE, showWarnings = FALSE)

panel <- readRDS(path_panel)
couce_rectangle <- readRDS(path_couce)
haul <- readRDS(path_haul)
residual_panel <- build_h2_rectangle_residuals(haul)

# Primary OLS + sensitivity models across min-haul thresholds
ols_by_threshold <- lapply(H2_MIN_HAULS_SENSITIVITY, function(min_h) {
  p <- if (min_h == H2_MIN_HAULS_DEFAULT) {
    panel
  } else {
    build_h2_analysis_panel(
      residual_panel,
      couce_rectangle,
      min_hauls = min_h,
      require_fishing = TRUE
    )
  }
  fit_h2_model_set(p, min_hauls_label = paste0("min", min_h))
})
ols_results <- bind_rows(ols_by_threshold)

primary_fit <- lm(
  mean_abs_residual ~ mean_annual_hours_total,
  data = panel
)

# Spatial analysis
rectangles_sf <- load_ices_rectangles_sf(project_root)
weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)

spatial_diagnostics <- h2_global_spatial_tests(
  residuals(primary_fit),
  weights$listw
)
spatial_diagnostics$variable <- "ols_primary_abs_residuals"

moran_dv <- h2_global_spatial_tests(panel$mean_abs_residual, weights$listw)
moran_dv$variable <- "mean_abs_residual"
moran_iv <- h2_global_spatial_tests(
  log(panel$mean_annual_hours_total + 1),
  weights$listw
)
moran_iv$variable <- "log_mean_annual_hours_total"

spatial_diagnostics <- bind_rows(spatial_diagnostics, moran_dv, moran_iv)
spatial_diagnostics$n_rectangles <- nrow(panel)
spatial_diagnostics$n_isolated <- weights$n_isolated

sem_results <- fit_h2_sem(panel, weights$listw)

write_csv(ols_results, paths$ols_results)
write_csv(spatial_diagnostics, paths$spatial_diagnostics)
write_csv(sem_results, paths$sem_results)

model_summary <- list(
  panel_n = nrow(panel),
  primary_ols = ols_results %>%
    filter(model_id == "primary_abs_min10") %>%
    slice(1),
  sem = sem_results,
  spatial = spatial_diagnostics,
  min_hauls = H2_MIN_HAULS_DEFAULT
)
saveRDS(model_summary, paths$model_summary)

save_h2_figures(
  panel,
  weights$panel_sf,
  weights,
  primary_fit,
  paths$fig_dir
)

cat("Analysis rectangles:", nrow(panel), "\n")
cat("Primary OLS beta:", round(coef(primary_fit)[2], 5),
    "p =", signif(summary(primary_fit)$coefficients[2, 4], 3), "\n")
cat("Moran's I p-value:",
    signif(spatial_diagnostics$p_value[spatial_diagnostics$test == "morans_i"][1], 3), "\n")
cat("SEM beta:", round(sem_results$estimate, 5),
    "p =", signif(sem_results$p_value, 3), "\n")
cat("Saved", paths$ols_results, "\n")
cat("Saved", paths$spatial_diagnostics, "\n")
cat("Saved", paths$sem_results, "\n")
