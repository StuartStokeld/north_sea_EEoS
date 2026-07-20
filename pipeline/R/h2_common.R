# Shared constants and paths for Hypothesis 2 (StatRec × fishing pressure)

H2_YEAR_MIN <- 1985L
H2_YEAR_MAX <- 2015L
H2_MIN_HAULS_DEFAULT <- 10L
H2_MIN_HAULS_SENSITIVITY <- c(5L, 10L, 20L)

#' Normalise ICES statistical rectangle codes to FishGlob format (e.g. "38F5").
normalize_stat_rec <- function(x) {
  x <- as.character(x)
  x <- gsub('"', "", x, fixed = TRUE)
  trimws(x)
}

h2_couce_raw_path <- function(project_root) {
  file.path(
    project_root,
    "data",
    "external",
    "couce_trawling_effort",
    "NorthSea_trawling_effort_1985to2015_REVIEW_v2.csv"
  )
}

h2_ices_shapefile_path <- function(project_root) {
  file.path(
    project_root,
    "gis",
    "ICES_rectangles",
    "ICES_Statistical_Rectangles_Eco.shp"
  )
}

h2_output_paths <- function(project_root) {
  out <- file.path(project_root, "outputs")
  fig <- file.path(out, "figures")
  list(
    couce_year = file.path(out, "h2_couce_year_effort.rds"),
    couce_rectangle = file.path(out, "h2_couce_rectangle_effort.rds"),
    couce_diagnostics = file.path(out, "h2_couce_import_diagnostics.csv"),
    panel = file.path(out, "h2_rectangle_panel.rds"),
    panel_csv = file.path(out, "h2_rectangle_panel.csv"),
    ols_results = file.path(out, "h2_ols_results.csv"),
    spatial_diagnostics = file.path(out, "h2_spatial_diagnostics.csv"),
    sem_results = file.path(out, "h2_sem_results.csv"),
    model_summary = file.path(out, "h2_model_summary.rds"),
    sem_biomass_covariate = file.path(out, "h2_sem_biomass_covariate.rds"),
    spatial_sensitivity = file.path(out, "h2_spatial_sensitivity_min5_min20.rds"),
    leverage_diagnostics = file.path(out, "h2_leverage_diagnostics.rds"),
    lm_tests = file.path(out, "h2_lm_tests.rds"),
    provenance_stamp = file.path(out, "h2_provenance_stamp.rds"),
    sar_lag_results = file.path(out, "h2_sar_lag_results.rds"),
    fig_dir = fig
  )
}
