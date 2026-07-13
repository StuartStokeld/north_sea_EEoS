# Build rectangle-level H2 analysis panel (EEoS residuals + fishing pressure)
# Run after import_couce_fishing_effort.R and build_eeos_predictions.R

suppressPackageStartupMessages({
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/build_h2_rectangle_panel.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_panel_helpers.R"))

paths <- h2_output_paths(project_root)
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_couce <- paths$couce_rectangle

stopifnot(
  file.exists(path_haul),
  file.exists(path_couce)
)

dir.create(dirname(paths$panel), recursive = TRUE, showWarnings = FALSE)

haul <- readRDS(path_haul)
couce_rectangle <- readRDS(path_couce)

residual_panel <- build_h2_rectangle_residuals(haul)
panel <- build_h2_analysis_panel(
  residual_panel,
  couce_rectangle,
  min_hauls = H2_MIN_HAULS_DEFAULT,
  require_fishing = TRUE
)

saveRDS(panel, paths$panel)
write_csv(panel, paths$panel_csv)

cat("Rectangles with residuals:", nrow(residual_panel), "\n")
cat("Analysis panel (min_hauls =", H2_MIN_HAULS_DEFAULT, "):", nrow(panel), "\n")
cat("Saved", paths$panel, "\n")
cat("Saved", paths$panel_csv, "\n")
