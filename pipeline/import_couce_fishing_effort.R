# Import Couce et al. (2020) North Sea trawling effort for H2
# Run from workspace root: Rscript pipeline/import_couce_fishing_effort.R

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
  stop("Run from pipeline/ or Rscript pipeline/import_couce_fishing_effort.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_couce_helpers.R"))

paths <- h2_output_paths(project_root)
raw_path <- h2_couce_raw_path(project_root)

stopifnot(file.exists(raw_path))

path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
stopifnot(file.exists(path_haul))

dir.create(dirname(paths$couce_year), recursive = TRUE, showWarnings = FALSE)

year_effort <- read_couce_year_effort(raw_path)
rectangle_effort <- aggregate_couce_rectangle(year_effort)

haul <- readRDS(path_haul)
diag <- couce_import_diagnostics(
  year_effort,
  rectangle_effort,
  haul$stat_rec
)

saveRDS(year_effort, paths$couce_year)
saveRDS(rectangle_effort, paths$couce_rectangle)
write_csv(diag, paths$couce_diagnostics)

cat("Couce year-level rows:", nrow(year_effort), "\n")
cat("Couce rectangles:", nrow(rectangle_effort), "\n")
cat("Saved", paths$couce_year, "\n")
cat("Saved", paths$couce_rectangle, "\n")
cat("Saved", paths$couce_diagnostics, "\n")
