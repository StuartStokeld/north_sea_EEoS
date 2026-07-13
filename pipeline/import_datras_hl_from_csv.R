# Standalone: import ICES unaggregated HL CSV -> outputs/datras_hl_raw.rds
# Normally called automatically by build_datras_state_variables.R.

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

cli_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cli_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/import_datras_hl_from_csv.R")
}

r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_hl_helpers.R"))
source(file.path(r_dir, "datras_csv_import.R"))

path_hl_raw <- file.path(project_root, "outputs", "datras_hl_raw.rds")
ensure_datras_hl_raw(project_root, path_hl_raw, force = force)

hl <- readRDS(path_hl_raw)
coverage <- assess_hl_coverage(hl)
cat("\nCoverage:", coverage$message, "\n")
