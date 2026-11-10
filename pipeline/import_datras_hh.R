# Import DATRAS HH (haul) records -> outputs/datras_hh_raw.rds
#
# Source, in order:
#   1. HH block of the ICES "Unaggregated trawl and biological information_*" CSV
#      (gitignored, ~200 MB; same download as the HL records)
#   2. icesDatras::getDATRAS("HH", "NS-IBTS", 1985:2015, 1)
#
# Run from repo root:
#   Rscript pipeline/import_datras_hh.R            # skip if RDS exists
#   Rscript pipeline/import_datras_hh.R --force    # re-import
#
# Fields kept (standardised names): Survey, Year, Quarter, Country, Ship, HaulNo,
# Gear, HaulDur, HaulVal, DataType, StatRec, Month, haul_key, hh_key_n, hh_source.

suppressPackageStartupMessages({
  library(dplyr)
})

force <- "--force" %in% commandArgs(trailingOnly = TRUE)

cli_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cli_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/import_datras_hh.R")
}

r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_csv_import.R"))
source(file.path(r_dir, "datras_hh_helpers.R"))

path_hh_raw <- file.path(project_root, "outputs", "datras_hh_raw.rds")
hh <- ensure_datras_hh_raw(project_root, path_hh_raw, force = force)

cat("\n=== HH summary ===\n")
cat("Source:", unique(hh$hh_source), "\n")
cat("Records:", nrow(hh), " | years:", min(hh$Year), "-", max(hh$Year), "\n")
cat("Haul keys shared by >1 record (ambiguous, dropped at alignment):",
    n_distinct(hh$haul_key[hh$hh_key_n > 1L]), "\n\n")
print(as.data.frame(count(hh, HaulVal, DataType)))
cat("\nHaulDur (min):\n")
print(summary(hh$HaulDur))
