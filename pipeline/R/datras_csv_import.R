# Import ICES unaggregated HL CSV export -> icesDatras-compatible HL data frame.
# Source R/datras_constants.R before this file.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Resolve path to ICES unaggregated HL CSV in project data folder.
find_ices_hl_csv <- function(project_root) {
  pattern_dir <- "^Unaggregated trawl and biological information_"
  candidates <- list.dirs(project_root, full.names = TRUE, recursive = FALSE)
  candidates <- candidates[grepl(pattern_dir, basename(candidates), ignore.case = TRUE)]

  if (length(candidates) == 0L) {
    return(NULL)
  }

  dir_path <- candidates[[1L]]
  csv_name <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  csv_name <- csv_name[!grepl("^Disclaimer", basename(csv_name), ignore.case = TRUE)]

  if (length(csv_name) == 0L) {
    return(NULL)
  }

  csv_name[[1L]]
}

#' Read and filter ICES HL CSV; return icesDatras-style columns for clean_hl_raw().
#'
#' ICES export uses Platform/HaulNumber/LengthCode/NumberAtLength; pipeline helpers
#' expect Ship/HaulNo/LngtCode/HLNoAtLngt. NumberAtLength in the export is the
#' raised count at length (SubsamplingFactor is usually 1).
import_hl_from_ices_csv <- function(
    csv_path,
    year_min = ANALYSIS_YEAR_MIN,
    year_max = ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER,
    survey = "NS-IBTS") {
  if (!file.exists(csv_path)) {
    stop("ICES HL CSV not found: ", csv_path)
  }

  cat("Reading ICES HL CSV:", csv_path, "\n")
  raw <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("  Raw rows:", nrow(raw), "\n")

  required <- c(
    "RecordHeader", "Survey", "Quarter", "Country", "Platform", "HaulNumber",
    "Year", "LengthCode", "LengthClass", "NumberAtLength", "SubsamplingFactor",
    "AphiaID"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("CSV missing columns: ", paste(missing, collapse = ", "))
  }

  hl <- raw %>%
    filter(
      RecordHeader == "HL",
      Survey == survey,
      Quarter == quarter,
      Year >= year_min,
      Year <= year_max,
      LengthCode != "-9",
      !is.na(AphiaID),
      !is.na(NumberAtLength),
      NumberAtLength > 0,
      !is.na(LengthClass)
    ) %>%
    transmute(
      Survey = Survey,
      Year = as.integer(Year),
      Quarter = as.integer(Quarter),
      Country = Country,
      Ship = Platform,
      HaulNo = as.integer(HaulNumber),
      Valid_Aphia = as.integer(AphiaID),
      LngtCode = as.character(LengthCode),
      LngtClass = as.numeric(LengthClass),
      HLNoAtLngt = as.numeric(NumberAtLength),
      SubFactor = as.numeric(SubsamplingFactor)
    )

  # ICES export NumberAtLength is already raised; SubFactor is usually 1.
  hl$SubFactor[is.na(hl$SubFactor) | hl$SubFactor <= 0] <- 1

  cat("  Filtered HL rows (Q", quarter, ", ", year_min, "-", year_max, "): ",
      nrow(hl), "\n", sep = "")
  cat("  Unique hauls:", n_distinct(
    paste(hl$Survey, hl$Year, hl$Quarter, hl$Country, hl$Ship, hl$HaulNo, sep = "_")
  ), "\n")

  hl
}

#' Import CSV and save outputs/datras_hl_raw.rds if needed.
ensure_datras_hl_raw <- function(
    project_root,
    path_hl_raw,
    force = FALSE) {
  csv_path <- find_ices_hl_csv(project_root)

  needs_import <- force || !file.exists(path_hl_raw)
  if (!needs_import && file.exists(path_hl_raw)) {
    hl_existing <- readRDS(path_hl_raw)
    coverage <- assess_hl_coverage(hl_existing)
    needs_import <- !coverage$complete
    if (!needs_import) {
      cat("datras_hl_raw.rds already complete — skipping CSV import.\n")
      return(invisible(hl_existing))
    }
    cat("datras_hl_raw.rds incomplete (", coverage$message, ") — importing from CSV.\n", sep = "")
  }

  if (is.null(csv_path)) {
    if (file.exists(path_hl_raw)) {
      return(invisible(readRDS(path_hl_raw)))
    }
    stop(
      "No ICES HL CSV folder found under project root and no datras_hl_raw.rds present.\n",
      "Expected: Unaggregated trawl and biological information_*/ *.csv"
    )
  }

  hl <- import_hl_from_ices_csv(csv_path)
  dir.create(dirname(path_hl_raw), recursive = TRUE, showWarnings = FALSE)
  saveRDS(hl, path_hl_raw)
  cat("Saved:", path_hl_raw, "\n")
  invisible(hl)
}
