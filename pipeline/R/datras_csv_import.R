# Import ICES unaggregated HL CSV export -> icesDatras-compatible HL data frame.
# Source R/datras_constants.R before this file.
#
# The ICES "Unaggregated trawl and biological information" CSV holds several
# record types (HH haul records, then HL length records), each block with its
# OWN header line. A single read.csv()/fread() would parse HL rows under the HH
# header (or stop at the second header), so read_ices_csv_blocks() splits the
# file on header lines and parses each block separately.
#
# In current exports the `Survey` column holds a code (e.g. "G1022"); the survey
# name ("NS-IBTS") is in `DATRASSurvey`. Haul keys use the survey name so they
# match outputs/datras_hl_raw.rds and FishGlob join keys.

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

#' Split an ICES unaggregated CSV into one data frame per record type.
#'
#' Every line starting with `RecordHeader` opens a new block; each block is parsed
#' with its own header. Returns a named list (e.g. `$HH`, `$HL`). Errors if a
#' block mixes record types or a record type appears in more than one block.
#' Code-like columns are read as character so codes such as "06NI" or "000056"
#' are not coerced.
read_ices_csv_blocks <- function(csv_path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required to read the ICES CSV: renv::install('data.table')")
  }
  if (!file.exists(csv_path)) {
    stop("ICES CSV not found: ", csv_path)
  }

  cat("Reading ICES CSV (all record blocks):", csv_path, "\n")
  lines <- readLines(csv_path, encoding = "UTF-8", warn = FALSE)
  lines[1L] <- sub("^﻿", "", lines[1L])

  hdr <- which(grepl('^"?RecordHeader"?,', lines))
  if (length(hdr) == 0L) {
    stop("No 'RecordHeader' header line found in ", csv_path)
  }
  ends <- c(hdr[-1L] - 1L, length(lines))

  char_cols <- c(
    "RecordHeader", "Survey", "DATRASSurvey", "Country", "Platform", "Gear",
    "StationName", "DataType", "HaulValidity", "StatisticalRectangle",
    "LengthCode", "SpeciesSex", "SpeciesCodeType"
  )

  blocks <- list()
  for (i in seq_along(hdr)) {
    if (ends[[i]] <= hdr[[i]]) {
      next
    }
    header_names <- gsub('"', "", strsplit(lines[[hdr[[i]]]], ",", fixed = TRUE)[[1L]])
    cc <- intersect(char_cols, header_names)
    block <- data.table::fread(
      text = lines[hdr[[i]]:ends[[i]]],
      colClasses = list(character = cc),
      na.strings = c("", "NA"),
      showProgress = FALSE
    )
    types <- unique(block$RecordHeader)
    if (length(types) != 1L) {
      stop(
        "Block starting at line ", hdr[[i]], " mixes record types: ",
        paste(types, collapse = ", ")
      )
    }
    if (types %in% names(blocks)) {
      stop("Record type ", types, " appears in more than one block of ", csv_path)
    }
    blocks[[types]] <- as.data.frame(block)
    cat(sprintf("  %s block: %s rows (lines %d-%d)\n", types, nrow(block), hdr[[i]], ends[[i]]))
  }
  blocks
}

#' Survey name column: DATRASSurvey in current exports, Survey in older ones.
ices_survey_name <- function(df) {
  if ("DATRASSurvey" %in% names(df)) {
    as.character(df$DATRASSurvey)
  } else {
    as.character(df$Survey)
  }
}

#' Read and filter ICES HL rows; return icesDatras-style columns for clean_hl_raw().
#'
#' ICES export uses Platform/HaulNumber/LengthCode/NumberAtLength; pipeline helpers
#' expect Ship/HaulNo/LngtCode/HLNoAtLngt. Counts are NOT effort-aligned here:
#' DataType C rows are per hour. Align with align_hl_to_hh() before building state.
#'
#' @param blocks Optional output of read_ices_csv_blocks() (avoids re-reading).
import_hl_from_ices_csv <- function(
    csv_path,
    year_min = ANALYSIS_YEAR_MIN,
    year_max = ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER,
    survey = "NS-IBTS",
    blocks = NULL) {
  if (is.null(blocks)) {
    blocks <- read_ices_csv_blocks(csv_path)
  }
  if (!"HL" %in% names(blocks)) {
    stop("No HL block found in ICES CSV: ", csv_path)
  }
  raw <- blocks$HL
  cat("  Raw HL rows:", nrow(raw), "\n")

  required <- c(
    "RecordHeader", "Quarter", "Country", "Platform", "HaulNumber",
    "Year", "LengthCode", "LengthClass", "NumberAtLength", "SubsamplingFactor",
    "AphiaID"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("HL block missing columns: ", paste(missing, collapse = ", "))
  }
  raw$survey_name <- ices_survey_name(raw)
  if (!any(raw$survey_name == survey, na.rm = TRUE)) {
    stop(
      "No HL rows with survey '", survey, "' (values found: ",
      paste(head(unique(raw$survey_name), 5L), collapse = ", "), ")"
    )
  }

  hl <- raw %>%
    filter(
      survey_name == survey,
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
      Survey = survey_name,
      Year = as.integer(Year),
      Quarter = as.integer(Quarter),
      Country = as.character(Country),
      Ship = as.character(Platform),
      HaulNo = as.integer(HaulNumber),
      Valid_Aphia = as.integer(AphiaID),
      LngtCode = as.character(LengthCode),
      LngtClass = as.numeric(LengthClass),
      HLNoAtLngt = as.numeric(NumberAtLength),
      SubFactor = as.numeric(SubsamplingFactor)
    )

  # Missing / -9 SubsamplingFactor means no subsampling.
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
