# DATRAS HH (haul) records: import from the ICES CSV (fallback icesDatras) and standardise.
# Source R/datras_constants.R and R/datras_csv_import.R before this file.
#
# Why: DATRAS HL counts for DataType C are catch per HOUR (DATRAS FAQ; SubFactor = 1),
# while R/P/S counts are per haul (raise with SubFactor). FishGlob B_obs is the raw
# haul weight (FishGlob get_datras.R converts C back with HaulDur / 60). S, N, E must
# be built from raw per-haul counts of the same tow, so HL is aligned to HH first:
#
#   n_raw = HLNoAtLngt * SubFactor * f,   f = HaulDur / 60 (C), 1 (R, P, S)
#
# Hauls kept: HaulVal == "V", DataType in HH_DATATYPES_KEEP, and
# dur_min <= HaulDur <= dur_max (default exactly 30 min = primary; 26-34 min is the
# sensitivity run). No duration rescaling.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Build the 6-field haul key used by HL, HH and the FishGlob join.
make_haul_key <- function(Survey, Year, Quarter, Country, Ship, HaulNo) {
  paste(Survey, Year, Quarter, Country, Ship, HaulNo, sep = "_")
}

#' Standardise HH records from either the ICES CSV (new header) or icesDatras.
#'
#' New-header names are mapped Platform -> Ship, HaulNumber -> HaulNo,
#' HaulDuration -> HaulDur, HaulValidity -> HaulVal; old icesDatras names pass
#' through. Survey is the survey NAME (DATRASSurvey in current CSV exports).
#' HaulDur <= 0 (DATRAS -9 = missing) becomes NA. `hh_key_n` counts records per
#' haul_key: > 1 means the 6-field key is ambiguous (e.g. two gears sharing a
#' haul number) and cannot be aligned to HL, which carries no Gear.
standardise_hh <- function(
    hh_raw,
    year_min = ANALYSIS_YEAR_MIN,
    year_max = ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER,
    survey = "NS-IBTS") {
  pick <- function(options, required = TRUE) {
    hit <- intersect(options, names(hh_raw))
    if (length(hit) == 0L) {
      if (required) {
        stop("HH data missing column (any of): ", paste(options, collapse = ", "))
      }
      return(rep(NA, nrow(hh_raw)))
    }
    hh_raw[[hit[[1L]]]]
  }
  if ("RecordHeader" %in% names(hh_raw)) {
    hh_raw <- hh_raw[hh_raw$RecordHeader == "HH", , drop = FALSE]
  }

  hh <- tibble(
    Survey = ices_survey_name(hh_raw),
    Year = as.integer(pick("Year")),
    Quarter = as.integer(pick("Quarter")),
    Country = trimws(as.character(pick("Country"))),
    Ship = trimws(as.character(pick(c("Ship", "Platform")))),
    HaulNo = as.integer(pick(c("HaulNo", "HaulNumber"))),
    Gear = trimws(as.character(pick("Gear", required = FALSE))),
    HaulDur = as.numeric(pick(c("HaulDur", "HaulDuration"))),
    HaulVal = trimws(as.character(pick(c("HaulVal", "HaulValidity")))),
    DataType = trimws(as.character(pick("DataType"))),
    StatRec = trimws(as.character(pick(c("StatRec", "StatisticalRectangle"), required = FALSE))),
    Month = as.integer(pick("Month", required = FALSE))
  ) %>%
    filter(
      Survey == survey,
      Quarter == quarter,
      Year >= year_min,
      Year <= year_max
    ) %>%
    mutate(
      HaulDur = if_else(!is.na(HaulDur) & HaulDur > 0, HaulDur, NA_real_),
      DataType = if_else(DataType %in% c("", "-9"), NA_character_, DataType),
      haul_key = make_haul_key(Survey, Year, Quarter, Country, Ship, HaulNo)
    ) %>%
    distinct() %>%
    add_count(haul_key, name = "hh_key_n")

  if (nrow(hh) == 0L) {
    stop("No HH records left for survey '", survey, "' Q", quarter, " ", year_min, "-", year_max)
  }
  hh
}

#' HH records from the ICES unaggregated CSV (HH block).
import_hh_from_ices_csv <- function(csv_path, blocks = NULL, ...) {
  if (is.null(blocks)) {
    blocks <- read_ices_csv_blocks(csv_path)
  }
  if (!"HH" %in% names(blocks)) {
    stop("No HH block found in ICES CSV: ", csv_path)
  }
  hh <- standardise_hh(blocks$HH, ...)
  hh$hh_source <- "ICES CSV"
  hh
}

#' Fallback: HH records downloaded from DATRAS.
fetch_hh_from_icesdatras <- function(
    years = ANALYSIS_YEAR_MIN:ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER,
    ...) {
  if (!requireNamespace("icesDatras", quietly = TRUE)) {
    stop("Package icesDatras not installed: renv::install('icesDatras')")
  }
  cat("Downloading HH from DATRAS (NS-IBTS Q", quarter, " ", min(years), "-", max(years), ")...\n", sep = "")
  hh_raw <- icesDatras::getDATRAS("HH", "NS-IBTS", years = years, quarters = quarter)
  if (is.null(hh_raw) || nrow(hh_raw) == 0L) {
    stop("icesDatras::getDATRAS('HH', ...) returned no rows")
  }
  if (!"Survey" %in% names(hh_raw)) {
    hh_raw$Survey <- "NS-IBTS"
  }
  hh <- standardise_hh(hh_raw, ...)
  hh$hh_source <- "icesDatras download"
  hh
}

#' Import HH (CSV first, icesDatras fallback) and save outputs/datras_hh_raw.rds.
ensure_datras_hh_raw <- function(project_root, path_hh_raw, force = FALSE) {
  if (!force && file.exists(path_hh_raw)) {
    cat("datras_hh_raw.rds present — skipping HH import.\n")
    return(invisible(readRDS(path_hh_raw)))
  }
  csv_path <- find_ices_hl_csv(project_root)
  hh <- if (!is.null(csv_path)) {
    import_hh_from_ices_csv(csv_path)
  } else {
    cat("No ICES CSV found — falling back to icesDatras download.\n")
    fetch_hh_from_icesdatras()
  }
  dir.create(dirname(path_hh_raw), recursive = TRUE, showWarnings = FALSE)
  saveRDS(hh, path_hh_raw)
  cat("Saved:", path_hh_raw, "(", nrow(hh), "HH records from", hh$hh_source[[1L]], ")\n")
  invisible(hh)
}
