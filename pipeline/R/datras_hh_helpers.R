# DATRAS HH (haul) records and HL <-> HH effort alignment.
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

#' Effort factor converting a DATRAS count to raw fish caught in the haul.
effort_factor <- function(DataType, HaulDur) {
  if_else(DataType == "C", HaulDur / 60, 1)
}

#' Align HL length records to HH haul records and convert counts to raw per haul.
#'
#' Each HL haul gets one effort_reason, assigned in this order:
#'   no_hh            no HH record with this haul key
#'   ambiguous_hh_key >1 HH record share the key (HL cannot be attributed)
#'   invalid          HaulVal not in `haulval`
#'   datatype         DataType missing or not in `datatypes`
#'   duration         HaulDur missing or outside [dur_min, dur_max]
#'   kept
#'
#' Returned `hl` holds kept hauls only, with HLNoAtLngt replaced by the raw count
#' (HLNoAtLngt * SubFactor * f) and SubFactor set to 1, so clean_hl_raw() yields
#' NumberAtLength = raw count. Each row is converted with its own SubFactor, so
#' species split into categories (DataType P, several SubFactors) sum correctly.
#' Originals kept as HLNoAtLngt_orig / SubFactor_orig.
#'
#' @return list(hl, hauls, summary, raw_count_check, match_rate, window)
align_hl_to_hh <- function(
    hl,
    hh,
    dur_min = HAUL_DUR_MIN,
    dur_max = HAUL_DUR_MAX,
    datatypes = HH_DATATYPES_KEEP,
    haulval = HH_HAULVAL_KEEP,
    min_match_rate = HH_HL_MIN_MATCH_RATE,
    year_min = ANALYSIS_YEAR_MIN,
    year_max = ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER) {
  req_hl <- c("Survey", "Year", "Quarter", "Country", "Ship", "HaulNo", "HLNoAtLngt", "SubFactor")
  req_hh <- c("haul_key", "HaulDur", "HaulVal", "DataType", "hh_key_n")
  miss <- c(setdiff(req_hl, names(hl)), setdiff(req_hh, names(hh)))
  if (length(miss) > 0L) {
    stop("align_hl_to_hh(): missing columns: ", paste(miss, collapse = ", "))
  }

  hl <- hl %>%
    filter(Year >= year_min, Year <= year_max, Quarter == quarter) %>%
    mutate(haul_key = make_haul_key(Survey, Year, Quarter, Country, Ship, HaulNo))

  hh_one <- hh %>%
    group_by(haul_key) %>%
    summarise(
      hh_key_n = first(hh_key_n),
      HaulDur = if (n() == 1L) first(HaulDur) else NA_real_,
      HaulVal = if (n() == 1L) first(HaulVal) else NA_character_,
      DataType = if (n() == 1L) first(DataType) else NA_character_,
      .groups = "drop"
    )

  hauls <- hl %>%
    distinct(haul_key, Survey, Year, Quarter, Country, Ship, HaulNo) %>%
    left_join(hh_one, by = "haul_key") %>%
    mutate(
      effort_reason = case_when(
        is.na(hh_key_n) ~ "no_hh",
        hh_key_n > 1L ~ "ambiguous_hh_key",
        is.na(HaulVal) | !HaulVal %in% haulval ~ "invalid",
        is.na(DataType) | !DataType %in% datatypes ~ "datatype",
        is.na(HaulDur) | HaulDur < dur_min | HaulDur > dur_max ~ "duration",
        TRUE ~ "kept"
      ),
      effort_f = effort_factor(DataType, HaulDur)
    )

  match_rate <- mean(!is.na(hauls$hh_key_n))
  cat(sprintf(
    "HL -> HH match: %.2f%% of %d HL hauls (window %g-%g min)\n",
    100 * match_rate, nrow(hauls), dur_min, dur_max
  ))
  if (match_rate < min_match_rate) {
    stop(sprintf(
      "Only %.1f%% of HL hauls match an HH record (< %.0f%%): check Survey/Country/Ship/HaulNo codes.",
      100 * match_rate, 100 * min_match_rate
    ))
  }

  kept_info <- hauls %>%
    filter(effort_reason == "kept") %>%
    select(haul_key, DataType, HaulDur, effort_f)

  hl_kept <- hl %>%
    inner_join(kept_info, by = "haul_key") %>%
    mutate(
      HLNoAtLngt_orig = HLNoAtLngt,
      SubFactor_orig = SubFactor,
      HLNoAtLngt = HLNoAtLngt * SubFactor * effort_f,
      SubFactor = 1
    )

  summary <- hauls %>%
    mutate(DataType = coalesce(DataType, "(none)")) %>%
    count(effort_reason, DataType, name = "n_hauls") %>%
    mutate(
      dur_min = dur_min,
      dur_max = dur_max,
      pct_of_hl_hauls = round(100 * n_hauls / nrow(hauls), 2)
    ) %>%
    arrange(factor(effort_reason, levels = c(
      "kept", "no_hh", "ambiguous_hh_key", "invalid", "datatype", "duration"
    )), DataType)

  raw_count_check <- summarise_raw_counts(hl_kept)

  list(
    hl = hl_kept,
    hauls = hauls,
    summary = summary,
    raw_count_check = raw_count_check,
    match_rate = match_rate,
    window = c(min = dur_min, max = dur_max)
  )
}

#' Share of non-integer raw counts by DataType (and C duration), after alignment.
summarise_raw_counts <- function(hl_aligned, tol = 1e-6) {
  is_nonint <- function(x) abs(x - round(x)) > tol
  hl_aligned %>%
    mutate(
      group = if_else(
        DataType == "C",
        paste0("C (", if_else(HaulDur == 30, "30 min", "not 30 min"), ")"),
        DataType
      ),
      nonint_raw = is_nonint(HLNoAtLngt),
      nonint_reported = is_nonint(HLNoAtLngt_orig),
      subfactor_ne_1 = abs(SubFactor_orig - 1) > tol
    ) %>%
    group_by(DataType, group) %>%
    summarise(
      n_hauls = n_distinct(haul_key),
      n_rows = n(),
      n_rows_subfactor_ne_1 = sum(subfactor_ne_1),
      pct_rows_nonint_reported = round(100 * mean(nonint_reported), 3),
      pct_rows_nonint_raw = round(100 * mean(nonint_raw), 3),
      n_hauls_any_nonint_raw = n_distinct(haul_key[nonint_raw]),
      .groups = "drop"
    ) %>%
    arrange(DataType, group)
}

#' Load HL + HH RDS and return the effort-aligned HL (see align_hl_to_hh()).
#'
#' Every live consumer that computes counts from HL must use this (or
#' align_hl_to_hh()) rather than outputs/datras_hl_raw.rds directly.
load_aligned_hl <- function(project_root, window = haul_dur_window(), ...) {
  path_hl <- file.path(project_root, "outputs", "datras_hl_raw.rds")
  path_hh <- file.path(project_root, "outputs", "datras_hh_raw.rds")
  if (!file.exists(path_hl)) stop("Missing ", path_hl)
  if (!file.exists(path_hh)) stop("Missing ", path_hh, " — run Rscript pipeline/import_datras_hh.R")
  align_hl_to_hh(
    readRDS(path_hl), readRDS(path_hh),
    dur_min = window[["min"]], dur_max = window[["max"]], ...
  )
}
