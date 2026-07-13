# DATRAS HL length-bin processing helpers
# Source R/datras_constants.R before this file.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Convert LngtClass to bin midpoint length in centimetres (ICES LngtCode vocabulary).
#'
#' @param lngt_code Character or numeric LngtCode from DATRAS HL.
#' @param lngt_class Numeric LngtClass (lower edge of bin).
#' @return Numeric midpoint in cm, or NA for unknown codes.
lngt_midpoint_cm <- function(lngt_code, lngt_class) {
  code <- as.character(lngt_code)
  cls <- as.numeric(lngt_class)
  out <- rep(NA_real_, length(cls))

  is1 <- code == "1"
  out[is1] <- cls[is1] + 0.5

  is_mm <- code %in% c("0", ".", "5")
  out[is_mm] <- (cls[is_mm] + 0.5) / 10

  out
}

#' Summarise HL raw coverage for audit diagnostics.
assess_hl_coverage <- function(hl, year_min = ANALYSIS_YEAR_MIN, year_max = ANALYSIS_YEAR_MAX) {
  hl_f <- hl %>%
    filter(
      Year >= year_min,
      Year <= year_max,
      Quarter == ANALYSIS_QUARTER
    )

  n_hauls <- n_distinct(
    paste(
      hl_f$Survey, hl_f$Year, hl_f$Quarter, hl_f$Country,
      hl_f$Ship, hl_f$HaulNo, sep = "_"
    )
  )
  years_present <- sort(unique(hl_f$Year))
  year_span_ok <- length(years_present) >= (year_max - year_min + 1L) * 0.9

  list(
    n_rows = nrow(hl_f),
    n_hauls = n_hauls,
    years_present = years_present,
    year_min = if (length(years_present)) min(years_present) else NA_integer_,
    year_max = if (length(years_present)) max(years_present) else NA_integer_,
    year_span_ok = year_span_ok,
    complete = year_span_ok && n_hauls >= 10000L,
    message = sprintf(
      "HL rows=%s hauls=%s years=%s-%s (complete=%s)",
      nrow(hl_f), n_hauls,
      if (length(years_present)) min(years_present) else NA,
      if (length(years_present)) max(years_present) else NA,
      year_span_ok && n_hauls >= 10000L
    )
  )
}

#' Standardise raw DATRAS HL to length-bin rows with raised counts.
#'
#' Expects columns: Survey, Year, Quarter, Country, Ship, HaulNo, Valid_Aphia,
#' LngtCode, LngtClass, HLNoAtLngt, SubFactor.
clean_hl_raw <- function(
    hl,
    non_fish = NON_FISH_APHIAIDS,
    year_min = ANALYSIS_YEAR_MIN,
    year_max = ANALYSIS_YEAR_MAX,
    quarter = ANALYSIS_QUARTER) {
  required <- c(
    "Survey", "Year", "Quarter", "Country", "Ship", "HaulNo",
    "Valid_Aphia", "LngtCode", "LngtClass", "HLNoAtLngt", "SubFactor"
  )
  missing <- setdiff(required, names(hl))
  if (length(missing) > 0L) {
    stop("HL data missing columns: ", paste(missing, collapse = ", "))
  }

  hl %>%
    filter(
      Year >= year_min,
      Year <= year_max,
      Quarter == quarter,
      !is.na(Valid_Aphia),
      !Valid_Aphia %in% non_fish,
      !is.na(HLNoAtLngt),
      HLNoAtLngt > 0,
      !is.na(LngtClass)
    ) %>%
    mutate(
      AphiaID = as.integer(Valid_Aphia),
      Platform = Ship,
      HaulNumber = HaulNo,
      haul_key = paste(Survey, Year, Quarter, Country, Ship, HaulNo, sep = "_"),
      length_cm = lngt_midpoint_cm(LngtCode, LngtClass),
      length_mm = length_cm * 10,
      NumberAtLength = HLNoAtLngt * SubFactor,
      lngt_code = as.character(LngtCode)
    ) %>%
    filter(
      !is.na(length_cm),
      length_mm > 0,
      NumberAtLength > 0
    )
}

#' Join LW parameters and compute mass (g) per length bin row.
#'
#' Drops bins with missing/non-finite LW parameters or mass so that a single
#' bad bin cannot make haul-level sum()/min() return NA.
add_lw_mass_to_bins <- function(hl_bins, lw_lookup) {
  lw <- lw_lookup %>%
    select(aphia_id, a_final, b_final) %>%
    filter(
      !is.na(a_final),
      !is.na(b_final),
      is.finite(a_final),
      is.finite(b_final)
    )

  hl_bins %>%
    inner_join(lw, by = c("AphiaID" = "aphia_id")) %>%
    mutate(
      mass_g = a_final * (length_cm ^ b_final),
      epsilon = mass_g ^ 0.75
    ) %>%
    filter(is.finite(mass_g), mass_g > 0, is.finite(epsilon), epsilon > 0)
}

#' Species x haul aggregates from length-bin rows (for mean length table).
aggregate_species_haul_from_bins <- function(hl_bins_with_mass) {
  hl_bins_with_mass %>%
    group_by(
      haul_key, Survey, Year, Quarter, Country, Platform, HaulNumber,
      AphiaID
    ) %>%
    summarise(
      mean_length_mm = weighted.mean(length_mm, w = NumberAtLength, na.rm = TRUE),
      min_length_mm = min(length_mm, na.rm = TRUE),
      total_n_measured = sum(NumberAtLength, na.rm = TRUE),
      n_length_bins = n(),
      .groups = "drop"
    ) %>%
    filter(total_n_measured > 0, is.finite(mean_length_mm), mean_length_mm > 0)
}

#' Haul-level state variables from length-bin rows (audit-compliant S, N, E_raw, m_min).
aggregate_haul_state_from_bins <- function(hl_bins_with_mass) {
  hl_bins_with_mass %>%
    group_by(
      haul_key, Survey, Year, Quarter, Country, Platform, HaulNumber
    ) %>%
    summarise(
      S = n_distinct(AphiaID),
      N = sum(NumberAtLength, na.rm = TRUE),
      E_raw = sum(NumberAtLength * mass_g ^ 0.75, na.rm = TRUE),
      m_min = min(mass_g, na.rm = TRUE),
      n_species_with_lw = n_distinct(AphiaID),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(E_raw), E_raw > 0,
      is.finite(m_min), m_min > 0,
      N >= S
    ) %>%
    mutate(m_min_method = "min_length_bin")
}

#' Minimum m_min per haul from HL bins only (for hybrid mode).
derive_m_min_from_hl <- function(hl, lw_lookup, non_fish = NON_FISH_APHIAIDS) {
  hl_bins <- clean_hl_raw(hl, non_fish = non_fish)
  if (nrow(hl_bins) == 0L) {
    return(tibble(haul_key = character(), m_min = numeric(), m_min_method = character()))
  }

  hl_bins %>%
    add_lw_mass_to_bins(lw_lookup) %>%
    group_by(haul_key) %>%
    summarise(
      m_min = min(mass_g, na.rm = TRUE),
      m_min_method = "min_length_bin",
      .groups = "drop"
    ) %>%
    filter(is.finite(m_min), m_min > 0)
}

#' LngtCode audit table for diagnostics.
lngt_code_audit_table <- function(hl) {
  hl %>%
    mutate(lngt_code = as.character(LngtCode)) %>%
    group_by(lngt_code) %>%
    summarise(
      n_rows = n(),
      example_midpoint_cm = lngt_midpoint_cm(lngt_code[1L], LngtClass[1L]),
      .groups = "drop"
    ) %>%
    mutate(known_code = lngt_code %in% KNOWN_LNGT_CODES)
}
