# Couce et al. (2020) fishing effort import helpers for H2
# Requires h2_common.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

#' Read raw Couce CSV and return year-level effort (Reconstructed gear hours only).
read_couce_year_effort <- function(path,
                                   year_min = H2_YEAR_MIN,
                                   year_max = H2_YEAR_MAX) {
  stopifnot(file.exists(path))

  raw <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  required <- c(
    "Gear", "Year", "ICES_rect", "Data_type", "Hours_trawling"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Couce CSV missing columns: ", paste(missing, collapse = ", "))
  }

  raw %>%
    mutate(
      stat_rec = normalize_stat_rec(ICES_rect),
      year = as.integer(Year),
      gear = Gear,
      data_type = Data_type,
      hours = as.numeric(Hours_trawling)
    ) %>%
    filter(
      data_type == "Reconstructed",
      gear %in% c("Otter", "Beam"),
      year >= year_min,
      year <= year_max,
      !is.na(stat_rec),
      stat_rec != "",
      is.finite(hours)
    ) %>%
    select(stat_rec, year, gear, hours, lon, lat) %>%
    pivot_wider(
      names_from = gear,
      values_from = hours,
      values_fn = list(hours = sum)
    ) %>%
    mutate(
      hours_otter = coalesce(Otter, 0),
      hours_beam = coalesce(Beam, 0),
      hours_total = hours_otter + hours_beam
    ) %>%
    select(stat_rec, year, hours_otter, hours_beam, hours_total, lon, lat)
}

#' Aggregate year-level Couce effort to one row per rectangle.
aggregate_couce_rectangle <- function(year_effort) {
  year_effort %>%
    group_by(stat_rec) %>%
    summarise(
      n_years_with_effort = n_distinct(year),
      total_hours = sum(hours_total, na.rm = TRUE),
      mean_annual_hours_total = mean(hours_total, na.rm = TRUE),
      mean_annual_hours_otter = mean(hours_otter, na.rm = TRUE),
      mean_annual_hours_beam = mean(hours_beam, na.rm = TRUE),
      log_mean_annual_hours_total = log(mean_annual_hours_total + 1),
      rect_lon = mean(lon, na.rm = TRUE),
      rect_lat = mean(lat, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Build import diagnostics comparing Couce rectangles to EEoS haul coverage.
couce_import_diagnostics <- function(year_effort, rectangle_effort, haul_stat_rec) {
  couce_recs <- unique(year_effort$stat_rec)
  haul_recs <- unique(normalize_stat_rec(haul_stat_rec))
  overlap <- intersect(couce_recs, haul_recs)

  data.frame(
    metric = c(
      "n_couce_year_rows",
      "n_couce_rectangles",
      "n_haul_stat_rec",
      "n_overlap_stat_rec",
      "pct_haul_stat_rec_with_couce",
      "year_min",
      "year_max"
    ),
    value = c(
      nrow(year_effort),
      length(couce_recs),
      length(haul_recs),
      length(overlap),
      round(100 * length(overlap) / max(length(haul_recs), 1L), 1),
      min(year_effort$year, na.rm = TRUE),
      max(year_effort$year, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}
