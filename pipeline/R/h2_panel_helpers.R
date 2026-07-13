# Rectangle-level EEoS panel helpers for H2
# Requires h2_common.R

suppressPackageStartupMessages({
  library(dplyr)
})

#' Aggregate haul-level EEoS residuals to one row per ICES rectangle.
build_h2_rectangle_residuals <- function(haul,
                                         year_min = H2_YEAR_MIN,
                                         year_max = H2_YEAR_MAX) {
  required <- c("stat_rec", "year", "residual", "abs_residual", "ln_B_obs")
  missing <- setdiff(required, names(haul))
  if (length(missing) > 0L) {
    stop("haul predictions missing columns: ", paste(missing, collapse = ", "))
  }

  haul %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(
      !is.na(stat_rec),
      stat_rec != "",
      year >= year_min,
      year <= year_max,
      is.finite(residual),
      is.finite(abs_residual)
    ) %>%
    group_by(stat_rec) %>%
    summarise(
      n_hauls = n(),
      mean_abs_residual = mean(abs_residual, na.rm = TRUE),
      mean_residual = mean(residual, na.rm = TRUE),
      median_abs_residual = median(abs_residual, na.rm = TRUE),
      sd_abs_residual = sd(abs_residual, na.rm = TRUE),
      mean_ln_B_obs = mean(ln_B_obs, na.rm = TRUE),
      year_min = min(year, na.rm = TRUE),
      year_max = max(year, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Join residuals to Couce fishing pressure and apply inclusion rules.
build_h2_analysis_panel <- function(residual_panel,
                                    couce_rectangle,
                                    min_hauls = H2_MIN_HAULS_DEFAULT,
                                    require_fishing = TRUE) {
  panel <- residual_panel %>%
    left_join(couce_rectangle, by = "stat_rec")

  if (require_fishing) {
    panel <- panel %>% filter(!is.na(mean_annual_hours_total))
  }

  panel %>%
    filter(n_hauls >= min_hauls) %>%
    mutate(
      log_mean_annual_hours_total = log(mean_annual_hours_total + 1)
    )
}
