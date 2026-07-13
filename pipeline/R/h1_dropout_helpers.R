# H1 pipeline dropout diagnostics
# Source R/h1_common.R before this file when using log metrics.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Summarise haul exclusion counts by year for joined FishGlob panel.
dropout_by_year <- function(haul_state_variables, haul_predictions) {
  all_years <- haul_state_variables %>%
    count(year, name = "n_joined")

  pred_years <- haul_predictions %>%
    count(year, name = "n_predicted")

  all_years %>%
    left_join(pred_years, by = "year") %>%
    mutate(
      n_predicted = coalesce(n_predicted, 0L),
      n_dropped = n_joined - n_predicted,
      pct_dropped = round(100 * n_dropped / n_joined, 1)
    ) %>%
    arrange(year)
}

#' Summarise haul exclusion by ICES rectangle.
dropout_by_stat_rec <- function(haul_state_variables, haul_predictions, min_joined = 5L) {
  joined <- haul_state_variables %>%
    count(stat_rec, name = "n_joined")

  pred <- haul_predictions %>%
    count(stat_rec, name = "n_predicted")

  joined %>%
    left_join(pred, by = "stat_rec") %>%
    mutate(
      n_predicted = coalesce(n_predicted, 0L),
      n_dropped = n_joined - n_predicted,
      pct_dropped = round(100 * n_dropped / pmax(n_joined, 1L), 1)
    ) %>%
    filter(n_joined >= min_joined) %>%
    arrange(desc(pct_dropped), desc(n_dropped))
}

#' Per-haul exclusion reasons after join.
dropout_reason_table <- function(haul_joined, haul_predictions, classify_fn = NULL) {
  if (is.null(classify_fn)) {
    classify_fn <- function(df) {
      # Inline fallback if join helpers not sourced
      dplyr::case_when(
        is.na(df$E) | !is.finite(df$E) ~ "na_E",
        is.na(df$E_raw) | !is.finite(df$E_raw) ~ "na_E_raw",
        is.na(df$m_min) | df$m_min <= 0 ~ "bad_m_min",
        is.na(df$S) | df$S < 2 ~ "S_lt_2",
        df$E <= df$N ~ "E_le_N",
        df$N < df$S ~ "N_lt_S",
        is.na(df$B_obs) | df$B_obs <= 0 ~ "bad_B_obs",
        TRUE ~ "passed"
      )
    }
  }

  joined <- haul_joined
  joined$exclusion_reason <- classify_fn(joined)

  excluded <- joined %>%
    anti_join(haul_predictions %>% select(haul_id), by = "haul_id") %>%
    count(year, exclusion_reason, name = "n_hauls") %>%
    arrange(year, desc(n_hauls))

  passed <- joined %>%
    semi_join(haul_predictions %>% select(haul_id), by = "haul_id") %>%
    count(year, name = "n_passed")

  list(
    excluded = excluded,
    passed = passed,
    summary = excluded %>%
      group_by(exclusion_reason) %>%
      summarise(n_hauls = sum(n_hauls), .groups = "drop") %>%
      arrange(desc(n_hauls))
  )
}

#' Year-level funnel: HL raw hauls -> DATRAS state -> FishGlob join -> EEoS predictions.
dropout_funnel_by_year <- function(
    hl_raw,
    datras_state,
    haul_state_variables,
    haul_predictions,
    year_min = 1985L,
    year_max = 2015L,
    quarter = 1L) {
  hl_years <- hl_raw %>%
    filter(Year >= year_min, Year <= year_max, Quarter == quarter) %>%
    mutate(haul_key = paste(Survey, Year, Quarter, Country, Ship, HaulNo, sep = "_")) %>%
    distinct(haul_key, Year) %>%
    count(Year, name = "n_hl_hauls")

  datras_years <- datras_state %>%
    count(Year, name = "n_datras_state")

  joined_years <- haul_state_variables %>%
    count(year, name = "n_fishglob_joined")

  pred_years <- haul_predictions %>%
    count(year, name = "n_eeos_predictions")

  hl_years %>%
    full_join(datras_years, by = c("Year" = "Year")) %>%
    full_join(joined_years, by = c("Year" = "year")) %>%
    full_join(pred_years, by = c("Year" = "year")) %>%
    mutate(
      n_hl_hauls = coalesce(n_hl_hauls, 0L),
      n_datras_state = coalesce(n_datras_state, 0L),
      n_fishglob_joined = coalesce(n_fishglob_joined, 0L),
      n_eeos_predictions = coalesce(n_eeos_predictions, 0L),
      hl_to_datras_drop = n_hl_hauls - n_datras_state,
      datras_to_join_drop = n_datras_state - n_fishglob_joined,
      join_to_pred_drop = n_fishglob_joined - n_eeos_predictions,
      pct_pred_of_hl = round(100 * n_eeos_predictions / pmax(n_hl_hauls, 1L), 1)
    ) %>%
    arrange(Year)
}

#' Focused diagnosis for years that showed spurious dropout spikes in earlier drafts.
diagnose_spike_years <- function(funnel_by_year, spike_years = c(1998L, 2013L, 2014L)) {
  funnel_by_year %>%
    filter(Year %in% spike_years) %>%
    mutate(
      spike_status = case_when(
        hl_to_datras_drop > 0 ~ "HL hauls lost at LW/state build",
        datras_to_join_drop > 0 ~ "DATRAS hauls without FishGlob B_obs",
        join_to_pred_drop > 0 ~ "Joined hauls failed EEoS filters",
        TRUE ~ "no_dropout_spike"
      )
    )
}

#' One-row verdict for audit documentation.
dropout_diagnosis_verdict <- function(funnel_by_year, spike_years = c(1998L, 2013L, 2014L)) {
  spike <- diagnose_spike_years(funnel_by_year, spike_years)
  tibble(
    check = c(
      "spike_years_hl_retention",
      "spike_years_join_retention",
      "spike_years_prediction_retention"
    ),
    status = c(
      if (all(spike$hl_to_datras_drop == 0L)) "PASS" else "WARN",
      if (all(spike$datras_to_join_drop <= 2L)) "PASS" else "WARN",
      if (all(spike$join_to_pred_drop <= 1L)) "PASS" else "WARN"
    ),
    detail = c(
      paste(spike_years, collapse = ", "),
      sprintf(
        "max datras->join drop = %d (%s)",
        max(spike$datras_to_join_drop, na.rm = TRUE),
        paste(spike$Year[spike$datras_to_join_drop == max(spike$datras_to_join_drop, na.rm = TRUE)], collapse = ", ")
      ),
      sprintf(
        "max join->pred drop = %d (%s)",
        max(spike$join_to_pred_drop, na.rm = TRUE),
        paste(spike$Year[spike$join_to_pred_drop == max(spike$join_to_pred_drop, na.rm = TRUE)], collapse = ", ")
      )
    )
  )
}
