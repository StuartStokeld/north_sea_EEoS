# DATRAS haul-level state variable assembly (S, N, E_raw, m_min)
# Source R/datras_constants.R and R/datras_hl_helpers.R before this file.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Derive S, N, E_raw from species x haul mean-length table + LW lookup.
#'
#' Uses the same species filter and LW join as the original E pipeline, but
#' computes S and N alongside E_raw so all three share one source.
derive_state_from_mean_length <- function(
    haul_mean_length,
    lw_lookup,
    non_fish = NON_FISH_APHIAIDS) {
  lw <- lw_lookup %>%
    select(aphia_id, a_final, b_final)

  species_haul <- haul_mean_length %>%
    filter(!AphiaID %in% non_fish) %>%
    inner_join(lw, by = c("AphiaID" = "aphia_id")) %>%
    filter(
      !is.na(a_final),
      !is.na(b_final),
      is.finite(a_final),
      is.finite(b_final)
    ) %>%
    mutate(
      length_cm = mean_length_mm / 10,
      mass_g = a_final * (length_cm ^ b_final)
    ) %>%
    filter(is.finite(mass_g), mass_g > 0)

  species_haul %>%
    group_by(
      haul_key, Survey, Year, Quarter, Country, Platform, HaulNumber
    ) %>%
    summarise(
      S = n_distinct(AphiaID),
      N = sum(total_n_measured, na.rm = TRUE),
      E_raw = sum(total_n_measured * mass_g ^ 0.75, na.rm = TRUE),
      n_species_with_lw = n_distinct(AphiaID),
      m_min_mean_length = min(mass_g, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(E_raw), E_raw > 0,
      is.finite(m_min_mean_length), m_min_mean_length > 0,
      N >= S
    )
}

#' Full haul state from HL length bins (preferred when HL covers analysis period).
build_state_from_hl <- function(hl, lw_lookup, non_fish = NON_FISH_APHIAIDS) {
  hl_bins <- clean_hl_raw(hl, non_fish = non_fish)
  hl_mass <- add_lw_mass_to_bins(hl_bins, lw_lookup)

  list(
    haul_mean_length = aggregate_species_haul_from_bins(hl_mass),
    haul_state = aggregate_haul_state_from_bins(hl_mass)
  )
}

#' Merge HL-based m_min onto mean-length state; fallback where HL missing.
apply_m_min_with_fallback <- function(haul_state, m_min_hl) {
  haul_state %>%
    left_join(
      m_min_hl %>% select(haul_key, m_min_hl = m_min),
      by = "haul_key"
    ) %>%
    mutate(
      m_min = coalesce(m_min_hl, m_min_mean_length),
      m_min_method = if_else(
        !is.na(m_min_hl),
        "min_length_bin",
        "mean_length_fallback"
      )
    ) %>%
    select(-m_min_hl, -m_min_mean_length)
}

#' Keep hauls with complete, finite state variables for EEoS.
filter_valid_haul_state <- function(haul_state) {
  haul_state %>%
    filter(
      !is.na(S),
      !is.na(N),
      N >= S,
      !is.na(E_raw),
      is.finite(E_raw),
      E_raw > 0,
      !is.na(m_min),
      is.finite(m_min),
      m_min > 0
    )
}

#' Add E normalisation columns required by EEoS (min individual metabolic rate = 1).
normalize_haul_E <- function(haul_state) {
  haul_state %>%
    mutate(
      m_min_epsilon = m_min ^ 0.75,
      E = E_raw / m_min_epsilon,
      min_epsilon = 1
    )
}

#' Validate internal consistency of haul state (audit checks).
validate_haul_state <- function(haul_state, tol = 1e-8) {
  checks <- list(
    n_species_le_s = all(haul_state$n_species_with_lw <= haul_state$S + tol, na.rm = TRUE),
    n_species_eq_s = all(haul_state$n_species_with_lw == haul_state$S, na.rm = TRUE),
    n_ge_s = all(haul_state$N >= haul_state$S, na.rm = TRUE),
    e_gt_n = all(haul_state$E > haul_state$N, na.rm = TRUE),
    e_norm_ratio = all(
      abs(haul_state$E * haul_state$m_min_epsilon - haul_state$E_raw) < tol * pmax(1, haul_state$E_raw),
      na.rm = TRUE
    ),
    min_epsilon_one = all(abs(haul_state$min_epsilon - 1) < tol, na.rm = TRUE)
  )

  tibble(
    check = names(checks),
    pass = vapply(checks, isTRUE, logical(1)),
    detail = c(
      "n_species_with_lw <= S",
      "n_species_with_lw == S",
      "N >= S",
      "E > N (after normalisation)",
      "E == E_raw / m_min^0.75",
      "min_epsilon == 1"
    )
  )
}
