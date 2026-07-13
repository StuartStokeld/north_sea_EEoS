# FishGlob <-> DATRAS join helpers (H1 pipeline)
# Ensures only hauls with valid normalised E enter the join — prevents null-E artifacts.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
})

#' Build FishGlob Q1 haul table with B_obs (grams) and 6-field join key.
#'
#' Join key matches DATRAS haul_key: Survey_Year_Quarter_Country_Ship_HaulNo.
#' Gear and StNo from FishGlob haul_id are dropped (not in DATRAS HL download).
build_fishglob_haul_table <- function(
    path_fishglob,
    year_min = 1985L,
    year_max = 2015L,
    analysis_quarter = 1L) {
  env <- new.env(parent = emptyenv())
  load(path_fishglob, envir = env)
  if (!exists("data", envir = env)) {
    stop("FishGlob RData must contain object 'data'")
  }
  fg <- env$data

  fg %>%
    filter(
      quarter == analysis_quarter,
      year >= year_min,
      year <= year_max
    ) %>%
    group_by(haul_id) %>%
    summarise(
      B_obs = sum(wgt, na.rm = TRUE) * 1000,
      stat_rec = first(stat_rec),
      year = first(year),
      .groups = "drop"
    ) %>%
    mutate(
      haul_parts = str_split(haul_id, " ", simplify = FALSE),
      join_key = map_chr(haul_parts, function(parts) {
        if (length(parts) < 8L) {
          return(NA_character_)
        }
        paste(parts[1], parts[2], parts[3], parts[4], parts[5], parts[8], sep = "_")
      })
    ) %>%
    filter(!is.na(join_key)) %>%
    select(-haul_parts)
}

#' Restrict DATRAS state to hauls with complete, finite normalised E before join.
prepare_datras_for_join <- function(datras_state) {
  required <- c("haul_key", "S", "N", "E", "E_raw", "m_min")
  missing <- setdiff(required, names(datras_state))
  if (length(missing) > 0L) {
    stop("datras_state missing columns: ", paste(missing, collapse = ", "))
  }

  datras_state %>%
    filter(
      !is.na(E),
      is.finite(E),
      !is.na(E_raw),
      is.finite(E_raw),
      E_raw > 0,
      !is.na(m_min),
      is.finite(m_min),
      m_min > 0,
      !is.na(S),
      !is.na(N),
      N >= S
    )
}

#' Inner join FishGlob B_obs to DATRAS state; validate no null E after join.
join_fishglob_datras <- function(fishglob_haul, datras_ready) {
  joined <- fishglob_haul %>%
    inner_join(datras_ready, by = c("join_key" = "haul_key"))

  n_null_e <- sum(is.na(joined$E))
  if (n_null_e > 0L) {
    stop(
      "Join produced ", n_null_e,
      " hauls with NA E — datras_ready must be pre-filtered; check prepare_datras_for_join()"
    )
  }

  joined
}

#' Document FishGlob and DATRAS hauls that fail to match.
summarise_join_gaps <- function(fishglob_haul, datras_ready) {
  fg_only <- fishglob_haul %>%
    anti_join(datras_ready, by = c("join_key" = "haul_key")) %>%
    mutate(gap_side = "fishglob_no_datras")

  datras_only <- datras_ready %>%
    anti_join(fishglob_haul, by = c("haul_key" = "join_key")) %>%
    mutate(gap_side = "datras_no_fishglob", join_key = haul_key)

  bind_rows(
    fg_only %>% select(gap_side, join_key, year, B_obs, stat_rec),
    datras_only %>%
      select(gap_side, join_key, Year, S, N, E, E_raw, m_min) %>%
      rename(year = Year)
  )
}

#' Assign a single exclusion reason per haul for EEoS filter step.
classify_eeos_exclusion <- function(haul) {
  case_when(
    is.na(haul$E) | !is.finite(haul$E) ~ "na_E",
    is.na(haul$E_raw) | !is.finite(haul$E_raw) ~ "na_E_raw",
    is.na(haul$m_min) | haul$m_min <= 0 ~ "bad_m_min",
    is.na(haul$S) | haul$S < 2 ~ "S_lt_2",
    haul$E <= haul$N ~ "E_le_N",
    haul$N < haul$S ~ "N_lt_S",
    is.na(haul$B_obs) | haul$B_obs <= 0 ~ "bad_B_obs",
    TRUE ~ "passed"
  )
}

#' Apply EEoS input filters (same rules as build_eeos_predictions.R).
filter_eeos_inputs <- function(haul_joined) {
  haul_joined %>%
    filter(
      !is.na(E),
      is.finite(E),
      !is.na(E_raw),
      !is.na(m_min),
      m_min > 0,
      S >= 2,
      E > N,
      N >= S,
      !is.na(B_obs),
      B_obs > 0
    )
}

#' Standard haul_state_variables columns for downstream scripts.
select_haul_state_columns <- function(haul_joined) {
  haul_joined %>%
    select(
      haul_id,
      haul_key = join_key,
      year,
      stat_rec,
      S,
      N,
      E,
      E_raw,
      B_obs,
      m_min,
      m_min_epsilon,
      min_epsilon,
      n_species_with_lw,
      m_min_method
    )
}
