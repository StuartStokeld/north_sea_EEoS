# Pre-H3 exploratory visualisation & variance decomposition helpers.
# See CURSOR_BRIEFING "Pre-H3 Exploratory Visualisation & Variance
# Decomposition" (chat-supplied, not a repo file) for the full spec these
# helpers implement.
#
# DESCRIPTIVE ONLY: no temporal model is fit here (no panel regression, no
# mediation test, no fixed effects). Rectangle-level "slopes" (fit_rectangle_slopes)
# and variance decompositions (variance_components_anova) are plotting/reporting
# aids only — p-values are deliberately not computed or reported for slopes.
#
# Requires h1_common.R (project root), h2_common.R (normalize_stat_rec,
# H2_YEAR_MIN/MAX).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Section A helpers — haul temporal/spatial coverage
# ---------------------------------------------------------------------------

#' One row per stat_rec x year haul count, from the full Q1 haul table
#' (build_fishglob_haul_table() output) — independent of DATRAS join / EEoS
#' filter success, so this reflects raw survey coverage, not modelling dropout.
build_rect_year_haul_counts <- function(haul_full, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX) {
  required <- c("stat_rec", "year")
  missing <- setdiff(required, names(haul_full))
  if (length(missing) > 0L) {
    stop("haul_full missing columns: ", paste(missing, collapse = ", "))
  }
  haul_full %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(!is.na(stat_rec), stat_rec != "", year >= year_min, year <= year_max) %>%
    count(stat_rec, year, name = "n_hauls")
}

#' Complete stat_rec x year grid (all rectangles ever observed x all years in
#' range), zero-filled — for heatmap display only (A.2). Zero is a genuine
#' "no haul that year", not a missing value.
complete_rect_year_grid <- function(rect_year_hauls, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX) {
  rect_year_hauls %>%
    tidyr::complete(stat_rec, year = year_min:year_max, fill = list(n_hauls = 0L))
}

#' Rectangle usability flags for within-rectangle temporal analysis (A.3).
#' "Sparse" years (present but below sparse_threshold hauls) are counted, not
#' dropped. `usable_temporal` is the sole A.3 threshold (haul-coverage based);
#' Couce fishing-pressure coverage is layered on separately (see F / add_couce_coverage_flag()).
build_rectangle_usability_flags <- function(rect_year_hauls,
                                            sparse_threshold,
                                            min_years) {
  rect_year_hauls %>%
    group_by(stat_rec) %>%
    summarise(
      n_years_present = dplyr::n_distinct(year[n_hauls > 0]),
      n_years_qualifying = sum(n_hauls >= sparse_threshold),
      n_years_sparse = sum(n_hauls > 0 & n_hauls < sparse_threshold),
      total_hauls = sum(n_hauls),
      year_min = suppressWarnings(min(year[n_hauls > 0], na.rm = TRUE)),
      year_max = suppressWarnings(max(year[n_hauls > 0], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      sparse_haul_threshold = sparse_threshold,
      min_years_per_rect = min_years,
      usable_temporal = n_years_qualifying >= min_years
    )
}

#' Add a binary Couce fishing-pressure coverage flag (any year, any
#' rectangle-level record in couce_year) and the combined usability flag used
#' by D/E. Rectangles with no Couce record at all are the "structural
#' exclusion" documented in the briefing (F) — logged by the caller, not here.
add_couce_coverage_flag <- function(rect_flags, couce_year) {
  covered <- unique(normalize_stat_rec(couce_year$stat_rec))
  rect_flags %>%
    mutate(
      has_couce_coverage = stat_rec %in% covered,
      usable_for_fishing_analysis = usable_temporal & has_couce_coverage
    )
}

# ---------------------------------------------------------------------------
# Section B/C helpers — EEoS haul-level prep and year-level summaries
# ---------------------------------------------------------------------------

#' Add the briefing-defined signed residual to haul_eeos_predictions.rds.
#' NOTE ON SIGN CONVENTION: this is the *negative* of the H1/H2 pipeline's
#' primary `residual` column (residual = log(B_obs) - log(B_pred); see
#' pipeline/README.md "Key conventions"). The briefing specifies
#' resid_signed = log(B_pred) - log(B_obs). Magnitude (abs_residual) is
#' identical under either sign convention and is reused as-is.
add_resid_signed <- function(haul_eeos) {
  required <- c("stat_rec", "year", "ln_B_pred", "ln_B_obs", "abs_residual")
  missing <- setdiff(required, names(haul_eeos))
  if (length(missing) > 0L) {
    stop("haul_eeos missing columns: ", paste(missing, collapse = ", "))
  }
  haul_eeos %>%
    mutate(
      stat_rec = normalize_stat_rec(stat_rec),
      resid_signed = ln_B_pred - ln_B_obs,
      resid_magnitude = abs_residual
    )
}

#' Generic per-year summary: center (mean/median) +/- spread (95% CI / IQR),
#' used for A.1-style series in B.1, C.1, D.1. Returns year, n, center, lo, hi.
summarise_year_stat <- function(df, value_col, year_min, year_max,
                                center = c("mean", "median"),
                                spread = c("ci", "iqr")) {
  center <- match.arg(center)
  spread <- match.arg(spread)
  stopifnot(value_col %in% names(df))

  df <- df %>%
    filter(year >= year_min, year <= year_max, is.finite(.data[[value_col]]))

  df %>%
    group_by(year) %>%
    summarise(
      n = dplyr::n(),
      center = if (center == "mean") mean(.data[[value_col]]) else median(.data[[value_col]]),
      sd_val = sd(.data[[value_col]]),
      lo = if (spread == "ci") {
        center - 1.96 * sd_val / sqrt(n)
      } else {
        stats::quantile(.data[[value_col]], 0.25, names = FALSE)
      },
      hi = if (spread == "ci") {
        center + 1.96 * sd_val / sqrt(n)
      } else {
        stats::quantile(.data[[value_col]], 0.75, names = FALSE)
      },
      .groups = "drop"
    ) %>%
    select(-sd_val) %>%
    mutate(center_type = center, spread_type = spread)
}

# ---------------------------------------------------------------------------
# Section C/D helpers — decade binning, spatial maps, rectangle-level slopes
# ---------------------------------------------------------------------------

#' Provisional decade-bin label for a year vector. decade_bins is a list of
#' c(start_year, end_year) pairs (DECADE_BINS). Years outside all bins -> NA
#' (dropped by the caller, not silently coerced into a bin).
assign_decade_bin <- function(year, decade_bins) {
  labels <- vapply(decade_bins, function(b) sprintf("%d-%d", b[1], b[2]), character(1))
  out <- rep(NA_character_, length(year))
  for (i in seq_along(decade_bins)) {
    b <- decade_bins[[i]]
    idx <- year >= b[1] & year <= b[2]
    out[idx] <- labels[i]
  }
  factor(out, levels = labels)
}

#' Rectangle x decade mean of a value column, for decade-faceted spatial maps
#' (C.2 / D.2). Rows with year outside all DECADE_BINS are dropped and should
#' be counted by the caller.
build_rect_decade_summary <- function(df, value_col, decade_bins) {
  stopifnot(value_col %in% names(df))
  df %>%
    mutate(decade = assign_decade_bin(year, decade_bins)) %>%
    filter(!is.na(decade)) %>%
    group_by(stat_rec, decade) %>%
    summarise(
      n = dplyr::n(),
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

#' Simple OLS slope of value_col ~ year, one fit per rectangle, restricted to
#' `usable_rects`. Descriptive plotting aid only (C.3 / D.3) — returns slope
#' magnitude and direction; does NOT return or report a p-value by design.
fit_rectangle_slopes <- function(rect_year_values, usable_rects, value_col, min_points = 2L) {
  stopifnot(value_col %in% names(rect_year_values))

  df <- rect_year_values %>%
    filter(stat_rec %in% usable_rects) %>%
    rename(value = all_of(value_col))

  if (nrow(df) == 0L) {
    return(tibble(
      stat_rec = character(0), n_years = integer(0), year_min = integer(0),
      year_max = integer(0), slope = numeric(0), direction = character(0)
    ))
  }

  df %>%
    group_by(stat_rec) %>%
    group_modify(function(sub, key) {
      if (nrow(sub) < min_points) {
        return(tibble(
          n_years = nrow(sub), year_min = suppressWarnings(min(sub$year)),
          year_max = suppressWarnings(max(sub$year)), slope = NA_real_,
          direction = NA_character_
        ))
      }
      fit <- stats::lm(value ~ year, data = sub)
      slope_val <- unname(stats::coef(fit)["year"])
      tibble(
        n_years = nrow(sub),
        year_min = min(sub$year),
        year_max = max(sub$year),
        slope = slope_val,
        direction = dplyr::case_when(
          slope_val > 0 ~ "increasing",
          slope_val < 0 ~ "decreasing",
          TRUE ~ "flat"
        )
      )
    }) %>%
    ungroup()
}

# ---------------------------------------------------------------------------
# Section E helpers — one-way variance-components decomposition
# ---------------------------------------------------------------------------

#' One-way random-effects variance decomposition by the classical ANOVA
#' method-of-moments estimator for unbalanced designs (Searle, Casella &
#' McCulloch, 1992, "Variance Components", sec. 3.4 / Henderson's Method 1).
#' Statistically equivalent to lme4::VarCorr() on an intercept-only mixed
#' model with `group` as a random effect ("... or equivalent" per the
#' briefing) — implemented directly here rather than adding lme4 as a project
#' dependency.
#'
#' Returns NA fields if there are fewer than 2 groups or no within-group
#' degrees of freedom (cannot decompose).
variance_components_anova <- function(value, group) {
  keep <- is.finite(value) & !is.na(group)
  value <- value[keep]
  group <- as.character(group[keep])

  empty <- list(
    n_groups = dplyr::n_distinct(group), n_obs = length(value),
    mean_n_per_group = NA_real_, min_n_per_group = NA_real_, max_n_per_group = NA_real_,
    var_between = NA_real_, var_within = NA_real_,
    sd_between = NA_real_, sd_within = NA_real_,
    icc = NA_real_, sd_ratio_within_over_between = NA_real_
  )

  agg <- tapply(value, group, function(x) list(n = length(x), mean = mean(x)))
  k <- length(agg)
  if (k < 2L) {
    return(empty)
  }

  n_i <- vapply(agg, function(a) a$n, numeric(1))
  mean_i <- vapply(agg, function(a) a$mean, numeric(1))
  N <- sum(n_i)
  df_between <- k - 1L
  df_within <- N - k
  if (df_within <= 0L) {
    return(empty)
  }

  grand_mean <- mean(value)
  ss_between <- sum(n_i * (mean_i - grand_mean)^2)
  group_mean_lookup <- setNames(mean_i, names(agg))
  ss_within <- sum((value - group_mean_lookup[group])^2)

  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within
  n0 <- (N - sum(n_i^2) / N) / df_between

  var_within <- ms_within
  var_between <- max(0, (ms_between - ms_within) / n0)
  sd_between <- sqrt(var_between)
  sd_within <- sqrt(var_within)
  icc <- if ((var_between + var_within) > 0) var_between / (var_between + var_within) else NA_real_

  list(
    n_groups = k,
    n_obs = N,
    mean_n_per_group = mean(n_i),
    min_n_per_group = min(n_i),
    max_n_per_group = max(n_i),
    var_between = var_between,
    var_within = var_within,
    sd_between = sd_between,
    sd_within = sd_within,
    icc = icc,
    sd_ratio_within_over_between = if (sd_between > 0) sd_within / sd_between else NA_real_
  )
}

#' Convenience wrapper: run variance_components_anova() and return a one-row
#' tibble tagged with a `variable` label, for stacking into a summary table.
variance_components_row <- function(value, group, variable_label) {
  vc <- variance_components_anova(value, group)
  tibble(variable = variable_label) %>%
    bind_cols(as_tibble(vc))
}
