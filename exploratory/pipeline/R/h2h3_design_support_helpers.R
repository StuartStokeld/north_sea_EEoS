# H2/H3 design-decision quantitative support helpers.
# See CURSOR_BRIEFING "Quantitative Support for H2/H3 Design Decision"
# (chat-supplied, not a repo file) for the full spec these helpers implement.
#
# REPORTING ONLY: no design recommendation is made or implied by any function
# here. Phase boundaries, decade bins, tercile cutoffs, and the relaxed
# rectangle-inclusion rule are all provisional inputs, surfaced by the caller
# alongside every output, not committed choices.
#
# Requires h1_common.R, h2_common.R (normalize_stat_rec, H2_YEAR_MIN/MAX),
# h3_pre_exploration_helpers.R (assign_decade_bin, variance_components_anova,
# variance_components_row), h3_policy_zones_helpers.R (assign_period,
# required_n_for_reliability, build_scheme_a_blocks, build_contiguous_zones,
# summarise_unit_period_hauls, summarise_unit_period_fishing,
# pivot_unit_period_hauls_wide).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Section A helpers — whole-study-area temporal trend in fishing pressure
# ---------------------------------------------------------------------------

#' One row per year: mean/median fishing hours across all Couce-covered
#' rectangles that have a record that year, plus the coverage denominator
#' (n_rect_contributing) so a reader can see the mean means something
#' different in a year with many rectangles vs few.
build_year_fishing_summary <- function(couce_year, year_min, year_max) {
  stopifnot(all(c("stat_rec", "year", "hours_total") %in% names(couce_year)))
  couce_year %>%
    filter(year >= year_min, year <= year_max, is.finite(hours_total)) %>%
    group_by(year) %>%
    summarise(
      n_rect_contributing = dplyr::n_distinct(stat_rec),
      mean_hours = mean(hours_total),
      median_hours = median(hours_total),
      sd_hours = sd(hours_total),
      .groups = "drop"
    ) %>%
    arrange(year)
}

#' Absolute + percent change in mean fishing hours between two specific years
#' (looked up directly from year_summary, e.g. output of
#' build_year_fishing_summary()). Returns NA if either year is missing.
phase_boundary_change <- function(year_summary, year_from, year_to, label) {
  row_from <- year_summary %>% filter(year == year_from)
  row_to <- year_summary %>% filter(year == year_to)
  mean_from <- if (nrow(row_from) == 1L) row_from$mean_hours else NA_real_
  mean_to <- if (nrow(row_to) == 1L) row_to$mean_hours else NA_real_
  abs_change <- mean_to - mean_from
  tibble(
    transition = label,
    year_from = year_from,
    year_to = year_to,
    mean_hours_from = mean_from,
    mean_hours_to = mean_to,
    abs_change = abs_change,
    pct_change = ifelse(is.finite(mean_from) && mean_from != 0, 100 * abs_change / mean_from, NA_real_)
  )
}

#' Simple OLS trend of mean_hours ~ year, fit separately WITHIN one phase
#' (a year range, inclusive of both endpoints — boundary years shared between
#' adjacent phases are therefore counted in both, per the briefing's
#' provisional phase definition; flagged by the caller, not hidden here).
#' Reports slope, direction, R-squared and slope as a percent of the phase's
#' own mean level (for readability across phases with very different
#' absolute fishing-hours levels). No structural-break detection is
#' performed — see the run log for the strucchange follow-up flag.
fit_phase_trend <- function(year_summary, year_from, year_to, phase_label) {
  sub <- year_summary %>% filter(year >= year_from, year <= year_to)
  n_years <- nrow(sub)
  if (n_years < 2L) {
    return(tibble(
      phase = phase_label, year_from = year_from, year_to = year_to, n_years = n_years,
      slope_hours_per_year = NA_real_, direction = NA_character_, r_squared = NA_real_,
      phase_mean_hours = suppressWarnings(mean(sub$mean_hours)), slope_pct_of_mean_per_year = NA_real_
    ))
  }
  fit <- stats::lm(mean_hours ~ year, data = sub)
  slope_val <- unname(stats::coef(fit)["year"])
  r2 <- summary(fit)$r.squared
  phase_mean <- mean(sub$mean_hours)
  tibble(
    phase = phase_label,
    year_from = year_from,
    year_to = year_to,
    n_years = n_years,
    slope_hours_per_year = slope_val,
    direction = dplyr::case_when(
      slope_val > 0 ~ "increasing",
      slope_val < 0 ~ "decreasing",
      TRUE ~ "flat"
    ),
    r_squared = r2,
    phase_mean_hours = phase_mean,
    slope_pct_of_mean_per_year = ifelse(phase_mean != 0, 100 * slope_val / phase_mean, NA_real_)
  )
}

# ---------------------------------------------------------------------------
# Section B helpers — spatial heterogeneity by decade, persistence
# ---------------------------------------------------------------------------

#' Wide stat_rec x decade table of mean fishing hours (one column per
#' DECADE_BINS entry), from build_rect_decade_summary()-style long input.
pivot_rect_decade_wide <- function(rect_decade_long, value_col = "mean_val") {
  rect_decade_long %>%
    select(stat_rec, decade, all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = decade, values_from = all_of(value_col))
}

#' Per-rectangle absolute + percent change between two named decade columns.
#' Only rectangles with non-NA values in BOTH columns are included; n_used
#' reports how many that was (attached as an attribute, not a column, so the
#' table stays tidy for direct CSV export).
rect_decade_change <- function(rect_decade_wide, col_from, col_to, label) {
  stopifnot(all(c(col_from, col_to) %in% names(rect_decade_wide)))
  out <- rect_decade_wide %>%
    filter(!is.na(.data[[col_from]]), !is.na(.data[[col_to]])) %>%
    transmute(
      stat_rec,
      comparison = label,
      decade_from = col_from,
      decade_to = col_to,
      value_from = .data[[col_from]],
      value_to = .data[[col_to]],
      abs_change = .data[[col_to]] - .data[[col_from]],
      pct_change = ifelse(.data[[col_from]] != 0, 100 * (.data[[col_to]] - .data[[col_from]]) / .data[[col_from]], NA_real_)
    )
  attr(out, "n_used") <- nrow(out)
  out
}

#' Tercile classification of per-rectangle |abs_change| into
#' "stable" (bottom tercile) / "moderate" (middle) / "large_change" (top
#' tercile), via dplyr::ntile — a transparent, provisional threshold, NOT a
#' formal cutoff. Returns the input with a `change_class` column added and
#' the two tercile boundary values as attributes (for run-log reporting).
classify_change_tercile <- function(change_df, abs_change_col = "abs_change") {
  x <- abs(change_df[[abs_change_col]])
  tier <- dplyr::ntile(x, 3L)
  out <- change_df %>%
    mutate(
      abs_change_magnitude = x,
      change_tercile = tier,
      change_class = dplyr::case_when(
        tier == 1L ~ "stable",
        tier == 3L ~ "large_change",
        TRUE ~ "moderate"
      )
    )
  ord <- sort(x)
  n <- length(ord)
  boundary_low <- if (n > 0L) ord[max(1L, round(n / 3))] else NA_real_
  boundary_high <- if (n > 0L) ord[max(1L, round(2 * n / 3))] else NA_real_
  attr(out, "tercile_boundary_low") <- boundary_low
  attr(out, "tercile_boundary_high") <- boundary_high
  out
}

#' Pearson + Spearman correlation of rectangle-level fishing pressure between
#' two decade columns (spatial persistence: does the high/low pressure
#' pattern hold its shape over time, or reshuffle?). Only complete-pair
#' rectangles are used; n reported alongside r for transparency.
decade_persistence_correlation <- function(rect_decade_wide, col_a, col_b, label) {
  stopifnot(all(c(col_a, col_b) %in% names(rect_decade_wide)))
  sub <- rect_decade_wide %>% filter(!is.na(.data[[col_a]]), !is.na(.data[[col_b]]))
  n <- nrow(sub)
  if (n < 3L) {
    return(tibble(
      comparison = label, decade_a = col_a, decade_b = col_b, n_rectangles = n,
      pearson_r = NA_real_, spearman_rho = NA_real_
    ))
  }
  tibble(
    comparison = label,
    decade_a = col_a,
    decade_b = col_b,
    n_rectangles = n,
    pearson_r = stats::cor(sub[[col_a]], sub[[col_b]], method = "pearson"),
    spearman_rho = stats::cor(sub[[col_a]], sub[[col_b]], method = "spearman")
  )
}

# ---------------------------------------------------------------------------
# Section C helpers — ICC by spatial-unit definition (5 resolutions)
# ---------------------------------------------------------------------------

#' Assign each rectangle-year (or rectangle-year-mean) row a spatial unit_id
#' via unit_map (stat_rec -> unit_id), then run the existing one-way
#' variance-components decomposition (variance_components_anova(), from
#' h3_pre_exploration_helpers.R) grouped by unit_id instead of raw stat_rec.
#' This is the whole "extension" — same estimator, coarser group variable.
icc_by_unit <- function(value_df, value_col, unit_map) {
  stopifnot(value_col %in% names(value_df), all(c("stat_rec", "unit_id") %in% names(unit_map)))
  joined <- value_df %>% inner_join(unit_map, by = "stat_rec")
  variance_components_anova(joined[[value_col]], joined$unit_id)
}

#' Full comparison row for one spatial-unit definition: n spatial units
#' (both "constructed" and "with data", which can differ once input data is
#' restricted to the relaxed rectangle universe), ICC for all three metrics,
#' and Spearman-Brown required-n at each reliability target, for each metric.
#' fishing_df/resid_df are rectangle-year(-mean) long tables already
#' restricted to the shared relaxed rectangle universe (see run script);
#' unit_map may span a larger construction universe (e.g. Scheme A's 197) —
#' natural attrition when joined is expected and reported, not corrected.
icc_comparison_row <- function(resolution_label, resolution_order,
                               n_units_constructed, unit_map,
                               fishing_df, resid_df, reliability_targets) {
  vc_fish <- icc_by_unit(fishing_df, "hours_total", unit_map)
  vc_mag <- icc_by_unit(resid_df, "mean_mag", unit_map)
  vc_signed <- icc_by_unit(resid_df, "mean_signed", unit_map)

  req_fish <- required_n_for_reliability(vc_fish$icc, reliability_targets)
  req_mag <- required_n_for_reliability(vc_mag$icc, reliability_targets)
  req_signed <- required_n_for_reliability(vc_signed$icc, reliability_targets)

  out <- tibble(
    resolution_order = resolution_order,
    resolution = resolution_label,
    n_units_constructed = n_units_constructed,
    n_units_with_fishing_data = vc_fish$n_groups,
    n_units_with_residual_data = vc_mag$n_groups,
    icc_fishing_pressure = vc_fish$icc,
    icc_residual_magnitude = vc_mag$icc,
    icc_residual_signed = vc_signed$icc
  )
  for (i in seq_along(reliability_targets)) {
    out[[sprintf("required_n_fishing_%.1f", reliability_targets[i])]] <- req_fish[i]
    out[[sprintf("required_n_residual_magnitude_%.1f", reliability_targets[i])]] <- req_mag[i]
    out[[sprintf("required_n_residual_signed_%.1f", reliability_targets[i])]] <- req_signed[i]
  }
  out
}

# ---------------------------------------------------------------------------
# Section D helpers — haul-effort unevenness vs fishing pressure
# ---------------------------------------------------------------------------

#' Gini coefficient of a non-negative numeric vector (standard sorted-values
#' formula: G = (2*sum(i*x_i))/(n*sum(x_i)) - (n+1)/n, i = 1..n after
#' ascending sort). NA if fewer than 2 finite values or the vector sums to 0.
gini_coefficient <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L || sum(x) == 0) {
    return(NA_real_)
  }
  x <- sort(x)
  idx <- seq_len(n)
  (2 * sum(idx * x)) / (n * sum(x)) - (n + 1) / n
}

#' Coefficient of variation (%) of a numeric vector: 100 * sd / mean.
cv_percent <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || mean(x) == 0) {
    return(NA_real_)
  }
  100 * sd(x) / mean(x)
}

#' One row of haul-effort-unevenness summary stats for one scheme/parameter
#' combination x period: CV and Gini of haul count across zones, plus
#' Pearson/Spearman correlation of zone mean fishing pressure with zone haul
#' count, plus the max/min haul-count ratio (two variants: strict, which can
#' be Inf if any zone has zero hauls that period, and "excluding zero-haul
#' zones", both reported so neither hides the other).
haul_unevenness_row <- function(scheme_label, param_label, period_label,
                                hauls_wide_col, fish_period_df) {
  n_hauls <- hauls_wide_col
  cv_val <- cv_percent(n_hauls)
  gini_val <- gini_coefficient(n_hauls)

  nonzero <- n_hauls[n_hauls > 0 & is.finite(n_hauls)]
  ratio_strict <- if (length(n_hauls) > 0 && min(n_hauls, na.rm = TRUE) > 0) {
    max(n_hauls, na.rm = TRUE) / min(n_hauls, na.rm = TRUE)
  } else {
    Inf
  }
  ratio_excl_zero <- if (length(nonzero) >= 2L) {
    max(nonzero) / min(nonzero)
  } else {
    NA_real_
  }

  join_df <- fish_period_df
  n_units_for_cor <- nrow(join_df)
  pearson_r <- NA_real_
  spearman_rho <- NA_real_
  if (n_units_for_cor >= 3L) {
    pearson_r <- suppressWarnings(stats::cor(join_df$mean_fishing_hours, join_df$n_hauls, method = "pearson"))
    spearman_rho <- suppressWarnings(stats::cor(join_df$mean_fishing_hours, join_df$n_hauls, method = "spearman"))
  }

  tibble(
    scheme = scheme_label,
    parameter = param_label,
    period = period_label,
    n_units = length(n_hauls),
    cv_pct_haul_count = cv_val,
    gini_haul_count = gini_val,
    max_min_ratio_strict = ratio_strict,
    max_min_ratio_excl_zero_units = ratio_excl_zero,
    n_units_zero_hauls = sum(n_hauls == 0, na.rm = TRUE),
    n_units_used_for_correlation = n_units_for_cor,
    pearson_r_fishing_vs_haulcount = pearson_r,
    spearman_rho_fishing_vs_haulcount = spearman_rho
  )
}
