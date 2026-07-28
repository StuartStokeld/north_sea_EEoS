# Step 0 diagnostic: is D (Berger-Parker dominance) / size_CV associated with
# Couce et al. (2020) fishing pressure, independent of the primary H2
# residual? Load-bearing diagnostic, not a robustness check — see
# CURSOR_BRIEFING referenced in run_h2_dominance_fishing_pressure_diagnostic.R.
#
# Requires h2_common.R (normalize_stat_rec, H2_YEAR_MIN, H2_YEAR_MAX).
#
# Two thresholds used throughout this diagnostic are PROVISIONAL — the final
# H2/H3 haul-inclusion threshold is pending supervisor review, and the H3
# "sufficient temporal coverage" definition has not been finalised either.
# Both are kept as named constants (not hard-coded inline) so this diagnostic
# is easy to re-run once either is decided.
STEP0_MIN_HAULS_PROVISIONAL <- 5L
STEP0_MIN_YEARS_PROVISIONAL <- 3L

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Step 1: cross-sectional rectangle aggregation (H2 case)
# ---------------------------------------------------------------------------

#' Unweighted mean-of-hauls rectangle aggregation of D and size_CV, pooled
#' across the full year_min:year_max period.
#'
#' Mirrors build_h2_rectangle_residuals()'s mean-of-hauls approach exactly —
#' no alternative weighting scheme is introduced here, per the briefing.
build_rectangle_dominance_panel <- function(haul_dom,
                                            year_min = H2_YEAR_MIN,
                                            year_max = H2_YEAR_MAX) {
  required <- c("stat_rec", "year", "D", "size_CV")
  missing <- setdiff(required, names(haul_dom))
  if (length(missing) > 0L) {
    stop("haul_dom missing columns: ", paste(missing, collapse = ", "))
  }

  haul_dom %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(
      !is.na(stat_rec), stat_rec != "",
      year >= year_min, year <= year_max,
      is.finite(D), is.finite(size_CV)
    ) %>%
    group_by(stat_rec) %>%
    summarise(
      n_hauls = n(),
      mean_D = mean(D, na.rm = TRUE),
      mean_size_CV = mean(size_CV, na.rm = TRUE),
      year_min_seen = min(year, na.rm = TRUE),
      year_max_seen = max(year, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Apply the provisional haul-count threshold, reporting the drop explicitly
#' rather than dropping rectangles silently.
apply_min_hauls_threshold <- function(rectangle_panel, min_hauls = STEP0_MIN_HAULS_PROVISIONAL) {
  n_before <- nrow(rectangle_panel)
  kept <- rectangle_panel %>% filter(n_hauls >= min_hauls)
  list(
    panel = kept,
    n_before = n_before,
    n_after = nrow(kept),
    n_dropped = n_before - nrow(kept),
    min_hauls = min_hauls
  )
}

# ---------------------------------------------------------------------------
# Step 2: cross-sectional association test (H2 case)
# ---------------------------------------------------------------------------

#' Pearson correlation + simple OLS (metric ~ fishing_hours) for one
#' dominance metric against rectangle-mean fishing hours.
#'
#' Regression direction follows the briefing spec exactly:
#' `mean_D ~ fishing_hours` / `mean_size_CV ~ fishing_hours` (metric is the
#' response). Pearson r is symmetric so direction doesn't affect it.
cross_sectional_association <- function(panel, x_col, y_col = "mean_annual_hours_total") {
  df <- panel %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  n <- nrow(df)
  if (n < 3L) {
    return(tibble(
      variable = x_col, n = n, correlation = NA_real_, r_squared = NA_real_,
      slope = NA_real_, p_value = NA_real_
    ))
  }
  test <- suppressWarnings(stats::cor.test(df[[y_col]], df[[x_col]], method = "pearson"))
  fit <- stats::lm(stats::reformulate(y_col, response = x_col), data = df)
  tibble(
    variable = x_col,
    n = n,
    correlation = unname(test$estimate),
    r_squared = summary(fit)$r.squared,
    slope = unname(stats::coef(fit)[2L]),
    p_value = unname(test$p.value)
  )
}

# ---------------------------------------------------------------------------
# Step 3: temporal within-rectangle association test (H3 case)
# ---------------------------------------------------------------------------

#' Annual (stat_rec x year) mean D / mean size_CV, unweighted mean-of-hauls
#' (same convention as the cross-sectional aggregation above).
build_annual_rectangle_dominance <- function(haul_dom, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX) {
  required <- c("stat_rec", "year", "D", "size_CV")
  missing <- setdiff(required, names(haul_dom))
  if (length(missing) > 0L) {
    stop("haul_dom missing columns: ", paste(missing, collapse = ", "))
  }

  haul_dom %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(
      !is.na(stat_rec), stat_rec != "",
      year >= year_min, year <= year_max,
      is.finite(D), is.finite(size_CV)
    ) %>%
    group_by(stat_rec, year) %>%
    summarise(
      n_hauls = n(),
      mean_D = mean(D, na.rm = TRUE),
      mean_size_CV = mean(size_CV, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Rectangles with hauls present in >= min_years distinct years — provisional
#' "sufficient temporal coverage" definition for the H3 case.
apply_min_years_threshold <- function(annual_panel, min_years = STEP0_MIN_YEARS_PROVISIONAL) {
  coverage <- annual_panel %>%
    group_by(stat_rec) %>%
    summarise(n_years = dplyr::n_distinct(year), .groups = "drop")
  kept <- coverage %>% filter(n_years >= min_years)
  list(
    kept_rectangles = kept$stat_rec,
    coverage = coverage,
    n_before = nrow(coverage),
    n_after = nrow(kept),
    n_dropped = nrow(coverage) - nrow(kept),
    min_years = min_years
  )
}

#' Per-rectangle within-rectangle correlation of an annual dominance metric
#' vs annual fishing hours (Couce year-level effort joined on stat_rec x year).
#'
#' Only rectangles in `rectangles` are considered. Within each, only years
#' with BOTH a dominance-metric observation and a Couce fishing-hours
#' observation are used; `n_years_paired` is reported per rectangle so any
#' loss from the join (missing Couce coverage in a given year) is visible
#' rather than silently reducing the effective sample.
within_rectangle_temporal_correlation <- function(annual_panel, couce_year, rectangles,
                                                   metric_col, min_paired_years = 3L) {
  joined <- annual_panel %>%
    filter(stat_rec %in% rectangles) %>%
    inner_join(couce_year %>% select(stat_rec, year, hours_total), by = c("stat_rec", "year")) %>%
    filter(is.finite(.data[[metric_col]]), is.finite(hours_total))

  if (nrow(joined) == 0L) {
    return(tibble(
      stat_rec = character(), n_years_paired = integer(),
      correlation = numeric(), p_value = numeric()
    ))
  }

  joined %>%
    group_by(stat_rec) %>%
    group_modify(function(chunk, key) {
      n_paired <- nrow(chunk)
      if (n_paired < min_paired_years) {
        return(tibble(n_years_paired = n_paired, correlation = NA_real_, p_value = NA_real_))
      }
      test <- suppressWarnings(
        stats::cor.test(chunk$hours_total, chunk[[metric_col]], method = "pearson")
      )
      tibble(
        n_years_paired = n_paired,
        correlation = unname(test$estimate),
        p_value = unname(test$p.value)
      )
    }) %>%
    ungroup()
}

#' Summarise the distribution of per-rectangle within-rectangle correlations
#' (median / IQR / % positive vs negative), rather than pooling across
#' rectangles — per the briefing, pooling would need a mixed-model structure
#' that hasn't been decided yet.
summarise_within_rectangle_correlations <- function(per_rectangle) {
  valid <- per_rectangle %>% filter(is.finite(correlation))
  n_valid <- nrow(valid)
  if (n_valid == 0L) {
    return(list(
      n_rectangles = nrow(per_rectangle), n_with_valid_corr = 0L,
      median_corr = NA_real_, iqr_low = NA_real_, iqr_high = NA_real_,
      pct_positive = NA_real_, pct_negative = NA_real_, pct_p_lt_05 = NA_real_
    ))
  }
  list(
    n_rectangles = nrow(per_rectangle),
    n_with_valid_corr = n_valid,
    median_corr = median(valid$correlation),
    iqr_low = unname(quantile(valid$correlation, 0.25, na.rm = TRUE)),
    iqr_high = unname(quantile(valid$correlation, 0.75, na.rm = TRUE)),
    pct_positive = round(100 * mean(valid$correlation > 0), 1),
    pct_negative = round(100 * mean(valid$correlation < 0), 1),
    pct_p_lt_05 = round(100 * mean(valid$p_value < 0.05, na.rm = TRUE), 1)
  )
}

# ---------------------------------------------------------------------------
# Step 4: assemble the single summary table (deliverable 1)
# ---------------------------------------------------------------------------

#' Assemble the Step 0 summary table (one row per test) per the briefing spec:
#' test | variable | n | correlation | R2 | slope | p_value | notes
build_step0_summary_table <- function(cross_D, cross_cv, temporal_D_summary, temporal_cv_summary) {
  temporal_note <- function(s) {
    if (s$n_with_valid_corr == 0L) {
      return(sprintf(
        "median/IQR of per-rectangle correlations: no rectangle had >= %d paired years (n_rectangles=%d)",
        STEP0_MIN_YEARS_PROVISIONAL, s$n_rectangles
      ))
    }
    sprintf(
      paste0(
        "median/IQR of per-rectangle correlations: median r=%.3f, IQR [%.3f, %.3f]; ",
        "%s%% positive / %s%% negative (n=%d of %d rectangles with >=3 paired years; ",
        "%s%% of per-rectangle tests p<0.05); provisional temporal-coverage threshold = %d years; ",
        "policy-change dates not yet finalised, so before/after split was not computed."
      ),
      s$median_corr, s$iqr_low, s$iqr_high, s$pct_positive, s$pct_negative,
      s$n_with_valid_corr, s$n_rectangles, s$pct_p_lt_05, STEP0_MIN_YEARS_PROVISIONAL
    )
  }

  bind_rows(
    tibble(
      test = "cross-sectional", variable = "D vs fishing_hours",
      n = cross_D$n, correlation = cross_D$correlation, r_squared = cross_D$r_squared,
      slope = cross_D$slope, p_value = cross_D$p_value,
      notes = NA_character_
    ),
    tibble(
      test = "cross-sectional", variable = "size_CV vs fishing_hours",
      n = cross_cv$n, correlation = cross_cv$correlation, r_squared = cross_cv$r_squared,
      slope = cross_cv$slope, p_value = cross_cv$p_value,
      notes = NA_character_
    ),
    tibble(
      test = "temporal (within-rectangle, summarised)", variable = "D vs fishing_hours",
      n = temporal_D_summary$n_with_valid_corr,
      correlation = temporal_D_summary$median_corr,
      r_squared = temporal_D_summary$median_corr^2,
      slope = NA_real_, p_value = NA_real_,
      notes = temporal_note(temporal_D_summary)
    ),
    tibble(
      test = "temporal (within-rectangle, summarised)", variable = "size_CV vs fishing_hours",
      n = temporal_cv_summary$n_with_valid_corr,
      correlation = temporal_cv_summary$median_corr,
      r_squared = temporal_cv_summary$median_corr^2,
      slope = NA_real_, p_value = NA_real_,
      notes = temporal_note(temporal_cv_summary)
    )
  )
}
