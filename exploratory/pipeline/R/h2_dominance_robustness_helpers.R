# Step 0 ROBUSTNESS CHECK helpers — is the D-vs-fishing-pressure cross-sectional
# correlation a stable signal or an artefact of a few high-leverage rectangles
# / biased missingness? Does NOT modify the Step 0 pipeline, thresholds, or
# outputs (h2_dominance_diagnostic_helpers.R, run_h2_dominance_fishing_pressure_
# diagnostic.R, step0_*.csv/.rds). These are additive checks only.
#
# Requires h2_dominance_diagnostic_helpers.R (reuses cross_sectional_association()
# unchanged, so the Pearson/OLS re-fits below use exactly the Step 0 convention:
# metric ~ fishing_hours, Pearson r symmetric).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Check 1: Pearson vs Spearman
# ---------------------------------------------------------------------------

#' Pearson and Spearman correlation for one dominance metric vs fishing hours,
#' on the same panel/rows used for the original Step 0 Pearson result.
pearson_spearman_comparison <- function(panel, x_col, y_col = "mean_annual_hours_total") {
  df <- panel %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  n <- nrow(df)
  if (n < 3L) {
    return(tibble(
      variable = x_col, method = c("pearson", "spearman"),
      statistic = NA_real_, p_value = NA_real_, n = n
    ))
  }
  pearson <- suppressWarnings(stats::cor.test(df[[y_col]], df[[x_col]], method = "pearson"))
  spearman <- suppressWarnings(
    stats::cor.test(df[[y_col]], df[[x_col]], method = "spearman", exact = FALSE)
  )
  tibble(
    variable = x_col,
    method = c("pearson", "spearman"),
    statistic = c(unname(pearson$estimate), unname(spearman$estimate)),
    p_value = c(pearson$p.value, spearman$p.value),
    n = n
  )
}

# ---------------------------------------------------------------------------
# Check 2: leverage / influence diagnostics
# ---------------------------------------------------------------------------

#' Cook's distance for every rectangle in the metric ~ fishing_hours OLS fit
#' (same regression direction as cross_sectional_association()), sorted
#' descending so the top-N is just head().
cooks_distance_table <- function(panel, x_col, y_col = "mean_annual_hours_total", id_col = "stat_rec") {
  df <- panel %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  fit <- stats::lm(stats::reformulate(y_col, response = x_col), data = df)
  df %>%
    mutate(cooks_distance = stats::cooks.distance(fit), cooks_threshold_4_over_n = 4 / n()) %>%
    arrange(desc(cooks_distance)) %>%
    select(all_of(id_col), all_of(x_col), all_of(y_col), cooks_distance, cooks_threshold_4_over_n)
}

#' Recompute the Step 0 Pearson/OLS association after dropping the top-N
#' rectangles by Cook's distance (reuses cross_sectional_association()
#' unchanged; source h2_dominance_diagnostic_helpers.R before this file).
refit_excluding_top_cooks <- function(panel, x_col, y_col = "mean_annual_hours_total", n_remove) {
  cooks_tbl <- cooks_distance_table(panel, x_col, y_col)
  drop_ids <- cooks_tbl$stat_rec[seq_len(min(n_remove, nrow(cooks_tbl)))]
  trimmed <- panel %>% filter(!.data$stat_rec %in% drop_ids)
  out <- cross_sectional_association(trimmed, x_col, y_col)
  out$n_removed <- length(drop_ids)
  out$removed_stat_rec <- paste(drop_ids, collapse = ";")
  out
}

# ---------------------------------------------------------------------------
# Check 3: linear vs loess fit comparison
# ---------------------------------------------------------------------------

#' Maximum absolute deviation between the linear fit and a loess fit
#' (default span) across the observed range of fishing_hours, and where it
#' occurs. Formalises the visual linear-vs-loess overlay already shown on
#' the Step 0 scatterplot.
linear_vs_loess_deviation <- function(panel, x_col, y_col = "mean_annual_hours_total",
                                      span = 0.75, n_points = 200L) {
  df <- panel %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  fit_lm <- stats::lm(stats::reformulate(y_col, response = x_col), data = df)
  fit_loess <- stats::loess(stats::reformulate(y_col, response = x_col), data = df, span = span)

  grid <- seq(min(df[[y_col]]), max(df[[y_col]]), length.out = n_points)
  grid_df <- setNames(data.frame(grid), y_col)
  pred_lm <- as.numeric(stats::predict(fit_lm, newdata = grid_df))
  pred_loess <- as.numeric(stats::predict(fit_loess, newdata = grid_df))

  diffs <- abs(pred_lm - pred_loess)
  idx <- which.max(diffs)
  tibble(
    variable = x_col,
    max_abs_diff = diffs[idx],
    at_fishing_hours = grid[idx],
    linear_fitted_at_max = pred_lm[idx],
    loess_fitted_at_max = pred_loess[idx],
    loess_span = span,
    n_grid_points = n_points,
    n = nrow(df)
  )
}

# ---------------------------------------------------------------------------
# Check 4: missingness check on rectangles dropped for no Couce coverage
# ---------------------------------------------------------------------------

#' Compare a rectangle-level metric (already aggregated from the same
#' haul-level source as Step 0, e.g. mean_D / mean_size_CV in
#' step0_rectangle_panel.csv) between the dropped (no Couce coverage) and
#' retained rectangle groups. Reports both a t-test and a Wilcoxon test —
#' n is expected to be small (26 dropped) and possibly non-normal, so the
#' choice between them is left open rather than prejudged.
missingness_comparison <- function(dropped, retained, metric_col) {
  x <- dropped[[metric_col]]
  y <- retained[[metric_col]]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  t_test <- if (length(x) >= 2L && length(y) >= 2L) {
    suppressWarnings(stats::t.test(x, y))
  } else {
    NULL
  }
  wilcox_test <- if (length(x) >= 1L && length(y) >= 1L) {
    suppressWarnings(stats::wilcox.test(x, y))
  } else {
    NULL
  }

  tibble(
    variable = metric_col,
    n_dropped = length(x),
    n_retained = length(y),
    mean_dropped = mean(x),
    mean_retained = mean(y),
    t_statistic = if (!is.null(t_test)) unname(t_test$statistic) else NA_real_,
    t_p_value = if (!is.null(t_test)) t_test$p.value else NA_real_,
    wilcox_statistic = if (!is.null(wilcox_test)) unname(wilcox_test$statistic) else NA_real_,
    wilcox_p_value = if (!is.null(wilcox_test)) wilcox_test$p.value else NA_real_
  )
}
