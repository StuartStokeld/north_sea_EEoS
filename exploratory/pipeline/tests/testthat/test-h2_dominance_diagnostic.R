source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h2_dominance_diagnostic_helpers.R", local = TRUE)

# --- Synthetic haul-level dominance fixture --------------------------------
# 3 rectangles, "38F5" has 6 hauls (>= provisional threshold of 5), "39F3" has
# 2 hauls (below threshold), "40F2" has 5 hauls across 4 years (temporal case).
haul_dom <- data.frame(
  stat_rec = c(
    rep("38F5", 6),
    rep("39F3", 2),
    rep("40F2", 5)
  ),
  year = c(
    1985, 1986, 1987, 1988, 1989, 1990,
    1990, 1991,
    2000, 2001, 2002, 2003, 2003
  ),
  D = c(
    0.5, 0.6, 0.55, 0.62, 0.58, 0.6,
    0.4, 0.42,
    0.3, 0.5, 0.7, 0.9, 0.85
  ),
  size_CV = c(
    0.2, 0.25, 0.22, 0.28, 0.24, 0.26,
    0.5, 0.52,
    0.1, 0.2, 0.4, 0.6, 0.55
  ),
  stringsAsFactors = FALSE
)

test_that("build_rectangle_dominance_panel aggregates unweighted mean-of-hauls per rectangle", {
  panel <- build_rectangle_dominance_panel(haul_dom, year_min = 1985L, year_max = 2015L)
  expect_equal(nrow(panel), 3L)
  row_38F5 <- panel[panel$stat_rec == "38F5", ]
  expect_equal(row_38F5$n_hauls, 6L)
  expect_equal(row_38F5$mean_D, mean(c(0.5, 0.6, 0.55, 0.62, 0.58, 0.6)))
})

test_that("apply_min_hauls_threshold drops rectangles below the provisional threshold and reports the drop", {
  panel <- build_rectangle_dominance_panel(haul_dom, year_min = 1985L, year_max = 2015L)
  out <- apply_min_hauls_threshold(panel, min_hauls = 5L)
  expect_equal(out$n_before, 3L)
  expect_equal(out$n_after, 2L)
  expect_equal(out$n_dropped, 1L)
  expect_false("39F3" %in% out$panel$stat_rec)
})

test_that("cross_sectional_association reports r, R2, slope, p_value, n and matches lm/cor.test", {
  panel <- data.frame(
    mean_D = c(0.1, 0.2, 0.3, 0.4, 0.5),
    mean_annual_hours_total = c(10, 20, 30, 40, 50),
    stringsAsFactors = FALSE
  )
  out <- cross_sectional_association(panel, "mean_D")
  expect_equal(out$n, 5L)
  expect_equal(out$correlation, 1, tolerance = 1e-8)
  expect_equal(out$r_squared, 1, tolerance = 1e-8)
  expect_equal(out$slope, 0.01, tolerance = 1e-8)
})

test_that("cross_sectional_association returns NA row when n < 3", {
  panel <- data.frame(mean_D = c(0.1, 0.2), mean_annual_hours_total = c(10, 20))
  out <- cross_sectional_association(panel, "mean_D")
  expect_equal(out$n, 2L)
  expect_true(is.na(out$correlation))
})

test_that("build_annual_rectangle_dominance produces one row per rectangle x year", {
  annual <- build_annual_rectangle_dominance(haul_dom, year_min = 1985L, year_max = 2015L)
  n_years_40F2 <- annual[annual$stat_rec == "40F2", ]
  expect_equal(nrow(n_years_40F2), 4L) # 2003 has 2 hauls -> collapses to 1 row
  expect_equal(n_years_40F2$n_hauls[n_years_40F2$year == 2003], 2L)
})

test_that("apply_min_years_threshold keeps only rectangles with sufficient distinct years", {
  annual <- build_annual_rectangle_dominance(haul_dom, year_min = 1985L, year_max = 2015L)
  out <- apply_min_years_threshold(annual, min_years = 3L)
  expect_true("40F2" %in% out$kept_rectangles) # 4 distinct years
  expect_true("38F5" %in% out$kept_rectangles) # 6 distinct years
  expect_false("39F3" %in% out$kept_rectangles) # 2 distinct years
})

test_that("within_rectangle_temporal_correlation computes per-rectangle Pearson correlation vs joined fishing hours", {
  annual <- build_annual_rectangle_dominance(haul_dom, year_min = 1985L, year_max = 2015L)
  couce_year <- data.frame(
    stat_rec = rep("40F2", 4),
    year = c(2000, 2001, 2002, 2003),
    hours_total = c(100, 200, 300, 400),
    stringsAsFactors = FALSE
  )
  out <- within_rectangle_temporal_correlation(
    annual, couce_year, rectangles = "40F2", metric_col = "mean_D", min_paired_years = 3L
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_years_paired, 4L)
  expect_true(out$correlation > 0.9) # mean_D rises monotonically with year/hours in the fixture
})

test_that("within_rectangle_temporal_correlation returns NA correlation below min_paired_years", {
  annual <- build_annual_rectangle_dominance(haul_dom, year_min = 1985L, year_max = 2015L)
  couce_year <- data.frame(
    stat_rec = "39F3", year = c(1990, 1991), hours_total = c(50, 60),
    stringsAsFactors = FALSE
  )
  out <- within_rectangle_temporal_correlation(
    annual, couce_year, rectangles = "39F3", metric_col = "mean_D", min_paired_years = 3L
  )
  expect_equal(out$n_years_paired, 2L)
  expect_true(is.na(out$correlation))
})

test_that("summarise_within_rectangle_correlations reports median/IQR/pct positive-negative", {
  per_rect <- data.frame(
    stat_rec = c("a", "b", "c", "d"),
    n_years_paired = c(4L, 4L, 4L, 2L),
    correlation = c(0.8, 0.6, -0.2, NA_real_),
    p_value = c(0.01, 0.04, 0.5, NA_real_),
    stringsAsFactors = FALSE
  )
  out <- summarise_within_rectangle_correlations(per_rect)
  expect_equal(out$n_rectangles, 4L)
  expect_equal(out$n_with_valid_corr, 3L)
  expect_equal(out$median_corr, 0.6)
  expect_equal(out$pct_positive, round(100 * 2 / 3, 1))
  expect_equal(out$pct_negative, round(100 * 1 / 3, 1))
})

test_that("build_step0_summary_table assembles the 4-row spec table with correct column names", {
  cross_D <- tibble::tibble(variable = "mean_D", n = 10L, correlation = 0.5, r_squared = 0.25, slope = 1.2, p_value = 0.03)
  cross_cv <- tibble::tibble(variable = "mean_size_CV", n = 10L, correlation = -0.1, r_squared = 0.01, slope = -0.4, p_value = 0.6)
  temporal_D_summary <- list(
    n_rectangles = 8L, n_with_valid_corr = 6L, median_corr = 0.3,
    iqr_low = 0.1, iqr_high = 0.5, pct_positive = 66.7, pct_negative = 33.3, pct_p_lt_05 = 20
  )
  temporal_cv_summary <- list(
    n_rectangles = 8L, n_with_valid_corr = 0L, median_corr = NA_real_,
    iqr_low = NA_real_, iqr_high = NA_real_, pct_positive = NA_real_, pct_negative = NA_real_, pct_p_lt_05 = NA_real_
  )
  out <- build_step0_summary_table(cross_D, cross_cv, temporal_D_summary, temporal_cv_summary)
  expect_equal(nrow(out), 4L)
  expect_equal(names(out), c("test", "variable", "n", "correlation", "r_squared", "slope", "p_value", "notes"))
  expect_equal(out$test, c(
    "cross-sectional", "cross-sectional",
    "temporal (within-rectangle, summarised)", "temporal (within-rectangle, summarised)"
  ))
  expect_false(is.na(out$notes[3]))
  expect_true(grepl("no rectangle had", out$notes[4]))
})
