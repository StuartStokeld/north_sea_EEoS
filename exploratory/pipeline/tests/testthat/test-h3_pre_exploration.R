source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h3_pre_exploration_helpers.R", local = TRUE)

test_that("build_rect_year_haul_counts counts per rectangle-year and normalises stat_rec", {
  haul_full <- data.frame(
    stat_rec = c('"38F5"', "38F5", "39F3", "38F5"),
    year = c(1990L, 1990L, 1990L, 1984L),
    stringsAsFactors = FALSE
  )
  out <- build_rect_year_haul_counts(haul_full, year_min = 1985L, year_max = 2015L)
  expect_equal(nrow(out), 2L)
  expect_equal(out$n_hauls[out$stat_rec == "38F5" & out$year == 1990], 2L)
  expect_false(1984L %in% out$year)
})

test_that("complete_rect_year_grid zero-fills missing rectangle-year combinations", {
  rect_year_hauls <- data.frame(stat_rec = "38F5", year = 1985L, n_hauls = 3L, stringsAsFactors = FALSE)
  grid <- complete_rect_year_grid(rect_year_hauls, year_min = 1985L, year_max = 1987L)
  expect_equal(nrow(grid), 3L)
  expect_equal(grid$n_hauls[grid$year == 1986L], 0L)
})

test_that("build_rectangle_usability_flags applies sparse/min-years thresholds correctly", {
  rect_year_hauls <- data.frame(
    stat_rec = c(rep("38F5", 12), rep("39F3", 2)),
    year = c(1985:1996, 1985:1986),
    n_hauls = c(rep(5L, 12), 1L, 1L),
    stringsAsFactors = FALSE
  )
  flags <- build_rectangle_usability_flags(rect_year_hauls, sparse_threshold = 5L, min_years = 10L)
  expect_true(flags$usable_temporal[flags$stat_rec == "38F5"])
  expect_false(flags$usable_temporal[flags$stat_rec == "39F3"])
  expect_equal(flags$n_years_qualifying[flags$stat_rec == "38F5"], 12L)
  expect_equal(flags$n_years_sparse[flags$stat_rec == "39F3"], 2L)
})

test_that("add_couce_coverage_flag flags rectangles absent from Couce data", {
  rect_flags <- data.frame(stat_rec = c("38F5", "39F3"), usable_temporal = c(TRUE, TRUE), stringsAsFactors = FALSE)
  couce_year <- data.frame(stat_rec = "38F5", year = 1985L, hours_total = 100, stringsAsFactors = FALSE)
  out <- add_couce_coverage_flag(rect_flags, couce_year)
  expect_true(out$has_couce_coverage[out$stat_rec == "38F5"])
  expect_false(out$has_couce_coverage[out$stat_rec == "39F3"])
  expect_false(out$usable_for_fishing_analysis[out$stat_rec == "39F3"])
})

test_that("add_resid_signed is the negative of the pipeline's primary residual convention", {
  haul_eeos <- data.frame(
    stat_rec = "38F5", year = 1990L,
    ln_B_pred = 10, ln_B_obs = 9,
    residual = -1, abs_residual = 1,
    stringsAsFactors = FALSE
  )
  out <- add_resid_signed(haul_eeos)
  expect_equal(out$resid_signed, 1)
  expect_equal(out$resid_signed, -out$residual)
  expect_equal(out$resid_magnitude, out$abs_residual)
})

test_that("summarise_year_stat computes mean+CI and median+IQR correctly", {
  df <- data.frame(year = rep(1990L, 5), x = c(1, 2, 3, 4, 5))
  ci_out <- summarise_year_stat(df, "x", 1985L, 2015L, center = "mean", spread = "ci")
  expect_equal(ci_out$center, 3)
  expect_equal(ci_out$n, 5L)

  iqr_out <- summarise_year_stat(df, "x", 1985L, 2015L, center = "median", spread = "iqr")
  expect_equal(iqr_out$center, 3)
  expect_equal(iqr_out$lo, 2)
  expect_equal(iqr_out$hi, 4)
})

test_that("assign_decade_bin labels years correctly and returns NA outside bins", {
  bins <- list(c(1985L, 1994L), c(1995L, 2004L), c(2005L, 2015L))
  out <- assign_decade_bin(c(1990L, 2000L, 2010L, 1980L), bins)
  expect_equal(as.character(out), c("1985-1994", "1995-2004", "2005-2015", NA))
})

test_that("fit_rectangle_slopes returns NA slope for rectangles below min_points", {
  rect_year_values <- data.frame(
    stat_rec = c("38F5", "38F5", "38F5", "39F3"),
    year = c(1985L, 1990L, 1995L, 1985L),
    mean_val = c(1, 2, 3, 5),
    stringsAsFactors = FALSE
  )
  out <- fit_rectangle_slopes(rect_year_values, c("38F5", "39F3"), "mean_val", min_points = 2L)
  expect_true(out$slope[out$stat_rec == "38F5"] > 0)
  expect_true(is.na(out$slope[out$stat_rec == "39F3"]))
})

test_that("variance_components_anova recovers zero between-group variance for identical group means", {
  value <- c(1, 1, 1, 1, 1, 1)
  group <- rep(c("A", "B"), each = 3)
  vc <- variance_components_anova(value, group)
  expect_equal(vc$var_between, 0)
  expect_equal(vc$var_within, 0)
})

test_that("variance_components_anova gives high ICC when between-group variance dominates", {
  set.seed(1)
  group <- rep(c("A", "B", "C"), each = 50)
  group_means <- c(A = 0, B = 20, C = 40)
  value <- group_means[group] + rnorm(150, sd = 0.1)
  vc <- variance_components_anova(value, group)
  expect_gt(vc$icc, 0.99)
})

test_that("variance_components_anova gives low ICC when within-group variance dominates", {
  set.seed(1)
  group <- rep(c("A", "B", "C"), each = 50)
  value <- rnorm(150, mean = 0, sd = 10)
  vc <- variance_components_anova(value, group)
  expect_lt(vc$icc, 0.2)
})

test_that("variance_components_anova returns NA fields for fewer than 2 groups", {
  vc <- variance_components_anova(c(1, 2, 3), rep("A", 3))
  expect_true(is.na(vc$icc))
  expect_equal(vc$n_groups, 1L)
})
