source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h2_dominance_diagnostic_helpers.R", local = TRUE)
source("../../R/h2_dominance_robustness_helpers.R", local = TRUE)

# --- Synthetic rectangle panel fixture -------------------------------------
# A near-linear D-vs-fishing_hours relationship, plus one deliberate
# high-leverage outlier rectangle ("hi_lev") to exercise Cook's distance.
set.seed(1)
n_normal <- 20L
panel <- data.frame(
  stat_rec = c(sprintf("R%02d", seq_len(n_normal)), "hi_lev"),
  mean_D = c(0.5 + 0.001 * seq_len(n_normal) + rnorm(n_normal, sd = 0.01), 0.95),
  mean_size_CV = c(0.4 + rnorm(n_normal, sd = 0.01), 0.41),
  mean_annual_hours_total = c(seq(100, 2000, length.out = n_normal), 15000),
  stringsAsFactors = FALSE
)

test_that("pearson_spearman_comparison reports both methods with matching n", {
  out <- pearson_spearman_comparison(panel, "mean_D")
  expect_equal(nrow(out), 2L)
  expect_setequal(out$method, c("pearson", "spearman"))
  expect_true(all(out$n == nrow(panel)))
  expect_true(all(is.finite(out$statistic)))
})

test_that("pearson_spearman_comparison returns NA row when n < 3", {
  small <- panel[1:2, ]
  out <- pearson_spearman_comparison(small, "mean_D")
  expect_true(all(is.na(out$statistic)))
})

test_that("cooks_distance_table flags the deliberate high-leverage rectangle at the top", {
  cooks_tbl <- cooks_distance_table(panel, "mean_D")
  expect_equal(nrow(cooks_tbl), nrow(panel))
  expect_equal(cooks_tbl$stat_rec[1], "hi_lev")
  expect_true(cooks_tbl$cooks_distance[1] > cooks_tbl$cooks_distance[2])
  expect_true(all(cooks_tbl$cooks_threshold_4_over_n == 4 / nrow(panel)))
})

test_that("refit_excluding_top_cooks drops the requested rectangles and recomputes association", {
  refit1 <- refit_excluding_top_cooks(panel, "mean_D", n_remove = 1L)
  expect_equal(refit1$n_removed, 1L)
  expect_equal(refit1$removed_stat_rec, "hi_lev")
  expect_equal(refit1$n, nrow(panel) - 1L)

  refit3 <- refit_excluding_top_cooks(panel, "mean_D", n_remove = 3L)
  expect_equal(refit3$n_removed, 3L)
  expect_equal(refit3$n, nrow(panel) - 3L)
})

test_that("linear_vs_loess_deviation reports a max deviation and where it occurs, within the observed range", {
  out <- linear_vs_loess_deviation(panel, "mean_D", n_points = 50L)
  expect_equal(nrow(out), 1L)
  expect_true(out$max_abs_diff >= 0)
  expect_true(out$at_fishing_hours >= min(panel$mean_annual_hours_total))
  expect_true(out$at_fishing_hours <= max(panel$mean_annual_hours_total))
})

test_that("missingness_comparison reports group means and both t-test and Wilcoxon statistics", {
  dropped <- data.frame(mean_D = c(0.4, 0.42, 0.45, 0.5, 0.48))
  retained <- data.frame(mean_D = c(0.6, 0.62, 0.65, 0.58, 0.63, 0.61))
  out <- missingness_comparison(dropped, retained, "mean_D")
  expect_equal(out$n_dropped, 5L)
  expect_equal(out$n_retained, 6L)
  expect_equal(out$mean_dropped, mean(dropped$mean_D))
  expect_equal(out$mean_retained, mean(retained$mean_D))
  expect_true(is.finite(out$t_p_value))
  expect_true(is.finite(out$wilcox_p_value))
  expect_true(out$t_p_value < 0.05) # groups are clearly separated by construction
})

test_that("missingness_comparison handles empty dropped group gracefully", {
  dropped <- data.frame(mean_D = numeric(0))
  retained <- data.frame(mean_D = c(0.6, 0.62, 0.65))
  out <- missingness_comparison(dropped, retained, "mean_D")
  expect_equal(out$n_dropped, 0L)
  expect_true(is.na(out$t_p_value))
  expect_true(is.na(out$wilcox_p_value))
})
