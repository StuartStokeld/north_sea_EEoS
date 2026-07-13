source("../../R/h1_common.R", local = TRUE)

test_that("log_r2 is 1 for perfect log prediction", {
  x <- log(c(100, 1000, 10000))
  expect_equal(log_r2(x, x), 1)
})

test_that("log_r2 is negative for constant multiplicative error", {
  x <- log(c(100, 1000, 10000))
  y <- x + log(4)
  expect_lt(log_r2(x, y), 0)
})

test_that("log_cor2 is 1 for scaled log prediction", {
  x <- log(c(100, 1000, 10000))
  y <- x + log(4)
  expect_equal(log_cor2(x, y), 1)
})

test_that("log_rmse matches constant offset", {
  x <- log(c(100, 1000, 10000))
  y <- x + log(4)
  expect_equal(log_rmse(x, y), log(4), tolerance = 1e-10)
})

test_that("productivity_ratio scales with B", {
  E <- 1000
  B <- c(100, 1000)
  r <- productivity_ratio(E, B)
  expect_equal(r[1] / r[2], (1000 / 100)^(3 / 4), tolerance = 1e-10)
})

test_that("harte_fig2_metrics Pearson R² matches cor² on raw ratios", {
  x <- c(1, 2, 4, 8)
  y <- c(1.1, 2.2, 3.5, 7)
  m <- harte_fig2_metrics(x, y)
  expect_equal(m$fig2_r2_pearson, cor(x, y)^2, tolerance = 1e-10)
  expect_equal(m$fig2_r2_pearson_all, m$fig2_r2_pearson)
  expect_equal(m$cor2, m$fig2_r2_pearson)
  expect_false(is.na(m$fig2_r2_cod_extended))
  expect_equal(m$fig2_r2_cod_extended, m$log_r2)
})

test_that("harte_fig2_metrics trimmed R² drops after excluding high-leverage points", {
  # Bulk near 1:1; two extreme points inflate Pearson R² on raw axes
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 700, 1000)
  y <- c(1.2, 1.8, 3.5, 3.2, 6, 5, 8, 7, 3400, 5300)
  m_all <- harte_fig2_metrics(x, y, n_exclude = 0L)
  m_trim <- harte_fig2_metrics(x, y, n_exclude = 2L)
  expect_equal(m_all$fig2_r2_pearson_all, cor(x, y)^2, tolerance = 1e-10)
  expect_equal(m_trim$fig2_r2_pearson_n_excluded, 2L)
  expect_lt(m_trim$fig2_r2_pearson_trimmed, m_trim$fig2_r2_pearson_all)
  expect_equal(
    m_trim$fig2_r2_pearson_trimmed,
    cor(x[1:8], y[1:8])^2,
    tolerance = 1e-10
  )
})

test_that("evaluate_prediction returns expected columns", {
  obs <- c(100, 1000, 10000)
  pred <- obs * 4
  m <- evaluate_prediction(obs, pred)
  expect_named(
    m,
    c(
      "log_r2", "cor2", "log_rmse", "raw_rmse", "median_ratio",
      "median_abs_log_resid", "ss_res_log", "n"
    )
  )
  expect_equal(m$n, 3L)
  expect_equal(m$median_ratio, 4)
})
