source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h2_couce_helpers.R", local = TRUE)
source("../../R/h2_panel_helpers.R", local = TRUE)

test_that("normalize_stat_rec strips quotes and whitespace", {
  expect_equal(normalize_stat_rec('"""38F5"""'), "38F5")
  expect_equal(normalize_stat_rec(" 39F3 "), "39F3")
})

test_that("read_couce_year_effort filters reconstructed gear hours", {
  path <- testthat::test_path("fixtures", "couce_sample.csv")
  year_effort <- read_couce_year_effort(path, year_min = 1985L, year_max = 2015L)
  expect_equal(nrow(year_effort), 3L)
  expect_true(all(year_effort$hours_total > 0))
  expect_equal(year_effort$hours_total[year_effort$stat_rec == "38F5" & year_effort$year == 1985], 300)
})

test_that("build_h2_rectangle_residuals enforces year filter and abs_residual >= 0", {
  haul <- data.frame(
    stat_rec = c("38F5", "38F5", "39F3", "39F3"),
    year = c(1984L, 1990L, 1990L, 2016L),
    residual = c(-0.5, -0.3, 0.2, 0.1),
    abs_residual = abs(c(-0.5, -0.3, 0.2, 0.1)),
    ln_B_obs = c(8, 8.1, 7.9, 8.2),
    stringsAsFactors = FALSE
  )
  panel <- build_h2_rectangle_residuals(haul)
  expect_equal(nrow(panel), 2L)
  expect_true(all(panel$mean_abs_residual >= 0))
  expect_equal(panel$n_hauls[panel$stat_rec == "38F5"], 1L)
})

test_that("build_h2_analysis_panel joins couce effort and applies min hauls", {
  residual_panel <- data.frame(
    stat_rec = c("38F5", "39F3"),
    n_hauls = c(12L, 8L),
    mean_abs_residual = c(1.1, 0.9),
    mean_residual = c(-0.2, 0.1),
    median_abs_residual = c(1.0, 0.8),
    sd_abs_residual = c(0.2, 0.1),
    mean_ln_B_obs = c(8, 7.8),
    year_min = 1985L,
    year_max = 2015L,
    stringsAsFactors = FALSE
  )
  couce_rectangle <- data.frame(
    stat_rec = c("38F5", "39F3"),
    mean_annual_hours_total = c(1000, 500),
    mean_annual_hours_otter = c(600, 300),
    mean_annual_hours_beam = c(400, 200),
    stringsAsFactors = FALSE
  )
  panel <- build_h2_analysis_panel(
    residual_panel,
    couce_rectangle,
    min_hauls = 10L,
    require_fishing = TRUE
  )
  expect_equal(nrow(panel), 1L)
  expect_equal(panel$stat_rec, "38F5")
})

test_that("couce import diagnostics reports overlap when outputs exist", {
  project_root <- normalizePath(
    file.path(testthat::test_path(), "..", "..", ".."),
    mustWork = TRUE
  )
  path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
  path_couce <- h2_couce_raw_path(project_root)
  skip_if_not(file.exists(path_haul))
  skip_if_not(file.exists(path_couce))

  year_effort <- read_couce_year_effort(path_couce)
  rectangle_effort <- aggregate_couce_rectangle(year_effort)
  haul <- readRDS(path_haul)
  diag <- couce_import_diagnostics(year_effort, rectangle_effort, haul$stat_rec)

  overlap <- diag$value[diag$metric == "n_overlap_stat_rec"]
  expect_gt(overlap, 100)
})
