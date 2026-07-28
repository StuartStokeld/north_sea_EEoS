source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h2_dominance_missingness_helpers.R", local = TRUE)

test_that("summarise_haul_counts reports min/median/max and flags the threshold check", {
  df <- data.frame(n_hauls = c(6L, 30L, 8L, 150L, 59L))
  out <- summarise_haul_counts(df, min_hauls = 5L)
  expect_equal(out$n_rectangles, 5L)
  expect_equal(out$min_n_hauls, 6L)
  expect_equal(out$median_n_hauls, 30L)
  expect_equal(out$max_n_hauls, 150L)
  expect_true(out$all_meet_threshold)
})

test_that("summarise_haul_counts flags when not all rectangles meet the threshold", {
  df <- data.frame(n_hauls = c(2L, 10L, 20L))
  out <- summarise_haul_counts(df, min_hauls = 5L)
  expect_false(out$all_meet_threshold)
  expect_equal(out$min_n_hauls, 2L)
})

test_that("ecoregion_breakdown_table reports counts, percentages and the group label", {
  df <- data.frame(Ecoregion = c("Greater North Sea", "Greater North Sea", "Baltic Sea"))
  out <- ecoregion_breakdown_table(df, "dropped (n=3)")
  expect_equal(sum(out$n_rectangles), 3L)
  expect_true(all(out$group == "dropped (n=3)"))
  expect_equal(out$n_rectangles[out$Ecoregion == "Greater North Sea"], 2L)
  expect_equal(out$pct[out$Ecoregion == "Greater North Sea"], round(100 * 2 / 3, 1))
})

test_that("geo_bbox_summary reports the correct lat/lon bounding box", {
  df <- data.frame(
    SOUTH = c(49.5, 55.5), NORTH = c(50, 56),
    WEST = c(-1, 8), EAST = c(0, 9)
  )
  out <- geo_bbox_summary(df, "dropped")
  expect_equal(out$lat_min, 49.5)
  expect_equal(out$lat_max, 56)
  expect_equal(out$lon_min, -1)
  expect_equal(out$lon_max, 9)
  expect_equal(out$n_rectangles, 2L)
})

test_that("join_rectangle_geo left-joins geometry onto a panel by stat_rec", {
  panel <- data.frame(stat_rec = c("28E9", "99Z9"), mean_D = c(0.6, 0.5))
  geo <- data.frame(stat_rec = "28E9", SOUTH = 49.5, NORTH = 50, WEST = -1, EAST = 0, Ecoregion = "Greater North Sea")
  out <- join_rectangle_geo(panel, geo)
  expect_equal(nrow(out), 2L)
  expect_equal(out$Ecoregion[out$stat_rec == "28E9"], "Greater North Sea")
  expect_true(is.na(out$Ecoregion[out$stat_rec == "99Z9"]))
})
