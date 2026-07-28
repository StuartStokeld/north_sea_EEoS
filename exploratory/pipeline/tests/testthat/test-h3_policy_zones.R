source("../../R/h1_common.R", local = TRUE)
source("../../R/h2_common.R", local = TRUE)
source("../../R/h3_pre_exploration_helpers.R", local = TRUE) # variance_components_anova()
source("../../R/h3_policy_zones_helpers.R", local = TRUE)

fake_ices_sf <- data.frame(
  stat_rec = c("40F2", "41F2", "42F2", "41F3", "40F4"),
  SOUTH = c(55.5, 56.0, 56.5, 56.0, 55.5),
  WEST = c(2, 2, 2, 3, 4),
  stringsAsFactors = FALSE
)

test_that("build_rect_grid_index derives correct row/col indices from SOUTH/WEST", {
  gi <- build_rect_grid_index(fake_ices_sf)
  expect_equal(gi$row_idx[gi$stat_rec == "40F2"], 111L)
  expect_equal(gi$row_idx[gi$stat_rec == "41F2"], 112L)
  expect_equal(gi$col_idx[gi$stat_rec == "41F3"], 3L)
})

test_that("build_scheme_a_blocks groups rectangles into contiguous blocks", {
  gi <- build_rect_grid_index(fake_ices_sf)
  blocks <- build_scheme_a_blocks(gi, fake_ices_sf$stat_rec, block_size = 2L)
  expect_equal(nrow(blocks), 5L)
  # row_idx: 40F2=111, 41F2=112, 42F2=113 -> with block_size=2, 111 %/% 2 = 55
  # but 112, 113 %/% 2 = 56, so 41F2 and 42F2 share a block; 40F2 falls in the block below.
  expect_equal(blocks$unit_id[blocks$stat_rec == "41F2"], blocks$unit_id[blocks$stat_rec == "42F2"])
  expect_false(blocks$unit_id[blocks$stat_rec == "40F2"] == blocks$unit_id[blocks$stat_rec == "41F2"])
})

test_that("build_rook_adjacency finds edge neighbours only (not diagonal)", {
  gi <- build_rect_grid_index(fake_ices_sf)
  adj <- build_rook_adjacency(gi, fake_ices_sf$stat_rec)
  expect_true("41F2" %in% adj[["40F2"]])
  expect_true("41F3" %in% adj[["41F2"]])
  expect_false("40F4" %in% adj[["40F2"]])
})

test_that("assign_pressure_tiers produces near-equal-sized quantile groups", {
  rect_pressure <- data.frame(
    stat_rec = paste0("R", 1:9),
    mean_annual_hours_total = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    stringsAsFactors = FALSE
  )
  tiers <- assign_pressure_tiers(rect_pressure, n_tiers = 3L)
  expect_equal(as.integer(table(tiers$tier)), c(3L, 3L, 3L))
  expect_equal(tiers$tier[tiers$stat_rec == "R1"], 1L)
  expect_equal(tiers$tier[tiers$stat_rec == "R9"], 3L)
})

test_that("build_contiguous_zones merges connected same-tier rectangles and isolates disconnected ones", {
  # Grid: A-B-C in a row (all tier 1), D isolated (tier 1, but not adjacent to A/B/C)
  adj <- list(A = c("B"), B = c("A", "C"), C = c("B"), D = character(0))
  zones <- build_contiguous_zones(c("A", "B", "C", "D"), c(1L, 1L, 1L, 1L), adj)
  expect_equal(zones$unit_id[zones$stat_rec == "A"], zones$unit_id[zones$stat_rec == "B"])
  expect_equal(zones$unit_id[zones$stat_rec == "B"], zones$unit_id[zones$stat_rec == "C"])
  expect_false(zones$unit_id[zones$stat_rec == "D"] == zones$unit_id[zones$stat_rec == "A"])
})

test_that("build_contiguous_zones does not merge adjacent rectangles of different tiers", {
  adj <- list(A = c("B"), B = c("A"))
  zones <- build_contiguous_zones(c("A", "B"), c(1L, 2L), adj)
  expect_false(zones$unit_id[zones$stat_rec == "A"] == zones$unit_id[zones$stat_rec == "B"])
})

test_that("assign_period splits at break_year with pre strictly before", {
  years <- c(2000L, 2002L, 2003L, 2010L)
  periods <- assign_period(years, break_year = 2003L)
  expect_equal(periods, c("pre", "pre", "post", "post"))
})

test_that("feasibility_summary_row counts units usable in both periods and reports distribution", {
  unit_period_hauls <- data.frame(
    unit_id = c("u1", "u1", "u2", "u2", "u3"),
    period = c("pre", "post", "pre", "post", "pre"),
    n_hauls = c(10L, 8L, 3L, 20L, 50L),
    stringsAsFactors = FALSE
  )
  out <- feasibility_summary_row(unit_period_hauls, "SchemeX", "paramY", min_hauls = 5L, n_units_total = 3L)
  expect_equal(out$n_units, 3L)
  # u1: pre=10,post=8 both >=5 -> usable. u2: pre=3 (<5) -> not usable. u3: post=0 (<5, missing) -> not usable.
  expect_equal(out$n_units_usable_both_periods, 1L)
  expect_equal(out$post_min_hauls, 0L) # u3 has no post row -> filled to 0
})

test_that("summarise_unit_period_hauls counts hauls per unit x period correctly", {
  haul_full <- data.frame(
    stat_rec = c("40F2", "40F2", "41F2", "41F2"),
    year = c(1990L, 2005L, 1990L, 2005L),
    stringsAsFactors = FALSE
  )
  unit_map <- data.frame(stat_rec = c("40F2", "41F2"), unit_id = c("blockA", "blockA"), stringsAsFactors = FALSE)
  out <- summarise_unit_period_hauls(haul_full, unit_map, break_year = 2003L, year_min = 1985L, year_max = 2015L)
  expect_equal(out$n_hauls[out$period == "pre"], 2L)
  expect_equal(out$n_hauls[out$period == "post"], 2L)
})

test_that("feasibility_across_thresholds reproduces feasibility_summary_row at each threshold", {
  unit_period_hauls <- data.frame(
    unit_id = c("u1", "u1", "u2", "u2", "u3"),
    period = c("pre", "post", "pre", "post", "pre"),
    n_hauls = c(10L, 8L, 3L, 20L, 50L),
    stringsAsFactors = FALSE
  )
  out <- feasibility_across_thresholds(
    unit_period_hauls, "SchemeX", "paramY",
    thresholds = c(1L, 5L, 10L), threshold_labels = c("fixed_1", "fixed_5", "fixed_10"),
    n_units_total = 3L
  )
  expect_equal(nrow(out), 3L)
  # threshold=1: u1 (10,8) and u2 (3,20) both usable, u3 has no post row (0) -> not usable
  expect_equal(out$n_units_usable_both_periods[out$threshold_value == 1L], 2L)
  # threshold=5: matches feasibility_summary_row's min_hauls=5 case
  ref <- feasibility_summary_row(unit_period_hauls, "SchemeX", "paramY", min_hauls = 5L, n_units_total = 3L)
  expect_equal(out$n_units_usable_both_periods[out$threshold_value == 5L], ref$n_units_usable_both_periods)
  # threshold=10: only u2 has both pre>=10 and post>=10 (pre=3 fails) -> none usable at 10 for u1/u2/u3
  expect_equal(out$n_units_usable_both_periods[out$threshold_value == 10L], 0L)
  expect_equal(out$pct_units_usable_both_periods[out$threshold_value == 1L], round(100 * 2 / 3, 1))
})

test_that("required_n_for_reliability inverts the Spearman-Brown formula correctly", {
  # icc = 0.2: at n=4, reliability = 4*0.2/(1+3*0.2) = 0.8/1.6 = 0.5
  expect_equal(4 * 0.2 / (1 + 3 * 0.2), 0.5)
  n <- required_n_for_reliability(icc = 0.2, target_reliability = 0.5)
  expect_equal(n, 4)
  # higher target reliability requires more replicates
  n_targets <- required_n_for_reliability(icc = 0.2, target_reliability = c(0.5, 0.7, 0.9))
  expect_true(all(diff(n_targets) >= 0))
  # icc <= 0 -> no amount of replication reaches any positive target reliability
  expect_equal(required_n_for_reliability(icc = 0, target_reliability = 0.8), Inf)
  expect_equal(required_n_for_reliability(icc = -0.01, target_reliability = 0.8), Inf)
})

test_that("label_haul_cells tags hauls with unit_id x period cell_id and filters to universe/year window", {
  haul_eeos <- data.frame(
    stat_rec = c("40F2", "40F2", "41F2", "99Z9"),
    year = c(1990L, 2005L, 1990L, 1990L),
    residual = c(0.1, 0.2, -0.1, 0.5),
    stringsAsFactors = FALSE
  )
  unit_map <- data.frame(stat_rec = c("40F2", "41F2"), unit_id = c("blockA", "blockA"), stringsAsFactors = FALSE)
  out <- label_haul_cells(haul_eeos, unit_map, break_year = 2003L, year_min = 1985L, year_max = 2015L)
  expect_equal(nrow(out), 3L) # 99Z9 dropped (not in unit_map)
  expect_true(all(c("blockA__pre", "blockA__post") %in% out$cell_id))
})

test_that("sample_size_adequacy_row surfaces variance components and required-n columns", {
  set.seed(1)
  haul_cells <- data.frame(
    cell_id = rep(c("A", "B", "C"), each = 20),
    residual = c(rnorm(20, mean = 0), rnorm(20, mean = 1), rnorm(20, mean = 2))
  )
  out <- sample_size_adequacy_row(haul_cells, "SchemeX", "paramY", reliability_targets = c(0.7, 0.8))
  expect_equal(out$n_cells, 3L)
  expect_equal(out$n_hauls, 60L)
  expect_true(out$icc > 0 && out$icc <= 1)
  expect_true("required_n_reliability_0.7" %in% names(out))
  expect_true("required_n_reliability_0.8" %in% names(out))
  expect_true(out$required_n_reliability_0.8 >= out$required_n_reliability_0.7)
})
