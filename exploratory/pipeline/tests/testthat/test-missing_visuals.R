source("../../R/missing_visuals_helpers.R", local = TRUE)

test_that("build_composition_table joins S onto the dominance table without loss", {
  dominance <- data.frame(haul_id = c("a", "b", "c"), D = c(0.5, 0.6, 0.9), size_CV = c(0.3, 0.4, 0.1),
                           stringsAsFactors = FALSE)
  haul_eeos <- data.frame(haul_id = c("c", "a", "b"), S = c(10L, 12L, 8L), other = 1:3,
                           stringsAsFactors = FALSE)
  out <- build_composition_table(dominance, haul_eeos)
  expect_equal(nrow(out), 3L)
  expect_equal(out$S[out$haul_id == "a"], 12L)
  expect_equal(out$S[out$haul_id == "c"], 10L)
})

test_that("build_composition_table errors on an unmatched haul_id", {
  dominance <- data.frame(haul_id = c("a", "b"), D = c(0.5, 0.6), size_CV = c(0.3, 0.4), stringsAsFactors = FALSE)
  haul_eeos <- data.frame(haul_id = "a", S = 10L, stringsAsFactors = FALSE)
  expect_error(build_composition_table(dominance, haul_eeos), "unmatched")
})

test_that("read_D_top_decile_reference reads the pre-existing top-decile boundary, not a fresh quantile", {
  tmp <- tempfile(fileext = ".csv")
  bins <- data.frame(
    metric = c("D", "D", "size_CV"),
    n_bins = c(10L, 10L, 10L),
    bin = c(9L, 10L, 10L),
    n = c(1207L, 1206L, 1206L),
    metric_min = c(0.82, 0.9005684, 0.05),
    metric_max = c(0.90, 0.9997592, 3.5),
    stringsAsFactors = FALSE
  )
  write.csv(bins, tmp, row.names = FALSE)
  ref <- read_D_top_decile_reference(tmp)
  expect_equal(ref$threshold, 0.9005684)
  expect_equal(ref$n_in_decile, 1206L)
  unlink(tmp)
})

test_that("read_D_top_decile_reference errors if the expected row is missing", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(metric = "D", n_bins = 4L, bin = 4L, n = 1L, metric_min = 0.5), tmp, row.names = FALSE)
  expect_error(read_D_top_decile_reference(tmp), "Could not locate")
  unlink(tmp)
})
