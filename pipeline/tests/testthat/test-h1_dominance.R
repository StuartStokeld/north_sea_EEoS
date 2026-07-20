suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
source("../../R/datras_constants.R")
source("../../R/datras_hl_helpers.R")
source("../../R/h1_dominance_helpers.R")

# --- Synthetic hl_mass fixture -----------------------------------------
# haul_1: species A dominant (80/100), 3 length bins spanning a wide mass range
# haul_2: species B dominant (60/100), single length bin (n_bins_dominant_species edge case)
hl_mass <- tibble(
  haul_key = c(
    "h1", "h1", "h1", "h1", "h1",
    "h2", "h2"
  ),
  AphiaID = c(
    "A", "A", "A", "B", "B",
    "B", "C"
  ),
  length_cm = c(10, 20, 30, 15, 15, 25, 25),
  NumberAtLength = c(40, 20, 20, 10, 10, 60, 40),
  mass_g = c(10, 40, 90, 20, 20, 50, 30)
)

test_that("species_haul_abundance sums raised counts per haul x species", {
  out <- species_haul_abundance(hl_mass)
  h1_a <- out %>% filter(haul_key == "h1", AphiaID == "A") %>% pull(n_species)
  expect_equal(h1_a, 80)
  h1_b <- out %>% filter(haul_key == "h1", AphiaID == "B") %>% pull(n_species)
  expect_equal(h1_b, 20)
})

test_that("dominant_species_per_haul picks the numerically dominant species", {
  abund <- species_haul_abundance(hl_mass)
  dom <- dominant_species_per_haul(abund)
  expect_equal(dom$dominant_aphia_id[dom$haul_key == "h1"], "A")
  expect_equal(dom$n_max_species[dom$haul_key == "h1"], 80)
  expect_equal(dom$dominant_aphia_id[dom$haul_key == "h2"], "B")
  expect_equal(dom$n_max_species[dom$haul_key == "h2"], 60)
})

test_that("dominant_species_per_haul breaks ties by lowest AphiaID", {
  tie <- tibble(
    haul_key = c("h3", "h3"),
    AphiaID = c("Z", "A"),
    n_species = c(50, 50)
  )
  dom <- dominant_species_per_haul(tie)
  expect_equal(dom$dominant_aphia_id, "A")
})

test_that("compute_haul_dominance uses supplied production N, not a recomputed sum", {
  haul_N <- tibble(haul_key = c("h1", "h2"), N = c(100, 100))
  out <- compute_haul_dominance(hl_mass, haul_N)

  expect_equal(out$D[out$haul_key == "h1"], 80 / 100)
  expect_equal(out$D[out$haul_key == "h2"], 60 / 100)
  # recomputed check should match supplied N when the bin set is complete
  expect_equal(out$n_recomputed_check, out$N)
})

test_that("dominant_species_size_cv matches hand-computed weighted CV", {
  haul_N <- tibble(haul_key = c("h1", "h2"), N = c(100, 100))
  dom <- compute_haul_dominance(hl_mass, haul_N)
  cv <- dominant_species_size_cv(hl_mass, dom)

  # h1 dominant species A: counts 40,20,20; mass 10,40,90
  mean_mass <- sum(c(40, 20, 20) * c(10, 40, 90)) / sum(c(40, 20, 20))
  var_mass <- sum(c(40, 20, 20) * (c(10, 40, 90) - mean_mass)^2) / sum(c(40, 20, 20))
  expected_cv <- sqrt(var_mass) / mean_mass

  h1_row <- cv %>% filter(haul_key == "h1")
  expect_equal(h1_row$size_CV, expected_cv)
  expect_equal(h1_row$n_bins_dominant_species, 3L)

  # h2 dominant species B: single bin -> CV mechanically zero, flagged via n_bins == 1
  h2_row <- cv %>% filter(haul_key == "h2")
  expect_equal(h2_row$size_CV, 0)
  expect_equal(h2_row$n_bins_dominant_species, 1L)
})

test_that("build_haul_dominance_table assembles per-haul outputs with outcome columns", {
  haul_outcomes <- tibble(
    haul_key = c("h1", "h2"),
    haul_id = c("id1", "id2"),
    year = c(2000L, 2001L),
    stat_rec = c("39F1", "39F2"),
    N = c(100, 100),
    B_obs = c(1000, 2000),
    B_pred = c(3000, 4000),
    residual = c(log(1000) - log(3000), log(2000) - log(4000))
  )
  lw_lookup <- tibble(
    aphia_id = c("A", "B", "C"),
    accepted_name = c("Species A", "Species B", "Species C")
  )

  out <- build_haul_dominance_table(hl_mass, lw_lookup, haul_outcomes)

  expect_equal(nrow(out), 2L)
  expect_true(all(c("D", "size_CV", "n_bins_dominant_species", "dominant_species",
                     "pred_obs_ratio", "ln_ratio", "B_obs_quartile") %in% names(out)))
  expect_equal(out$dominant_species[out$haul_key == "h1"], "Species A")
  expect_equal(out$pred_obs_ratio[out$haul_key == "h1"], 3)
})

test_that("bin_ratio_table reports median and IQR of ratio per bin", {
  df <- tibble(
    metric = c(1, 2, 3, 4, 5, 6, 7, 8),
    pred_obs_ratio = c(1, 1, 2, 2, 3, 3, 4, 4)
  )
  out <- bin_ratio_table(df, "metric", n_bins = 4L)
  expect_equal(nrow(out), 4L)
  expect_equal(sum(out$n), 8L)
  expect_true(all(c("median_ratio", "iqr_ratio_low", "iqr_ratio_high") %in% names(out)))
})

test_that("dominance_correlation_matrix flags strong pairwise correlations", {
  df <- tibble(
    x = 1:20,
    y = 1:20 * 2,
    z = rev(1:20)
  )
  out <- dominance_correlation_matrix(df, c("x", "y", "z"))
  expect_equal(out$pearson_r[out$var1 == "x" & out$var2 == "y"], 1)
  expect_true(out$flag_gt_0.4[out$var1 == "x" & out$var2 == "y"])
})

test_that("dataquality_crosstab reports flagged-year fraction among low-CV hauls", {
  df <- tibble(
    year = c(1998L, 1998L, 2000L, 2000L, 2000L, 2000L, 2000L, 2000L, 2000L, 2000L),
    size_CV = c(0.01, 0.02, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)
  )
  out <- dataquality_crosstab(df, decile = 1L)
  expect_true(out$n_low_cv_total >= 1L)
  expect_true("flagged" %in% names(out$flagged_summary))
})

test_that("taxonomic_breakdown is descriptive and does not filter the input", {
  df <- tibble(
    metric_extreme = 1:10,
    dominant_species = c(rep("Herring", 5), rep("Cod", 5))
  )
  out <- taxonomic_breakdown(df, "metric_extreme", extreme = "top")
  expect_true(all(c("dominant_species", "n_hauls", "pct_of_decile") %in% names(out)))
  expect_equal(sum(out$n_hauls), 1L) # top decile of 10 rows = 1 row
})
