suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("../../R/datras_constants.R")
source("../../R/datras_hl_helpers.R")
source("../../R/datras_state_helpers.R")

test_that("derive_state_from_mean_length keeps S == n_species_with_lw", {
  lw <- tibble(
    aphia_id = c(126417L, 126425L),
    a_final = c(0.004, 0.005),
    b_final = c(3.1, 3.1)
  )
  hml <- tibble(
    haul_key = rep("NS-IBTS_1985_1_DE_06DA_1", 2L),
    Survey = "NS-IBTS",
    Year = 1985L,
    Quarter = 1L,
    Country = "DE",
    Platform = "06DA",
    HaulNumber = 1L,
    AphiaID = c(126417L, 126425L),
    mean_length_mm = c(100, 80),
    total_n_measured = c(50, 30)
  )

  state <- derive_state_from_mean_length(hml, lw)
  expect_equal(state$S, state$n_species_with_lw)
  expect_equal(state$N, 80)
  expect_gt(state$E_raw, 0)
})

test_that("aggregate_haul_state_from_bins uses minimum bin mass for m_min", {
  lw <- tibble(aphia_id = 126417L, a_final = 0.004, b_final = 3.1)
  hl <- tibble(
    Survey = "NS-IBTS",
    Year = 1985L,
    Quarter = 1L,
    Country = "DE",
    Ship = "06DA",
    HaulNo = 1L,
    Valid_Aphia = 126417L,
    LngtCode = c("1", "1"),
    LngtClass = c(5, 20),
    HLNoAtLngt = c(10L, 5L),
    SubFactor = 1
  )

  bins <- clean_hl_raw(hl) %>% add_lw_mass_to_bins(lw)
  state <- aggregate_haul_state_from_bins(bins)

  small_mass <- bins$mass_g[which.min(bins$length_cm)]
  expect_equal(state$m_min, small_mass)
  expect_equal(state$S, 1L)
  expect_equal(state$N, 15L)
})

test_that("normalize_haul_E satisfies E = E_raw / m_min^0.75", {
  state <- tibble(
    haul_key = "k1",
    Survey = "NS-IBTS",
    Year = 1985L,
    Quarter = 1L,
    Country = "DE",
    Platform = "06DA",
    HaulNumber = 1L,
    S = 2L,
    N = 100,
    E_raw = 500,
    n_species_with_lw = 2L,
    m_min = 0.5,
    m_min_method = "min_length_bin"
  )

  norm <- normalize_haul_E(state)
  expect_equal(norm$E, norm$E_raw / (norm$m_min ^ 0.75))
  expect_equal(norm$min_epsilon, 1)
})

test_that("aggregate_haul_state_from_bins ignores NA mass bins", {
  lw <- tibble(
    aphia_id = c(126417L, 126425L),
    a_final = c(0.004, NA),
    b_final = c(3.1, 3.1)
  )
  hl <- tibble(
    Survey = rep("NS-IBTS", 2L),
    Year = rep(1985L, 2L),
    Quarter = rep(1L, 2L),
    Country = rep("DE", 2L),
    Ship = rep("06DA", 2L),
    HaulNo = rep(1L, 2L),
    Valid_Aphia = c(126417L, 126425L),
    LngtCode = rep("1", 2L),
    LngtClass = c(10, 15),
    HLNoAtLngt = c(5L, 5L),
    SubFactor = 1
  )

  bins <- clean_hl_raw(hl) %>% add_lw_mass_to_bins(lw)
  state <- aggregate_haul_state_from_bins(bins)

  expect_equal(nrow(state), 1L)
  expect_true(is.finite(state$E_raw))
  expect_true(is.finite(state$m_min))
})

test_that("validate_haul_state passes on consistent synthetic data", {
  state <- tibble(
    haul_key = "k1",
    S = 2L,
    N = 100,
    E_raw = 500,
    n_species_with_lw = 2L,
    m_min = 0.5,
    m_min_method = "min_length_bin"
  )
  state <- normalize_haul_E(state)
  v <- validate_haul_state(state)
  expect_true(all(v$pass))
})
