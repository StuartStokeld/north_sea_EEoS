suppressPackageStartupMessages(library(tibble))
source("../../R/datras_constants.R")
source("../../R/datras_hl_helpers.R")

test_that("LngtCode conversions match ICES vocabulary", {
  expect_equal(lngt_midpoint_cm("1", 10), 10.5)
  expect_equal(lngt_midpoint_cm("0", 50), 5.05)
  expect_equal(lngt_midpoint_cm(".", 50), 5.05)
  expect_equal(lngt_midpoint_cm("5", 50), 5.05)
  expect_true(is.na(lngt_midpoint_cm("9", 10)))
})

test_that("clean_hl_raw applies Q1 year filter and raises counts", {
  hl <- tibble(
    Survey = "NS-IBTS",
    Year = c(1985L, 1985L, 2016L),
    Quarter = c(1L, 2L, 1L),
    Country = "DE",
    Ship = "06DA",
    HaulNo = c(1L, 1L, 1L),
    Valid_Aphia = c(126417L, 126417L, 126417L),
    LngtCode = c("1", "1", "1"),
    LngtClass = c(10, 10, 10),
    HLNoAtLngt = c(5, 5, 5),
    SubFactor = c(1, 1, 1)
  )

  out <- clean_hl_raw(hl)
  expect_equal(nrow(out), 1L)
  expect_equal(out$NumberAtLength, 5)
  expect_equal(out$length_cm, 10.5)
  expect_equal(out$haul_key, "NS-IBTS_1985_1_DE_06DA_1")
})

test_that("non-fish AphiaIDs are excluded", {
  hl <- tibble(
    Survey = "NS-IBTS",
    Year = 1985L,
    Quarter = 1L,
    Country = "DE",
    Ship = "06DA",
    HaulNo = 1L,
    Valid_Aphia = NON_FISH_APHIAIDS[1L],
    LngtCode = "1",
    LngtClass = 10,
    HLNoAtLngt = 5L,
    SubFactor = 1
  )
  expect_equal(nrow(clean_hl_raw(hl)), 0L)
})
