suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("../../R/datras_constants.R")
source("../../R/datras_hl_helpers.R")
source("../../R/datras_csv_import.R")
source("../../R/datras_hh_helpers.R")

# Two-block ICES CSV: HH block then HL block, each with its own header.
write_two_block_csv <- function(path) {
  writeLines(c(
    "RecordHeader,Quarter,Country,Platform,Gear,HaulNumber,Year,Month,HaulDuration,HaulValidity,DataType,StatisticalRectangle,Survey,DATRASSurvey",
    "HH,1,DK,26D4,GOV,1,1990,2,30,V,C,38F0,G1022,NS-IBTS",
    "HH,1,GB-SCT,74E9,GOV,2,1990,2,60,V,R,44E8,G1022,NS-IBTS",
    "RecordHeader,Quarter,Country,Platform,Gear,HaulNumber,Year,SpeciesCategory,SubsamplingFactor,LengthCode,LengthClass,NumberAtLength,Survey,AphiaID,DATRASSurvey",
    "HL,1,DK,26D4,GOV,1,1990,1,1,1,20,4,G1022,126436,NS-IBTS",
    "HL,1,DK,26D4,GOV,1,1990,1,1,1,21,6,G1022,126436,NS-IBTS",
    "HL,1,GB-SCT,74E9,GOV,2,1990,1,2.5,1,30,2,G1022,126437,NS-IBTS"
  ), path)
}

test_that("read_ices_csv_blocks parses each block with its own header", {
  path <- tempfile(fileext = ".csv")
  write_two_block_csv(path)
  blocks <- read_ices_csv_blocks(path)
  expect_setequal(names(blocks), c("HH", "HL"))
  expect_equal(nrow(blocks$HH), 2L)
  expect_equal(nrow(blocks$HL), 3L)
  expect_true("HaulDuration" %in% names(blocks$HH))
  expect_true("NumberAtLength" %in% names(blocks$HL))
  expect_equal(blocks$HL$NumberAtLength, c(4, 6, 2))
})

test_that("HH and HL imports use DATRASSurvey for the haul key", {
  path <- tempfile(fileext = ".csv")
  write_two_block_csv(path)
  blocks <- read_ices_csv_blocks(path)
  hh <- import_hh_from_ices_csv(path, blocks = blocks, year_min = 1985L, year_max = 2015L)
  hl <- import_hl_from_ices_csv(path, blocks = blocks)
  expect_equal(unique(hh$Survey), "NS-IBTS")
  expect_equal(unique(hl$Survey), "NS-IBTS")
  expect_true(all(
    make_haul_key(hl$Survey, hl$Year, hl$Quarter, hl$Country, hl$Ship, hl$HaulNo) %in% hh$haul_key
  ))
  expect_equal(sort(hh$HaulDur), c(30, 60))
})

make_hh <- function(...) {
  tibble(...) %>%
    mutate(haul_key = make_haul_key(Survey, Year, Quarter, Country, Ship, HaulNo)) %>%
    add_count(haul_key, name = "hh_key_n")
}

hl_fixture <- tibble(
  Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX", Ship = "S1",
  HaulNo = c(1L, 1L, 2L, 3L, 4L, 5L, 6L, 6L, 7L),
  Valid_Aphia = 126436L, LngtCode = "1", LngtClass = 20,
  HLNoAtLngt = c(4, 6, 10, 10, 10, 10, 3, 5, 8),
  SubFactor = c(1, 1, 1, 1, 1, 1, 2, 4, 1)
)
hh_fixture <- make_hh(
  Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX",
  Ship = c("S1", "S1", "S1", "S1", "S1", "S1", "S1"),
  HaulNo = c(1L, 2L, 3L, 4L, 5L, 5L, 6L),
  HaulDur = c(30, 60, 30, 30, 30, 30, 26),
  HaulVal = c("V", "V", "I", "V", "V", "V", "V"),
  DataType = c("C", "C", "C", "X", "R", "R", "P")
)

test_that("align_hl_to_hh assigns one reason per haul in the documented order", {
  # 26-34 min sensitivity window, so the 26-min P haul sits on the window edge
  out <- align_hl_to_hh(hl_fixture, hh_fixture, dur_min = 26, dur_max = 34, min_match_rate = 0)
  reasons <- setNames(out$hauls$effort_reason, out$hauls$HaulNo)
  expect_equal(reasons[["1"]], "kept")            # C, 30 min
  expect_equal(reasons[["2"]], "duration")        # 60 min
  expect_equal(reasons[["3"]], "invalid")         # HaulVal I
  expect_equal(reasons[["4"]], "datatype")        # unknown DataType
  expect_equal(reasons[["5"]], "ambiguous_hh_key")
  expect_equal(reasons[["6"]], "kept")            # P, 26 min (window edge)
  expect_equal(reasons[["7"]], "no_hh")
  expect_equal(sum(out$summary$n_hauls), 7L)
})

test_that("default (primary) window keeps exactly 30-min hauls only", {
  expect_equal(c(HAUL_DUR_MIN, HAUL_DUR_MAX), c(30, 30))
  out <- align_hl_to_hh(hl_fixture, hh_fixture, min_match_rate = 0)
  reasons <- setNames(out$hauls$effort_reason, out$hauls$HaulNo)
  expect_equal(reasons[["1"]], "kept")            # C, 30 min
  expect_equal(reasons[["6"]], "duration")        # P, 26 min: outside exact-30 window
  expect_true(all(out$hl$HaulDur == 30))
})

test_that("C counts are converted per haul with HaulDur/60; R/P raised by SubFactor", {
  out <- align_hl_to_hh(hl_fixture, hh_fixture, dur_min = 26, dur_max = 34, min_match_rate = 0)
  c_rows <- out$hl %>% filter(HaulNo == 1L)
  expect_equal(c_rows$HLNoAtLngt, c(4, 6) * 30 / 60)
  # P haul with two subcategories / SubFactors: each row raised by its own factor
  p_rows <- out$hl %>% filter(HaulNo == 6L)
  expect_equal(sum(p_rows$HLNoAtLngt), 3 * 2 + 5 * 4)
  expect_true(all(out$hl$SubFactor == 1))
  # clean_hl_raw() then yields NumberAtLength = raw count
  bins <- clean_hl_raw(out$hl)
  expect_equal(sum(bins$NumberAtLength[bins$HaulNo == 6L]), 26)
})

test_that("align_hl_to_hh stops when HL hauls do not match HH", {
  expect_error(align_hl_to_hh(hl_fixture, hh_fixture[0, ]), "match an HH record")
})

test_that("haul_dur_window reads env overrides and rejects bad values", {
  withr_set <- function(vars) {
    old <- Sys.getenv(names(vars), unset = NA)
    do.call(Sys.setenv, as.list(vars))
    old
  }
  restore <- function(old) {
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
    }
  }
  old <- withr_set(c(EEOS_HAUL_DUR_MIN = "26", EEOS_HAUL_DUR_MAX = "34"))
  on.exit(restore(old), add = TRUE)
  expect_equal(haul_dur_window(), c(min = 26, max = 34))
  Sys.setenv(EEOS_HAUL_DUR_MIN = "abc")
  expect_error(haul_dur_window(), "non-negative number")
  Sys.setenv(EEOS_HAUL_DUR_MIN = "40")
  expect_error(haul_dur_window(), "empty")
  Sys.unsetenv(c("EEOS_HAUL_DUR_MIN", "EEOS_HAUL_DUR_MAX"))
  expect_equal(haul_dur_window(), c(min = HAUL_DUR_MIN, max = HAUL_DUR_MAX))
})
