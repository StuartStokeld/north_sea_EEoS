# Layer 4.2: input data contracts for DATRAS HL and HH (committed RDS inputs).
# Code vocabularies: ICES vocab LngtCode (ref 18), DataType (ref 9), HaulVal (ref 1); see audit/sources.md.

suppressPackageStartupMessages(library(dplyr))

hl_key <- function(hl) paste(hl$Survey, hl$Year, hl$Quarter, hl$Country, hl$Ship, hl$HaulNo, sep = "_")

test_that("HL schema: required columns with the expected types", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  types <- c(Survey = "character", Year = "integer", Quarter = "integer", Country = "character",
             Ship = "character", HaulNo = "integer", Valid_Aphia = "integer", LngtCode = "character",
             LngtClass = "numeric", HLNoAtLngt = "numeric", SubFactor = "numeric")
  expect_true(all(names(types) %in% names(hl)))
  expect_identical(vapply(hl[names(types)], function(x) class(x)[1], ""), types)
})

test_that("HL scope: NS-IBTS, quarter 1, 1985-2015", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  expect_identical(unique(hl$Survey), "NS-IBTS")
  expect_true(all(hl$Quarter == 1L))
  expect_true(all(hl$Year >= 1985L & hl$Year <= 2015L))
})

test_that("HL LngtCode takes only ICES vocabulary values (ref 18)", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  ices_codes <- c(".", "0", "1", "2", "5", "9", "-9")
  expect_true(all(hl$LngtCode %in% ices_codes))
})

test_that("HL LngtClass is consistent with the class width of its code", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  # "0": 0.5 cm classes reported in mm -> LngtClass multiple of 5 (one documented exception)
  expect_lte(sum(hl$LngtCode == "0" & hl$LngtClass %% 5 != 0), 1L)
  # "5": 5 cm classes reported in cm -> LngtClass multiple of 5
  expect_equal(sum(hl$LngtCode == "5" & hl$LngtClass %% 5 != 0), 0L)
  # under the ICES reading, every fish is between 0 and 300 cm
  len_cm <- ifelse(hl$LngtCode %in% c(".", "0"), hl$LngtClass / 10, hl$LngtClass)
  expect_true(all(len_cm >= 0 & len_cm <= 300))
})

test_that("HL missing-value sentinels are absent (-9 removed or NA handled at import)", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  for (col in c("LngtClass", "HLNoAtLngt", "SubFactor", "Valid_Aphia")) {
    expect_false(any(hl[[col]] == -9, na.rm = TRUE), info = col)
  }
  expect_true(all(hl$HLNoAtLngt > 0))
  expect_false(anyNA(hl$Valid_Aphia))
})

test_that("HL SubFactor >= 1", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  expect_true(all(hl$SubFactor >= 1))
})

test_that("[F-B11] HL keeps the fields needed to tell sex/category splits from duplicated rows", {
  # The committed HL has 4,097 exact duplicate rows and 13,360 duplicated haul x species x length
  # keys. In the raw ICES CSV these are distinct SpeciesSex / SpeciesCategory rows (0 exact
  # duplicates), but the import drops those columns, so duplicate rows cannot be detected here.
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  expect_true(all(c("Sex", "CatIdentifier") %in% names(hl)) ||
                all(c("SpeciesSex", "SpeciesCategory") %in% names(hl)))
})

test_that("HH schema and vocabularies (DataType ref 9, HaulVal ref 1)", {
  hh <- cached_rds("outputs", "datras_hh_raw.rds")
  expect_true(all(c("haul_key", "HaulDur", "HaulVal", "DataType", "StatRec", "Gear", "hh_key_n") %in% names(hh)))
  expect_true(all(hh$DataType %in% c("C", "R", "P", "S")))
  expect_true(all(hh$HaulVal %in% c("V", "I", "A", "N")))
  expect_true(all(hh$Quarter == 1L & hh$Year >= 1985L & hh$Year <= 2015L))
  expect_true(all(is.finite(hh$HaulDur) & hh$HaulDur > 0 & hh$HaulDur <= 360))
  expect_true(all(grepl("^[0-9]{2}[A-M][0-9]$", hh$StatRec)))
})

test_that("HH haul keys are unique except the documented 1985 NO 58EJ GOV/HOB collision", {
  hh <- cached_rds("outputs", "datras_hh_raw.rds")
  dup <- hh$haul_key[duplicated(hh$haul_key)]
  expect_true(all(grepl("^NS-IBTS_1985_1_NO_58EJ_(4[0-8])$", dup)))
  expect_equal(length(unique(dup)), 9L)
})

test_that("HL hauls are a subset of HH hauls", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  hh <- cached_rds("outputs", "datras_hh_raw.rds")
  expect_true(all(unique(hl_key(hl)) %in% hh$haul_key))
})

test_that("DataType C hauls have SubFactor = 1 on every HL row (DATRAS FAQ v3)", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  hh <- cached_rds("outputs", "datras_hh_raw.rds")
  dt <- hh$DataType[match(hl_key(hl), hh$haul_key)]
  expect_equal(sum(dt == "C" & hl$SubFactor != 1, na.rm = TRUE), 0L)
})

test_that("raw ICES CSV: duplicate HL keys only in the documented collision (slow; CSV is local only)", {
  skip_if_not_slow()
  csv <- list.files(audit_repo_root(), pattern = "^Unaggregated trawl and biological information.*\\.csv$",
                    recursive = TRUE, full.names = TRUE)
  skip_if(length(csv) == 0L, "raw ICES CSV not present")
  source_pipeline("datras_csv_import.R")
  b <- read_ices_csv_blocks(csv[1])
  hl <- b$HL %>% filter(Quarter == 1, Year >= 1985, Year <= 2015)
  expect_equal(sum(duplicated(hl)), 0L)
  k <- hl %>%
    count(Year, Country, Platform, HaulNumber, AphiaID, SpeciesSex, SpeciesCategory,
          DevelopmentStage, LengthCode, LengthClass) %>%
    filter(n > 1)
  expect_true(all(k$Year == 1985 & k$Country == "NO" & k$Platform == "58EJ" & k$HaulNumber %in% 40:48))
})
