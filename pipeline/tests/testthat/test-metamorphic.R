# Layer 4.6 (scoped): metamorphic tests - C <-> R equivalence and row-order invariance.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source_pipeline("datras_constants.R")
source_pipeline("datras_hl_helpers.R")
source_pipeline("datras_state_helpers.R")
source_pipeline("datras_csv_import.R")
source_pipeline("datras_hh_helpers.R")

lw_m <- tibble(aphia_id = c(1L, 2L, 3L), a_final = c(0.01, 0.02, 0.005), b_final = c(3, 2.5, 3.2))

state_of <- function(hl, hh) {
  al <- align_hl_to_hh(hl, hh, min_match_rate = 0)
  build_state_from_hl(al$hl, lw_m)$haul_state %>% normalize_haul_E() %>% arrange(haul_key)
}

hh_one <- function(haul_no, dt) {
  tibble(Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX", Ship = "S1", HaulNo = haul_no,
         HaulDur = 30, HaulVal = "V", DataType = dt) %>%
    mutate(haul_key = make_haul_key(Survey, Year, Quarter, Country, Ship, HaulNo), hh_key_n = 1L)
}

test_that("a 30-min C haul and the equivalent R haul (counts x 60/HaulDur) give identical state variables", {
  raw <- tibble(Valid_Aphia = c(1L, 1L, 2L, 3L), LngtCode = c("1", "1", "0", "."),
                LngtClass = c(10, 12, 150, 85), n_raw = c(2, 3, 5, 1))
  base <- tibble(Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX", Ship = "S1")
  # R haul: raw counts with SubFactor 1; C haul: same catch reported per hour (x 60/30), SubFactor 1
  hl_r <- bind_cols(base, raw) %>% mutate(HaulNo = 1L, HLNoAtLngt = n_raw, SubFactor = 1) %>% select(-n_raw)
  hl_c <- bind_cols(base, raw) %>% mutate(HaulNo = 2L, HLNoAtLngt = n_raw * 60 / 30, SubFactor = 1) %>% select(-n_raw)
  st <- state_of(bind_rows(hl_r, hl_c), bind_rows(hh_one(1L, "R"), hh_one(2L, "C")))
  expect_equal(nrow(st), 2L)
  r <- st[st$HaulNumber == 1L, ]
  c <- st[st$HaulNumber == 2L, ]
  for (v in c("S", "N", "E_raw", "m_min", "E")) expect_equal(c[[v]], r[[v]], info = v)
})

test_that("an R haul with SubFactor 2 equals the same catch reported unsubsampled", {
  base <- tibble(Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX", Ship = "S1",
                 Valid_Aphia = c(1L, 2L), LngtCode = "1", LngtClass = c(10, 20))
  hl <- bind_rows(base %>% mutate(HaulNo = 1L, HLNoAtLngt = c(3, 4), SubFactor = 2),
                  base %>% mutate(HaulNo = 2L, HLNoAtLngt = c(6, 8), SubFactor = 1))
  st <- state_of(hl, bind_rows(hh_one(1L, "R"), hh_one(2L, "R")))
  for (v in c("S", "N", "E_raw", "m_min")) expect_equal(st[[v]][1], st[[v]][2], info = v)
})

test_that("shuffling HL row order changes nothing (real data, first 200 kept hauls)", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  hh <- cached_rds("outputs", "datras_hh_raw.rds")
  keys <- paste(hl$Survey, hl$Year, hl$Quarter, hl$Country, hl$Ship, hl$HaulNo, sep = "_")
  pick <- head(unique(keys[keys %in% hh$haul_key[hh$HaulDur == 30 & hh$HaulVal == "V"]]), 200)
  sub <- hl[keys %in% pick, ]
  lw <- read.csv(repo_path("outputs", "fishbase_lw_lookup_v2.csv"))
  run <- function(x) {
    al <- align_hl_to_hh(x, hh, min_match_rate = 0)
    build_state_from_hl(al$hl, lw)$haul_state %>% arrange(haul_key) %>% as.data.frame()
  }
  a <- run(sub)
  set.seed(11)
  b <- run(sub[sample(nrow(sub)), ])
  expect_equal(nrow(a), length(pick))
  expect_equal(a, b)
})
