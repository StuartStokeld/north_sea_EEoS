# Layer 4.3: unit tests with hand-computed expected values (arithmetic in the comments).
# Tests tagged [F-..] encode an audit finding and are expected to fail until it is fixed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source_pipeline("h1_common.R")
source_pipeline("datras_constants.R")
source_pipeline("datras_hl_helpers.R")
source_pipeline("datras_state_helpers.R")
source_pipeline("h1_join_helpers.R")
source_pipeline("h2_common.R")
source_pipeline("h2h3_feasibility_helpers.R")
source_pipeline("h2h3_within_between_helpers.R")

# --- Length midpoints (ICES vocab ref 18) -------------------------------------

test_that("LngtCode '.' (1 mm classes, mm) and '1' (1 cm classes, cm) midpoints", {
  expect_equal(lngt_midpoint_cm(".", 50), 5.05) # (50 + 0.5) mm = 5.05 cm
  expect_equal(lngt_midpoint_cm("1", 10), 10.5) # 10 cm + 0.5 cm
})

test_that("[F-A1] LngtCode '0' = 0.5 cm classes reported in mm: midpoint = class + 2.5 mm", {
  expect_equal(lngt_midpoint_cm("0", 50), 5.25) # (50 + 2.5) mm = 5.25 cm
  expect_equal(lngt_midpoint_cm("0", 195), 19.75)
})

test_that("[F-A2] LngtCode '5' = 5 cm classes reported in cm: midpoint = class + 2.5 cm", {
  expect_equal(lngt_midpoint_cm("5", 60), 62.5)
  expect_equal(lngt_midpoint_cm("5", 145), 147.5)
})

test_that("LngtCode '9' (+ group) and '-9' (missing) have no midpoint", {
  expect_true(is.na(lngt_midpoint_cm("9", 10)))
  expect_true(is.na(lngt_midpoint_cm("-9", 10)))
})

# --- LW mass, S, N, E, m_min --------------------------------------------------

lw_fix <- tibble(aphia_id = c(1L, 2L), a_final = c(0.01, 0.02), b_final = c(3, 2.5))
hl_fix <- tibble(
  Survey = "NS-IBTS", Year = 1990L, Quarter = 1L, Country = "XX", Ship = "S1", HaulNo = 1L,
  Valid_Aphia = c(1L, 2L), LngtCode = "1", LngtClass = c(9, 19),
  HLNoAtLngt = c(2, 3), SubFactor = c(1, 2)
)

test_that("LW mass in grams: W = a L^b with L in cm", {
  bins <- clean_hl_raw(hl_fix)
  m <- add_lw_mass_to_bins(bins, lw_fix)
  # species 1: L = 9.5 cm, W = 0.01 * 9.5^3 = 8.57375 g
  # species 2: L = 19.5 cm, W = 0.02 * 19.5^2.5 = 33.5827706942 g
  expect_equal(m$mass_g[m$AphiaID == 1L], 8.57375)
  expect_equal(m$mass_g[m$AphiaID == 2L], 33.5827706942, tolerance = 1e-10)
})

test_that("S, N, E_raw, m_min and E on a two-species haul", {
  st <- build_state_from_hl(hl_fix, lw_fix)$haul_state %>% normalize_haul_E()
  # S = 2 species; N = 2*1 + 3*2 = 8 individuals (SubFactor raises the second row)
  expect_equal(st$S, 2L)
  expect_equal(st$N, 8)
  # E_raw = 2 * 8.57375^0.75 + 6 * 33.5827706942^0.75 = 93.7234730471
  expect_equal(st$E_raw, 93.7234730471, tolerance = 1e-10)
  # m_min = lightest bin = 8.57375 g; E = E_raw / m_min^0.75 = 18.7055404976
  expect_equal(st$m_min, 8.57375)
  expect_equal(st$E, 18.7055404976, tolerance = 1e-10)
})

test_that("non-fish taxa, missing AphiaID, zero counts and missing lengths are dropped", {
  hl <- bind_rows(
    hl_fix,
    hl_fix[1, ] %>% mutate(Valid_Aphia = NON_FISH_APHIAIDS[1]),
    hl_fix[1, ] %>% mutate(Valid_Aphia = NA_integer_),
    hl_fix[1, ] %>% mutate(HLNoAtLngt = 0),
    hl_fix[1, ] %>% mutate(LngtClass = NA_real_)
  )
  expect_equal(nrow(clean_hl_raw(hl)), 2L)
})

test_that("haul key is Survey_Year_Quarter_Country_Ship_HaulNo", {
  expect_identical(clean_hl_raw(hl_fix)$haul_key[1], "NS-IBTS_1990_1_XX_S1_1")
})

# --- FishGlob join: haul_id parsing, kg -> g, NA weights ----------------------

write_fishglob_fixture <- function(wgt) {
  data <- tibble(
    haul_id = c("NS-IBTS 1990 1 XX S1 GOV 7 1", "NS-IBTS 1990 1 XX S1 GOV 7 1", "NS-IBTS 1990 1 XX S1 GOV 8 2"),
    quarter = 1L, year = 1990L, stat_rec = "38F0", wgt = wgt
  )
  p <- tempfile(fileext = ".RData")
  save(data, file = p)
  p
}

test_that("FishGlob haul_id -> 6-field join key and B_obs kg -> g", {
  fg <- build_fishglob_haul_table(write_fishglob_fixture(c(1.5, 0.25, 2)))
  # haul_id fields 1-5 and 8 (gear and station dropped)
  expect_setequal(fg$join_key, c("NS-IBTS_1990_1_XX_S1_1", "NS-IBTS_1990_1_XX_S1_2"))
  # haul 1: (1.5 + 0.25) kg = 1750 g; haul 2: 2 kg = 2000 g
  expect_equal(fg$B_obs[fg$join_key == "NS-IBTS_1990_1_XX_S1_1"], 1750)
  expect_equal(fg$B_obs[fg$join_key == "NS-IBTS_1990_1_XX_S1_2"], 2000)
})

test_that("[F-B9] a haul with an NA species weight does not get a silently reduced B_obs", {
  fg <- build_fishglob_haul_table(write_fishglob_fixture(c(1.5, NA, 2)))
  b1 <- fg$B_obs[fg$join_key == "NS-IBTS_1990_1_XX_S1_1"]
  # correct behaviour: NA (or the haul flagged / excluded), never 1500 g
  expect_true(is.na(b1) || length(b1) == 0L)
})

# --- Fishing pressure: ln(hours + 1), within-between decomposition, phases ------

fp_fixture <- function() {
  # rectangle A: 1990 (2 hauls, 0 h), 1991 (1 haul, e^2 - 1 h); rectangle B: 1990 (1 haul, e - 1 h)
  haul <- tibble(stat_rec = c("38F0", "38F0", "38F0", "39F0"), year = c(1990L, 1990L, 1991L, 1990L),
                 residual = c(-0.5, -0.6, -0.7, -0.8))
  panel <- tibble(stat_rec = c("38F0", "39F0"), mean_ln_B_obs = 10, rect_lon = c(0.5, 0.5),
                  rect_lat = c(54.75, 55.25))
  couce <- tibble(stat_rec = c("38F0", "38F0", "39F0"), year = c(1990L, 1991L, 1990L),
                  hours_total = c(0, exp(2) - 1, exp(1) - 1))
  suppressPackageStartupMessages(library(glmmTMB)) # numFactor() used in build_feasibility_data
  build_feasibility_data(haul, panel, couce, 1985L, 2015L)$data
}

test_that("fishing pressure is the natural log of (hours + 1)", {
  d <- fp_fixture()
  # ln(0 + 1) = 0; ln(e^2 - 1 + 1) = 2; ln(e - 1 + 1) = 1
  expect_equal(sort(unique(d$log_hours_total)), c(0, 1, 2))
})

test_that("FP_between is the haul-weighted rectangle mean; FP_within averages to 0 over each rectangle's hauls", {
  d <- add_fp_within_between(fp_fixture())
  a <- d[d$stat_rec == "38F0", ]
  # rectangle A: hauls carry 0, 0, 2 -> FP_between = 2/3 (haul-weighted, per decision 2)
  expect_equal(unique(a$FP_between), 2 / 3)
  expect_equal(sort(a$FP_within), c(-2 / 3, -2 / 3, 4 / 3))
  by_rect <- tapply(d$FP_within, d$stat_rec, mean)
  expect_equal(as.numeric(by_rect), c(0, 0))
  expect_equal(d$FP_between + d$FP_within, d$log_hours_total)
  expect_true(all(tapply(d$FP_between, d$stat_rec, function(x) length(unique(x))) == 1L))
})

test_that("phase_v2 boundaries: 1985-91, 1992-2001, 2002-07, 2008-15 inclusive", {
  script <- pipeline_script("run_h2h3_phase_v2_refit.R")
  exprs <- as.list(parse(script, keep.source = FALSE))
  hit <- Filter(function(e) is.call(e) && identical(as.character(e[[1]]), "<-") &&
                  identical(deparse(e[[2]]), "dat$phase_v2"), exprs)
  expect_length(hit, 1L)
  env <- new.env()
  env$dat <- data.frame(year = 1985:2015)
  eval(hit[[1]], env)
  ph <- setNames(as.character(env$dat$phase_v2), 1985:2015)
  expect_equal(unname(ph[c("1985", "1991")]), rep("1985-1991", 2))
  expect_equal(unname(ph[c("1992", "2001")]), rep("1992-2001", 2))
  expect_equal(unname(ph[c("2002", "2007")]), rep("2002-2007", 2))
  expect_equal(unname(ph[c("2008", "2015")]), rep("2008-2015", 2))
})

# --- Metrics -------------------------------------------------------------------

test_that("log_r2 hand value for a constant x4 offset (corrects test-h1_metrics.R:11)", {
  x <- log(c(100, 1000, 10000))
  # SS_tot = 2 (ln 10)^2 = 10.6038; SS_res = 3 (ln 4)^2 = 5.7654; R2 = 1 - 5.7654/10.6038 = 0.4563
  expect_equal(log_r2(x, x + log(4)), 1 - 3 * log(4)^2 / (2 * log(10)^2))
  expect_equal(log_r2(x, x + log(4)), 0.456285650263, tolerance = 1e-10)
  # a x100 offset is worse than predicting the mean: 1 - 3 (ln 100)^2 / (2 (ln 10)^2) = -5
  expect_equal(log_r2(x, x + log(100)), -5)
})

test_that("Nakagawa R2: intercept-only model has marginal R2 = 0; marginal <= conditional <= 1", {
  set.seed(1)
  d <- data.frame(stat_rec = factor(rep(1:20, each = 10)), x = rnorm(200))
  d$y <- rep(rnorm(20), each = 10) + 0.5 * d$x + rnorm(200, sd = 0.5)
  source_pipeline("h2h3_results_helpers.R")
  f0 <- glmmTMB::glmmTMB(y ~ 1 + (1 | stat_rec), data = d)
  f1 <- glmmTMB::glmmTMB(y ~ x + (1 | stat_rec), data = d)
  r0 <- nakagawa_r2_glmmtmb(f0)
  r1 <- nakagawa_r2_glmmtmb(f1)
  expect_equal(r0$r2_marginal, 0)
  expect_true(r1$r2_marginal > 0 && r1$r2_marginal <= r1$r2_conditional && r1$r2_conditional <= 1)
  # components add up: var_f / (var_f + var_r + var_e)
  expect_equal(r1$r2_marginal, r1$var_fixed / (r1$var_fixed + r1$var_random + r1$var_residual))
})

test_that("Bonferroni: m = 4 per family, alpha = 0.05 (alpha_adj = 0.0125)", {
  script <- pipeline_script("run_h2h3_h2_multiplicity_subsampling.R")
  env <- new.env()
  extract_definition(script, "M_FAMILY", env)
  extract_definition(script, "ALPHA", env)
  expect_equal(env$M_FAMILY, 4L)
  expect_equal(env$ALPHA, 0.05)
  extract_definition(script, "build_multiplicity_table", env)
  body_txt <- paste(deparse(body(env$build_multiplicity_table)), collapse = " ")
  expect_match(body_txt, "pmin(1, raw_p * M_FAMILY)", fixed = TRUE)
})

test_that("permutation p = (1 + #{|null| >= |obs|}) / (n + 1), non-finite nulls dropped", {
  env <- new.env()
  extract_definition(pipeline_script("permutation_bootstrap_FP_between_CAR.R"), "empirical_p", env)
  # |null| = 3, 1, 0, 1, 2 vs |obs| = 2 -> 2 exceed or tie -> (2 + 1) / (5 + 1) = 0.5
  expect_equal(env$empirical_p(c(-3, -1, 0, 1, 2), 2), 0.5)
  # NA and Inf replicates (failed fits) are dropped from both numerator and n
  expect_equal(env$empirical_p(c(-3, -1, 0, 1, 2, NA, Inf), 2), 0.5)
  expect_equal(env$empirical_p(c(0.1, 0.2), 5), 1 / 3)
})

test_that("permutation reassigns FP_between among rectangles, leaving everything else unchanged", {
  env <- new.env()
  extract_definition(pipeline_script("permutation_bootstrap_FP_between_CAR.R"), "permute_fp_between_data", env)
  d <- data.frame(stat_rec = rep(c("A", "B", "C"), times = c(3, 2, 4)),
                  FP_between = rep(c(1, 2, 3), times = c(3, 2, 4)), FP_within = seq(-1, 1, length.out = 9),
                  residual = 1:9)
  set.seed(3)
  p <- env$permute_fp_between_data(d)
  expect_true(all(tapply(p$FP_between, p$stat_rec, function(x) length(unique(x))) == 1L))
  expect_setequal(tapply(p$FP_between, p$stat_rec, unique), c(1, 2, 3))
  expect_identical(p[, c("stat_rec", "FP_within", "residual")], d[, c("stat_rec", "FP_within", "residual")])
})

# --- Residual specification (decision 1: absolute residual) --------------------

test_that("[F-A5] H2/H3 primary model response is |ln B_obs - ln B_pred| (model object)", {
  v2 <- cached_rds("outputs", "primary_model_v2.rds")
  d <- v2$data
  y <- d[[all.vars(formula(v2$primary_model_v2))[1]]]
  expect_equal(y, abs(log(d$B_obs) - log(d$B_pred)))
})

test_that("[F-A5] H2 CAR and H3 RE formulas in the live scripts use the absolute residual", {
  lhs_of <- function(script, var) {
    env <- new.env()
    f <- extract_definition(pipeline_script(script), var, env)
    all.vars(f)[1]
  }
  expect_identical(lhs_of("run_h2h3_phase_v2_reporting.R", "formula_car"), "abs_residual")
  expect_identical(lhs_of("run_h2h3_phase_v2_refit.R", "formula_v2"), "abs_residual")
})
