# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H2 corollary: rectangle-level biomass vs. fishing-pressure correlation test
#
# Standalone check of whether mean rectangle-level log observed biomass is
# spatially correlated with Couce et al. (2020) fishing pressure, independent
# of the H2 EEoS-residual analysis. This is a bivariate test only (no
# covariates). It does NOT touch EEoS residuals and does NOT modify or
# overwrite any existing H2 model/panel output.
#
# Run after build_h2_rectangle_panel.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_biomass_fishing_pressure_test.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_panel_helpers.R"))
source(file.path(script_dir, "R", "h2_model_helpers.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h2_robustness_helpers.R"))
source(file.path(script_dir, "R", "h2_provenance_helpers.R"))

paths <- h2_output_paths(project_root)
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_couce <- paths$couce_rectangle

stopifnot(
  file.exists(path_haul),
  file.exists(path_couce),
  file.exists(paths$panel)
)

out_panel_rds <- file.path(project_root, "outputs", "h2_biomass_fishing_panel.rds")
out_panel_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_panel.csv")
out_ols_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_ols.csv")
out_lm_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_lm_tests.csv")
out_spatial_csv <- file.path(project_root, "outputs", "h2_biomass_fishing_spatial_model.csv")
out_summary_rds <- file.path(project_root, "outputs", "h2_biomass_fishing_test_summary.rds")
out_notes_txt <- file.path(project_root, "outputs", "h2_biomass_fishing_notes.txt")

notes <- character(0)
add_note <- function(x) notes <<- c(notes, x)

haul <- readRDS(path_haul)
couce_rectangle <- readRDS(path_couce)
h2_main_panel <- readRDS(paths$panel)

# ---------------------------------------------------------------------------
# Step 1 — rectangle-level biomass aggregation (reuse, not duplicate, H2 logic)
# ---------------------------------------------------------------------------
# build_h2_rectangle_residuals() is NOT parameterised by a single input column:
# it hardcodes one summarise() call per rectangle that jointly produces
# mean_residual, mean_abs_residual, AND mean_ln_B_obs (mean of log-transformed
# haul-level biomass) using identical row filters (stat_rec normalisation,
# H2_YEAR_MIN/MAX window, finite residual/abs_residual). Because mean_ln_B_obs
# is already a byproduct of that exact call, this task reuses that column
# directly rather than adding a column argument (the function accepts none)
# or duplicating the aggregation logic in a new function. No refactor was
# performed or required.
add_note(paste(
  "Aggregation: build_h2_rectangle_residuals() is not parameterised by input",
  "column (it hardcodes a single summarise() computing mean_residual,",
  "mean_abs_residual, AND mean_ln_B_obs together in one pass). mean_ln_B_obs is",
  "already a byproduct of that exact call, so it was reused directly rather than",
  "calling the function with a column argument or duplicating aggregation code.",
  "No refactor was performed or was necessary; flagged explicitly per task",
  "instructions rather than silently reusing it without comment."
))

residual_panel <- build_h2_rectangle_residuals(haul) # identical call to H2 main pipeline
n_rect_hauls_window <- nrow(residual_panel)

joined <- residual_panel %>% left_join(couce_rectangle, by = "stat_rec")
n_missing_couce <- sum(is.na(joined$mean_annual_hours_total))
n_missing_couce_min_hauls <- sum(is.na(joined$mean_annual_hours_total) & joined$n_hauls >= H2_MIN_HAULS_DEFAULT)
n_zero_couce <- sum(!is.na(joined$mean_annual_hours_total) & joined$mean_annual_hours_total == 0)

add_note(sprintf(
  paste0(
    "Structural exclusions: of %d rectangles with >=1 NS-IBTS Q1 haul in %d-%d, ",
    "%d lack Couce et al. (2020) coverage (NA mean_annual_hours_total after the ",
    "rectangle-level join) unconditionally; %d of those also meet the H2 primary ",
    "min_hauls >= %d threshold (this %d figure matches the '~23 rectangles' ",
    "Couce-coverage exclusion already documented for the main H2 analysis in ",
    "display_discussion/H2_methods_draft.md). RESOLVED (2026-07-21): the task ",
    "briefing's '26 rectangles lacking Couce coverage' figure was traced to ",
    "outputs/step0_exclusion_counts.csv / step0_robustness_check.md / ",
    "step0_robustness_followups.md, i.e. the separate Step 0 (Berger-Parker ",
    "dominance D vs fishing pressure) cross-sectional analysis, which uses a ",
    "different, looser provisional haul threshold (STEP0_MIN_HAULS_PROVISIONAL = 5, ",
    "not H2's default min_hauls = 10). Under that 5-haul denominator (187 ",
    "qualifying rectangles), 26 lack Couce coverage; under H2's 10-haul denominator ",
    "(181 qualifying rectangles), 23 lack Couce coverage. The two counts are both ",
    "internally correct and differ by exactly 3 rectangles (28E9, 40G2, 46G1; 6-8 ",
    "hauls each) that clear the Step 0 5-haul bar but not H2's 10-haul bar and also ",
    "lack Couce coverage. The briefing's '26' and its Skagerrak/Kattegat + eastern ",
    "English Channel geographic description match step0_robustness_followups.md's ",
    "descriptive spatial-clustering table verbatim, confirming it was carried over ",
    "from that Step 0 follow-up doc rather than computed from the H2 primary ",
    "(min_hauls = 10) pipeline used here. This does not change the panel used in ",
    "this test: the final rectangle count below is the quantity that was gated ",
    "against the main H2 panel, and it matches exactly. ",
    "%d rectangles had non-missing but exactly zero mean_annual_hours_total after the ",
    "join; these were retained as real zero-effort observations, not imputed or excluded."
  ),
  n_rect_hauls_window, H2_YEAR_MIN, H2_YEAR_MAX,
  n_missing_couce, n_missing_couce_min_hauls, H2_MIN_HAULS_DEFAULT, n_missing_couce_min_hauls,
  n_zero_couce
))

# ---------------------------------------------------------------------------
# Step 2 — apply identical H2 inclusion rule (min_hauls, require_fishing)
# ---------------------------------------------------------------------------
panel <- build_h2_analysis_panel(
  residual_panel,
  couce_rectangle,
  min_hauls = H2_MIN_HAULS_DEFAULT,
  require_fishing = TRUE
)

n_match <- nrow(panel) == nrow(h2_main_panel) && setequal(panel$stat_rec, h2_main_panel$stat_rec)
if (!n_match) {
  stop(
    "STOP: rectangle set for the biomass-vs-fishing panel (n=", nrow(panel),
    ") does not match the saved H2 main-analysis panel (n=", nrow(h2_main_panel),
    "). Flagging discrepancy rather than proceeding, per task instructions."
  )
}
add_note(sprintf(
  "Final rectangle count = %d, identical rectangle set to the H2 main-analysis panel (min_hauls = %d, %s).",
  nrow(panel), H2_MIN_HAULS_DEFAULT, paths$panel
))

if (any(!is.finite(panel$mean_ln_B_obs))) {
  stop("STOP: non-finite mean_ln_B_obs values in analysis panel; stopping rather than imputing.")
}

# ---------------------------------------------------------------------------
# Cleaned rectangle-level dataset for independent verification
# ---------------------------------------------------------------------------
biomass_panel <- panel %>%
  select(stat_rec, mean_ln_B_obs, mean_annual_hours_total, n_hauls)

skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  mean((x - m)^3) / s^3
}
log_biomass_skew <- skewness(biomass_panel$mean_ln_B_obs)
raw_biomass_skew <- skewness(exp(biomass_panel$mean_ln_B_obs))
add_note(sprintf(
  paste0(
    "Distribution: skewness of mean_ln_B_obs (log scale, as specified) = %.3f; ",
    "skewness of the back-transformed raw scale (exp(mean_ln_B_obs)) = %.3f. ",
    "Reported for information only; no transformation beyond the log already ",
    "specified was applied."
  ),
  log_biomass_skew, raw_biomass_skew
))

saveRDS(biomass_panel, out_panel_rds)
write_csv(biomass_panel, out_panel_csv)

# ---------------------------------------------------------------------------
# 4a — Simple OLS: mean_log_biomass ~ fishing_hours (bivariate, no covariates)
# ---------------------------------------------------------------------------
BIOMASS_FISHING_FORMULA <- mean_ln_B_obs ~ mean_annual_hours_total

ols_fit <- lm(BIOMASS_FISHING_FORMULA, data = panel)
ols_row <- fit_h2_ols(
  panel,
  BIOMASS_FISHING_FORMULA,
  "biomass_fishing_ols",
  term = "mean_annual_hours_total"
)
write_csv(ols_row, out_ols_csv)

# ---------------------------------------------------------------------------
# 4b — Spatial weights + LM specification tests
# ---------------------------------------------------------------------------
# No standalone spatial-weights object is persisted anywhere in outputs/ for
# H2: run_h2_models.R, run_h2_sar_lag_models.R, and run_h2_robustness_checks.R
# each independently call build_h2_spatial_weights() on the min10 panel and
# the ICES rectangle shapefile, rather than loading a cached listw object.
# This task follows that same convention. Because the biomass-fishing panel
# has the identical rectangle set (n and stat_rec) as the H2 main-analysis
# panel (verified above), build_h2_spatial_weights() on that panel and the
# same shapefile reconstructs a weights matrix that is deterministic and
# therefore equivalent to the one underlying the primary H2 SEM (queen
# contiguity, row-standardised style "W", zero.policy = TRUE). No missingness
# mismatch arose: mean_ln_B_obs is non-missing for exactly the 158 rectangles
# with non-missing residual data, because both are byproducts of the same
# build_h2_rectangle_residuals() call.
rectangles_sf <- load_ices_rectangles_sf(project_root)
weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)

add_note(sprintf(
  paste0(
    "Spatial weights: queen contiguity, row-standardised (style = 'W'), ",
    "zero.policy = TRUE, n = %d, n_isolated = %d. No cached weights object exists ",
    "in outputs/ for H2 (each H2 script rebuilds it via build_h2_spatial_weights()); ",
    "this task reconstructed it the same way, on the identical panel/rectangle set ",
    "used for the primary H2 SEM, so it is equivalent rather than loaded verbatim ",
    "from a saved object."
  ),
  nrow(panel), weights$n_isolated
))

lm_tests <- h2_lm_spatial_tests(ols_fit, weights$listw)
write_csv(lm_tests, out_lm_csv)

# ---------------------------------------------------------------------------
# 4c — Spatial model matching the LM-test result (SEM or SAR, not assumed)
# ---------------------------------------------------------------------------
rlm_err_p <- lm_tests$p_value[lm_tests$test == "RLM_error"]
rlm_lag_p <- lm_tests$p_value[lm_tests$test == "RLM_lag"]
rlm_err_stat <- lm_tests$statistic[lm_tests$test == "RLM_error"]
rlm_lag_stat <- lm_tests$statistic[lm_tests$test == "RLM_lag"]

alpha <- 0.05
spatial_spec <- if (length(rlm_err_p) == 1L && length(rlm_lag_p) == 1L) {
  if (rlm_err_p < alpha && rlm_lag_p >= alpha) {
    "sem"
  } else if (rlm_lag_p < alpha && rlm_err_p >= alpha) {
    "sar_lag"
  } else if (rlm_err_p < alpha && rlm_lag_p < alpha) {
    if (rlm_err_stat >= rlm_lag_stat) "sem" else "sar_lag"
  } else {
    "sem" # neither robust LM test significant: default to H2 SEM convention
  }
} else {
  "sem"
}

add_note(sprintf(
  paste0(
    "LM decision rule: robust LM-error p = %s, robust LM-lag p = %s (alpha = 0.05). ",
    "Decision: fit %s."
  ),
  ifelse(length(rlm_err_p) == 1L, signif(rlm_err_p, 4), "NA"),
  ifelse(length(rlm_lag_p) == 1L, signif(rlm_lag_p, 4), "NA"),
  toupper(spatial_spec)
))

if (spatial_spec == "sem") {
  spatial_results <- fit_h2_sem(
    panel, weights$listw,
    formula = BIOMASS_FISHING_FORMULA,
    model_id = "biomass_fishing_sem"
  )
  add_note(paste(
    "Spatial model: LM tests favoured (or did not clearly reject) a spatial ERROR",
    "process, consistent with the H2 main-analysis convention. A spatial error model",
    "(SEM, spatialreg::errorsarlm()) was fit with the same specification convention",
    "used in H2."
  ))
} else {
  spatial_results <- fit_h2_sar_lag(
    panel, weights$listw,
    formula = BIOMASS_FISHING_FORMULA,
    model_id = "biomass_fishing_sar_lag"
  )
  add_note(paste(
    "Spatial model: LM tests favoured a spatial LAG process (robust LM-lag significant,",
    "robust LM-error not significant, or with the larger robust statistic) for this",
    "bivariate biomass-vs-fishing specification. A spatial lag model (SAR/SLM,",
    "spatialreg::lagsarlm()) was fit instead of a spatial error model. This DEPARTS",
    "from the H2 main-analysis convention (SEM fit a priori as the primary model",
    "there); flagged explicitly rather than silently substituting SEM."
  ))
}
write_csv(spatial_results, out_spatial_csv)

# ---------------------------------------------------------------------------
# Save summary bundle + provenance + human-readable notes
# ---------------------------------------------------------------------------
provenance <- h2_collect_provenance(project_root)
summary_out <- h2_stamp_result(
  list(
    n_rectangles = nrow(panel),
    n_rectangles_h2_main = nrow(h2_main_panel),
    n_missing_couce = n_missing_couce,
    n_missing_couce_min_hauls = n_missing_couce_min_hauls,
    n_zero_couce = n_zero_couce,
    log_biomass_skew = log_biomass_skew,
    raw_biomass_skew = raw_biomass_skew,
    biomass_panel = biomass_panel,
    ols_results = ols_row,
    lm_tests = lm_tests,
    spatial_spec = spatial_spec,
    spatial_results = spatial_results,
    weights_info = data.frame(
      n_rectangles = nrow(panel),
      n_isolated = weights$n_isolated,
      nb_type = "queen",
      nb_style = "W",
      zero_policy = TRUE,
      stringsAsFactors = FALSE
    ),
    notes = notes
  ),
  provenance
)
saveRDS(summary_out, out_summary_rds)
writeLines(notes, out_notes_txt)

cat("=== H2 corollary: biomass vs fishing-pressure test ===\n")
cat("Final rectangle count:", nrow(panel), "(H2 main panel:", nrow(h2_main_panel), ")\n")
cat("Rectangles excluded for missing Couce coverage (unconditional):", n_missing_couce, "\n")
cat("Rectangles excluded for missing Couce coverage (n_hauls >= min threshold):", n_missing_couce_min_hauls, "\n")
cat("Rectangles with zero (non-missing) Couce hours:", n_zero_couce, "\n")
cat("\n--- OLS: mean_ln_B_obs ~ mean_annual_hours_total ---\n")
print(ols_row)
cat("\n--- LM specification tests ---\n")
print(lm_tests)
cat("\nSpatial model fit:", toupper(spatial_spec), "\n")
print(spatial_results)
cat("\nSaved", out_panel_csv, "\n")
cat("Saved", out_ols_csv, "\n")
cat("Saved", out_lm_csv, "\n")
cat("Saved", out_spatial_csv, "\n")
cat("Saved", out_summary_rds, "\n")
cat("Saved", out_notes_txt, "\n")
