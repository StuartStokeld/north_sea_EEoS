# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H2/H3 shared hierarchical model — RESULTS RUN (finalised design)
#
# PURPOSE: produce the substantive H2/H3 results under the corrected model
# design. THIS IS A RESULTS RUN, NOT A FEASIBILITY CHECK. Interpretation of
# fishing-pressure and fishing-pressure x phase coefficients as the actual
# H2/H3 findings is IN SCOPE.
#
# PRIMARY SPECIFICATION (corrected — no biomass covariate):
#   residual ~ log_hours_total * phase + (1 | stat_rec)
#   [glmmTMB, REML; plain random intercept — no spatial correlation]
# Outcome = canonical pipeline residual (= log(B_obs) - log(B_pred)).
# phase = 4-level factor from breaks 1989, 2001, 2008 (marginal 1997
# excluded) — see build_phase_factor() in R/h2h3_feasibility_helpers.R.
#
# SENSITIVITY SPECIFICATIONS (report alongside, not as headline):
#   1. CAR adjacency: residual ~ log_hours_total * phase +
#      adjacency(1 | stat_rec)  [spaMM REML]
#   2. Varying-coefficient GAM: residual ~ s(year, k=8) +
#      s(year, by = log_hours_total, k=8) + s(stat_rec, bs = "re")
#      [mgcv REML]
#
# REUSE / REFIT: feasibility / temporal-robustness RDS fits included
# mean_ln_B_obs and therefore CANNOT be reused for this corrected spec.
# This script always refits all three models without biomass. The Round 2
# adjacency matrix is still reused (neighbour definition unchanged).
# Data prep still uses build_feasibility_data() (same 158-rectangle panel).
#
# MUST be run with `Rscript --vanilla` (renv not activated) — same
# environment note as the feasibility scripts (glmmTMB / spaMM ad hoc).
#
# OUT OF SCOPE: no further model-specification changes beyond the biomass
# removal requested here. If a problem requiring a further specification
# change is found, flag it back rather than resolving it unilaterally.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_shared_model_results.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' is required and is not on the library path. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_shared_model_results.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))
source(file.path(script_dir, "R", "h2h3_temporal_robustness_helpers.R"))
source(file.path(script_dir, "R", "h2h3_results_helpers.R"))

has_spamm <- requireNamespace("spaMM", quietly = TRUE)
has_mgcv <- requireNamespace("mgcv", quietly = TRUE)
has_dharma <- requireNamespace("DHARMa", quietly = TRUE)
if (!has_spamm) {
  stop("Package 'spaMM' is required to load/report the CAR sensitivity model. Run with --vanilla.")
}
if (!has_mgcv) {
  stop("Package 'mgcv' is required to load/report the GAM sensitivity model.")
}
suppressPackageStartupMessages(library(spaMM))
suppressPackageStartupMessages(library(mgcv))

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_panel <- file.path(project_root, "outputs", "h2_rectangle_panel.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
path_round1 <- file.path(project_root, "outputs", "h2h3_feasibility_model_objects.rds")
path_round2 <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
path_temporal <- file.path(project_root, "outputs", "h2h3_temporal_robustness_model_objects.rds")
stopifnot(
  file.exists(path_haul), file.exists(path_panel), file.exists(path_couce_year),
  file.exists(path_round1), file.exists(path_round2), file.exists(path_temporal)
)

path_out_primary_fe <- file.path(project_root, "outputs", "h2h3_results_primary_fixed_effects.csv")
path_out_sensitivity <- file.path(project_root, "outputs", "h2h3_results_sensitivity_comparison.csv")
path_out_phase_slopes <- file.path(project_root, "outputs", "h2h3_results_fp_slopes_by_phase.csv")
path_out_year_slopes <- file.path(project_root, "outputs", "h2h3_results_fp_slopes_by_year.csv")
path_out_pooling <- file.path(project_root, "outputs", "h2h3_results_partial_pooling.csv")
path_out_r2 <- file.path(project_root, "outputs", "h2h3_results_model_fit.csv")
path_out_influence <- file.path(project_root, "outputs", "h2h3_results_influence_flags.csv")
path_out_wald <- file.path(project_root, "outputs", "h2h3_results_wald_tests.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_results_model_objects.rds")
path_out_fig_effect <- file.path(fig_dir, "h2h3_results_fp_effect_by_phase.png")
path_out_fig_pooling <- file.path(fig_dir, "h2h3_results_partial_pooling.png")
path_out_fig_resid <- file.path(fig_dir, "h2h3_results_residual_diagnostics.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_results_run_log.md")
path_out_session <- file.path(project_root, "outputs", "h2h3_results_sessionInfo.txt")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 shared hierarchical model — RESULTS run log")
logmsg("")
logmsg(
  "RESULTS RUN (not a feasibility check). Fishing-pressure and fishing-pressure x phase ",
  "coefficients ARE interpreted below as the substantive H2 and H3 findings, within the ",
  "limits of what the coefficients directly support."
)
logmsg("")
logmsg(
  "CORRECTION / RE-RUN: biomass (`mean_ln_B_obs`) has been REMOVED from the primary model and ",
  "both sensitivity models. Corrected fixed effects = `log_hours_total * phase` only. This ",
  "re-run overwrites the previous with-biomass results outputs; a before/after comparison of ",
  "headline numbers is logged below."
)

# Snapshot previous (with-biomass) headline numbers before overwrite, if present.
prev_slopes <- if (file.exists(path_out_phase_slopes)) {
  tryCatch(read_csv(path_out_phase_slopes, show_col_types = FALSE), error = function(e) NULL)
} else {
  NULL
}
prev_fe <- if (file.exists(path_out_primary_fe)) {
  tryCatch(read_csv(path_out_primary_fe, show_col_types = FALSE), error = function(e) NULL)
} else {
  NULL
}
prev_fit <- if (file.exists(path_out_r2)) {
  tryCatch(read_csv(path_out_r2, show_col_types = FALSE), error = function(e) NULL)
} else {
  NULL
}
prev_had_biomass <- !is.null(prev_fe) && "mean_ln_B_obs" %in% prev_fe$term

# Hardcoded snapshot from the completed with-biomass results run (2026-07-26),
# used if on-disk previous outputs were already overwritten by a no-biomass re-run.
prev_biomass_fallback <- list(
  slopes = data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    fp_slope_prev = c(-0.014865463181557114, -0.025053167770713196,
                      0.072357533894124, 0.0429907719229194),
    stringsAsFactors = FALSE
  ),
  interactions = data.frame(
    term = c(
      "log_hours_total:phase1989-2000",
      "log_hours_total:phase2001-2007",
      "log_hours_total:phase2008-2015"
    ),
    estimate_prev = c(-0.010187704589156082, 0.08722299707568111, 0.05785623510447651),
    stringsAsFactors = FALSE
  ),
  r2_marginal = 0.10420402314295504,
  r2_conditional = 0.1630348389364165,
  cor_primary_car = 1.000,
  cor_primary_gam = 0.989
)

# ---------------------------------------------------------------------------
# Session / package versions
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session and package versions")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("Full sessionInfo() written to: ", path_out_session)
logmsg(sprintf(
  "Key packages: R %s; glmmTMB %s; spaMM %s; mgcv %s; DHARMa %s (available = %s).",
  getRversion(),
  as.character(utils::packageVersion("glmmTMB")),
  as.character(utils::packageVersion("spaMM")),
  as.character(utils::packageVersion("mgcv")),
  if (has_dharma) as.character(utils::packageVersion("DHARMa")) else "NA",
  has_dharma
))
logmsg(
  "Environment note (unchanged from feasibility): glmmTMB and spaMM were installed ad hoc into ",
  "the ambient/user library, NOT added to renv.lock. This script MUST be run with ",
  "`Rscript --vanilla pipeline/run_h2h3_shared_model_results.R`."
)

# ---------------------------------------------------------------------------
# Data: rebuild via the same feasibility pipeline; validate against stored
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data and analysis universe")

haul <- readRDS(path_haul)
panel <- readRDS(path_panel)
couce_year <- readRDS(path_couce_year)
built <- build_feasibility_data(haul, panel, couce_year, H2_YEAR_MIN, H2_YEAR_MAX)
dat <- built$data
dat$stat_rec <- factor(dat$stat_rec, levels = panel$stat_rec)

logmsg(sprintf(
  paste0(
    "Rebuilt analysis data via build_feasibility_data() (identical pipeline to feasibility checks): ",
    "%d hauls across %d rectangles, years %d–%d. Universe = established H2 panel ",
    "(>= %d hauls AND Couce coverage) = %d rectangles. Dropped %d hauls lacking rectangle-year Couce."
  ),
  built$n_hauls_final, dplyr::n_distinct(dat$stat_rec), min(dat$year), max(dat$year),
  H2_MIN_HAULS_DEFAULT, nrow(panel), built$n_hauls_dropped_no_rect_year_couce
))

phase_table <- dat %>%
  group_by(phase) %>%
  summarise(
    year_min = min(year), year_max = max(year),
    n_hauls = dplyr::n(), n_rectangles = dplyr::n_distinct(stat_rec),
    .groups = "drop"
  )
logmsg("Phase definition (break year = first year of new phase; 1997 excluded):")
for (i in seq_len(nrow(phase_table))) {
  r <- phase_table[i, ]
  logmsg(sprintf(
    "  - %s (years %d–%d): n_hauls = %d, n_rectangles = %d",
    as.character(r$phase), r$year_min, r$year_max, r$n_hauls, r$n_rectangles
  ))
}

data_checksum <- function(d) {
  c(
    n = nrow(d),
    n_rect = dplyr::n_distinct(d$stat_rec),
    sum_resid = sum(d$residual),
    sum_logh = sum(d$log_hours_total),
    sum_biom = sum(d$mean_ln_B_obs),
    sum_year = sum(as.numeric(d$year))
  )
}
checksum_new <- data_checksum(dat)

validate_stored_data <- function(stored, label) {
  if (is.null(stored)) {
    return(list(ok = FALSE, detail = paste0(label, ": no stored data")))
  }
  cs <- data_checksum(stored)
  ok <- isTRUE(all.equal(cs, checksum_new, tolerance = 1e-8))
  list(
    ok = ok,
    detail = sprintf(
      "%s stored vs rebuilt checksum: n=%s/%s, n_rect=%s/%s, sum_resid match=%s",
      label, cs["n"], checksum_new["n"], cs["n_rect"], checksum_new["n_rect"],
      isTRUE(all.equal(cs["sum_resid"], checksum_new["sum_resid"], tolerance = 1e-8))
    )
  )
}

# ---------------------------------------------------------------------------
# Fit models (corrected: no biomass). Prior feasibility fits cannot be reused.
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fitting (corrected specification — no biomass)")

round2 <- readRDS(path_round2)
v_round2 <- validate_stored_data(round2$data, "Round 2 data")
logmsg(v_round2$detail)
logmsg(
  "Prior feasibility / temporal-robustness fitted objects included mean_ln_B_obs and are ",
  "therefore not reused. All three models are REFIT under the corrected formula. Round 2 ",
  "adjMatrix is still reused (same queen-adjacency neighbour definition)."
)

formula_primary <- residual ~ log_hours_total * phase + (1 | stat_rec)
formula_car <- residual ~ log_hours_total * phase + adjacency(1 | stat_rec)
formula_gam <- residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) +
  s(stat_rec, bs = "re")

logmsg("Corrected formulas used for this re-run:")
logmsg("  - Primary: residual ~ log_hours_total * phase + (1 | stat_rec)  [glmmTMB REML]")
logmsg("  - CAR:     residual ~ log_hours_total * phase + adjacency(1 | stat_rec)  [spaMM REML]")
logmsg("  - GAM:     residual ~ s(year, k=8) + s(year, by=log_hours_total, k=8) + s(stat_rec, bs=\"re\")  [mgcv REML]")

logmsg("Fitting primary model ...")
time_primary <- system.time({
  fit_primary <- glmmTMB(formula_primary, data = dat, REML = TRUE)
})
logmsg(sprintf("Primary fit time: %.2f sec.", time_primary["elapsed"]))

if (!is.null(round2$adjMatrix)) {
  adjMatrix <- round2$adjMatrix
  logmsg("CAR adjacency matrix: REUSED from Round 2 RDS.")
} else {
  logmsg("CAR adjacency matrix missing from Round 2 RDS — rebuilding via established H2 weights path.")
  source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
  suppressPackageStartupMessages(library(spdep))
  rectangles_sf <- load_ices_rectangles_sf(project_root)
  weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)
  adjMatrix <- build_car_adjacency_matrix(weights$nb, panel)
}
dat_car <- dat
dat_car$stat_rec <- factor(dat_car$stat_rec, levels = rownames(adjMatrix))

logmsg("Fitting CAR sensitivity model ...")
time_car <- system.time({
  fit_car <- spaMM::fitme(formula_car, data = dat_car, adjMatrix = adjMatrix, method = "REML")
})
logmsg(sprintf("CAR fit time: %.2f sec.", time_car["elapsed"]))

logmsg("Fitting GAM sensitivity model (expect ~1–2 min) ...")
time_gam <- system.time({
  fit_gam <- mgcv::gam(formula_gam, data = dat, method = "REML")
})
logmsg(sprintf("GAM fit time: %.2f sec.", time_gam["elapsed"]))

# Confirm no biomass term slipped into any fitted formula
fml_check <- c(
  primary = paste(deparse(stats::formula(fit_primary), width.cutoff = 500L), collapse = " "),
  car = paste(deparse(stats::formula(fit_car), width.cutoff = 500L), collapse = " "),
  gam = paste(deparse(stats::formula(fit_gam), width.cutoff = 500L), collapse = " ")
)
biomass_leaked <- any(grepl("mean_ln_B_obs", fml_check))
logmsg("Fitted formulas:")
for (nm in names(fml_check)) logmsg("  - ", nm, ": ", fml_check[[nm]])
if (biomass_leaked) {
  stop("BIOMASS LEAK: mean_ln_B_obs present in at least one fitted formula after correction.")
}
logmsg("Confirmed: mean_ln_B_obs absent from all three fitted formulas.")

logmsg("")
logmsg("### Reuse summary")
logmsg("  - Primary (1|stat_rec): REFIT (biomass removed)")
logmsg("  - CAR adjacency model:  REFIT (biomass removed); adjMatrix REUSED from Round 2")
logmsg("  - GAM smooth year:      REFIT (biomass removed)")

logmsg("")
logmsg("## Specification note")
logmsg(
  "Biomass removed from primary and both sensitivities relative to the immediately preceding ",
  "results run. Otherwise unchanged: plain RE primary; CAR and varying-coefficient GAM ",
  "sensitivities; phase breaks 1989/2001/2008; canonical residual; 158-rectangle H2 panel."
)

saveRDS(
  list(
    fit_primary = fit_primary,
    fit_car = fit_car,
    fit_gam = fit_gam,
    data = dat,
    adjMatrix = adjMatrix,
    formula_primary = formula_primary,
    formula_car = formula_car,
    formula_gam = formula_gam,
    biomass_removed = TRUE,
    primary_reused = FALSE,
    car_reused = FALSE,
    gam_reused = FALSE,
    adjMatrix_reused = TRUE
  ),
  path_out_models
)
logmsg("Saved model objects: ", path_out_models)

# ---------------------------------------------------------------------------
# 1. Primary fixed-effects table (full, with H2/H3 labels)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Primary model — full fixed-effects table")

primary_fe <- tidy_fixed_effects(fit_primary, "primary_plain_re") %>%
  annotate_fixed_effects() %>%
  select(model_id, term, term_plain, hypothesis, estimate, std_error, statistic, p_value)
write_csv(primary_fe, path_out_primary_fe)

logmsg("Primary formula: residual ~ log_hours_total * phase + (1 | stat_rec) [REML; no biomass]")
logmsg("Fixed effects (estimate, SE, z, p) with H2/H3 labels:")
for (i in seq_len(nrow(primary_fe))) {
  r <- primary_fe[i, ]
  logmsg(sprintf(
    "  - %-40s  est = %+.4f  SE = %.4f  z = %+.3f  p = %.3e  [%s]",
    r$term, r$estimate, r$std_error, r$statistic, r$p_value, r$hypothesis
  ))
}
logmsg("Saved: ", path_out_primary_fe)

# ---------------------------------------------------------------------------
# Phase-specific FP slopes (all three models) + Wald tests
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Phase-specific fishing-pressure slopes (derived; used for H2 reporting)")

phase_slopes_primary <- extract_phase_fp_slopes(fit_primary) %>%
  mutate(model_id = "primary_plain_re")
phase_slopes_car <- extract_phase_fp_slopes_spamm(fit_car)

years_grid <- seq.int(min(dat$year), max(dat$year))
gam_year_slopes <- extract_gam_fp_slopes(fit_gam, years_grid, dat)
phase_slopes_gam <- summarise_gam_slopes_by_phase(gam_year_slopes)

phase_slopes_all <- bind_rows(
  phase_slopes_primary %>% select(model_id, phase, year_start, year_end, fp_slope, fp_slope_se, fp_slope_lo, fp_slope_hi) %>%
    mutate(note = "parametric phase-model slope (delta method)"),
  phase_slopes_car %>% select(model_id, phase, year_start, year_end, fp_slope, fp_slope_se, fp_slope_lo, fp_slope_hi) %>%
    mutate(note = "parametric CAR phase-model slope (delta method)"),
  phase_slopes_gam
)
write_csv(phase_slopes_all, path_out_phase_slopes)

logmsg("Per-phase FP slopes (primary):")
for (i in seq_len(nrow(phase_slopes_primary))) {
  r <- phase_slopes_primary[i, ]
  logmsg(sprintf(
    "  - %s: slope = %+.4f, SE = %.4f, 95%% CI = [%+.4f, %+.4f]",
    r$phase, r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi
  ))
}
logmsg("Saved: ", path_out_phase_slopes)

# Wald joint tests on primary model
b_pri <- glmmTMB::fixef(fit_primary)$cond
V_pri <- stats::vcov(fit_primary)$cond
h3_terms <- grep("^log_hours_total:phase", names(b_pri), value = TRUE)
h2_terms_all_slopes <- c("log_hours_total", h3_terms)
# Joint H3: all interactions = 0 (FP effect identical across phases)
wald_h3 <- wald_joint_zero(b_pri, V_pri, h3_terms)
# Joint H2-style: all phase-specific slopes = 0 is equivalent to
# log_hours_total = 0 AND all interactions = 0
wald_h2 <- wald_joint_zero(b_pri, V_pri, h2_terms_all_slopes)

wald_table <- data.frame(
  test_id = c("H3_interactions_joint_zero", "H2_all_phase_fp_slopes_joint_zero"),
  description = c(
    "Joint Wald: all log_hours_total:phase interactions = 0 (no phase change in FP effect)",
    "Joint Wald: reference FP slope and all interactions = 0 (FP effect zero in every phase)"
  ),
  statistic_chisq = c(wald_h3$statistic, wald_h2$statistic),
  df = c(wald_h3$df, wald_h2$df),
  p_value = c(wald_h3$p_value, wald_h2$p_value),
  terms = c(paste(wald_h3$terms, collapse = "; "), paste(wald_h2$terms, collapse = "; ")),
  stringsAsFactors = FALSE
)
write_csv(wald_table, path_out_wald)
logmsg(sprintf(
  "Joint Wald H3 (all FP x phase interactions = 0): chi2(%.0f) = %.3f, p = %.3e",
  wald_h3$df, wald_h3$statistic, wald_h3$p_value
))
logmsg(sprintf(
  "Joint Wald H2-style (FP slope = 0 in every phase): chi2(%.0f) = %.3f, p = %.3e",
  wald_h2$df, wald_h2$statistic, wald_h2$p_value
))
logmsg("Saved: ", path_out_wald)

# ---------------------------------------------------------------------------
# 2. Side-by-side sensitivity comparison
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Side-by-side sensitivity comparison")

fe_primary <- tidy_fixed_effects(fit_primary, "primary_plain_re")
fe_car <- tidy_fixed_effects_spamm(fit_car, "car_adjacency")
# GAM: parametric terms are intercept-only under the no-biomass varying-
# coefficient form; phase terms do not exist — handled via phase-slope table.
fe_gam_param <- tidy_gam_effects(fit_gam, "smooth_year_gam") %>%
  filter(term_type == "parametric") %>%
  select(model_id, term, estimate, std_error, statistic, p_value)

# Wide comparison for terms shared by primary and CAR
shared_terms <- intersect(fe_primary$term, fe_car$term)
cmp <- fe_primary %>%
  filter(term %in% shared_terms) %>%
  select(term, estimate, std_error, statistic, p_value) %>%
  rename(
    estimate_primary = estimate, se_primary = std_error,
    z_primary = statistic, p_primary = p_value
  ) %>%
  left_join(
    fe_car %>%
      filter(term %in% shared_terms) %>%
      select(term, estimate, std_error, statistic, p_value) %>%
      rename(
        estimate_car = estimate, se_car = std_error,
        z_car = statistic, p_car = p_value
      ),
    by = "term"
  ) %>%
  left_join(
    fe_gam_param %>%
      select(term, estimate, std_error, statistic, p_value) %>%
      rename(
        estimate_gam = estimate, se_gam = std_error,
        z_gam = statistic, p_gam = p_value
      ),
    by = "term"
  ) %>%
  annotate_fixed_effects() %>%
  select(
    term, term_plain, hypothesis,
    estimate_primary, se_primary, z_primary, p_primary,
    estimate_car, se_car, z_car, p_car,
    estimate_gam, se_gam, z_gam, p_gam
  )

# Append a block of phase-slope comparisons as extra rows with term = "fp_slope:<phase>"
slope_cmp <- phase_slopes_all %>%
  select(model_id, phase, fp_slope, fp_slope_se) %>%
  mutate(
    term = paste0("fp_slope:", phase),
    model_col = dplyr::recode(
      model_id,
      primary_plain_re = "primary",
      car_adjacency = "car",
      smooth_year_gam = "gam"
    )
  )
slope_wide <- slope_cmp %>%
  select(term, phase, model_col, fp_slope, fp_slope_se) %>%
  tidyr::pivot_wider(
    names_from = model_col,
    values_from = c(fp_slope, fp_slope_se)
  ) %>%
  mutate(
    term_plain = paste0("Fishing-pressure slope in phase ", phase),
    hypothesis = "H2 (phase-specific FP effect; GAM = mean of year slopes within phase)",
    estimate_primary = fp_slope_primary,
    se_primary = fp_slope_se_primary,
    z_primary = estimate_primary / se_primary,
    p_primary = 2 * stats::pnorm(-abs(z_primary)),
    estimate_car = fp_slope_car,
    se_car = fp_slope_se_car,
    z_car = estimate_car / se_car,
    p_car = 2 * stats::pnorm(-abs(z_car)),
    estimate_gam = fp_slope_gam,
    se_gam = fp_slope_se_gam,
    z_gam = estimate_gam / se_gam,
    p_gam = 2 * stats::pnorm(-abs(z_gam))
  ) %>%
  select(
    term, term_plain, hypothesis,
    estimate_primary, se_primary, z_primary, p_primary,
    estimate_car, se_car, z_car, p_car,
    estimate_gam, se_gam, z_gam, p_gam
  )

sensitivity_table <- bind_rows(
  cmp %>% mutate(block = "parametric_fixed_effects"),
  slope_wide %>% mutate(block = "phase_fp_slopes")
) %>%
  select(block, everything())

write_csv(sensitivity_table, path_out_sensitivity)
logmsg(
  "Comparison table has two blocks: (1) parametric fixed effects shared by primary and CAR ",
  "(GAM columns filled only where the GAM has a matching parametric term — typically ",
  "Intercept only under the no-biomass GAM; phase terms are NA for GAM by design); ",
  "(2) phase-specific FP slopes for all three models (GAM = mean of year-by-year by-smooth ",
  "within each phase)."
)
logmsg("Parametric block — primary vs CAR estimate shift:")
for (i in seq_len(nrow(cmp))) {
  r <- cmp[i, ]
  pct <- 100 * (r$estimate_car - r$estimate_primary) / abs(r$estimate_primary)
  logmsg(sprintf(
    "  - %s: primary %+.4f -> CAR %+.4f (%+.1f%%)",
    r$term, r$estimate_primary, r$estimate_car, pct
  ))
}
logmsg("Phase-slope block:")
for (i in seq_len(nrow(slope_wide))) {
  r <- slope_wide[i, ]
  logmsg(sprintf(
    "  - %s: primary %+.4f | CAR %+.4f | GAM(mean) %+.4f",
    r$term, r$estimate_primary, r$estimate_car, r$estimate_gam
  ))
}
logmsg("Saved: ", path_out_sensitivity)

# ---------------------------------------------------------------------------
# 3. Effect visualisation: phase steps + GAM smooth
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Effect visualisation: FP effect by phase (primary) + GAM overlay")

phase_year <- expand_phase_slopes_to_years(
  phase_slopes_primary %>% mutate(model_id = "categorical_phase"),
  years_grid
)
slopes_plot <- bind_rows(
  phase_year,
  gam_year_slopes
)
write_csv(slopes_plot, path_out_year_slopes)

plot_phase_seg <- phase_slopes_primary %>%
  mutate(model_label = "Primary (categorical phase)")
plot_gam <- gam_year_slopes %>%
  mutate(model_label = "Sensitivity: smooth year (GAM)")

col_map <- c(
  "Primary (categorical phase)" = "#d73027",
  "Sensitivity: smooth year (GAM)" = "#1a9850"
)

p_effect <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_rect(
    data = plot_phase_seg,
    aes(
      xmin = year_start - 0.5, xmax = year_end + 0.5,
      ymin = fp_slope_lo, ymax = fp_slope_hi,
      fill = model_label
    ),
    alpha = 0.15, colour = NA
  ) +
  geom_segment(
    data = plot_phase_seg,
    aes(
      x = year_start - 0.5, xend = year_end + 0.5,
      y = fp_slope, yend = fp_slope,
      colour = model_label
    ),
    linewidth = 1.15
  ) +
  geom_ribbon(
    data = plot_gam,
    aes(x = year, ymin = fp_slope_lo, ymax = fp_slope_hi, fill = model_label),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = plot_gam,
    aes(x = year, y = fp_slope, colour = model_label),
    linewidth = 1.0
  ) +
  geom_vline(xintercept = c(1989, 2001, 2008), linetype = "dotted", colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_fill_manual(values = col_map, name = NULL) +
  labs(
    x = "Year",
    y = "Fishing-pressure effect\n(∂ residual / ∂ log(hours+1))",
    title = "H2/H3 results: fishing-pressure effect by phase",
    subtitle = paste0(
      "Primary model (categorical phase, red) with 95% CI bands; GAM varying-coefficient ",
      "smooth overlaid (green). Dotted lines = structural breaks 1989, 2001, 2008."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.5, colour = "grey30")
  )
ggsave(path_out_fig_effect, p_effect, width = 10, height = 5.5, dpi = 150)
logmsg("Saved figure: ", path_out_fig_effect)
logmsg("Saved year-grid slopes: ", path_out_year_slopes)

# ---------------------------------------------------------------------------
# 4. Random-effects / partial-pooling summary
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 4. Random-effects summary (partial pooling)")

pooling <- extract_partial_pooling(fit_primary, dat, "primary_plain_re", spatial = FALSE)
write_csv(pooling, path_out_pooling)

cor_pool <- suppressWarnings(cor(log(pooling$n_hauls), pooling$shrinkage_ratio, method = "spearman"))
logmsg(sprintf(
  paste0(
    "Rectangle random-intercept SD and partial pooling for the primary model: Spearman cor ",
    "log(n_hauls) vs shrinkage_ratio = %.3f across %d rectangles (positive = more hauls, less ",
    "shrinkage — expected pattern)."
  ),
  cor_pool, nrow(pooling)
))
vc_pri <- extract_variance_components(fit_primary, "primary_plain_re", spatial = FALSE)
logmsg(sprintf(
  "Rectangle intercept SD = %.4f (95%% CI [%.4f, %.4f]); residual SD = %.4f.",
  vc_pri$estimate, vc_pri$ci_low, vc_pri$ci_high, vc_pri$residual_sd
))

p_pool <- ggplot(pooling, aes(x = n_hauls, y = random_intercept)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.7, colour = "#4575b4") +
  scale_x_log10() +
  labs(
    x = "Haul count per rectangle (log scale)",
    y = "Estimated rectangle random intercept",
    title = "H2/H3 results: rectangle random intercepts vs haul count",
    subtitle = "Primary model (1 | stat_rec) — partial pooling evidence for the write-up"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_pooling, p_pool, width = 8, height = 5, dpi = 150)
logmsg("Saved: ", path_out_pooling)
logmsg("Saved figure: ", path_out_fig_pooling)

# ---------------------------------------------------------------------------
# 5. Model diagnostics (results-chapter appropriate)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 5. Model diagnostics")

# Fit statistics
r2 <- nakagawa_r2_glmmtmb(fit_primary)
fit_stats <- data.frame(
  model_id = "primary_plain_re",
  n_hauls = nrow(dat),
  n_rectangles = dplyr::n_distinct(dat$stat_rec),
  r2_marginal = r2$r2_marginal,
  r2_conditional = r2$r2_conditional,
  var_fixed = r2$var_fixed,
  var_random = r2$var_random,
  var_residual = r2$var_residual,
  residual_sd = as.numeric(stats::sigma(fit_primary)),
  rectangle_intercept_sd = vc_pri$estimate,
  logLik = as.numeric(stats::logLik(fit_primary)),
  aic = as.numeric(stats::AIC(fit_primary)),
  bic = as.numeric(stats::BIC(fit_primary)),
  stringsAsFactors = FALSE
)
write_csv(fit_stats, path_out_r2)
logmsg(sprintf(
  paste0(
    "Primary model fit: Nakagawa R2 marginal = %.4f, conditional = %.4f ",
    "(manual Nakagawa–Schielzeth implementation — performance/MuMIn not installed). ",
    "logLik = %.1f, AIC = %.1f, BIC = %.1f."
  ),
  r2$r2_marginal, r2$r2_conditional, fit_stats$logLik, fit_stats$aic, fit_stats$bic
))
logmsg("Saved: ", path_out_r2)

# Residual diagnostics figure
fitted_vals <- as.numeric(stats::predict(fit_primary, re.form = NULL))
resids <- dat$residual - fitted_vals
diag_df <- data.frame(fitted = fitted_vals, resid = resids)

p_rvf <- ggplot(diag_df, aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.15, size = 0.6, colour = "#4575b4") +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7, colour = "#d73027", formula = y ~ x) +
  labs(x = "Fitted values", y = "Response residual", title = "Residuals vs fitted") +
  theme_minimal(base_size = 10)

p_qq <- ggplot(diag_df, aes(sample = resid)) +
  stat_qq(alpha = 0.2, size = 0.6, colour = "#4575b4") +
  stat_qq_line(colour = "#d73027", linewidth = 0.7) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", title = "Normal QQ plot of residuals") +
  theme_minimal(base_size = 10)

# Combine via cowplot if available, else save a simple two-panel via patchwork/grid
p_diag <- NULL
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  p_diag <- p_rvf + p_qq + plot_annotation(
    title = "H2/H3 results: primary-model residual diagnostics",
    subtitle = "Response residuals from glmmTMB primary fit (including random intercepts)"
  )
  ggsave(path_out_fig_resid, p_diag, width = 10, height = 4.5, dpi = 150)
} else {
  # Fall back: save a single stacked image via gridExtra or two files — use png device + print
  grDevices::png(path_out_fig_resid, width = 10, height = 4.5, units = "in", res = 150)
  grid::grid.newpage()
  lay <- grid::grid.layout(1, 2)
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_rvf, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_qq, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}
logmsg("Saved residual diagnostics figure: ", path_out_fig_resid)

if (has_dharma) {
  logmsg("DHARMa simulated-residual screen (primary model; n = 250 simulations):")
  set.seed(123)
  dharma_ok <- tryCatch({
    sim <- DHARMa::simulateResiduals(fit_primary, n = 250, plot = FALSE)
    ut <- DHARMa::testUniformity(sim, plot = FALSE)
    dt <- DHARMa::testDispersion(sim, plot = FALSE)
    logmsg(sprintf(
      "  - Uniformity (KS): D = %.4f, p = %.3e",
      unname(ut$statistic), ut$p.value
    ))
    logmsg(sprintf(
      "  - Dispersion: ratio/stat = %.4f, p = %.3e",
      unname(dt$statistic), dt$p.value
    ))
    logmsg(
      "  - NOTE: with n = 10,464, the KS uniformity test is almost always significant even for ",
      "small departures; read alongside the QQ plot rather than as a hard reject of the model."
    )
    TRUE
  }, error = function(e) {
    logmsg("  - DHARMa diagnostics failed: ", conditionMessage(e))
    FALSE
  })
} else {
  logmsg("DHARMa not available — skipped simulated-residual tests (QQ / resid-vs-fitted still produced).")
}

# Influence / leverage flags.
# 4/n is a common textbook cut but flags ~10% of rows at this n; for a
# results-chapter screen we require BOTH elevated Cook's D (top 1% within
# the sample) AND |Pearson residual| > 3, OR hat > 3 * mean hat (high
# leverage). Still approximate (FE hat only).
infl <- approximate_cooks_glmmtmb_noforte(fit_primary, dat)
cook_p99 <- as.numeric(stats::quantile(infl$cooks_approx, 0.99, na.rm = TRUE))
hat_cut <- 3 * mean(infl$hat, na.rm = TRUE)
infl <- infl %>%
  mutate(
    flag_cooks_extreme = is.finite(cooks_approx) & cooks_approx >= cook_p99 &
      is.finite(pearson_resid) & abs(pearson_resid) > 3,
    flag_hat = is.finite(hat) & hat > hat_cut,
    flag_any = flag_cooks_extreme | flag_hat
  )
flags <- infl %>% filter(flag_any) %>% arrange(desc(cooks_approx))
write_csv(flags, path_out_influence)
logmsg(sprintf(
  paste0(
    "Influence screen (approximate Cook's D using FE hat matrix): flag if ",
    "(Cook's D >= sample 99th percentile = %.5f AND |Pearson| > 3) OR ",
    "(hat > 3 * mean hat = %.4f). Flagged %d / %d hauls. Top 5 by Cook's D among flags:"
  ),
  cook_p99, hat_cut, nrow(flags), nrow(dat)
))
top5 <- head(flags, 5)
if (nrow(top5) > 0L) {
  for (i in seq_len(nrow(top5))) {
    r <- top5[i, ]
    logmsg(sprintf(
      "  - row %d, %s, year %d: cooks_approx = %.4f, hat = %.4f, pearson = %+.3f",
      r$row_id, r$stat_rec, r$year, r$cooks_approx, r$hat, r$pearson_resid
    ))
  }
} else {
  logmsg("  - none flagged.")
}
logmsg("Saved flagged observations: ", path_out_influence)
logmsg(
  "NOTE: Cook's D here is approximate (ignores random-effect structure in the hat matrix). ",
  "Use as a screen for write-up caveats, not as a formal deletion diagnostic."
)

conv_pri <- glmmtmb_convergence_report(fit_primary, "primary_plain_re")
logmsg(sprintf(
  "Convergence (primary): code = %s, Hessian PD = %s, max|grad| = %.3e, any NA SE = %s.",
  conv_pri$optimizer_convergence_code, conv_pri$hessian_positive_definite,
  conv_pri$max_abs_gradient, conv_pri$any_na_std_error
))

# ---------------------------------------------------------------------------
# 6. Plain-language H2 / H3 findings
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 6. Plain-language summary of H2 and H3 findings")
logmsg("")
logmsg("### Framing")
logmsg(
  "Outcome = haul-level residual = log(B_obs) - log(B_pred). A positive fishing-pressure ",
  "slope means higher Couce fishing hours associate with larger (more positive) residuals ",
  "(observed biomass above EEOS prediction); a negative slope means the opposite. No ",
  "biomass covariate. Effect sizes are on the log-residual scale per unit log(hours+1)."
)
logmsg("")
logmsg("### H2 — Does fishing pressure predict residual?")

ref_row <- primary_fe %>% filter(term == "log_hours_total")
logmsg(sprintf(
  paste0(
    "In the reference phase (1985–1988), the fishing-pressure coefficient is %+.4f ",
    "(SE %.4f, z = %+.3f, p = %.3e). Per-phase slopes (primary model):"
  ),
  ref_row$estimate, ref_row$std_error, ref_row$statistic, ref_row$p_value
))
for (i in seq_len(nrow(phase_slopes_primary))) {
  r <- phase_slopes_primary[i, ]
  # two-sided Wald p for the derived slope
  z <- r$fp_slope / r$fp_slope_se
  p <- 2 * stats::pnorm(-abs(z))
  logmsg(sprintf(
    "  - %s: %+.4f (95%% CI [%+.4f, %+.4f]; z = %+.2f, p = %.3e)",
    r$phase, r$fp_slope, r$fp_slope_lo, r$fp_slope_hi, z, p
  ))
}
logmsg(sprintf(
  paste0(
    "Joint Wald test that the FP slope is zero in every phase: chi2(%d) = %.3f, p = %.3e. ",
    "This is the global H2-style test under the primary model; phase-specific slopes and ",
    "CIs above are the effect-size detail."
  ),
  wald_h2$df, wald_h2$statistic, wald_h2$p_value
))

# Factual statement based on coefficients
sig_phases <- phase_slopes_primary %>%
  mutate(z = fp_slope / fp_slope_se, p = 2 * stats::pnorm(-abs(z))) %>%
  filter(p < 0.05)
logmsg("Factual H2 statement (coefficients only):")
if (nrow(sig_phases) == 0L) {
  logmsg(paste0(
    "  No phase-specific fishing-pressure slope differs from zero at p < 0.05 under the primary ",
    "model. Joint test p = ", sprintf("%.3e", wald_h2$p_value),
    ". Fishing pressure is not a statistically detectable predictor of residual within individual ",
    "phases at conventional alpha = 0.05; see joint test and effect sizes above rather than ",
    "significance alone."
  ))
} else {
  logmsg(paste0(
    "  Fishing-pressure slopes differ from zero at p < 0.05 in: ",
    paste(sprintf("%s (%+.4f)", sig_phases$phase, sig_phases$fp_slope), collapse = "; "),
    ". Joint test that all phase slopes are zero: p = ", sprintf("%.3e", wald_h2$p_value),
    ". Magnitudes remain small on the residual scale (see slopes above); significance and ",
    "effect size should be read together."
  ))
}

logmsg("")
logmsg("### H3 — Does the fishing-pressure–residual relationship change across phases?")
h3_rows <- primary_fe %>% filter(grepl("^log_hours_total:phase", term))
for (i in seq_len(nrow(h3_rows))) {
  r <- h3_rows[i, ]
  logmsg(sprintf(
    "  - %s: %+.4f (SE %.4f, z = %+.3f, p = %.3e)",
    r$term_plain, r$estimate, r$std_error, r$statistic, r$p_value
  ))
}
logmsg(sprintf(
  "Joint Wald test that all FP x phase interactions are zero: chi2(%d) = %.3f, p = %.3e.",
  wald_h3$df, wald_h3$statistic, wald_h3$p_value
))
sig_int <- h3_rows %>% filter(p_value < 0.05)
logmsg("Factual H3 statement (coefficients only):")
if (wald_h3$p_value >= 0.05 && nrow(sig_int) == 0L) {
  logmsg(paste0(
    "  The joint test of FP x phase interactions is not significant at alpha = 0.05 ",
    "(p = ", sprintf("%.3e", wald_h3$p_value),
    "), and no individual interaction reaches p < 0.05. The data do not support a detectable ",
    "change in the fishing-pressure–residual slope across the 1989/2001/2008 phases under the ",
    "primary model, at conventional significance. Effect sizes of the interactions are reported ",
    "above for completeness."
  ))
} else if (nrow(sig_int) > 0L) {
  logmsg(paste0(
    "  At least one FP x phase interaction differs from zero at p < 0.05: ",
    paste(sprintf("%s (%+.4f, p = %.3e)", sig_int$term, sig_int$estimate, sig_int$p_value), collapse = "; "),
    ". Joint interaction test p = ", sprintf("%.3e", wald_h3$p_value),
    ". This supports a phase-dependent change in the fishing-pressure–residual relationship ",
    "(H3) under the primary model; read alongside the per-phase slopes and the sensitivity ",
    "overlay (GAM), which should agree if the design justification holds."
  ))
} else {
  logmsg(paste0(
    "  Joint interaction test p = ", sprintf("%.3e", wald_h3$p_value),
    " with no individual interaction at p < 0.05 (or the reverse pattern). See the coefficient ",
    "table and joint test rather than a binary yes/no overstatement."
  ))
}

logmsg("")
logmsg("### Sensitivity agreement (primary vs CAR vs GAM)")
# Correlation of phase slopes primary vs CAR vs GAM means
ps <- phase_slopes_all %>%
  select(model_id, phase, fp_slope) %>%
  tidyr::pivot_wider(names_from = model_id, values_from = fp_slope)
cor_pc <- suppressWarnings(cor(ps$primary_plain_re, ps$car_adjacency))
cor_pg <- suppressWarnings(cor(ps$primary_plain_re, ps$smooth_year_gam))
logmsg(sprintf(
  paste0(
    "Phase-slope correlation across the four phases: primary vs CAR = %.3f; primary vs ",
    "GAM-within-phase mean = %.3f. Primary vs CAR parametric FE shifts are logged in section 2; ",
    "the overlay figure shows primary step vs GAM smooth year-by-year."
  ),
  cor_pc, cor_pg
))
logmsg(
  "No conclusions beyond the coefficients: this summary does not claim mechanism, management ",
  "causality, or that non-significance equals evidence of no effect."
)

# ---------------------------------------------------------------------------
# Before/after vs previous with-biomass results run
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Before/after: with-biomass previous run vs this no-biomass re-run")
if (isTRUE(prev_had_biomass) && !is.null(prev_slopes) && !is.null(prev_fit)) {
  prev_pri <- prev_slopes %>%
    filter(model_id == "primary_plain_re") %>%
    select(phase, fp_slope_prev = fp_slope)
  prev_int <- prev_fe %>%
    filter(grepl("^log_hours_total:phase", term)) %>%
    select(term, estimate_prev = estimate)
  prev_r2_m <- prev_fit$r2_marginal[1]
  prev_r2_c <- prev_fit$r2_conditional[1]
  prev_cor_note <- "previous with-biomass run correlations were computed in that run's log"
  logmsg("Source of previous numbers: on-disk results CSVs (still contained mean_ln_B_obs).")
} else {
  prev_pri <- prev_biomass_fallback$slopes
  prev_int <- prev_biomass_fallback$interactions
  prev_r2_m <- prev_biomass_fallback$r2_marginal
  prev_r2_c <- prev_biomass_fallback$r2_conditional
  prev_cor_note <- sprintf(
    "previous with-biomass run: primary-CAR ≈ %.3f, primary-GAM ≈ %.3f",
    prev_biomass_fallback$cor_primary_car, prev_biomass_fallback$cor_primary_gam
  )
  logmsg(
    "Source of previous numbers: hardcoded snapshot from the completed with-biomass results ",
    "run (on-disk CSVs had already been overwritten by a no-biomass re-run)."
  )
}
new_pri <- phase_slopes_primary %>% select(phase, fp_slope_new = fp_slope)
cmp_slopes <- prev_pri %>%
  left_join(new_pri, by = "phase") %>%
  mutate(delta = fp_slope_new - fp_slope_prev)
logmsg("Primary phase slopes (previous with biomass -> this re-run without):")
for (i in seq_len(nrow(cmp_slopes))) {
  r <- cmp_slopes[i, ]
  logmsg(sprintf(
    "  - %s: %+.4f -> %+.4f (delta %+.4f)",
    r$phase, r$fp_slope_prev, r$fp_slope_new, r$delta
  ))
}
new_int <- primary_fe %>%
  filter(grepl("^log_hours_total:phase", term)) %>%
  select(term, estimate_new = estimate)
cmp_int <- prev_int %>%
  left_join(new_int, by = "term") %>%
  mutate(delta = estimate_new - estimate_prev)
logmsg("Primary H3 interaction terms (previous -> this re-run):")
for (i in seq_len(nrow(cmp_int))) {
  r <- cmp_int[i, ]
  logmsg(sprintf(
    "  - %s: %+.4f -> %+.4f (delta %+.4f)",
    r$term, r$estimate_prev, r$estimate_new, r$delta
  ))
}
logmsg(sprintf(
  "Nakagawa R2 (previous -> this): marginal %.4f -> %.4f; conditional %.4f -> %.4f.",
  prev_r2_m, r2$r2_marginal, prev_r2_c, r2$r2_conditional
))
logmsg(sprintf(
  "Sensitivity phase-slope correlations this re-run: primary-CAR = %.3f; primary-GAM = %.3f (%s).",
  cor_pc, cor_pg, prev_cor_note
))

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_primary_fe)
logmsg("- ", path_out_sensitivity)
logmsg("- ", path_out_phase_slopes)
logmsg("- ", path_out_year_slopes)
logmsg("- ", path_out_pooling)
logmsg("- ", path_out_r2)
logmsg("- ", path_out_influence)
logmsg("- ", path_out_wald)
logmsg("- ", path_out_models)
logmsg("- ", path_out_fig_effect)
logmsg("- ", path_out_fig_pooling)
logmsg("- ", path_out_fig_resid)
logmsg("- ", path_out_session)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 shared model RESULTS run complete. ===\n")
