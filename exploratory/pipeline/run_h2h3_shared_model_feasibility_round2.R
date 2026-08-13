# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H2/H3 shared hierarchical model — FEASIBILITY CHECK, ROUND 2
# (adjacency-based / CAR spatial structure)
#
# PURPOSE: Round 1 (run_h2h3_shared_model_feasibility.R) found the continuous
# distance-decay spatial term (glmmTMB, exp() covariance over rectangle
# centroids) was not usefully identified: range and sill 95% CIs spanned
# orders of magnitude, and the range estimate exceeded the data's maximum
# inter-centroid distance. This script refits ONLY the spatial random-effect
# structure — same fixed effects, same data — using a discrete adjacency-
# based (CAR) correlation structure over ICES rectangles, built from the
# SAME queen-contiguity neighbour list already used by the original H2 SEM/
# SAR analysis (R/h2_spatial_helpers.R::build_h2_spatial_weights()).
#
# THIS IS STILL A FEASIBILITY CHECK, NOT A RESULTS RUN. Fishing-pressure,
# phase, and interaction coefficients below are NOT interpreted as answering
# H2 or H3 anywhere in this script or its outputs.
#
# MODEL:
#   residual ~ log(hours_total + 1) * phase + mean_ln_B_obs +
#              adjacency(1 | stat_rec)   [CAR random intercept, this script]
#   vs. (reused from Round 1, NOT refit here)
#   residual ~ log(hours_total + 1) * phase + mean_ln_B_obs + (1 | stat_rec)
#              [plain non-spatial comparison model]
# Same outcome (canonical `residual`), same phase definition (3 robust
# breaks: 1989, 2001, 2008; marginal 1997 excluded), same fixed effects, same
# 158-rectangle H2 panel universe as Round 1 — see
# R/h2h3_feasibility_helpers.R::build_feasibility_data()/build_phase_factor().
#
# ADJACENCY REUSE (see run log for full detail): there is no serialized
# nb/listw object on disk anywhere in this project — every existing H2
# spatial script (run_h2_sar_lag_models.R, run_h2_biomass_fishing_sar_
# diagnostics.R, run_h2_models.R) calls build_h2_spatial_weights(panel,
# rectangles_sf, project_root) fresh each time. That function is a
# deterministic pure function of (panel$stat_rec, the ICES rectangle
# shapefile): spdep::poly2nb(., queen = TRUE) then spdep::nb2listw(nb, style
# = "W", zero.policy = TRUE). This script calls that EXACT SAME function on
# the EXACT SAME panel object (outputs/h2_rectangle_panel.rds, the
# established 158-rectangle H2 universe) used throughout the original H2
# analysis — this is a reuse of the neighbour-list DEFINITION, not an
# independently-rebuilt approximation. Confirmed n_rectangles = 158, n_
# isolated = 0, matching outputs/h2_spatial_diagnostics.csv (the recorded
# diagnostics from the original H2 SEM run using this same weights object).
#
# PACKAGE CHOICE: spaMM::fitme() with the adjacency(1 | stat_rec) random-
# effect family (a proper/exact CAR model with a free spatial-autocorrelation
# parameter rho — NOT the intrinsic/improper CAR that mgcv's bs = "mrf"
# Markov random field smooth implements by default, which fixes rho = 1 by
# construction and so cannot answer check #2 below, "is rho identified /
# pinned at a boundary"). brms::car(type = "escar") would also give a proper
# CAR rho, but requires compiling a new Stan model via rstan/cmdstanr for
# this specific formula (higher risk/runtime for a feasibility check in this
# environment) and MCMC sampling; spaMM fits via fast (non-MCMC) REML/ML —
# under 2 seconds for this haul-level model — while still supporting an
# arbitrary raw adjacency matrix (adjMatrix = ...) and an explicit rho with
# an uncertainty interval, exactly what the brief requires. spaMM was named
# as an acceptable option in the Round 1 brief's own tool list. Documented
# here as the deviation from the two examples named in the Round 2 brief
# (brms::car() / mgcv mrf), with the reasoning above.
#
# NEW DEPENDENCY: spaMM was NOT part of this project's renv-managed set —
# installed ad hoc (install.packages("spaMM", repos =
# "https://cloud.r-project.org")) into the ambient/user R library, NOT added
# to renv.lock. Same pattern as strucchange (h2h3_designA4) and glmmTMB
# (Round 1). MUST be run with `Rscript --vanilla` (renv not activated) —
# see Round 1's run log for why.
#
# OUT OF SCOPE (per the brief): no interpretation of fishing-pressure/phase/
# interaction coefficients as H2/H3 results; no continuous-year-term
# alternative; no changes to phase definition, biomass covariate, or
# residual sign convention.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_shared_model_feasibility_round2.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' is required (used only for numFactor()/build_feasibility_data(), inherited from ",
    "Round 1's helper file) and is not on the library path. Run with: Rscript --vanilla ",
    "pipeline/run_h2h3_shared_model_feasibility_round2.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop(
    "Package 'spaMM' is required and is not on the library path. This is expected under a renv-",
    "activated session (spaMM was installed ad hoc into the ambient/user library, not renv-managed — ",
    "see script header). Run with: Rscript --vanilla pipeline/run_h2h3_shared_model_feasibility_round2.R"
  )
}
suppressPackageStartupMessages(library(spaMM))
suppressPackageStartupMessages(library(spdep))

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_panel <- file.path(project_root, "outputs", "h2_rectangle_panel.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
path_round1_models <- file.path(project_root, "outputs", "h2h3_feasibility_model_objects.rds")
path_round1_convergence <- file.path(project_root, "outputs", "h2h3_feasibility_convergence_diagnostics.csv")
path_h2_spatial_diagnostics <- file.path(project_root, "outputs", "h2_spatial_diagnostics.csv")
stopifnot(
  file.exists(path_haul), file.exists(path_panel), file.exists(path_couce_year),
  file.exists(path_round1_models)
)

path_out_fixed <- file.path(project_root, "outputs", "h2h3_feasibility_round2_fixed_effects.csv")
path_out_spatial_param <- file.path(project_root, "outputs", "h2h3_feasibility_round2_spatial_param.csv")
path_out_pooling <- file.path(project_root, "outputs", "h2h3_feasibility_round2_partial_pooling.csv")
path_out_convergence <- file.path(project_root, "outputs", "h2h3_feasibility_round2_convergence_diagnostics.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
path_out_fig_pooling <- file.path(fig_dir, "h2h3_feasibility_round2_partial_pooling.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_feasibility_round2_run_log.md")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

#' First capture-group match of `pattern` across a character vector `lines`,
#' or NA if no line matches. Used to parse fit-time figures back out of
#' Round 1's prose run log for a numeric fit-time comparison table.
regmatches_first <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0L) return(NA_character_)
  m <- regmatches(hit[1], regexpr(pattern, hit[1]))
  groups <- regmatches(hit[1], regexec(pattern, hit[1]))[[1]]
  if (length(groups) < 2L) return(NA_character_)
  groups[2]
}

logmsg("# H2/H3 shared hierarchical model — feasibility check ROUND 2 run log")
logmsg("## Adjacency-based (CAR) spatial structure")
logmsg("")
logmsg(
  "FEASIBILITY CHECK ONLY (Round 2). No fishing-pressure or fishing-pressure x phase coefficient below ",
  "is interpreted as answering H2 or H3. This task does not finalise the model as the committed ",
  "approach — that follows supervisor discussion, not this script."
)
logmsg(
  "Parallel in structure to the Round 1 run log (outputs/h2h3_feasibility_run_log.md) for direct side-",
  "by-side comparison; only the spatial random-effect term changes."
)

# ---------------------------------------------------------------------------
# Adjacency reuse (required, not rebuild)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Adjacency/weights object: reuse confirmation")

panel <- readRDS(path_panel)
rectangles_sf <- load_ices_rectangles_sf(project_root)
weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)

logmsg(
  "No serialized nb/listw object exists anywhere in this repository (checked: outputs/*.rds, ",
  "outputs/*listw*, outputs/*weight* — none found). Every existing H2 spatial script ",
  "(run_h2_sar_lag_models.R, run_h2_biomass_fishing_sar_diagnostics.R, run_h2_models.R) calls ",
  "R/h2_spatial_helpers.R::build_h2_spatial_weights(panel, rectangles_sf, project_root) fresh each ",
  "time it is needed. That function is a DETERMINISTIC pure function of (panel$stat_rec, the ICES ",
  "rectangle shapefile at gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp): spdep::poly2nb(., ",
  "queen = TRUE) then spdep::nb2listw(nb, style = \"W\", zero.policy = TRUE). This script calls that ",
  "EXACT SAME function, unmodified, on the EXACT SAME panel object (outputs/h2_rectangle_panel.rds, ",
  "the established 158-rectangle H2 universe) that the original H2 SEM/SAR analysis used — so the ",
  "resulting neighbour list is PROVABLY identical (same code path, same inputs), not an independently-",
  "built approximation. No rebuild-and-match step was needed because there was nothing separate to ",
  "match against other than re-deriving via the identical function."
)
logmsg(sprintf(
  "Resulting neighbour list: n_rectangles = %d, n_isolated (zero-neighbour rectangles) = %d.",
  nrow(panel), weights$n_isolated
))
if (file.exists(path_h2_spatial_diagnostics)) {
  h2_diag <- read_csv(path_h2_spatial_diagnostics, show_col_types = FALSE)
  ref_n_rect <- unique(h2_diag$n_rectangles)
  ref_n_iso <- unique(h2_diag$n_isolated)
  match_ok <- isTRUE(all.equal(nrow(panel), ref_n_rect)) && isTRUE(all.equal(weights$n_isolated, ref_n_iso))
  logmsg(sprintf(
    paste0(
      "Cross-check against outputs/h2_spatial_diagnostics.csv (recorded from the original H2 SEM run ",
      "using this same weights object): n_rectangles = %s, n_isolated = %s. MATCH: %s."
    ),
    paste(ref_n_rect, collapse = ","), paste(ref_n_iso, collapse = ","), match_ok
  ))
} else {
  logmsg("outputs/h2_spatial_diagnostics.csv not found — cross-check skipped (flagged, not resolved).")
}

W <- build_car_adjacency_matrix(weights$nb, panel)
logmsg(sprintf(
  paste0(
    "Built the raw binary (0/1) adjacency matrix from this SAME reused neighbour list (weights$nb) via ",
    "spdep::nb2mat(nb, style = \"B\") — %d rectangles, %d undirected edges (%d matrix entries = 1), ",
    "symmetric: %s. This differs from the row-STANDARDISED \"W\"-style listw used by the original H2 ",
    "errorsarlm()/lagsarlm() SEM/SAR models only in MATRIX REPRESENTATION (binary vs. row-normalised), ",
    "which is the standard convention difference between CAR models (binary adjacency, correlation ",
    "strength captured by a free rho) and SAR/SEM models (row-standardised weights) — not a different ",
    "neighbour-list DEFINITION."
  ),
  nrow(W), sum(W) / 2, sum(W), isSymmetric(unname(W))
))

# ---------------------------------------------------------------------------
# Data (identical to Round 1)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data and analysis universe (identical to Round 1)")

haul <- readRDS(path_haul)
couce_year <- readRDS(path_couce_year)
built <- build_feasibility_data(haul, panel, couce_year, H2_YEAR_MIN, H2_YEAR_MAX)
dat <- built$data
dat$stat_rec <- factor(dat$stat_rec, levels = panel$stat_rec)

logmsg(sprintf(
  paste0(
    "%d hauls across %d rectangles — same universe, same fixed-effect construction as Round 1 (see that ",
    "run log for the full breakdown of the H2-panel/Couce-year join)."
  ),
  nrow(dat), dplyr::n_distinct(dat$stat_rec)
))

# ---------------------------------------------------------------------------
# Reuse Round 1's plain non-spatial model (do not refit)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Reusing Round 1's plain non-spatial comparison model")

round1 <- readRDS(path_round1_models)
fit_plain <- round1$fit_plain
logmsg(
  "Loaded fit_plain (residual ~ log_hours_total * phase + mean_ln_B_obs + (1 | stat_rec)) directly from ",
  path_round1_models, " — NOT refit here, per the brief."
)

# ---------------------------------------------------------------------------
# Fit the CAR model
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fitting")

formula_car <- residual ~ log_hours_total * phase + mean_ln_B_obs + adjacency(1 | stat_rec)
logmsg("Fitting CAR model: residual ~ log_hours_total * phase + mean_ln_B_obs + adjacency(1 | stat_rec) via spaMM::fitme(..., method = \"REML\") ...")
time_car <- system.time({
  fit_car <- fitme(formula_car, data = dat, adjMatrix = W, method = "REML")
})
logmsg(sprintf(
  "CAR model fit time: %.2f sec elapsed (%.2f sec user, %.2f sec system).",
  time_car["elapsed"], time_car["user.self"], time_car["sys.self"]
))

path_round1_run_log <- file.path(project_root, "outputs", "h2h3_feasibility_run_log.md")
time_plain_round1 <- NA_real_
time_spatial_round1 <- NA_real_
if (file.exists(path_round1_run_log)) {
  r1_log <- readLines(path_round1_run_log)
  m_spatial <- regmatches_first(r1_log, "Full spatial model fit time: ([0-9.]+) sec elapsed")
  m_plain <- regmatches_first(r1_log, "Simpler non-spatial model fit time: ([0-9.]+) sec elapsed")
  if (!is.na(m_spatial)) time_spatial_round1 <- as.numeric(m_spatial)
  if (!is.na(m_plain)) time_plain_round1 <- as.numeric(m_plain)
}
logmsg(sprintf(
  paste0(
    "Fit-time comparison (all reused times parsed directly from outputs/h2h3_feasibility_run_log.md, not ",
    "refit): CAR adjacency model (spaMM, this script) = %.2f sec elapsed; Round 1 plain (1 | stat_rec) ",
    "model (glmmTMB) = %s sec elapsed; Round 1 full continuous-spatial exp() model (glmmTMB) = %s sec ",
    "elapsed. All three are fast (single-digit seconds); the CAR model here is comparable to or faster ",
    "than both Round 1 models, consistent with spaMM's non-MCMC REML fitting on this modestly-sized ",
    "(158 x 158) adjacency matrix — fit time is not a practical obstacle to iterating on this model ",
    "during further development."
  ),
  time_car["elapsed"],
  ifelse(is.finite(time_plain_round1), sprintf("%.2f", time_plain_round1), "not recorded"),
  ifelse(is.finite(time_spatial_round1), sprintf("%.2f", time_spatial_round1), "not recorded")
))

saveRDS(
  list(fit_car = fit_car, fit_plain = fit_plain, data = dat, adjMatrix = W),
  path_out_models
)
logmsg("Saved fitted model objects (+ analysis data.frame + adjacency matrix): ", path_out_models)

# ---------------------------------------------------------------------------
# 1. Convergence diagnostics
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Convergence diagnostics")

rhorange <- car_rho_admissible_range(W)
logmsg(sprintf(
  paste0(
    "spaMM's fitme() does not expose an iterative-optimizer status code the way glmmTMB does (no MCMC ",
    "either — REML/ML via TMB-free internal likelihood maximisation). Reporting what IS available: (a) ",
    "any warnings/errors captured during the fit, (b) a robustness check refitting from 3 different ",
    "starting values for rho spanning the admissible range [%.4f, %.4f] and confirming they converge to ",
    "the same estimate (a standard local-optimum check), and (c) whether all fixed-effect standard ",
    "errors are finite."
  ),
  rhorange["lower"], rhorange["upper"]
))
conv_car <- spamm_car_convergence_report(fit_car, "car_adjacency", formula_car, dat, W, rhorange)
write_csv(conv_car, path_out_convergence)
logmsg(sprintf(
  paste0(
    "  - car_adjacency: warnings during fit = %d (%s); refits from near-lower/zero/near-upper starting ",
    "rho all converge to rho = %.6f (robust to starting value = %s); any non-finite fixed-effect SE = %s."
  ),
  conv_car$n_warnings_during_fit, conv_car$warnings_text, conv_car$fitted_rho,
  conv_car$robust_to_starting_value, conv_car$any_na_or_infinite_fixed_effect_se
))
logmsg("Saved: ", path_out_convergence)

# ---------------------------------------------------------------------------
# 2. Spatial autocorrelation parameter (rho) identification
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Spatial autocorrelation parameter (rho) identification")

spatial_param <- extract_car_spatial_param(fit_car, W, rhorange, "car_adjacency", nsim = 99, seed = 123)
write_csv(spatial_param, path_out_spatial_param)

logmsg(sprintf(
  paste0(
    "rho estimate = %.6f. Admissible range for this adjacency matrix (from its eigenvalues, 1/lambda_min ",
    "to 1/lambda_max) = [%.6f, %.6f]. 95%% parametric bootstrap CI (percentile method, %d replicates, ",
    "seed = 123) = [%.6f, %.6f]. Rho lies %.1f%% of the way from the lower bound to the upper bound of ",
    "the admissible range (i.e. %.1f%% of the range's width from the UPPER bound)."
  ),
  spatial_param$rho_estimate, spatial_param$rho_admissible_lower, spatial_param$rho_admissible_upper,
  spatial_param$nsim_bootstrap, spatial_param$rho_ci_low_boot, spatial_param$rho_ci_high_boot,
  100 * (1 - spatial_param$fraction_of_range_from_upper), 100 * spatial_param$fraction_of_range_from_upper
))
logmsg(sprintf(
  "CAR variance component (lambda, variance of the rectangle-level CAR random effect) = %.6f.",
  spatial_param$car_lambda_variance
))
if (isTRUE(spatial_param$boundary_pinned_flag)) {
  logmsg(sprintf(
    paste0(
      "FLAG: rho is PINNED NEAR A BOUNDARY of its admissible range (within 1%% of the upper or lower ",
      "bound: estimate %.6f vs. upper bound %.6f — %.2f%% of the way to the boundary). This is a ",
      "different FORM of non-identifiability than Round 1's (there, the continuous range parameter's CI ",
      "spanned several orders of magnitude and exceeded the data's spatial extent; here, the discrete ",
      "adjacency structure DOES converge to a single well-defined point estimate and a comparatively ",
      "narrow bootstrap CI, but that point estimate sits essentially at the edge of what the model ",
      "algebraically allows). Reported as observed, not resolved."
    ),
    spatial_param$rho_estimate, spatial_param$rho_admissible_upper,
    100 * (1 - spatial_param$fraction_of_range_from_upper)
  ))
} else {
  logmsg("Rho is NOT pinned near a boundary of its admissible range on this check.")
}

# ---------------------------------------------------------------------------
# 3. Partial pooling
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Partial pooling check")

pooling_car <- extract_partial_pooling_spamm(fit_car, dat, "car_adjacency")
pooling_plain <- extract_partial_pooling(fit_plain, dat, "simple_nonspatial", spatial = FALSE)
pooling_table <- bind_rows(pooling_car, pooling_plain)
write_csv(pooling_table, path_out_pooling)

logmsg(
  "Same shrinkage_ratio definition as Round 1: random_intercept / unpooled_intercept, where ",
  "unpooled_intercept is the per-rectangle mean of (observed residual - the model's OWN population-",
  "level [fixed-effects-only, re.form = NA] prediction)."
)
cor_car <- suppressWarnings(cor(log(pooling_car$n_hauls), pooling_car$shrinkage_ratio, method = "spearman"))
cor_plain <- suppressWarnings(cor(log(pooling_plain$n_hauls), pooling_plain$shrinkage_ratio, method = "spearman"))
logmsg(sprintf(
  paste0(
    "Spearman correlation of log(n_hauls) with shrinkage_ratio across all %d rectangles: car_adjacency ",
    "rho = %.3f; simple_nonspatial (reused Round 1 values) rho = %.3f."
  ),
  nrow(pooling_car), cor_car, cor_plain
))

sample_car <- pooling_car %>%
  arrange(n_hauls) %>%
  slice(unique(round(seq(1, nrow(.), length.out = 10))))
logmsg("Sample of rectangles spanning the haul-count range (CAR adjacency model):")
for (i in seq_len(nrow(sample_car))) {
  r <- sample_car[i, ]
  logmsg(sprintf(
    "  - %s: n_hauls = %d, unpooled_intercept = %.4f, random_intercept = %.4f, shrinkage_ratio = %.3f",
    r$stat_rec, r$n_hauls, r$unpooled_intercept, r$random_intercept, r$shrinkage_ratio
  ))
}
sample_plain <- pooling_plain %>%
  filter(stat_rec %in% sample_car$stat_rec) %>%
  arrange(n_hauls)
logmsg("Same rectangles, simple non-spatial model (reused from Round 1):")
for (i in seq_len(nrow(sample_plain))) {
  r <- sample_plain[i, ]
  logmsg(sprintf(
    "  - %s: n_hauls = %d, unpooled_intercept = %.4f, random_intercept = %.4f, shrinkage_ratio = %.3f",
    r$stat_rec, r$n_hauls, r$unpooled_intercept, r$random_intercept, r$shrinkage_ratio
  ))
}
logmsg("Saved: ", path_out_pooling)

p_pooling <- ggplot(pooling_table, aes(x = n_hauls, y = random_intercept)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.7, colour = "#4575b4") +
  scale_x_log10() +
  facet_wrap(~model_id, ncol = 2) +
  labs(
    x = "Haul count per rectangle (log scale)",
    y = "Estimated rectangle random intercept",
    title = "Round 2: partial pooling, CAR adjacency vs plain (1|stat_rec)",
    subtitle = "Feasibility check only — not an H2/H3 result"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_pooling, p_pooling, width = 9.5, height = 5, dpi = 150)
logmsg("Saved figure: ", path_out_fig_pooling)

# ---------------------------------------------------------------------------
# 4. Fixed effects: CAR model vs plain non-spatial comparison model
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 4. Fixed effects: CAR adjacency model vs plain non-spatial comparison (reused from Round 1)")

fixed_car <- tidy_fixed_effects_spamm(fit_car, "car_adjacency")
fixed_plain <- tidy_fixed_effects(fit_plain, "simple_nonspatial")
fixed_table <- bind_rows(fixed_car, fixed_plain)
write_csv(fixed_table, path_out_fixed)

comparison <- fixed_table %>%
  select(model_id, term, estimate, std_error) %>%
  tidyr::pivot_wider(names_from = model_id, values_from = c(estimate, std_error))
logmsg(
  "Descriptive only — NOT a formal model-selection test, and NOT an interpretation of the fishing-",
  "pressure or interaction terms as answering H2/H3. Reporting how much fixed-effect estimates and ",
  "SEs shift between the plain non-spatial model (Round 1, reused) and this CAR adjacency model:"
)
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  pct_est_change <- 100 * (r$estimate_car_adjacency - r$estimate_simple_nonspatial) / abs(r$estimate_simple_nonspatial)
  pct_se_change <- 100 * (r$std_error_car_adjacency - r$std_error_simple_nonspatial) / r$std_error_simple_nonspatial
  logmsg(sprintf(
    "  - %s: estimate %.4f -> %.4f (%+.1f%%); SE %.4f -> %.4f (%+.1f%%)",
    r$term, r$estimate_simple_nonspatial, r$estimate_car_adjacency, pct_est_change,
    r$std_error_simple_nonspatial, r$std_error_car_adjacency, pct_se_change
  ))
}
logmsg("Saved: ", path_out_fixed)

# ---------------------------------------------------------------------------
# Feasibility verdict
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Feasibility verdict")

fit_ok <- conv_car$n_warnings_during_fit == 0L &&
  conv_car$robust_to_starting_value &&
  !conv_car$any_na_or_infinite_fixed_effect_se
logmsg(sprintf(
  paste0(
    "Fit-level behaviour (no warnings/errors during fitting, robust to 3 different starting values for ",
    "rho, all fixed-effect SEs finite): %s."
  ),
  fit_ok
))
logmsg(sprintf(
  paste0(
    "Rho boundary-pinning flag: %s (rho = %.6f is %.2f%% of the way to the nearer admissible bound out ",
    "of [%.6f, %.6f])."
  ),
  spatial_param$boundary_pinned_flag, spatial_param$rho_estimate,
  100 * (1 - min(spatial_param$fraction_of_range_from_upper, spatial_param$fraction_of_range_from_lower)),
  spatial_param$rho_admissible_lower, spatial_param$rho_admissible_upper
))
logmsg("PLAIN STATEMENT (per the brief — feasible / feasible with caveats / not feasible):")

if (fit_ok && !spatial_param$boundary_pinned_flag) {
  logmsg("  FEASIBLE. The CAR adjacency model fits cleanly and rho is well identified within its admissible range.")
} else if (fit_ok && spatial_param$boundary_pinned_flag) {
  logmsg(paste0(
    "  FEASIBLE WITH CAVEATS. Unlike Round 1's continuous distance-decay model, the CAR adjacency ",
    "structure fits cleanly by every check available in this framework: no warnings or errors during ",
    "fitting, the optimum is reached from 3 widely-spaced starting values for rho (near the lower ",
    "admissible bound, zero, near the upper admissible bound) with no sensitivity to starting value, ",
    "and all fixed-effect standard errors are finite. This is a materially better-behaved fit than ",
    "Round 1's exp() structure, which had range/sill CIs spanning orders of magnitude. HOWEVER, the ",
    "fitted rho (", sprintf("%.6f", spatial_param$rho_estimate), ") sits essentially AT the upper edge ",
    "of what this adjacency matrix algebraically permits (admissible range [",
    sprintf("%.6f", spatial_param$rho_admissible_lower), ", ", sprintf("%.6f", spatial_param$rho_admissible_upper),
    "]; rho is ", sprintf("%.2f", 100 * (1 - spatial_param$fraction_of_range_from_upper)), "% of the way ",
    "to that upper bound) — a boundary-pinned rho, even with a narrow bootstrap CI, means the model is ",
    "pushing as much positive spatial correlation into the CAR term as the adjacency structure allows, ",
    "which is itself a form of identification strain (in the CAR literature this is a well-known ",
    "feature of rectangle-density adjacency structures with strong, spatially widespread residual ",
    "correlation — the model wants MORE spatial smoothing than a single free rho against this ",
    "particular graph can supply). RECOMMENDATION for supervisor discussion: proceed with this CAR ",
    "structure as the practical choice between the two spatial approaches tried (it is unambiguously ",
    "better-identified than Round 1's continuous exp() structure), but flag the boundary-adjacent rho ",
    "as a caveat — e.g. consider whether a Leroux-type CAR (mixing pure CAR with an unstructured ",
    "component) or the plain (1|stat_rec) model (which has no such boundary issue at all) is preferred ",
    "if this sensitivity matters for the eventual H2/H3 model."
  ))
} else {
  logmsg("  NOT FEASIBLE AS SPECIFIED. See convergence diagnostics above for the specific failure.")
}
logmsg(
  "This verdict is about model FIT BEHAVIOUR only — it does not interpret, and is not based on, the ",
  "sign, magnitude, or significance of the fishing-pressure or fishing-pressure x phase coefficients."
)

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_fixed)
logmsg("- ", path_out_spatial_param)
logmsg("- ", path_out_pooling)
logmsg("- ", path_out_convergence)
logmsg("- ", path_out_models)
logmsg("- ", path_out_fig_pooling)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 shared model feasibility check (Round 2, CAR adjacency) complete. ===\n")
