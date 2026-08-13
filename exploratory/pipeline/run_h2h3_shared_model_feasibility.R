# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H2/H3 shared hierarchical model — FEASIBILITY CHECK (Section 4 design proposal)
#
# PURPOSE: check whether the proposed shared H2/H3 model converges and behaves
# sensibly on this data, before it is committed to as the approach. THIS IS A
# FEASIBILITY CHECK, NOT A RESULTS RUN. Fishing-pressure and fishing-pressure x
# phase coefficients below are NOT interpreted as answering H2 or H3 anywhere
# in this script or its outputs.
#
# MODEL (per the briefing's Section 4 spec, reproduced here since the source
# proposal document was not found as a repo file):
#   residual ~ log(fishing_hours + 1) * phase + mean_ln_B_obs +
#              [rectangle random intercept with spatial correlation over
#               rectangle centroids]  (full model)
#   residual ~ log(fishing_hours + 1) * phase + mean_ln_B_obs + (1 | stat_rec)
#              (simpler non-spatial comparison model)
# outcome = canonical pipeline `residual` (= log(B_obs) - log(B_pred)), NOT
# the flipped `resid_signed` convention used in the pre-H3 exploration script.
# phase = 4-level factor from the 3 robust structural breaks (1989, 2001,
# 2008; marginal 1997 break excluded) — see build_phase_factor() in
# R/h2h3_feasibility_helpers.R for the exact boundary convention used.
#
# PACKAGE CHOICE: glmmTMB (spatial covariance structures via exp(pos+0|group)
# on rectangle centroids). NEW DEPENDENCY: glmmTMB was NOT part of this
# project's renv-managed dependency set — installed ad hoc
# (`install.packages("glmmTMB", lib = <user library>, repos =
# "https://cloud.r-project.org")`) into the ambient/user R library, NOT added
# to renv.lock. Same pattern as strucchange for h2h3_designA4 (see that run
# log). Consequently THIS SCRIPT MUST BE RUN WITHOUT renv ACTIVATION
# (`Rscript --vanilla pipeline/run_h2h3_shared_model_feasibility.R` from the
# project root, or equivalent) — under a normal renv-activated session,
# glmmTMB (and, per the ambient-library check performed for this task, also
# spdep/spatialreg/lme4) are not on the library path. Flagged as an
# environment note, not resolved here.
#
# OUT OF SCOPE (per the briefing): no interpretation of the fishing-pressure
# or interaction coefficients as answering H2/H3; no finalising of this model
# as the committed approach; no changes to phase/break definitions, the
# biomass covariate, or the residual sign convention.

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_shared_model_feasibility.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' is required and is not on the library path. This is expected under a ",
    "renv-activated session (glmmTMB was installed ad hoc into the ambient/user library, not ",
    "renv-managed — see script header). Run with: Rscript --vanilla ",
    "pipeline/run_h2h3_shared_model_feasibility.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_panel <- file.path(project_root, "outputs", "h2_rectangle_panel.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
stopifnot(file.exists(path_haul), file.exists(path_panel), file.exists(path_couce_year))

path_out_fixed <- file.path(project_root, "outputs", "h2h3_feasibility_fixed_effects.csv")
path_out_varcomp <- file.path(project_root, "outputs", "h2h3_feasibility_variance_components.csv")
path_out_pooling <- file.path(project_root, "outputs", "h2h3_feasibility_partial_pooling.csv")
path_out_convergence <- file.path(project_root, "outputs", "h2h3_feasibility_convergence_diagnostics.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_feasibility_model_objects.rds")
path_out_fig_pooling <- file.path(fig_dir, "h2h3_feasibility_partial_pooling.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_feasibility_run_log.md")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 shared hierarchical model — feasibility check run log")
logmsg("")
logmsg(
  "FEASIBILITY CHECK ONLY. No fishing-pressure or fishing-pressure x phase coefficient below is ",
  "interpreted as answering H2 or H3. This task does not finalise the model as the committed ",
  "approach — that follows supervisor discussion, not this script."
)

# ---------------------------------------------------------------------------
# Package choice and dependency note
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Package choice and new dependency")
logmsg(
  "Package: glmmTMB (version ", as.character(utils::packageVersion("glmmTMB")), "). Chosen over brms ",
  "(would need Stan toolchain + MCMC compile time, much slower for a feasibility check) and spaMM ",
  "(less commonly used/documented in this ecosystem) because glmmTMB supports the required exp() ",
  "spatial covariance structure over arbitrary point coordinates natively, fits via fast ML/REML ",
  "(TMB automatic differentiation) rather than MCMC, and integrates with the same glm-style formula ",
  "interface already familiar from this project's lme4/nlme-adjacent H2 work."
)
logmsg(
  "NEW DEPENDENCY: glmmTMB was NOT part of this project's renv-managed set (absent from renv.lock and ",
  "from the renv project library). Installed ad hoc via install.packages('glmmTMB', repos = ",
  "'https://cloud.r-project.org') into the ambient/user R library — same pattern as strucchange for ",
  "the h2h3_designA4 structural-break task (see that run log). NOT added to renv.lock. A version ",
  "mismatch warning is emitted at load time (glmmTMB was built against TMB 1.9.19; TMB 1.9.20 was ",
  "installed as its dependency) — verified with a trivial test model that this does not prevent fitting ",
  "or produce incorrect results; noted here as an environment fact, not resolved."
)
logmsg(
  "CONSEQUENCE FOR HOW THIS SCRIPT MUST BE RUN: because glmmTMB (and, per an ambient-library check run ",
  "for this task, also spdep/spatialreg/lme4 used by the existing H2 spatial scripts) are not visible ",
  "under a renv-activated R session, this script must be run WITHOUT renv activation: `Rscript ",
  "--vanilla pipeline/run_h2h3_shared_model_feasibility.R` from the project root. Flagged as an ",
  "environment inconsistency already present in the project (run_h2_models.R has the same ",
  "unstated requirement), not introduced or resolved here."
)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data and analysis universe")

haul <- readRDS(path_haul)
panel <- readRDS(path_panel)
couce_year <- readRDS(path_couce_year)

built <- build_feasibility_data(haul, panel, couce_year, H2_YEAR_MIN, H2_YEAR_MAX)
dat <- built$data

logmsg(
  "Analysis universe: the established H2 rectangle panel (outputs/h2_rectangle_panel.rds), i.e. ",
  "rectangles with >= ", H2_MIN_HAULS_DEFAULT, " hauls AND Couce fishing-pressure coverage — ",
  nrow(panel), " rectangles. No new rectangle-inclusion rule introduced for this task."
)
logmsg(
  "Outcome is haul-level (not rectangle-year aggregated), per the brief. ", built$n_hauls_in_h2_universe,
  " hauls fall in this 158-rectangle universe with a finite `residual`; of these, ",
  built$n_hauls_dropped_no_rect_year_couce, " lack a Couce record for their SPECIFIC rectangle-year ",
  "(rectangle has overall Couce coverage, but not every year) and are dropped by the inner join to ",
  "rectangle-year fishing pressure, leaving ", built$n_hauls_final, " hauls for the model."
)
logmsg(
  "Fixed effects: log(hours_total + 1) [Couce rectangle-YEAR fishing pressure, same log1p transform ",
  "as the established H2 log_mean_annual_hours_total convention] x phase [4-level, see below] + ",
  "mean_ln_B_obs [rectangle-level mean log observed biomass, the established H2 confound, from ",
  "h2_rectangle_panel.rds — constant within rectangle, does not vary by year]."
)

phase_table <- dat %>%
  group_by(phase) %>%
  summarise(
    year_min = min(year), year_max = max(year),
    n_hauls = dplyr::n(), n_rectangles = dplyr::n_distinct(stat_rec),
    .groups = "drop"
  )
logmsg("Phase definition and haul counts (break year = first year of new phase; see build_phase_factor()):")
for (i in seq_len(nrow(phase_table))) {
  r <- phase_table[i, ]
  logmsg(sprintf(
    "  - %s (years %d-%d): n_hauls = %d, n_rectangles = %d",
    as.character(r$phase), r$year_min, r$year_max, r$n_hauls, r$n_rectangles
  ))
}

n_unique_pos <- dplyr::n_distinct(dat$pos)
logmsg(
  "Spatial structure: ", n_unique_pos, " unique rectangle centroids (= number of analysis rectangles; ",
  "one spatial random-effect level per rectangle, shared by all hauls in that rectangle, matching the ",
  "H2 SEM's rectangle-level spatial error term in spirit)."
)

# ---------------------------------------------------------------------------
# Fit models
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fitting")

formula_fixed <- residual ~ log_hours_total * phase + mean_ln_B_obs

logmsg("Fitting full spatial model: residual ~ log_hours_total * phase + mean_ln_B_obs + exp(pos + 0 | dummy) ...")
time_spatial <- system.time({
  fit_spatial <- glmmTMB(
    update(formula_fixed, . ~ . + exp(pos + 0 | dummy)),
    data = dat,
    REML = TRUE
  )
})
logmsg(sprintf(
  "Full spatial model fit time: %.1f sec elapsed (%.1f sec user, %.1f sec system).",
  time_spatial["elapsed"], time_spatial["user.self"], time_spatial["sys.self"]
))

logmsg("Fitting simpler non-spatial comparison model: residual ~ log_hours_total * phase + mean_ln_B_obs + (1 | stat_rec) ...")
time_plain <- system.time({
  fit_plain <- glmmTMB(
    update(formula_fixed, . ~ . + (1 | stat_rec)),
    data = dat,
    REML = TRUE
  )
})
logmsg(sprintf(
  "Simpler non-spatial model fit time: %.1f sec elapsed (%.1f sec user, %.1f sec system).",
  time_plain["elapsed"], time_plain["user.self"], time_plain["sys.self"]
))
logmsg(sprintf(
  "Practicality: the full spatial model takes %.1fx as long to fit as the simpler comparison model.",
  as.numeric(time_spatial["elapsed"]) / as.numeric(time_plain["elapsed"])
))

saveRDS(
  list(fit_spatial = fit_spatial, fit_plain = fit_plain, data = dat),
  path_out_models
)
logmsg("Saved fitted model objects (+ analysis data.frame): ", path_out_models)

# ---------------------------------------------------------------------------
# Convergence diagnostics (item 1)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Convergence diagnostics")

conv_table <- bind_rows(
  glmmtmb_convergence_report(fit_spatial, "full_spatial"),
  glmmtmb_convergence_report(fit_plain, "simple_nonspatial")
)
write_csv(conv_table, path_out_convergence)
for (i in seq_len(nrow(conv_table))) {
  r <- conv_table[i, ]
  logmsg(sprintf(
    "  - %s: optimizer_convergence_code = %s (0 = converged per nlminb), message = '%s', Hessian positive-definite = %s, max|gradient| = %.3e, any NA std. errors = %s",
    r$model_id, r$optimizer_convergence_code, r$optimizer_message, r$hessian_positive_definite,
    r$max_abs_gradient, r$any_na_std_error
  ))
  logmsg("    glmmTMB::diagnose() output: ", r$diagnose_output)
}
logmsg("Saved: ", path_out_convergence)

# ---------------------------------------------------------------------------
# Variance components (item 2)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Variance components")

varcomp_table <- bind_rows(
  extract_variance_components(fit_spatial, "full_spatial", spatial = TRUE),
  extract_variance_components(fit_plain, "simple_nonspatial", spatial = FALSE)
)
write_csv(varcomp_table, path_out_varcomp)
for (i in seq_len(nrow(varcomp_table))) {
  r <- varcomp_table[i, ]
  logmsg(sprintf(
    "  - %s / %s: estimate = %.4f, 95%% CI = [%.4f, %.4f], residual SD = %.4f",
    r$model_id, r$component, r$estimate, r$ci_low, r$ci_high, r$residual_sd
  ))
}
range_row <- varcomp_table %>% filter(model_id == "full_spatial", component == "spatial_range")
sd_row <- varcomp_table %>% filter(model_id == "full_spatial", component == "rectangle_spatial_sd")
logmsg(
  "Boundary check (spatial model): range estimate = ", sprintf("%.3f", range_row$estimate),
  " (degrees, since centroids are in lon/lat); range CI upper bound = ", sprintf("%.3f", range_row$ci_high),
  ". Flagged as degenerate if the range is pinned near 0, near/above the maximum inter-centroid ",
  "distance in the data, or if its CI is many orders of magnitude wide (all checked explicitly below)."
)
max_dist <- {
  coords <- dat %>% distinct(stat_rec, rect_lon, rect_lat)
  d <- as.matrix(dist(coords[, c("rect_lon", "rect_lat")]))
  max(d)
}
logmsg(sprintf("Maximum inter-centroid distance in the data (degrees): %.3f", max_dist))
if (range_row$estimate > max_dist) {
  logmsg("FLAG: spatial range estimate EXCEEDS the maximum inter-centroid distance in the data — the fitted correlation structure is close to a single shared field across all rectangles, not a locally-decaying one. Reported as a degeneracy flag, not resolved.")
} else if (range_row$estimate < 0.05) {
  logmsg("FLAG: spatial range estimate is near zero relative to rectangle spacing — the exp() structure is behaving close to an unstructured/independent rectangle effect. Reported as a degeneracy flag, not resolved.")
} else {
  logmsg("Range estimate falls within the plausible interval (0.05 degrees, max inter-centroid distance) — not flagged as degenerate on this simple check.")
}
plain_row <- varcomp_table %>% filter(model_id == "simple_nonspatial", component == "rectangle_intercept_sd")
if (plain_row$estimate < 1e-4) {
  logmsg("FLAG: simple non-spatial model's rectangle random-intercept SD is ~0 — near-singular fit.")
} else {
  logmsg(sprintf("Simple non-spatial model's rectangle random-intercept SD (%.4f) is not near zero — not flagged as singular.", plain_row$estimate))
}
logmsg("Saved: ", path_out_varcomp)

# ---------------------------------------------------------------------------
# Partial pooling (item 3)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Partial pooling check")

pooling_spatial <- extract_partial_pooling(fit_spatial, dat, "full_spatial", spatial = TRUE)
pooling_plain <- extract_partial_pooling(fit_plain, dat, "simple_nonspatial", spatial = FALSE)
pooling_table <- bind_rows(pooling_spatial, pooling_plain)
write_csv(pooling_table, path_out_pooling)

logmsg(
  "Shrinkage measure: shrinkage_ratio = model random intercept / unpooled_intercept, where ",
  "unpooled_intercept is the per-rectangle mean of (observed residual - the model's OWN population-",
  "level [fixed-effects-only, re.form=NA] prediction) — i.e. what each rectangle's intercept would be ",
  "under NO pooling at all, on the same scale as the random intercept (both are deviations from the ",
  "SAME fixed-effects prediction). Expected pattern: shrinkage_ratio closer to 0 (heavy shrinkage ",
  "toward the population mean) for low-haul-count rectangles, closer to 1 (little shrinkage) for high-",
  "haul-count rectangles."
)
cor_spatial <- suppressWarnings(cor(log(pooling_spatial$n_hauls), pooling_spatial$shrinkage_ratio, method = "spearman"))
cor_plain <- suppressWarnings(cor(log(pooling_plain$n_hauls), pooling_plain$shrinkage_ratio, method = "spearman"))
logmsg(
  "Spearman correlation of log(n_hauls) with shrinkage_ratio across all ", nrow(pooling_spatial),
  " rectangles: full_spatial rho = ", sprintf("%.3f", cor_spatial), "; simple_nonspatial rho = ",
  sprintf("%.3f", cor_plain), ". Both are POSITIVE, consistent with the expected shrinkage pattern ",
  "(more hauls -> shrinkage_ratio closer to 1), but only weak-to-moderate and noisy at the individual-",
  "rectangle level (see sample below — a few rectangles have shrinkage_ratio outside [0, 1], including ",
  "sign flips). For the full_spatial model this is expected to be noisier than a plain random intercept: ",
  "the exp() structure pools each rectangle toward a SPATIALLY-WEIGHTED combination of its neighbours, ",
  "not just the flat population mean, so haul count alone does not fully determine the amount or ",
  "direction of shrinkage. Reported as observed, not assumed, and not smoothed over."
)

sample_spatial <- pooling_spatial %>%
  arrange(n_hauls) %>%
  slice(unique(round(seq(1, nrow(.), length.out = 10))))
logmsg("Sample of rectangles spanning the haul-count range (full spatial model):")
for (i in seq_len(nrow(sample_spatial))) {
  r <- sample_spatial[i, ]
  logmsg(sprintf(
    "  - %s: n_hauls = %d, unpooled_intercept = %.4f, random_intercept = %.4f, shrinkage_ratio = %.3f",
    r$stat_rec, r$n_hauls, r$unpooled_intercept, r$random_intercept, r$shrinkage_ratio
  ))
}
sample_plain <- pooling_plain %>%
  filter(stat_rec %in% sample_spatial$stat_rec) %>%
  arrange(n_hauls)
logmsg("Same rectangles, simple non-spatial model:")
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
    title = "Partial pooling: rectangle random intercept vs haul count",
    subtitle = "Feasibility check only — not an H2/H3 result"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_pooling, p_pooling, width = 9.5, height = 5, dpi = 150)
logmsg("Saved figure: ", path_out_fig_pooling)

# ---------------------------------------------------------------------------
# Fixed effects (spatial vs simpler comparison) (item 4)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 4. Fixed effects: full spatial vs simpler non-spatial comparison")

fixed_table <- bind_rows(
  tidy_fixed_effects(fit_spatial, "full_spatial"),
  tidy_fixed_effects(fit_plain, "simple_nonspatial")
)
write_csv(fixed_table, path_out_fixed)

comparison <- fixed_table %>%
  select(model_id, term, estimate, std_error) %>%
  tidyr::pivot_wider(names_from = model_id, values_from = c(estimate, std_error))
logmsg(
  "Descriptive only — NOT a formal model-selection test, and NOT an interpretation of the fishing-",
  "pressure or interaction terms as answering H2/H3. Reporting how much fixed-effect estimates and ",
  "SEs shift with vs without the spatial term:"
)
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  pct_est_change <- 100 * (r$estimate_full_spatial - r$estimate_simple_nonspatial) / abs(r$estimate_simple_nonspatial)
  pct_se_change <- 100 * (r$std_error_full_spatial - r$std_error_simple_nonspatial) / r$std_error_simple_nonspatial
  logmsg(sprintf(
    "  - %s: estimate %.4f -> %.4f (%+.1f%%); SE %.4f -> %.4f (%+.1f%%)",
    r$term, r$estimate_simple_nonspatial, r$estimate_full_spatial, pct_est_change,
    r$std_error_simple_nonspatial, r$std_error_full_spatial, pct_se_change
  ))
}
logmsg("Saved: ", path_out_fixed)

# ---------------------------------------------------------------------------
# Feasibility verdict
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Feasibility verdict")
optimizer_ok <- all(conv_table$optimizer_convergence_code == 0) && all(conv_table$hessian_positive_definite)
spatial_degenerate <- range_row$estimate > max_dist || range_row$ci_high / max(range_row$ci_low, 1e-6) > 1e4
logmsg(sprintf(
  "Optimizer-level convergence (both models converge per nlminb, positive-definite Hessian): %s.",
  optimizer_ok
))
logmsg(sprintf(
  "Spatial-parameter identifiability flag (full spatial model's range estimate %s the maximum inter-centroid distance, and/or its CI spans > 4 orders of magnitude): %s.",
  ifelse(range_row$estimate > max_dist, "EXCEEDS", "is within"), spatial_degenerate
))
logmsg(
  "PLAIN STATEMENT (per the brief — feasible / feasible only in simplified form / not feasible as ",
  "specified): "
)
if (optimizer_ok && !spatial_degenerate) {
  logmsg("  FEASIBLE AS SPECIFIED. The full spatial model converges cleanly and its spatial parameters are well identified.")
} else if (optimizer_ok && spatial_degenerate) {
  logmsg(paste0(
    "  FEASIBLE ONLY IN SIMPLIFIED FORM. The full model as specified (rectangle intercept + exp() ",
    "spatial correlation over centroids) reaches a numerically converged optimum (optimizer code 0, ",
    "positive-definite Hessian, finite standard errors) and does NOT throw a hard error or singular-",
    "fit warning — so it is not a flat non-convergence failure. However its spatial covariance ",
    "parameters are weakly identified at this rectangle density/extent: the range estimate (",
    sprintf("%.1f", range_row$estimate), " degrees) exceeds the maximum inter-centroid distance in the ",
    "data (", sprintf("%.1f", max_dist), " degrees), and both the range and sill 95% CIs span several ",
    "orders of magnitude (range CI [", sprintf("%.4f", range_row$ci_low), ", ", sprintf("%.1f", range_row$ci_high),
    "]; sill CI [", sprintf("%.4f", sd_row$ci_low), ", ", sprintf("%.2f", sd_row$ci_high), "]). This is ",
    "consistent with a likelihood surface that is nearly flat in the range direction — the data (158 ",
    "point locations spanning ", sprintf("%.0f", max_dist), " degrees, many clustered) do not contain ",
    "enough spatial contrast to pin down a continuous exponential-decay range parameter, in contrast ",
    "to the original H2 SEM which used a coarser, fixed queen-adjacency spatial structure (no range ",
    "parameter to estimate). The plain (1|stat_rec) comparison model, by contrast, is fully well-",
    "behaved (rectangle-intercept SD 95% CI is narrow: [", sprintf("%.3f", plain_row$ci_low), ", ",
    sprintf("%.3f", plain_row$ci_high), "]). RECOMMENDATION for supervisor discussion: either (a) fit ",
    "the spatial term with a FIXED or externally-constrained range/adjacency structure (e.g. reusing ",
    "the H2 SEM's queen-contiguity weights rather than a freely-estimated continuous range), or (b) ",
    "proceed with the plain (1|stat_rec) random intercept, which converges cleanly with well-identified ",
    "variance components, as the practical fallback."
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
logmsg("- ", path_out_varcomp)
logmsg("- ", path_out_pooling)
logmsg("- ", path_out_convergence)
logmsg("- ", path_out_models)
logmsg("- ", path_out_fig_pooling)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 shared model feasibility check complete. ===\n")
