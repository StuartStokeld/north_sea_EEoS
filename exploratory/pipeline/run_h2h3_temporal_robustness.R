# H2/H3 shared hierarchical model — TEMPORAL ROBUSTNESS CHECK
# (categorical phase vs continuous linear year vs smooth year)
#
# PURPOSE: the shared model so far has used a 4-level categorical phase
# variable (from the 3 structural breaks: 1989, 2001, 2008) to capture how
# fishing pressure's association with residual changes over time. This
# script checks whether that discrete structure is doing meaningful work,
# or whether a continuous year term tells the same story — the open
# question flagged in the research design proposal (Section 6, item 3).
#
# STILL A FEASIBILITY/ROBUSTNESS CHECK, NOT A RESULTS RUN. Do not treat
# either model's coefficients as the final H2/H3 answer.
#
# REUSE: same data-prep pipeline as the spatial feasibility checks
# (R/h2h3_feasibility_helpers.R::build_feasibility_data) — 158-rectangle,
# 10,464-haul panel; canonical residual; mean_ln_B_obs; log(hours_total+1).
# Only the time term changes. Base specification = the plain non-spatial
# (1 | stat_rec) random intercept recommended as the primary model after
# the spatial feasibility checks — NOT the CAR or continuous-exp spatial
# models. The categorical-phase simple_nonspatial fit is REUSED from
# Round 1 (outputs/h2h3_feasibility_model_objects.rds), not refit under
# REML; it IS re-fit under ML solely for a fair AIC/BIC comparison.
#
# MODELS:
#   1. categorical_phase (reused): residual ~ log_hours_total * phase +
#      mean_ln_B_obs + (1 | stat_rec)          [glmmTMB, REML from Round 1]
#   2. continuous_linear_year: residual ~ log_hours_total * year_c +
#      mean_ln_B_obs + (1 | stat_rec)          [glmmTMB, REML; year_c =
#      year - mean(year)]
#   3. smooth_year_gam: residual ~ s(year, k=8) +
#      s(year, by = log_hours_total, k=8) + mean_ln_B_obs +
#      s(stat_rec, bs = "re")                  [mgcv::gam, REML]
#      Varying-coefficient form chosen over the brief's ti(log_hours_total,
#      year) example so the fishing-pressure effect as a function of year
#      is directly recoverable as the by-smooth f(year) (predict type =
#      "terms" at log_hours_total = 1). Documented as a deliberate formula
#      choice, not a silent deviation.
#
# PACKAGE CHOICES / DEPENDENCIES:
#   - glmmTMB: already installed ad hoc for Round 1 (not in renv.lock).
#   - mgcv: ships with base R (here 1.9.1); NO new dependency. Chosen over
#     gamm4 because mgcv::gam with s(stat_rec, bs="re") already gives a
#     random-intercept GAMM equivalent, fits cleanly here (~90s), and
#     exposes the by-smooth / gam.check diagnostics needed for this check.
#   MUST be run with `Rscript --vanilla` (renv not activated) — same
#   environment note as Round 1 / Round 2.
#
# OUT OF SCOPE: no final H2/H3 interpretation; no refitting of CAR /
# continuous-exp spatial models with this new time structure; no changes
# to phase break years, biomass covariate, or residual sign convention.

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_temporal_robustness.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' is required and is not on the library path. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_temporal_robustness.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))
source(file.path(script_dir, "R", "h2h3_temporal_robustness_helpers.R"))

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required (ships with base R). Not found on the library path.")
}
suppressPackageStartupMessages(library(mgcv))

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_round1_models <- file.path(project_root, "outputs", "h2h3_feasibility_model_objects.rds")
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_panel <- file.path(project_root, "outputs", "h2_rectangle_panel.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
stopifnot(file.exists(path_round1_models), file.exists(path_haul), file.exists(path_panel), file.exists(path_couce_year))

path_out_effects <- file.path(project_root, "outputs", "h2h3_temporal_robustness_effects.csv")
path_out_ic <- file.path(project_root, "outputs", "h2h3_temporal_robustness_ic_comparison.csv")
path_out_slopes <- file.path(project_root, "outputs", "h2h3_temporal_robustness_fp_slope_by_year.csv")
path_out_conv <- file.path(project_root, "outputs", "h2h3_temporal_robustness_convergence.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_temporal_robustness_model_objects.rds")
path_out_fig <- file.path(fig_dir, "h2h3_temporal_robustness_fp_effect_overlay.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_temporal_robustness_run_log.md")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 shared hierarchical model — temporal robustness check run log")
logmsg("## Categorical phase vs continuous linear year vs smooth year")
logmsg("")
logmsg(
  "FEASIBILITY/ROBUSTNESS CHECK ONLY. No fishing-pressure or fishing-pressure x time coefficient ",
  "below is treated as the final H2/H3 answer. This task does not finalise the time structure as ",
  "the committed approach — that follows supervisor discussion, not this script."
)

# ---------------------------------------------------------------------------
# Tool / formula choices
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Tool and formula choices")
logmsg(paste0(
  "Base specification: the plain non-spatial (1 | stat_rec) random intercept recommended as the ",
  "primary model after the spatial feasibility checks (Rounds 1–2) — NOT the CAR or continuous-exp ",
  "spatial models. Categorical-phase simple_nonspatial fit is REUSED from Round 1 ",
  "(outputs/h2h3_feasibility_model_objects.rds), not refit under REML."
))
logmsg(paste0(
  "Continuous linear year: glmmTMB, residual ~ log_hours_total * year_c + mean_ln_B_obs + ",
  "(1 | stat_rec), REML = TRUE. year_c = year - mean(year) over the analysis hauls (mean ≈ 1999.5) ",
  "for numerical stability; the centring constant is recorded and used when recovering the ",
  "fishing-pressure slope as a function of calendar year."
))
logmsg(paste0(
  "Smooth year (GAM): mgcv::gam (version ", as.character(utils::packageVersion("mgcv")),
  "), residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) + mean_ln_B_obs + ",
  "s(stat_rec, bs = \"re\"), method = \"REML\". This is a VARYING-COEFFICIENT form: the by-smooth ",
  "contributes log_hours_total * f(year) to the linear predictor, so f(year) IS the ",
  "fishing-pressure effect as a function of year (recovered via predict(type = \"terms\") at ",
  "log_hours_total = 1). Chosen over the brief's ti(log_hours_total, year) example for that ",
  "direct recoverability; ti() would give a more general 2-D interaction surface that is harder ",
  "to read as \"the FP effect by year\" without additional slicing. mgcv chosen over gamm4 ",
  "because s(stat_rec, bs = \"re\") already implements the random intercept inside gam(), fits ",
  "cleanly here, and exposes gam.check() diagnostics. mgcv ships with base R — NO new ad-hoc ",
  "dependency (unlike glmmTMB / spaMM / strucchange)."
))
logmsg(
  "MUST be run with `Rscript --vanilla` (renv not activated) — same environment note as Round 1/2."
)

# ---------------------------------------------------------------------------
# Data (reuse Round 1 analysis frame where available; rebuild if needed)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data")

round1 <- readRDS(path_round1_models)
fit_phase_reml <- round1$fit_plain
dat <- round1$data

# Sanity: confirm the reused object is the expected simple_nonspatial fit
phase_fml <- deparse(stats::formula(fit_phase_reml), width.cutoff = 500L)
logmsg("Reused Round 1 simple_nonspatial formula: ", paste(phase_fml, collapse = " "))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d (identical to Round 1 / Round 2 universe).",
  nrow(dat), dplyr::n_distinct(dat$stat_rec), min(dat$year), max(dat$year)
))

yc <- build_year_centred(dat$year)
dat$year_c <- yc$year_c
year_centre <- yc$year_centre
dat$stat_rec <- factor(dat$stat_rec)
logmsg(sprintf("year_c centre (mean calendar year over analysis hauls): %.4f", year_centre))

years_grid <- seq.int(min(dat$year), max(dat$year))

# ---------------------------------------------------------------------------
# Fit models
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fitting")

logmsg("Fitting continuous_linear_year (glmmTMB, REML) ...")
time_linear <- system.time({
  fit_linear_reml <- glmmTMB(
    residual ~ log_hours_total * year_c + mean_ln_B_obs + (1 | stat_rec),
    data = dat,
    REML = TRUE
  )
})
logmsg(sprintf("continuous_linear_year REML fit time: %.2f sec.", time_linear["elapsed"]))

logmsg("Fitting smooth_year_gam (mgcv::gam, REML; expect ~1–2 min) ...")
time_gam <- system.time({
  fit_gam_reml <- gam(
    residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) +
      mean_ln_B_obs + s(stat_rec, bs = "re"),
    data = dat,
    method = "REML"
  )
})
logmsg(sprintf("smooth_year_gam REML fit time: %.2f sec.", time_gam["elapsed"]))

logmsg("Refitting all three models under ML for a fair AIC/BIC comparison ...")
time_ml <- system.time({
  fit_phase_ml <- glmmTMB(
    residual ~ log_hours_total * phase + mean_ln_B_obs + (1 | stat_rec),
    data = dat,
    REML = FALSE
  )
  fit_linear_ml <- glmmTMB(
    residual ~ log_hours_total * year_c + mean_ln_B_obs + (1 | stat_rec),
    data = dat,
    REML = FALSE
  )
  fit_gam_ml <- gam(
    residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) +
      mean_ln_B_obs + s(stat_rec, bs = "re"),
    data = dat,
    method = "ML"
  )
})
logmsg(sprintf("ML refits (all three) total time: %.2f sec.", time_ml["elapsed"]))

saveRDS(
  list(
    fit_phase_reml = fit_phase_reml,
    fit_linear_reml = fit_linear_reml,
    fit_gam_reml = fit_gam_reml,
    fit_phase_ml = fit_phase_ml,
    fit_linear_ml = fit_linear_ml,
    fit_gam_ml = fit_gam_ml,
    data = dat,
    year_centre = year_centre
  ),
  path_out_models
)
logmsg("Saved model objects: ", path_out_models)

# ---------------------------------------------------------------------------
# 1. Fit-quality / IC comparison
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Fit-quality comparison (AIC / BIC)")

ic_table <- dplyr::bind_rows(
  extract_ic(fit_phase_reml, "categorical_phase", "REML"),
  extract_ic(fit_linear_reml, "continuous_linear_year", "REML"),
  extract_ic(fit_gam_reml, "smooth_year_gam", "REML"),
  extract_ic(fit_phase_ml, "categorical_phase", "ML"),
  extract_ic(fit_linear_ml, "continuous_linear_year", "ML"),
  extract_ic(fit_gam_ml, "smooth_year_gam", "ML")
)
write_csv(ic_table, path_out_ic)

logmsg(paste0(
  "CAVEAT on information criteria: REML AIC/BIC are NOT strictly comparable across models with ",
  "different fixed-effect structures (the REML criterion profiles out fixed effects, so the ",
  "\"likelihood\" being compared is not the same quantity). The PRIMARY comparison below therefore ",
  "uses the ML refits. A further caveat specific to the GAM: mgcv's AIC uses an effective-degrees-",
  "of-freedom penalty based on the smooth's EDF, which is a DIFFERENT penalty construction from ",
  "glmmTMB's parametric AIC — so even the ML AIC comparison between the GAM and the two glmmTMB ",
  "models is approximate, not an exact nested-model likelihood-ratio comparison. Reported with ",
  "that caveat, not as a formal model-selection test."
))

ic_ml <- ic_table %>% filter(estimation_method == "ML") %>% arrange(aic)
logmsg("ML AIC/BIC (primary comparison):")
for (i in seq_len(nrow(ic_ml))) {
  r <- ic_ml[i, ]
  logmsg(sprintf(
    "  - %s: AIC = %.1f, BIC = %.1f, logLik = %.1f, n_params_reported = %s",
    r$model_id, r$aic, r$bic, r$logLik, r$n_params_reported
  ))
}
ic_reml <- ic_table %>% filter(estimation_method == "REML") %>% arrange(aic)
logmsg("REML AIC/BIC (reported for completeness; NOT the primary comparison — see caveat above):")
for (i in seq_len(nrow(ic_reml))) {
  r <- ic_reml[i, ]
  logmsg(sprintf(
    "  - %s: AIC = %.1f, BIC = %.1f, logLik = %.1f, n_params_reported = %s",
    r$model_id, r$aic, r$bic, r$logLik, r$n_params_reported
  ))
}
logmsg("Saved: ", path_out_ic)

# ---------------------------------------------------------------------------
# 2. Does the fishing-pressure × time interaction still show up?
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Does the fishing-pressure × time interaction still show up?")

effects_table <- dplyr::bind_rows(
  tidy_parametric_effects(fit_phase_reml, "categorical_phase"),
  tidy_parametric_effects(fit_linear_reml, "continuous_linear_year"),
  tidy_gam_effects(fit_gam_reml, "smooth_year_gam")
)
write_csv(effects_table, path_out_effects)

# Linear interaction
lin_int <- effects_table %>%
  filter(model_id == "continuous_linear_year", term == "log_hours_total:year_c")
logmsg(sprintf(
  paste0(
    "Continuous linear year — interaction log_hours_total:year_c: estimate = %.6f, SE = %.6f, ",
    "z = %.3f, p = %.3e. (Positive estimate means the fishing-pressure slope increases with year; ",
    "NOT interpreted here as an H2/H3 finding.)"
  ),
  lin_int$estimate, lin_int$std_error, lin_int$statistic, lin_int$p_value
))

# GAM smooth interaction
gam_int <- effects_table %>%
  filter(model_id == "smooth_year_gam", grepl("log_hours_total", term), term_type == "smooth")
logmsg(sprintf(
  paste0(
    "Smooth year GAM — by-smooth s(year):log_hours_total: edf = %.3f, Ref.df = %.3f, F = %.3f, ",
    "approx p = %.3e. (Highly significant approx p indicates a non-flat fishing-pressure effect ",
    "across years; NOT interpreted here as an H2/H3 finding.)"
  ),
  gam_int$edf, gam_int$ref_df, gam_int$statistic, gam_int$p_value
))

# Phase-model interaction terms (for context, not as H2/H3 results)
phase_ints <- effects_table %>%
  filter(model_id == "categorical_phase", grepl("^log_hours_total:phase", term))
logmsg("Categorical phase — fishing-pressure × phase interaction terms (reused Round 1 REML fit):")
for (i in seq_len(nrow(phase_ints))) {
  r <- phase_ints[i, ]
  logmsg(sprintf(
    "  - %s: estimate = %.4f, SE = %.4f, z = %.3f, p = %.3e",
    r$term, r$estimate, r$std_error, r$statistic, r$p_value
  ))
}
logmsg("Saved: ", path_out_effects)

# ---------------------------------------------------------------------------
# 3. Visual comparison — overlay figure
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Visual comparison: fishing-pressure effect vs year (overlay)")

phase_slopes <- extract_phase_fp_slopes(fit_phase_reml)
linear_slopes <- extract_linear_fp_slopes(fit_linear_reml, years_grid, year_centre)
gam_slopes <- extract_gam_fp_slopes(fit_gam_reml, years_grid, dat)
phase_year <- expand_phase_slopes_to_years(phase_slopes, years_grid)

slopes_table <- dplyr::bind_rows(phase_year, linear_slopes, gam_slopes)
write_csv(slopes_table, path_out_slopes)
logmsg("Per-phase fishing-pressure slopes (categorical phase model):")
for (i in seq_len(nrow(phase_slopes))) {
  r <- phase_slopes[i, ]
  logmsg(sprintf(
    "  - %s (%d–%d): slope = %.4f, SE = %.4f, 95%% CI = [%.4f, %.4f]",
    r$phase, r$year_start, r$year_end, r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi
  ))
}
logmsg("Saved: ", path_out_slopes)

# Build plot data: continuous curves + phase as step segments
plot_cont <- slopes_table %>%
  filter(model_id %in% c("continuous_linear_year", "smooth_year_gam")) %>%
  mutate(
    model_label = dplyr::recode(
      model_id,
      continuous_linear_year = "Continuous linear year",
      smooth_year_gam = "Smooth year (GAM)"
    )
  )
# Phase step: one horizontal segment per phase
plot_phase_seg <- phase_slopes %>%
  mutate(model_label = "Categorical phase (step)")

col_map <- c(
  "Categorical phase (step)" = "#d73027",
  "Continuous linear year" = "#4575b4",
  "Smooth year (GAM)" = "#1a9850"
)

p_overlay <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  # Phase step bands (CI) and segments
  geom_rect(
    data = plot_phase_seg,
    aes(
      xmin = year_start - 0.5, xmax = year_end + 0.5,
      ymin = fp_slope_lo, ymax = fp_slope_hi,
      fill = model_label
    ),
    alpha = 0.12, colour = NA
  ) +
  geom_segment(
    data = plot_phase_seg,
    aes(
      x = year_start - 0.5, xend = year_end + 0.5,
      y = fp_slope, yend = fp_slope,
      colour = model_label
    ),
    linewidth = 1.1
  ) +
  # Continuous / smooth ribbons and lines
  geom_ribbon(
    data = plot_cont,
    aes(x = year, ymin = fp_slope_lo, ymax = fp_slope_hi, fill = model_label),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = plot_cont,
    aes(x = year, y = fp_slope, colour = model_label),
    linewidth = 1.0
  ) +
  # Vertical lines at structural-break years (for visual reference only)
  geom_vline(xintercept = c(1989, 2001, 2008), linetype = "dotted", colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_fill_manual(values = col_map, name = NULL) +
  labs(
    x = "Year",
    y = "Fishing-pressure effect\n(∂ residual / ∂ log(hours+1))",
    title = "Temporal robustness: fishing-pressure effect vs year",
    subtitle = paste0(
      "Categorical phase (step) vs continuous linear vs smooth GAM — ",
      "feasibility/robustness check only; not an H2/H3 result. ",
      "Dotted lines = structural-break years (1989, 2001, 2008)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.5, colour = "grey30")
  )

ggsave(path_out_fig, p_overlay, width = 10, height = 5.5, dpi = 150)
logmsg("Saved figure: ", path_out_fig)
logmsg(paste0(
  "How to read the figure: agreement between the red step function and the green smooth curve ",
  "(within CIs, and especially whether the step jumps land near where the smooth bends) would ",
  "support the categorical phase structure as a reasonable simplification of a smoother underlying ",
  "pattern; material disagreement (e.g. the smooth staying flat across a phase boundary the step ",
  "treats as a level shift, or vice versa) would argue the discrete phases capture something a ",
  "continuous term misses — or overstate a smooth change as a break. Verdict below."
))

# ---------------------------------------------------------------------------
# 4. Convergence diagnostics
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 4. Convergence diagnostics")

conv_phase <- glmmtmb_convergence_brief(fit_phase_reml, "categorical_phase")
conv_linear <- glmmtmb_convergence_brief(fit_linear_reml, "continuous_linear_year")
conv_gam <- mgcv_convergence_brief(fit_gam_reml, "smooth_year_gam")
# Align columns (gam has gam_check_text; glmmTMB rows get NA for it)
conv_phase$gam_check_text <- NA_character_
conv_linear$gam_check_text <- NA_character_
conv_table <- dplyr::bind_rows(conv_phase, conv_linear, conv_gam)
write_csv(conv_table, path_out_conv)

for (i in seq_len(nrow(conv_table))) {
  r <- conv_table[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s (%s): converged = %s; hessian_positive_definite = %s; max|gradient| = %s; ",
      "optimizer message = '%s'"
    ),
    r$model_id, r$framework, r$converged, r$hessian_positive_definite,
    ifelse(is.finite(r$max_abs_gradient), sprintf("%.3e", r$max_abs_gradient), "NA"),
    r$optimizer_message
  ))
  if (!is.na(r$gam_check_text) && nzchar(r$gam_check_text)) {
    logmsg("    gam.check() (basis-dimension / k check): ", r$gam_check_text)
  }
}
logmsg("Saved: ", path_out_conv)

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Feasibility / robustness verdict")

# Simple quantitative agreement metrics between phase step and the two continuous forms
agree <- slopes_table %>%
  select(model_id, year, fp_slope) %>%
  tidyr::pivot_wider(names_from = model_id, values_from = fp_slope)
cor_phase_lin <- suppressWarnings(cor(agree$categorical_phase, agree$continuous_linear_year, use = "complete.obs"))
cor_phase_gam <- suppressWarnings(cor(agree$categorical_phase, agree$smooth_year_gam, use = "complete.obs"))
cor_lin_gam <- suppressWarnings(cor(agree$continuous_linear_year, agree$smooth_year_gam, use = "complete.obs"))
# Sign agreement: same sign of the FP effect in each year
sign_agree_phase_lin <- mean(sign(agree$categorical_phase) == sign(agree$continuous_linear_year), na.rm = TRUE)
sign_agree_phase_gam <- mean(sign(agree$categorical_phase) == sign(agree$smooth_year_gam), na.rm = TRUE)

logmsg(sprintf(
  paste0(
    "Year-by-year correlation of the fitted fishing-pressure effect: phase vs linear = %.3f; ",
    "phase vs GAM = %.3f; linear vs GAM = %.3f. Fraction of years with the same sign of the FP ",
    "effect: phase vs linear = %.1f%%; phase vs GAM = %.1f%%."
  ),
  cor_phase_lin, cor_phase_gam, cor_lin_gam,
  100 * sign_agree_phase_lin, 100 * sign_agree_phase_gam
))

ml_best <- ic_ml$model_id[1]
logmsg(sprintf("Lowest-ML-AIC model: %s (see caveat above — not a formal selection test).", ml_best))

# Interaction significance flags (descriptive)
lin_sig <- isTRUE(lin_int$p_value < 0.05)
gam_sig <- isTRUE(gam_int$p_value < 0.05)
phase_any_sig <- any(phase_ints$p_value < 0.05, na.rm = TRUE)

logmsg("PLAIN STATEMENT (per the brief — categorical phase is a reasonable simplification / ",
       "materially disagrees with the continuous alternatives):")

# Verdict logic: if the FP×time interaction shows up in all three AND the phase step
# broadly tracks the smooth (high correlation / sign agreement), phase is a reasonable
# simplification; if the smooth shows a qualitatively different pattern (e.g. opposite
# late-period slope, or no break-like bend where the step jumps), flag material disagreement.
if (lin_sig && gam_sig && phase_any_sig && cor_phase_gam >= 0.5 && sign_agree_phase_gam >= 0.7) {
  logmsg(paste0(
    "  CATEGORICAL PHASE IS A REASONABLE SIMPLIFICATION. The fishing-pressure × time interaction ",
    "shows up under all three time structures (phase interactions, linear interaction p = ",
    sprintf("%.2e", lin_int$p_value), ", GAM by-smooth approx p = ", sprintf("%.2e", gam_int$p_value),
    "), and the phase model's step-function FP effect tracks the continuous alternatives at the ",
    "year-by-year level (phase–GAM correlation = ", sprintf("%.2f", cor_phase_gam),
    "; same-sign years = ", sprintf("%.0f%%", 100 * sign_agree_phase_gam),
    "). The GAM is more flexible and preferred by ML AIC, but it does not tell a qualitatively ",
    "different story from the discrete phases — the step function is a coarse but recognisable ",
    "summary of the smoother curve. RECOMMENDATION for supervisor discussion: keep the ",
    "categorical phase structure as the primary (interpretable, break-aligned) specification; ",
    "treat the linear and GAM fits as supporting robustness checks, not replacements."
  ))
} else if (lin_sig && gam_sig && cor_phase_gam < 0.3) {
  logmsg(paste0(
    "  CATEGORICAL PHASE MATERIALLY DISAGREES WITH THE CONTINUOUS ALTERNATIVES. The FP × time ",
    "interaction is present under all three structures, but the phase model's step-function FP ",
    "effect correlates only weakly with the GAM curve (phase–GAM correlation = ",
    sprintf("%.2f", cor_phase_gam),
    "). The discrete breaks may be imposing level shifts the smooth data do not support, or the ",
    "smooth may be blurring genuine break-year jumps — either way, the two time structures are ",
    "not interchangeable. RECOMMENDATION for supervisor discussion: do not treat the categorical ",
    "phase model as a harmless simplification of a continuous year effect; choose deliberately ",
    "between them (or report both) based on whether the structural-break motivation is theoretically ",
    "primary."
  ))
} else {
  logmsg(paste0(
    "  MIXED / CAVEATED AGREEMENT. FP × time interaction present in linear model: ", lin_sig,
    "; in GAM by-smooth: ", gam_sig, "; in at least one phase interaction: ", phase_any_sig,
    ". Phase–GAM year-by-year correlation = ", sprintf("%.2f", cor_phase_gam),
    "; same-sign years = ", sprintf("%.0f%%", 100 * sign_agree_phase_gam),
    ". Read the overlay figure and the per-phase / per-year slope table alongside the IC comparison ",
    "(lowest ML AIC = ", ml_best, ") before treating the categorical phase structure as either a ",
    "safe simplification or a material disagreement — the quantitative agreement metrics sit in a ",
    "middle ground that a single binary verdict would overstate. RECOMMENDATION for supervisor ",
    "discussion: present the overlay figure as the primary evidence; keep phase as the ",
    "break-aligned primary specification only if the figure shows the step function landing on the ",
    "same broad pattern as the smooth, and otherwise report sensitivity to the time-structure choice."
  ))
}
logmsg(
  "This verdict is about TIME-STRUCTURE AGREEMENT only — it does not interpret, and is not based ",
  "on treating, the sign or magnitude of any fishing-pressure coefficient as an H2/H3 finding."
)

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_effects)
logmsg("- ", path_out_ic)
logmsg("- ", path_out_slopes)
logmsg("- ", path_out_conv)
logmsg("- ", path_out_models)
logmsg("- ", path_out_fig)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 temporal robustness check complete. ===\n")
