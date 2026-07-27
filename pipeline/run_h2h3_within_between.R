# H2/H3 within-between fishing-pressure decomposition
#
# PURPOSE: separate persistent between-rectangle fishing pressure (H2) from
# within-rectangle year-to-year deviations (H3), following the within-between
# panel approach (Bell & Jones 2015). The blended log_hours_total * phase
# primary model conflates these; this script answers each with its own
# coefficient.
#
# TRANSFORMATION#   FP_between_i  = mean_t log(hours_it + 1)   [time-invariant per rectangle]
#   FP_within_it  = log(hours_it + 1) - FP_between_i
#
# PRIMARY: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)
# Same data / phase / RE / no biomass as the biomass-free results run.
#
# SENSITIVITY: CAR with the same decomposed FE (if spaMM + Round 2 adjMatrix
# available). GAM continuous-year version DEFERRED (flagged in run log).
#
# Run: Rscript --vanilla pipeline/run_h2h3_within_between.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_within_between.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: Rscript --vanilla pipeline/run_h2h3_within_between.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))
source(file.path(script_dir, "R", "h2h3_results_helpers.R"))
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

has_spamm <- requireNamespace("spaMM", quietly = TRUE)
if (has_spamm) suppressPackageStartupMessages(library(spaMM))

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_results_models <- file.path(project_root, "outputs", "h2h3_results_model_objects.rds")
path_blended_slopes <- file.path(project_root, "outputs", "h2h3_results_fp_slopes_by_phase.csv")
path_round2 <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
stopifnot(file.exists(path_results_models), file.exists(path_blended_slopes))

path_out_fe <- file.path(project_root, "outputs", "h2h3_wb_primary_fixed_effects.csv")
path_out_slopes <- file.path(project_root, "outputs", "h2h3_wb_fp_slopes_by_phase.csv")
path_out_wald <- file.path(project_root, "outputs", "h2h3_wb_wald_tests.csv")
path_out_compare <- file.path(project_root, "outputs", "h2h3_wb_blended_comparison.csv")
path_out_pooling <- file.path(project_root, "outputs", "h2h3_wb_partial_pooling.csv")
path_out_fit <- file.path(project_root, "outputs", "h2h3_wb_model_fit.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_out_fig_pooling <- file.path(fig_dir, "h2h3_wb_partial_pooling.png")
path_out_fig_slopes <- file.path(fig_dir, "h2h3_wb_fp_slopes_by_phase.png")
path_out_car_fe <- file.path(project_root, "outputs", "h2h3_wb_car_fixed_effects.csv")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_wb_run_log.md")
path_out_session <- file.path(project_root, "outputs", "h2h3_wb_sessionInfo.txt")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 within-between fishing-pressure decomposition — run log")
logmsg("")
logmsg(
  "Separates between-rectangle (H2) and within-rectangle (H3) fishing-pressure ",
  "variation. Biomass excluded. Same panel / phase / (1|stat_rec) as the biomass-free ",
  "primary results run. No ecological interpretation beyond labelling which coefficients ",
  "answer which hypothesis."
)

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf(
  "glmmTMB %s; spaMM available = %s%s",
  as.character(utils::packageVersion("glmmTMB")),
  has_spamm,
  if (has_spamm) paste0(" (", utils::packageVersion("spaMM"), ")") else ""
))

# ---------------------------------------------------------------------------
# Data + decomposition
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data and FP within-between transformation")

res <- readRDS(path_results_models)
dat <- res$data
if (!"log_hours_total" %in% names(dat) || !"phase" %in% names(dat)) {
  stop("Results model objects data missing log_hours_total or phase.")
}
dat$stat_rec <- factor(dat$stat_rec)
dat <- add_fp_within_between(dat)

checks <- fp_within_between_checks(dat)
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d (reused from biomass-free results RDS).",
  nrow(dat), n_distinct(dat$stat_rec), min(dat$year), max(dat$year)
))
logmsg(sprintf(
  "FP_between: time-invariant rectangle mean of log(hours+1); var = %.4f; range = [%.3f, %.3f].",
  checks$var_between, min(dat$FP_between), max(dat$FP_between)
))
logmsg(sprintf(
  "FP_within: haul-level deviation from rectangle mean; var = %.4f; range = [%.3f, %.3f].",
  checks$var_within, min(dat$FP_within), max(dat$FP_within)
))
logmsg("")
logmsg("### Sanity check: FP_between vs FP_within correlation")
logmsg(sprintf(
  "Correlation(FP_between, FP_within) = %.6g (expected ≈ 0 by construction).",
  checks$cor_between_within
))
logmsg(sprintf(
  "Max |within-rectangle mean of FP_within| = %.3e; mean abs = %.3e (expected ≈ 0).",
  checks$max_abs_within_rect_mean, checks$mean_abs_within_rect_mean
))
if (abs(checks$cor_between_within) > 0.05) {
  logmsg("FLAG: between-within correlation exceeds 0.05 — unexpected; investigate.")
} else {
  logmsg("Sanity check passed: between and within components are effectively uncorrelated.")
}

# ---------------------------------------------------------------------------
# Fit primary
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Primary model fit")
formula_wb <- residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)
logmsg("Formula: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec) [REML]")
logmsg(
  "H2 answered by FP_between (+ phase interactions / phase-specific between slopes). ",
  "H3 answered by FP_within (+ phase interactions / phase-specific within slopes)."
)

time_pri <- system.time({
  fit_wb <- glmmTMB(formula_wb, data = dat, REML = TRUE)
})
logmsg(sprintf("Primary within-between fit time: %.2f sec.", time_pri["elapsed"]))

# ---------------------------------------------------------------------------
# CAR sensitivity (if possible)
# ---------------------------------------------------------------------------
fit_wb_car <- NULL
car_fitted <- FALSE
logmsg("")
logmsg("## Sensitivity: CAR with decomposed terms")

if (!has_spamm) {
  logmsg("spaMM not available — CAR sensitivity SKIPPED.")
} else if (!file.exists(path_round2)) {
  logmsg("Round 2 RDS missing — CAR sensitivity SKIPPED.")
} else {
  round2 <- readRDS(path_round2)
  if (is.null(round2$adjMatrix)) {
    logmsg("Round 2 adjMatrix missing — CAR sensitivity SKIPPED.")
  } else {
    adjMatrix <- round2$adjMatrix
    dat_car <- dat
    dat_car$stat_rec <- factor(as.character(dat_car$stat_rec), levels = rownames(adjMatrix))
    formula_car <- residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)
    logmsg("Fitting CAR: residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)")
    time_car <- system.time({
      fit_wb_car <- spaMM::fitme(
        formula_car, data = dat_car, adjMatrix = adjMatrix, method = "REML"
      )
    })
    logmsg(sprintf("CAR within-between fit time: %.2f sec.", time_car["elapsed"]))
    car_fitted <- TRUE
    fe_car <- tidy_fixed_effects_spamm(fit_wb_car, "wb_car") %>%
      annotate_wb_fixed_effects()
    write_csv(fe_car, path_out_car_fe)
    logmsg("Saved CAR FE table: ", path_out_car_fe)
  }
}

logmsg("")
logmsg("## Deferred sensitivity")
logmsg(
  "GAM / continuous-year version of the within-between decomposition is DEFERRED ",
  "(timeline). Not silently dropped — flag here for a follow-up if needed."
)

saveRDS(
  list(
    fit_wb = fit_wb,
    fit_wb_car = fit_wb_car,
    data = dat,
    formula_wb = formula_wb,
    car_fitted = car_fitted,
    checks = checks
  ),
  path_out_models
)
logmsg("Saved model objects: ", path_out_models)

# ---------------------------------------------------------------------------
# 1. Full FE table
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Primary fixed-effects table (H2/H3 labelled)")

fe <- tidy_fixed_effects(fit_wb, "wb_primary") %>%
  annotate_wb_fixed_effects() %>%
  mutate(
    ci_lo = estimate - 1.96 * std_error,
    ci_hi = estimate + 1.96 * std_error
  ) %>%
  select(
    model_id, term, term_plain, hypothesis, hypothesis_group,
    estimate, std_error, ci_lo, ci_hi, statistic, p_value
  )
write_csv(fe, path_out_fe)

logmsg("Fixed effects:")
for (i in seq_len(nrow(fe))) {
  r <- fe[i, ]
  logmsg(sprintf(
    "  - %-35s est=%+.4f SE=%.4f z=%+.2f p=%.3e  [%s]",
    r$term, r$estimate, r$std_error, r$statistic, r$p_value, r$hypothesis_group
  ))
}
logmsg("Saved: ", path_out_fe)

# Phase-specific slopes (8 = 2 components × 4 phases)
slopes_between <- extract_wb_phase_slopes(
  fit_wb, "FP_between", "wb_primary", "H2_spatial_between"
)
slopes_within <- extract_wb_phase_slopes(
  fit_wb, "FP_within", "wb_primary", "H3_temporal_within"
)
slopes_all <- bind_rows(slopes_between, slopes_within)

if (car_fitted) {
  slopes_all <- bind_rows(
    slopes_all,
    extract_wb_phase_slopes_spamm(fit_wb_car, "FP_between", "wb_car", "H2_spatial_between"),
    extract_wb_phase_slopes_spamm(fit_wb_car, "FP_within", "wb_car", "H3_temporal_within")
  )
}
write_csv(slopes_all, path_out_slopes)
logmsg("")
logmsg("Phase-specific slopes (primary):")
pri_slopes <- slopes_all %>% filter(model_id == "wb_primary")
for (i in seq_len(nrow(pri_slopes))) {
  r <- pri_slopes[i, ]
  logmsg(sprintf(
    "  - [%s] %s / %s: slope=%+.4f SE=%.4f CI=[%+.4f, %+.4f] p=%.3e",
    r$hypothesis_group, r$component, r$phase,
    r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi, r$p_value
  ))
}
logmsg("Saved: ", path_out_slopes)

# ---------------------------------------------------------------------------
# 2. Joint Wald tests
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Joint Wald tests")

b <- glmmTMB::fixef(fit_wb)$cond
V <- stats::vcov(fit_wb)$cond
terms_between_all <- wb_all_slope_terms(names(b), "FP_between")
terms_within_all <- wb_all_slope_terms(names(b), "FP_within")
terms_within_int <- wb_interaction_terms(names(b), "FP_within")
terms_between_int <- wb_interaction_terms(names(b), "FP_between")

wald_h2 <- wald_joint_zero(b, V, terms_between_all)
wald_h3_slopes <- wald_joint_zero(b, V, terms_within_all)
wald_h3_int <- wald_joint_zero(b, V, terms_within_int)
wald_h2_int <- wald_joint_zero(b, V, terms_between_int)

wald_table <- data.frame(
  test_id = c(
    "H2_all_between_phase_slopes_joint_zero",
    "H3_all_within_phase_slopes_joint_zero",
    "H3_within_phase_interactions_joint_zero",
    "H2_between_phase_interactions_joint_zero"
  ),
  hypothesis = c(
    "H2 — all FP_between phase slopes jointly zero",
    "H3 — all FP_within phase slopes jointly zero",
    "H3 — FP_within × phase interactions jointly zero (change vs reference)",
    "H2 — FP_between × phase interactions jointly zero (change vs reference; supplementary)"
  ),
  statistic_chisq = c(
    wald_h2$statistic, wald_h3_slopes$statistic,
    wald_h3_int$statistic, wald_h2_int$statistic
  ),
  df = c(wald_h2$df, wald_h3_slopes$df, wald_h3_int$df, wald_h2_int$df),
  p_value = c(wald_h2$p_value, wald_h3_slopes$p_value, wald_h3_int$p_value, wald_h2_int$p_value),
  terms = c(
    paste(wald_h2$terms, collapse = "; "),
    paste(wald_h3_slopes$terms, collapse = "; "),
    paste(wald_h3_int$terms, collapse = "; "),
    paste(wald_h2_int$terms, collapse = "; ")
  ),
  stringsAsFactors = FALSE
)
write_csv(wald_table, path_out_wald)
for (i in seq_len(nrow(wald_table))) {
  r <- wald_table[i, ]
  logmsg(sprintf(
    "  - %s: chi2(%d)=%.3f, p=%.3e",
    r$test_id, r$df, r$statistic_chisq, r$p_value
  ))
}
logmsg("Saved: ", path_out_wald)

# ---------------------------------------------------------------------------
# 3. Comparison to blended model
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Comparison to blended log_hours_total model")

blended <- read_csv(path_blended_slopes, show_col_types = FALSE) %>%
  filter(model_id == "primary_plain_re") %>%
  select(phase, blended_slope = fp_slope, blended_se = fp_slope_se,
         blended_lo = fp_slope_lo, blended_hi = fp_slope_hi)

cmp <- pri_slopes %>%
  select(hypothesis_group, component, phase, fp_slope, fp_slope_se, fp_slope_lo, fp_slope_hi, p_value) %>%
  tidyr::pivot_wider(
    id_cols = phase,
    names_from = component,
    values_from = c(fp_slope, fp_slope_se, fp_slope_lo, fp_slope_hi, p_value)
  ) %>%
  left_join(blended, by = "phase") %>%
  mutate(
    # Blended ≈ weighted mix; report both components beside it (do not sum as identity)
    note = "blended_slope is from residual ~ log_hours_total * phase + (1|stat_rec); not equal to between+within"
  ) %>%
  select(
    phase,
    blended_slope, blended_se, blended_lo, blended_hi,
    FP_between_slope = fp_slope_FP_between,
    FP_between_se = fp_slope_se_FP_between,
    FP_between_lo = fp_slope_lo_FP_between,
    FP_between_hi = fp_slope_hi_FP_between,
    FP_between_p = p_value_FP_between,
    FP_within_slope = fp_slope_FP_within,
    FP_within_se = fp_slope_se_FP_within,
    FP_within_lo = fp_slope_lo_FP_within,
    FP_within_hi = fp_slope_hi_FP_within,
    FP_within_p = p_value_FP_within,
    note
  )

write_csv(cmp, path_out_compare)
logmsg("Per-phase blended vs decomposed slopes:")
for (i in seq_len(nrow(cmp))) {
  r <- cmp[i, ]
  logmsg(sprintf(
    "  - %s: blended=%+.4f | between(H2)=%+.4f (p=%.3e) | within(H3)=%+.4f (p=%.3e)",
    r$phase, r$blended_slope, r$FP_between_slope, r$FP_between_p,
    r$FP_within_slope, r$FP_within_p
  ))
}
logmsg("Saved: ", path_out_compare)

# ---------------------------------------------------------------------------
# 5. Diagnostics (convergence, pooling, R2)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 5. Diagnostics (primary within-between model)")

conv <- glmmtmb_convergence_report(fit_wb, "wb_primary")
logmsg(sprintf(
  "Convergence: code=%s, Hessian PD=%s, max|grad|=%.3e, any NA SE=%s",
  conv$optimizer_convergence_code, conv$hessian_positive_definite,
  conv$max_abs_gradient, conv$any_na_std_error
))

r2 <- nakagawa_r2_glmmtmb(fit_wb)
vc <- extract_variance_components(fit_wb, "wb_primary", spatial = FALSE)
fit_stats <- data.frame(
  model_id = "wb_primary",
  n_hauls = nrow(dat),
  n_rectangles = n_distinct(dat$stat_rec),
  r2_marginal = r2$r2_marginal,
  r2_conditional = r2$r2_conditional,
  var_fixed = r2$var_fixed,
  var_random = r2$var_random,
  var_residual = r2$var_residual,
  residual_sd = vc$residual_sd,
  rectangle_intercept_sd = vc$estimate,
  logLik = as.numeric(stats::logLik(fit_wb)),
  aic = as.numeric(stats::AIC(fit_wb)),
  bic = as.numeric(stats::BIC(fit_wb)),
  cor_FP_between_FP_within = checks$cor_between_within,
  stringsAsFactors = FALSE
)
write_csv(fit_stats, path_out_fit)
logmsg(sprintf(
  "Nakagawa R2: marginal=%.4f, conditional=%.4f; rect SD=%.4f; resid SD=%.4f",
  r2$r2_marginal, r2$r2_conditional, vc$estimate, vc$residual_sd
))
logmsg("Saved: ", path_out_fit)

pooling <- extract_partial_pooling(fit_wb, dat, "wb_primary", spatial = FALSE)
write_csv(pooling, path_out_pooling)
cor_pool <- suppressWarnings(cor(log(pooling$n_hauls), pooling$shrinkage_ratio, method = "spearman"))
logmsg(sprintf(
  "Partial pooling: Spearman cor(log n_hauls, shrinkage_ratio)=%.3f across %d rectangles.",
  cor_pool, nrow(pooling)
))

p_pool <- ggplot(pooling, aes(x = n_hauls, y = random_intercept)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(alpha = 0.7, colour = "#4575b4") +
  scale_x_log10() +
  labs(
    x = "Haul count per rectangle (log scale)",
    y = "Estimated rectangle random intercept",
    title = "Within-between model: rectangle random intercepts vs haul count",
    subtitle = "Primary: residual ~ FP_between * phase + FP_within * phase + (1|stat_rec)"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_pooling, p_pool, width = 8, height = 5, dpi = 150)
logmsg("Saved: ", path_out_pooling)
logmsg("Saved figure: ", path_out_fig_pooling)

# Slope comparison figure: between vs within vs blended
plot_slopes <- bind_rows(
  blended %>%
    transmute(
      phase, component = "blended (old)",
      fp_slope = blended_slope, fp_slope_lo = blended_lo, fp_slope_hi = blended_hi
    ),
  pri_slopes %>%
    transmute(
      phase,
      component = dplyr::recode(
        component,
        FP_between = "FP_between (H2)",
        FP_within = "FP_within (H3)"
      ),
      fp_slope, fp_slope_lo, fp_slope_hi
    )
) %>%
  mutate(
    phase = factor(phase, levels = c("1985-1988", "1989-2000", "2001-2007", "2008-2015")),
    component = factor(
      component,
      levels = c("blended (old)", "FP_between (H2)", "FP_within (H3)")
    )
  )

p_slopes <- ggplot(plot_slopes, aes(x = phase, y = fp_slope, colour = component)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(ymin = fp_slope_lo, ymax = fp_slope_hi),
    position = position_dodge(width = 0.45), width = 0.18, linewidth = 0.6
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  scale_colour_manual(
    values = c(
      "blended (old)" = "#878787",
      "FP_between (H2)" = "#d73027",
      "FP_within (H3)" = "#4575b4"
    ),
    name = NULL
  ) +
  labs(
    x = "Phase",
    y = "Fishing-pressure slope\n(∂ residual / ∂ FP component)",
    title = "Within-between decomposition vs blended FP slopes",
    subtitle = paste0(
      "H2 = between-rectangle mean log(hours+1); H3 = within-rectangle deviation. ",
      "Blended = previous log_hours_total * phase model. Error bars = 95% CI."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.2, colour = "grey30"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
ggsave(path_out_fig_slopes, p_slopes, width = 9.5, height = 5.5, dpi = 150)
logmsg("Saved figure: ", path_out_fig_slopes)

# ---------------------------------------------------------------------------
# Coefficient-labelling reminder (no broader interpretation)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Coefficient → hypothesis mapping (no further interpretation)")
logmsg(
  "H2 (spatial / persistent pressure): FP_between main effect, FP_between × phase ",
  "interactions, and the four phase-specific FP_between slopes."
)
logmsg(
  "H3 (temporal / own-rectangle deviation): FP_within main effect, FP_within × phase ",
  "interactions, and the four phase-specific FP_within slopes. The joint test of ",
  "FP_within × phase interactions is the direct parallel to the previous blended ",
  "interaction test for change across phases."
)

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_fe)
logmsg("- ", path_out_slopes)
logmsg("- ", path_out_wald)
logmsg("- ", path_out_compare)
logmsg("- ", path_out_pooling)
logmsg("- ", path_out_fit)
logmsg("- ", path_out_models)
if (car_fitted) logmsg("- ", path_out_car_fe)
logmsg("- ", path_out_fig_pooling)
logmsg("- ", path_out_fig_slopes)
logmsg("- ", path_out_session)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 within-between decomposition complete. ===\n")
