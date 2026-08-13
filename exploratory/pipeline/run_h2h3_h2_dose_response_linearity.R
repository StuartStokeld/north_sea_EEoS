# ARCHIVED (exploratory): H2 dose-response GAM linearity check. Not part of
# current primary or sensitivity reporting. Production model unchanged.
# Re-run only intentionally:
#   Rscript --vanilla exploratory/pipeline/run_h2h3_h2_dose_response_linearity.R
#
# H2 dose-response linearity diagnostic (additive only — does not touch primary model)
#
# PURPOSE: test whether the linear FP_between × phase_v2 specification in the
# current primary within-between model is adequate, by fitting a phase-wise
# smooth of FP_between in mgcv and comparing fit / edf / marginal curves.
#
# Same data, response, phase_v2 cut points, and random-intercept structure as
# primary_model_v2. H3 (FP_within × phase_v2) left parametric and unchanged.
#
# Decision gate (report only — do NOT respecify production here):
#   - edf ≈ 1 and no AIC / nested-test improvement → linearity holds (reporting
#     addendum only)
#   - nonlinearity in one or more phases → stop; flag for design review before
#     any production-model change
#
# Run: Rscript --vanilla exploratory/pipeline/run_h2h3_h2_dose_response_linearity.R
#
# Prerequisite: outputs/primary_model_v2.rds (live). Writes under exploratory/outputs/.

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
  stop("Run: Rscript --vanilla exploratory/pipeline/run_h2h3_h2_dose_response_linearity.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' required (ships with base R).")
}
suppressPackageStartupMessages(library(mgcv))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
# Basis dimension per phase smooth; penalisation can shrink toward linear.
# k = 5 is enough to detect saturation / threshold / non-monotonicity without
# over-fitting the sparse FP_between tails (158 unique rectangle means).
K_SMOOTH <- 5L
# Soft edf gate used in the written verdict (brief: edf ≈ 1 → linear OK).
EDF_LINEAR_CEILING <- 1.5
ALPHA <- 0.05

# ---------------------------------------------------------------------------
# Paths (archived under exploratory/; live primary model stays at top-level)
# ---------------------------------------------------------------------------
out_root <- file.path(project_root, "exploratory", "outputs")
path_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_table <- file.path(out_root, "h2_dose_response_linearity_by_phase.csv")
path_out_ic <- file.path(out_root, "h2_dose_response_linearity_model_comparison.csv")
path_out_curves <- file.path(out_root, "h2_dose_response_linearity_curves.csv")
path_out_models <- file.path(out_root, "h2_dose_response_linearity_model_objects.rds")
path_out_summary <- file.path(out_root, "h2_dose_response_linearity_summary.md")
path_out_run_log <- file.path(out_root, "h2_dose_response_linearity_run_log.md")
path_out_session <- file.path(out_root, "h2_dose_response_linearity_sessionInfo.txt")
path_out_fig_panel <- file.path(fig_dir, "h2_dose_response_linearity_by_phase.png")

stopifnot(file.exists(path_v2))

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2 dose-response linearity diagnostic — run log")
logmsg("")
logmsg(
  "Additive diagnostic only. Does not modify `primary_model_v2` or any ",
  "downstream Moran / permutation / KNN pipeline. H3 left parametric."
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
logmsg(sprintf("mgcv %s", as.character(utils::packageVersion("mgcv"))))

# ---------------------------------------------------------------------------
# Data (same as primary_model_v2)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

v2 <- readRDS(path_v2)
dat <- v2$data
if (is.null(dat)) stop("primary_model_v2.rds missing data.")
stopifnot(
  all(c("residual", "FP_between", "FP_within", "phase_v2", "stat_rec") %in% names(dat))
)

dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
if (anyNA(dat$phase_v2)) stop("NA in phase_v2 after relevel.")
dat$stat_rec <- factor(as.character(dat$stat_rec))

n_haul <- nrow(dat)
n_rect <- nlevels(dat$stat_rec)
fp_rect <- dat %>%
  distinct(stat_rec, FP_between) %>%
  arrange(FP_between)

logmsg(sprintf(
  "Loaded %s — %d hauls, %d rectangles, years %d–%d.",
  path_v2, n_haul, n_rect, min(dat$year), max(dat$year)
))
logmsg(sprintf("phase_v2 levels: %s", paste(levels(dat$phase_v2), collapse = " | ")))
logmsg(sprintf(
  "FP_between (rectangle means): min=%.3f, p25=%.3f, median=%.3f, p75=%.3f, max=%.3f (n_unique=%d).",
  min(fp_rect$FP_between),
  quantile(fp_rect$FP_between, 0.25),
  median(fp_rect$FP_between),
  quantile(fp_rect$FP_between, 0.75),
  max(fp_rect$FP_between),
  nrow(fp_rect)
))
logmsg(sprintf(
  "Production formula (unchanged): %s",
  paste(deparse(v2$formula_v2), collapse = " ")
))
logmsg(sprintf(
  "Diagnostic linear (mgcv): residual ~ FP_between * phase_v2 + FP_within * phase_v2 + s(stat_rec, bs=\"re\")"
))
logmsg(sprintf(
  "Diagnostic smooth (mgcv): residual ~ phase_v2 + s(FP_between, by=phase_v2, k=%d) + FP_within * phase_v2 + s(stat_rec, bs=\"re\")",
  K_SMOOTH
))

# ---------------------------------------------------------------------------
# Fit linear vs smooth (ML for AIC / nested test; REML for curves)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fitting")
logmsg(
  "Both models fit in mgcv so AIC and the nested comparison share the same ",
  "likelihood / EDF penalty construction. Production glmmTMB AIC is recorded ",
  "for reference only (not used for the linearity gate)."
)

fml_linear <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 +
  s(stat_rec, bs = "re")
fml_smooth <- residual ~ phase_v2 + s(FP_between, by = phase_v2, k = K_SMOOTH) +
  FP_within * phase_v2 + s(stat_rec, bs = "re")

logmsg("Fitting linear GAM (ML) ...")
time_lin_ml <- system.time({
  fit_lin_ml <- gam(fml_linear, data = dat, method = "ML")
})
logmsg(sprintf("  elapsed: %.1f sec", time_lin_ml[["elapsed"]]))

logmsg("Fitting smooth GAM (ML) ...")
time_sm_ml <- system.time({
  fit_sm_ml <- gam(fml_smooth, data = dat, method = "ML")
})
logmsg(sprintf("  elapsed: %.1f sec", time_sm_ml[["elapsed"]]))

logmsg("Fitting REML copies for marginal-effect plots ...")
time_reml <- system.time({
  fit_lin_reml <- gam(fml_linear, data = dat, method = "REML")
  fit_sm_reml <- gam(fml_smooth, data = dat, method = "REML")
})
logmsg(sprintf("  elapsed: %.1f sec", time_reml[["elapsed"]]))

# ---------------------------------------------------------------------------
# AIC / nested nonlinearity test
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Fit comparison (ML)")

aic_lin <- as.numeric(AIC(fit_lin_ml))
aic_sm <- as.numeric(AIC(fit_sm_ml))
bic_lin <- as.numeric(BIC(fit_lin_ml))
bic_sm <- as.numeric(BIC(fit_sm_ml))
delta_aic <- aic_sm - aic_lin

logmsg(sprintf("Linear GAM ML AIC = %.2f; BIC = %.2f", aic_lin, bic_lin))
logmsg(sprintf("Smooth GAM ML AIC = %.2f; BIC = %.2f", aic_sm, bic_sm))
logmsg(sprintf(
  "Delta AIC (smooth − linear) = %+.2f  (negative favours smooth)",
  delta_aic
))

# Nested comparison: linear is nested in the smooth family when the smooth
# can collapse to a straight line (edf → 1). anova.gam uses a Chi-sq approx.
anova_cmp <- anova(fit_lin_ml, fit_sm_ml, test = "Chisq")
# anova.gam returns a data.frame-like object; second row is the smooth model.
anova_df <- as.data.frame(anova_cmp)
cn <- names(anova_df)
# Prefer the model-comparison Deviance column, not Resid. Dev.
dev_col <- if ("Deviance" %in% cn) {
  "Deviance"
} else {
  setdiff(grep("Deviance|Chi", cn, value = TRUE),
          grep("Resid|Res\\.", cn, value = TRUE))[1]
}
df_col <- if ("Df" %in% cn) "Df" else grep("^Df$", cn, value = TRUE)[1]
pval_col <- grep("^Pr\\(", cn, value = TRUE)[1]
if (length(pval_col) == 0L) pval_col <- grep("Pr\\(>|p.?value", cn, value = TRUE)[1]

nested_df <- if (length(df_col) && !is.na(df_col)) as.numeric(anova_df[[df_col]][2]) else NA_real_
nested_dev <- if (length(dev_col) && !is.na(dev_col)) as.numeric(anova_df[[dev_col]][2]) else NA_real_
nested_p <- if (length(pval_col) && !is.na(pval_col)) as.numeric(anova_df[[pval_col]][2]) else NA_real_

logmsg("Nested anova(linear, smooth, test=\"Chisq\"):")
logmsg(paste(capture.output(print(anova_cmp)), collapse = "\n"))
logmsg(sprintf(
  "Extracted: Df=%.3f, deviance/Chi=%.3f, p=%.4g",
  nested_df, nested_dev, nested_p
))

# Production glmmTMB AIC for reference (needs glmmTMB methods registered)
aic_prod <- tryCatch({
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    as.numeric(stats::AIC(v2$primary_model_v2))
  } else {
    NA_real_
  }
}, error = function(e) NA_real_)
logmsg(sprintf(
  "Production glmmTMB REML AIC (reference only) = %s",
  if (is.finite(aic_prod)) sprintf("%.2f", aic_prod) else "NA"
))

ic_table <- tibble::tibble(
  model_id = c("linear_gam_ml", "smooth_gam_ml", "production_glmmTMB_reml"),
  framework = c("mgcv::gam", "mgcv::gam", "glmmTMB"),
  estimation = c("ML", "ML", "REML"),
  aic = c(aic_lin, aic_sm, aic_prod),
  bic = c(bic_lin, bic_sm, NA_real_),
  role = c(
    "linearity null (same RE as primary)",
    "flexible H2 dose-response by phase",
    "production primary — AIC not used for gate"
  )
)
write_csv(ic_table, path_out_ic)
logmsg("Wrote: ", path_out_ic)

# ---------------------------------------------------------------------------
# Per-phase smooth edf
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Per-phase smooth diagnostics")

sm_table <- as.data.frame(summary(fit_sm_ml)$s.table)
sm_table$term <- rownames(sm_table)
rownames(sm_table) <- NULL

# Keep FP_between by-phase smooths only (drop random intercept row).
sm_fp <- sm_table %>%
  filter(grepl("FP_between", term, fixed = TRUE)) %>%
  mutate(
    phase_v2 = NA_character_
  )

# Map term labels like s(FP_between):phase_v21985-1991 to phase levels
for (ph in PHASE_V2) {
  hit <- grepl(ph, sm_fp$term, fixed = TRUE)
  sm_fp$phase_v2[hit] <- ph
}
if (anyNA(sm_fp$phase_v2) || nrow(sm_fp) != length(PHASE_V2)) {
  logmsg("Smooth-term table (raw):")
  logmsg(paste(capture.output(print(sm_table)), collapse = "\n"))
  stop("Could not map one smooth term per phase_v2 level.")
}

sm_fp <- sm_fp %>%
  mutate(phase_v2 = factor(phase_v2, levels = PHASE_V2)) %>%
  arrange(phase_v2)

phase_table <- sm_fp %>%
  transmute(
    phase_v2 = as.character(phase_v2),
    edf = edf,
    ref_df = Ref.df,
    F_approx = F,
    p_smooth_vs_zero = `p-value`,
    aic_linear = aic_lin,
    aic_smooth = aic_sm,
    delta_aic_smooth_minus_linear = delta_aic,
    nested_test = "anova(linear_gam_ml, smooth_gam_ml, test=\"Chisq\")",
    nested_df = nested_df,
    nested_statistic = nested_dev,
    nested_p = nested_p,
    edf_linear_ceiling = EDF_LINEAR_CEILING,
    phase_linearity_held = edf <= EDF_LINEAR_CEILING,
    term = term
  )

write_csv(phase_table, path_out_table)
logmsg("Wrote: ", path_out_table)
logmsg("Per-phase edf:")
for (i in seq_len(nrow(phase_table))) {
  r <- phase_table[i, ]
  logmsg(sprintf(
    "  %s: edf=%.3f (Ref.df=%.3f); F=%.3f; approx p(vs zero)=%.4g; phase_linearity_held=%s",
    r$phase_v2, r$edf, r$ref_df, r$F_approx, r$p_smooth_vs_zero, r$phase_linearity_held
  ))
}
logmsg(
  "Note: summary.gam p-values test the smooth against a zero function, not ",
  "against a linear null. The linearity gate uses edf + nested AIC/anova."
)

# ---------------------------------------------------------------------------
# Marginal-effect curves (REML fits; exclude random intercept)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Marginal-effect curves")

fp_grid <- seq(
  min(fp_rect$FP_between, na.rm = TRUE),
  max(fp_rect$FP_between, na.rm = TRUE),
  length.out = 120L
)
# Dummy rectangle level required by predict(); excluded via exclude=
dummy_rec <- levels(dat$stat_rec)[1]

build_newdata <- function(phase) {
  data.frame(
    FP_between = fp_grid,
    FP_within = 0,
    phase_v2 = factor(phase, levels = PHASE_V2),
    stat_rec = factor(dummy_rec, levels = levels(dat$stat_rec))
  )
}

# Contribution of FP_between terms only (phase intercepts removed) so the
# overlay compares dose-response shape, not phase main-effect level shifts.
fp_term_contribution <- function(fit, newdata, model_kind) {
  tr <- predict(fit, newdata = newdata, type = "terms", se.fit = TRUE)
  fit_mat <- as.matrix(tr$fit)
  se_mat <- as.matrix(tr$se.fit)
  cn <- colnames(fit_mat)
  if (model_kind == "smooth") {
    cols <- grepl("FP_between", cn, fixed = TRUE)
  } else {
    # parametric: FP_between main + FP_between:phase interactions
    cols <- grepl("FP_between", cn, fixed = TRUE) & !grepl("^s\\(", cn)
  }
  if (!any(cols)) {
    stop("No FP_between term columns found in predict(type=\"terms\") for ", model_kind)
  }
  list(
    fit = rowSums(fit_mat[, cols, drop = FALSE]),
    se = sqrt(rowSums(se_mat[, cols, drop = FALSE]^2))
  )
}

# Center both curves at median FP_between so overlays compare shape/slope,
# not arbitrary mgcv term-centering offsets between parametric and smooth.
fp_med <- stats::median(fp_rect$FP_between, na.rm = TRUE)
i_med <- which.min(abs(fp_grid - fp_med))

curve_rows <- list()
for (ph in PHASE_V2) {
  nd <- build_newdata(ph)
  sm <- fp_term_contribution(fit_sm_reml, nd, "smooth")
  ln <- fp_term_contribution(fit_lin_reml, nd, "linear")
  sm_c <- sm$fit - sm$fit[i_med]
  ln_c <- ln$fit - ln$fit[i_med]
  curve_rows[[ph]] <- tibble::tibble(
    phase_v2 = ph,
    FP_between = fp_grid,
    smooth_fit = sm_c,
    smooth_se = sm$se,
    smooth_lo = sm_c - 1.96 * sm$se,
    smooth_hi = sm_c + 1.96 * sm$se,
    linear_fit = ln_c,
    linear_se = ln$se,
    linear_lo = ln_c - 1.96 * ln$se,
    linear_hi = ln_c + 1.96 * ln$se,
    centered_at_FP_between = fp_med
  )
}
curves <- bind_rows(curve_rows)
write_csv(curves, path_out_curves)
logmsg("Wrote: ", path_out_curves)

# Four-panel figure
rug_df <- fp_rect %>% transmute(FP_between, phase_v2 = "rug")
edf_lab <- phase_table %>%
  transmute(
    phase_v2,
    label = sprintf("edf = %.2f", edf)
  )

p <- ggplot(curves, aes(x = FP_between)) +
  geom_ribbon(
    aes(ymin = smooth_lo, ymax = smooth_hi),
    fill = "#4393c3",
    alpha = 0.25
  ) +
  geom_line(aes(y = smooth_fit, colour = "Smooth (GAM)"), linewidth = 0.9) +
  geom_line(
    aes(y = linear_fit, colour = "Linear (GAM)"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_rug(
    data = fp_rect,
    aes(x = FP_between),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.35,
    length = unit(0.04, "npc")
  ) +
  geom_text(
    data = edf_lab,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5,
    size = 3.2,
    colour = "grey20"
  ) +
  facet_wrap(~phase_v2, ncol = 2, scales = "free_y") +
  scale_colour_manual(
    name = NULL,
    values = c("Smooth (GAM)" = "#2166ac", "Linear (GAM)" = "#b2182b")
  ) +
  labs(
    title = "H2 dose-response linearity diagnostic: FP_between by phase_v2",
    subtitle = paste0(
      "Partial effect of FP_between on EEoS residual (FP_within = 0; ",
      "RE excluded; curves centered at median FP_between). Rug = rectangles."
    ),
    x = expression(FP[between]~"(log fishing hours, rectangle long-run mean)"),
    y = "Partial effect on residual (centered)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", colour = NA),
    plot.subtitle = element_text(size = 9)
  )

ggsave(path_out_fig_panel, p, width = 10, height = 8, dpi = 150)
logmsg("Wrote: ", path_out_fig_panel)

# Also write one file per phase (brief asked for four plots)
for (ph in PHASE_V2) {
  sub <- curves %>% filter(phase_v2 == ph)
  edf_i <- phase_table$edf[phase_table$phase_v2 == ph]
  p1 <- ggplot(sub, aes(x = FP_between)) +
    geom_ribbon(
      aes(ymin = smooth_lo, ymax = smooth_hi),
      fill = "#4393c3",
      alpha = 0.25
    ) +
    geom_line(aes(y = smooth_fit, colour = "Smooth (GAM)"), linewidth = 1) +
    geom_line(
      aes(y = linear_fit, colour = "Linear (GAM)"),
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    geom_rug(
      data = fp_rect,
      aes(x = FP_between),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.4,
      length = unit(0.05, "npc")
    ) +
    scale_colour_manual(
      name = NULL,
      values = c("Smooth (GAM)" = "#2166ac", "Linear (GAM)" = "#b2182b")
    ) +
    labs(
      title = sprintf("H2 dose-response: %s (edf = %.2f)", ph, edf_i),
      subtitle = "Centered at median FP_between; 95% CI on smooth; rug = rectangles",
      x = expression(FP[between]),
      y = "Partial effect on residual (centered)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  fname <- file.path(
    fig_dir,
    sprintf("h2_dose_response_linearity_%s.png", gsub("-", "_", ph))
  )
  ggsave(fname, p1, width = 7, height = 5, dpi = 150)
  logmsg("Wrote: ", fname)
}

# ---------------------------------------------------------------------------
# Verdict / summary paragraph
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Decision gate")

any_edf_high <- any(phase_table$edf > EDF_LINEAR_CEILING)
nested_sig <- is.finite(nested_p) && nested_p < ALPHA
aic_favours_smooth <- delta_aic < -2  # conventional "positive" AIC evidence
nonlinearity_detected <- any_edf_high || nested_sig || aic_favours_smooth

phase_sentences <- vapply(seq_len(nrow(phase_table)), function(i) {
  r <- phase_table[i, ]
  if (r$edf <= EDF_LINEAR_CEILING) {
    sprintf(
      "%s: linear assumption held (edf = %.2f ≈ 1)",
      r$phase_v2, r$edf
    )
  } else {
    # Brief shape note from centered curve (for design-review handoff)
    sub <- curves %>% filter(phase_v2 == r$phase_v2)
    ds <- diff(sub$smooth_fit)
    n_sign_flip <- sum(diff(sign(ds + 1e-12)) != 0)
    # Mid vs upper tercile rise (data-dense region ~ p25–p90 of rectangles)
    q <- quantile(fp_rect$FP_between, c(0.25, 0.5, 0.9))
    near <- function(t) which.min(abs(sub$FP_between - t))
    y25 <- sub$smooth_fit[near(q[1])]
    y50 <- sub$smooth_fit[near(q[2])]
    y90 <- sub$smooth_fit[near(q[3])]
    shape <- if (n_sign_flip == 0L && all(ds >= -1e-8)) {
      "monotonic increasing"
    } else if (n_sign_flip == 0L && all(ds <= 1e-8)) {
      "monotonic decreasing"
    } else if (n_sign_flip >= 2L && (y50 < y25 || y50 < y90 - 0.05)) {
      "non-monotonic (mid-range dip / uneven rise)"
    } else if (abs(y90 - y50) > 2 * abs(y50 - y25) + 0.05) {
      "threshold-like (shallow then steeper rise at higher FP)"
    } else {
      "curved / non-linear"
    }
    sprintf(
      "%s: linearity NOT held (edf = %.2f > %.1f; shape ≈ %s)",
      r$phase_v2, r$edf, EDF_LINEAR_CEILING, shape
    )
  }
}, character(1))

if (!nonlinearity_detected) {
  decision <- "LINEARITY_HOLDS"
  next_step <- paste0(
    "No further action on the production model. Treat this as a reporting ",
    "addendum: keep the linear FP_between × phase_v2 specification; add ",
    "per-SD and wider percentile (e.g. 10th–90th) contrasts alongside the ",
    "existing IQR figure, and include the marginal-effect panel as supporting ",
    "evidence of linearity."
  )
} else {
  decision <- "NONLINEARITY_DETECTED"
  next_step <- paste0(
    "STOP — do not respecify the production model from this script. Review ",
    "the phase-wise curves (threshold vs saturation vs non-monotonic) before ",
    "any redesign of the fixed-effect form, phase interaction, or re-running ",
    "the Moran / permutation / KNN residual pipeline."
  )
}

summary_para <- paste0(
  "H2 dose-response linearity diagnostic (mgcv, same haul data and ",
  "(1|stat_rec) structure as primary_model_v2; H3 left linear): ML AIC ",
  sprintf("linear = %.1f, smooth = %.1f (Δ = %+.1f); ", aic_lin, aic_sm, delta_aic),
  sprintf(
    "nested Chi-sq test p = %s. ",
    if (is.finite(nested_p)) sprintf("%.4g", nested_p) else "NA"
  ),
  paste(phase_sentences, collapse = "; "),
  ". ",
  if (decision == "LINEARITY_HOLDS") {
    "Overall verdict: the linear FP_between specification is adequate across all phases."
  } else {
    "Overall verdict: nonlinearity is detected in at least one phase or by the nested fit comparison — production model unchanged pending design review."
  }
)

logmsg(sprintf("Decision flag: %s", decision))
logmsg(sprintf("any_edf_high=%s; nested_sig=%s; aic_favours_smooth=%s",
               any_edf_high, nested_sig, aic_favours_smooth))
logmsg("")
logmsg("### Summary paragraph")
logmsg(summary_para)
logmsg("")
logmsg("### Recommended next step (do not implement here)")
logmsg(next_step)

summary_md <- c(
  "# H2 dose-response linearity diagnostic — summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "## Decision flag",
  "",
  paste0("**`", decision, "`**"),
  "",
  "## One-paragraph result",
  "",
  summary_para,
  "",
  "## Results table (per phase)",
  "",
  paste0("See `exploratory/outputs/h2_dose_response_linearity_by_phase.csv`. ",
         "Gate: edf ≤ ", EDF_LINEAR_CEILING, " treated as linear-adequate; ",
         "nested test and ΔAIC are model-level."),
  "",
  "| phase_v2 | edf | Ref.df | F | p (vs zero) | phase_linearity_held |",
  "|----------|-----|--------|---|-------------|----------------------|",
  vapply(seq_len(nrow(phase_table)), function(i) {
    r <- phase_table[i, ]
    sprintf(
      "| %s | %.3f | %.3f | %.3f | %.4g | %s |",
      r$phase_v2, r$edf, r$ref_df, r$F_approx, r$p_smooth_vs_zero, r$phase_linearity_held
    )
  }, character(1)),
  "",
  "## Model comparison (ML)",
  "",
  sprintf("- Linear GAM AIC = %.2f", aic_lin),
  sprintf("- Smooth GAM AIC = %.2f", aic_sm),
  sprintf("- ΔAIC (smooth − linear) = %+.2f", delta_aic),
  sprintf(
    "- Nested anova Chi-sq: Df = %.3f, statistic = %.3f, p = %.4g",
    nested_df, nested_dev, nested_p
  ),
  "",
  "## Figures",
  "",
  "- Panel: `exploratory/outputs/figures/h2_dose_response_linearity_by_phase.png`",
  "- Per phase: `exploratory/outputs/figures/h2_dose_response_linearity_1985_1991.png` (and siblings)",
  "",
  "## What happens next (not implemented)",
  "",
  next_step,
  "",
  "## Explicitly unchanged",
  "",
  "- Production model `outputs/primary_model_v2.rds`",
  "- H3 (`FP_within × phase_v2`)",
  "- Random-intercept structure, phase cut points",
  "- Moran / permutation / KNN residual pipeline",
  ""
)
writeLines(summary_md, path_out_summary)
logmsg("Wrote: ", path_out_summary)

saveRDS(
  list(
    fit_lin_ml = fit_lin_ml,
    fit_sm_ml = fit_sm_ml,
    fit_lin_reml = fit_lin_reml,
    fit_sm_reml = fit_sm_reml,
    phase_table = phase_table,
    ic_table = ic_table,
    curves = curves,
    anova_comparison = anova_cmp,
    decision = decision,
    summary_paragraph = summary_para,
    k_smooth = K_SMOOTH,
    edf_linear_ceiling = EDF_LINEAR_CEILING,
    formula_linear = fml_linear,
    formula_smooth = fml_smooth,
    data_n = n_haul,
    n_rect = n_rect
  ),
  path_out_models
)
logmsg("Wrote: ", path_out_models)

writeLines(run_log, path_out_run_log)
logmsg("Wrote: ", path_out_run_log)
cat("=== H2 dose-response linearity diagnostic complete. Decision: ",
    decision, " ===\n", sep = "")
