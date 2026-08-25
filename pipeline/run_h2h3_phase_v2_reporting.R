# phase_v2 reporting — slopes, CAR, proportional effects, presentation figures
#
# PURPOSE: after run_h2h3_phase_v2_refit.R, produce the full reporting stack
# that the H2/H3 drafts need for policy-anchored phases:
#   - phase-specific FP_between / FP_within slopes (RE companion + CAR)
#   - CAR fit with phase_v2 (H2 and H3 reported from the same CAR object)
#   - pooled FP_between CAR contrast (H2 figure dashed overlay)
#   - proportional / gap-change tables (H2 and H3 from CAR)
#   - presentation figures (H2 + H3) and combined gap-change figure
#
# Supplementary CAR diagnostics (LOO influence, CAR R², signed vs unsigned
# residual audit) do NOT live here — they read the saved fit_car and must not
# change this signed primary fit:
#   Rscript --vanilla pipeline/run_h2h3_car_reporting_diagnostics.R
#
# Does NOT overwrite original-phase artifacts (h2h3_wb_*). All outputs use
# the phase_v2_ prefix.
#
# Run: Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R
#
# Prerequisite: outputs/primary_model_v2.rds from run_h2h3_phase_v2_refit.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_phase_v2_reporting.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

# Prefer project-local library (e.g. .R_libs/spaMM) when present
local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))
source(file.path(script_dir, "R", "h2h3_results_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("glmmTMB required. Run: Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R")
}
suppressPackageStartupMessages(library(glmmTMB))

has_spamm <- requireNamespace("spaMM", quietly = TRUE)
if (has_spamm) suppressPackageStartupMessages(library(spaMM))

PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
PHASE_V2_BREAKS <- c(1992, 2002, 2008)

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_round2 <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
stopifnot(file.exists(path_v2))

path_out_slopes <- file.path(project_root, "outputs", "phase_v2_fp_slopes_by_phase.csv")
path_out_fe <- file.path(project_root, "outputs", "phase_v2_primary_fixed_effects.csv")
path_out_car_fe <- file.path(project_root, "outputs", "phase_v2_car_fixed_effects.csv")
path_out_h2 <- file.path(project_root, "outputs", "phase_v2_proportional_effects_H2.csv")
path_out_h3 <- file.path(project_root, "outputs", "phase_v2_proportional_effects_H3.csv")
path_out_pooled <- file.path(project_root, "outputs", "phase_v2_pooled_between_coef.csv")
path_out_models <- file.path(project_root, "outputs", "phase_v2_reporting_model_objects.rds")
path_out_fig_combined <- file.path(fig_dir, "phase_v2_wb_gap_change_by_phase.png")
path_out_fig_h2 <- file.path(fig_dir, "phase_v2_presentation_H2_gap_change_by_phase.png")
path_out_fig_h3 <- file.path(fig_dir, "phase_v2_presentation_H3_gap_change_by_phase.png")
path_out_run_log <- file.path(project_root, "outputs", "phase_v2_reporting_run_log.md")
path_out_session <- file.path(project_root, "outputs", "phase_v2_reporting_sessionInfo.txt")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# phase_v2 reporting — slopes / CAR / proportional effects / figures — run log")
logmsg("")
logmsg(
  "Builds the reporting stack for the policy-anchored primary model ",
  "(`phase_v2`: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015). ",
  "Does not overwrite original-phase `h2h3_wb_*` artifacts."
)

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg(sprintf(
  "glmmTMB %s; spaMM available = %s%s",
  as.character(utils::packageVersion("glmmTMB")),
  has_spamm,
  if (has_spamm) paste0(" (", utils::packageVersion("spaMM"), ")") else ""
))

# ---------------------------------------------------------------------------
# Load phase_v2 primary
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

v2 <- readRDS(path_v2)
fit_wb <- v2$primary_model_v2
dat <- v2$data
stopifnot(!is.null(fit_wb), !is.null(dat))
stopifnot(all(c("phase_v2", "FP_between", "FP_within", "residual", "stat_rec") %in% names(dat)))
if (!identical(levels(dat$phase_v2), PHASE_V2)) {
  dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
}
logmsg(sprintf(
  "Loaded %s: %d hauls, %d rectangles; formula: %s",
  path_v2, nrow(dat), dplyr::n_distinct(dat$stat_rec),
  paste(deparse(formula(fit_wb)), collapse = " ")
))
logmsg(sprintf("phase_v2 levels: %s", paste(levels(dat$phase_v2), collapse = " | ")))

# ---------------------------------------------------------------------------
# Primary FE table + slopes
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Primary (RE) fixed effects and phase-specific slopes")

fe_pri <- tidy_fixed_effects(fit_wb, "wb_primary_v2") %>%
  annotate_wb_fixed_effects()
write_csv(fe_pri, path_out_fe)
logmsg("Saved: ", path_out_fe)

slopes_between_pri <- extract_wb_phase_slopes(
  fit_wb, "FP_between", "wb_primary_v2", "H2_spatial_between", phases = PHASE_V2
)
slopes_within_pri <- extract_wb_phase_slopes(
  fit_wb, "FP_within", "wb_primary_v2", "H3_temporal_within", phases = PHASE_V2
)

# ---------------------------------------------------------------------------
# CAR refit + pooled FP_between contrast
# ---------------------------------------------------------------------------
fit_car <- NULL
fit_pooled <- NULL
car_fitted <- FALSE
pooled_fitted <- FALSE
slopes_car <- NULL

logmsg("")
logmsg("## CAR sensitivity (phase_v2)")

if (!has_spamm) {
  logmsg("spaMM not available — CAR SKIPPED. H2 proportional table requires CAR between slopes.")
} else if (!file.exists(path_round2)) {
  logmsg("Round 2 RDS missing — CAR SKIPPED.")
} else {
  round2 <- readRDS(path_round2)
  if (is.null(round2$adjMatrix)) {
    logmsg("Round 2 adjMatrix missing — CAR SKIPPED.")
  } else {
    adjMatrix <- round2$adjMatrix
    dat_car <- dat
    dat_car$stat_rec <- factor(as.character(dat_car$stat_rec), levels = rownames(adjMatrix))
    formula_car <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
    if (!is.null(v2$fit_wb_car_v2)) {
      fit_car <- v2$fit_wb_car_v2
      car_fitted <- TRUE
      logmsg("Reusing existing fit_wb_car_v2 (no CAR refit).")
    } else {
      logmsg("Fitting CAR: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)")
      time_car <- system.time({
        fit_car <- spaMM::fitme(
          formula_car, data = dat_car, adjMatrix = adjMatrix, method = "REML"
        )
      })
      logmsg(sprintf("CAR fit time: %.2f sec.", time_car[["elapsed"]]))
      car_fitted <- TRUE
    }

    fe_car <- tidy_fixed_effects_spamm(fit_car, "wb_car_v2") %>%
      annotate_wb_fixed_effects()
    write_csv(fe_car, path_out_car_fe)
    logmsg("Saved: ", path_out_car_fe)

    slopes_car <- bind_rows(
      extract_wb_phase_slopes_spamm(
        fit_car, "FP_between", "wb_car_v2", "H2_spatial_between", phases = PHASE_V2
      ),
      extract_wb_phase_slopes_spamm(
        fit_car, "FP_within", "wb_car_v2", "H3_temporal_within", phases = PHASE_V2
      )
    )

    # Pooled FP_between contrast (for H2 figure overlay)
    logmsg("")
    logmsg("### Pooled FP_between CAR contrast (H2 figure overlay)")
    formula_pooled <- residual ~ FP_between + FP_within * phase_v2 + adjacency(1 | stat_rec)
    if (!is.null(v2$fit_pooled_between_car_v2)) {
      fit_pooled <- v2$fit_pooled_between_car_v2
      pooled_fitted <- TRUE
      logmsg("Reusing existing pooled CAR fit (no pooled refit).")
    } else {
      logmsg("Fitting: residual ~ FP_between + FP_within * phase_v2 + adjacency(1 | stat_rec)")
      time_pool <- system.time({
        fit_pooled <- spaMM::fitme(
          formula_pooled, data = dat_car, adjMatrix = adjMatrix, method = "REML"
        )
      })
      logmsg(sprintf("Pooled CAR fit time: %.2f sec.", time_pool[["elapsed"]]))
      pooled_fitted <- TRUE
    }

    fe_pooled <- tidy_fixed_effects_spamm(fit_pooled, "wb_pooled_between_car_v2")
    pooled_row <- fe_pooled %>% filter(term == "FP_between")
    if (nrow(pooled_row) != 1L) {
      stop("Expected one FP_between in pooled CAR; found: ", paste(fe_pooled$term, collapse = ", "))
    }
    pooled_out <- pooled_row %>%
      mutate(
        ci_lo = estimate - 1.96 * std_error,
        ci_hi = estimate + 1.96 * std_error,
        term_plain = "Pooled FP_between (no phase_v2 interaction)",
        hypothesis = "H2 contrast — time-stable between-rectangle FP slope (phase_v2)",
        model_note = "residual ~ FP_between + FP_within * phase_v2 + adjacency(1|stat_rec)"
      ) %>%
      select(
        model_id, term, term_plain, hypothesis,
        estimate, std_error, ci_lo, ci_hi, statistic, p_value, model_note
      )
    write_csv(pooled_out, path_out_pooled)
    logmsg(sprintf(
      "Pooled FP_between = %+.6f (SE %.6f; 95%% CI [%+.6f, %+.6f]; p = %.4g)",
      pooled_out$estimate, pooled_out$std_error,
      pooled_out$ci_lo, pooled_out$ci_hi, pooled_out$p_value
    ))
    logmsg("Saved: ", path_out_pooled)
  }
}

slopes_all <- bind_rows(slopes_between_pri, slopes_within_pri, slopes_car) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(model_id, component, phase)
write_csv(slopes_all, path_out_slopes)
logmsg("")
logmsg("## Phase-specific slopes")
logmsg("Saved: ", path_out_slopes)
for (mid in unique(slopes_all$model_id)) {
  logmsg(sprintf("### %s", mid))
  sub <- slopes_all %>% filter(model_id == mid)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    logmsg(sprintf(
      "  - [%s] %s / %s: slope=%+.4f SE=%.4f CI=[%+.4f, %+.4f] p=%.3g",
      r$hypothesis_group, r$component, as.character(r$phase),
      r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi, r$p_value
    ))
  }
}

h2_slope_model <- if (car_fitted) "wb_car_v2" else "wb_primary_v2"
h3_slope_model <- if (car_fitted) "wb_car_v2" else "wb_primary_v2"
if (!car_fitted) {
  logmsg("")
  logmsg(
    "FLAG: CAR unavailable — H2 and H3 proportional tables / figures will use ",
    "primary RE (`wb_primary_v2`) slopes. Install spaMM (+ GSL) and re-run for ",
    "CAR-based H2/H3 matching the presented shared-model convention."
  )
}

# ---------------------------------------------------------------------------
# Proportional / gap-change transforms
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Proportional effect sizes")

h2_slopes <- slopes_all %>%
  filter(model_id == h2_slope_model, component == "FP_between")
h3_slopes <- slopes_all %>%
  filter(model_id == h3_slope_model, component == "FP_within")
h3_re <- slopes_all %>%
  filter(model_id == "wb_primary_v2", component == "FP_within")

if (nrow(h2_slopes) != 4L || nrow(h3_slopes) != 4L) {
  stop("Expected 4 H2 between and 4 H3 within phase_v2 slopes.")
}

logmsg(sprintf("H2 source: `%s` FP_between phase slopes.", h2_slope_model))
logmsg(sprintf("H3 source: `%s` FP_within phase slopes.", h3_slope_model))
if (identical(h3_slope_model, "wb_car_v2") && nrow(h3_re) == 4L) {
  cmp_h3 <- h3_slopes %>%
    select(phase, car = fp_slope) %>%
    left_join(h3_re %>% select(phase, re = fp_slope), by = "phase") %>%
    mutate(abs_diff = abs(car - re))
  logmsg(sprintf(
    "Max |CAR − RE| FP_within phase slope = %.5f (RE retained as companion).",
    max(cmp_h3$abs_diff, na.rm = TRUE)
  ))
}

phase_base <- dat %>%
  group_by(phase_v2) %>%
  summarise(
    n_hauls = dplyr::n(),
    median_residual = median(residual, na.rm = TRUE),
    baseline_ratio = exp(median_residual),
    iqr_FP_within = IQR(FP_within, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(phase = as.character(phase_v2))

rect_between <- dat %>% distinct(stat_rec, FP_between)
iqr_FP_between <- IQR(rect_between$FP_between, na.rm = TRUE)

logmsg("Baselines (phase_v2 median residual → ratio):")
for (i in seq_len(nrow(phase_base))) {
  r <- phase_base[i, ]
  logmsg(sprintf(
    "  - %s: median residual=%.4f; ratio=%.4f (%.1f%%); IQR(FP_within)=%.3f",
    r$phase, r$median_residual, r$baseline_ratio, 100 * r$baseline_ratio, r$iqr_FP_within
  ))
}
logmsg(sprintf(
  "H2 IQR(FP_between) across %d rectangles = %.4f.",
  nrow(rect_between), iqr_FP_between
))

prop_metrics <- function(slope, slope_lo, slope_hi, dx, r0) {
  pct_ratio <- function(b) (exp(b * dx) - 1) * 100
  r1 <- function(b) r0 * exp(b * dx)
  pp_shift <- function(b) (r1(b) - r0) * 100
  gap0 <- 1 - r0
  pct_gap_change <- function(b) {
    g1 <- 1 - r1(b)
    ((gap0 - g1) / gap0) * 100
  }
  gap_dir <- function(b) {
    g <- pct_gap_change(b)
    if (!is.finite(g)) return(NA_character_)
    if (g > 0) "closed" else if (g < 0) "widened" else "unchanged"
  }
  list(
    pct_ratio_change = pct_ratio(slope),
    pct_ratio_change_lo = pct_ratio(slope_lo),
    pct_ratio_change_hi = pct_ratio(slope_hi),
    ratio_start = r0,
    ratio_end = r1(slope),
    pp_ratio_shift = pp_shift(slope),
    pp_ratio_shift_lo = pp_shift(slope_lo),
    pp_ratio_shift_hi = pp_shift(slope_hi),
    pct_gap_change = pct_gap_change(slope),
    pct_gap_change_lo = pct_gap_change(slope_lo),
    pct_gap_change_hi = pct_gap_change(slope_hi),
    gap_direction = gap_dir(slope)
  )
}

h2_rows <- lapply(seq_len(nrow(h2_slopes)), function(i) {
  s <- h2_slopes[i, ]
  ph <- as.character(s$phase)
  base <- phase_base %>% filter(phase == ph)
  r0 <- base$baseline_ratio
  d_dbl <- prop_metrics(s$fp_slope, s$fp_slope_lo, s$fp_slope_hi, log(2), r0)
  d_iqr <- prop_metrics(s$fp_slope, s$fp_slope_lo, s$fp_slope_hi, iqr_FP_between, r0)
  data.frame(
    hypothesis = "H2",
    component = "FP_between",
    model_id = h2_slope_model,
    phase = ph,
    year_start = s$year_start,
    year_end = s$year_end,
    fp_slope = s$fp_slope,
    fp_slope_lo = s$fp_slope_lo,
    fp_slope_hi = s$fp_slope_hi,
    baseline = "phase_median_residual",
    baseline_residual = base$median_residual,
    baseline_ratio = r0,
    ratio_start_pct = r0 * 100,
    delta_x_doubling = log(2),
    pct_ratio_change_doubling = d_dbl$pct_ratio_change,
    pct_ratio_change_doubling_lo = d_dbl$pct_ratio_change_lo,
    pct_ratio_change_doubling_hi = d_dbl$pct_ratio_change_hi,
    ratio_end_doubling_pct = d_dbl$ratio_end * 100,
    pp_ratio_shift_doubling = d_dbl$pp_ratio_shift,
    pp_ratio_shift_doubling_lo = d_dbl$pp_ratio_shift_lo,
    pp_ratio_shift_doubling_hi = d_dbl$pp_ratio_shift_hi,
    pct_gap_change_doubling = d_dbl$pct_gap_change,
    pct_gap_change_doubling_lo = d_dbl$pct_gap_change_lo,
    pct_gap_change_doubling_hi = d_dbl$pct_gap_change_hi,
    gap_direction_doubling = d_dbl$gap_direction,
    delta_x_iqr_FP_between = iqr_FP_between,
    pct_ratio_change_iqr = d_iqr$pct_ratio_change,
    pct_ratio_change_iqr_lo = d_iqr$pct_ratio_change_lo,
    pct_ratio_change_iqr_hi = d_iqr$pct_ratio_change_hi,
    ratio_end_iqr_pct = d_iqr$ratio_end * 100,
    pp_ratio_shift_iqr = d_iqr$pp_ratio_shift,
    pp_ratio_shift_iqr_lo = d_iqr$pp_ratio_shift_lo,
    pp_ratio_shift_iqr_hi = d_iqr$pp_ratio_shift_hi,
    pct_gap_change_iqr = d_iqr$pct_gap_change,
    pct_gap_change_iqr_lo = d_iqr$pct_gap_change_lo,
    pct_gap_change_iqr_hi = d_iqr$pct_gap_change_hi,
    gap_direction_iqr = d_iqr$gap_direction,
    stringsAsFactors = FALSE
  )
})
h2_out <- bind_rows(h2_rows) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase)
write_csv(h2_out, path_out_h2)
logmsg("Saved: ", path_out_h2)
for (i in seq_len(nrow(h2_out))) {
  r <- h2_out[i, ]
  logmsg(sprintf(
    "  H2 %s: slope=%+.4f; IQR gap %s by %.2f%% [%.2f, %.2f]",
    as.character(r$phase), r$fp_slope, r$gap_direction_iqr,
    abs(r$pct_gap_change_iqr), r$pct_gap_change_iqr_lo, r$pct_gap_change_iqr_hi
  ))
}

h3_rows <- lapply(seq_len(nrow(h3_slopes)), function(i) {
  s <- h3_slopes[i, ]
  ph <- as.character(s$phase)
  base <- phase_base %>% filter(phase == ph)
  r0 <- base$baseline_ratio
  dx <- base$iqr_FP_within
  d <- prop_metrics(s$fp_slope, s$fp_slope_lo, s$fp_slope_hi, dx, r0)
  data.frame(
    hypothesis = "H3",
    component = "FP_within",
    model_id = h3_slope_model,
    phase = ph,
    year_start = s$year_start,
    year_end = s$year_end,
    fp_slope = s$fp_slope,
    fp_slope_lo = s$fp_slope_lo,
    fp_slope_hi = s$fp_slope_hi,
    baseline = "phase_median_residual",
    baseline_residual = base$median_residual,
    baseline_ratio = r0,
    ratio_start_pct = r0 * 100,
    doubling_convention = "not_applicable",
    doubling_note = paste0(
      "FP_within is a mean-zero within-rectangle deviation (can be negative); ",
      "doubling of fishing hours is not a meaningful Δx on this scale."
    ),
    delta_x_iqr_FP_within = dx,
    pct_ratio_change_iqr = d$pct_ratio_change,
    pct_ratio_change_iqr_lo = d$pct_ratio_change_lo,
    pct_ratio_change_iqr_hi = d$pct_ratio_change_hi,
    ratio_end_iqr_pct = d$ratio_end * 100,
    pp_ratio_shift_iqr = d$pp_ratio_shift,
    pp_ratio_shift_iqr_lo = d$pp_ratio_shift_lo,
    pp_ratio_shift_iqr_hi = d$pp_ratio_shift_hi,
    pct_gap_change_iqr = d$pct_gap_change,
    pct_gap_change_iqr_lo = d$pct_gap_change_lo,
    pct_gap_change_iqr_hi = d$pct_gap_change_hi,
    gap_direction_iqr = d$gap_direction,
    stringsAsFactors = FALSE
  )
})
h3_out <- bind_rows(h3_rows) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase)
write_csv(h3_out, path_out_h3)
logmsg("Saved: ", path_out_h3)
for (i in seq_len(nrow(h3_out))) {
  r <- h3_out[i, ]
  logmsg(sprintf(
    "  H3 %s: slope=%+.4f; IQR gap %s by %.2f%% [%.2f, %.2f]",
    as.character(r$phase), r$fp_slope, r$gap_direction_iqr,
    abs(r$pct_gap_change_iqr), r$pct_gap_change_iqr_lo, r$pct_gap_change_iqr_hi
  ))
}

# ---------------------------------------------------------------------------
# Combined figure (like wb_proportional_effects)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Figures")

plot_df <- bind_rows(
  h2_out %>%
    transmute(
      phase, year_mid = (year_start + year_end) / 2,
      series = sprintf("H2: FP_between (%s), per doubling", h2_slope_model),
      pct_gap = pct_gap_change_doubling,
      pct_gap_lo = pct_gap_change_doubling_lo,
      pct_gap_hi = pct_gap_change_doubling_hi
    ),
  h3_out %>%
    transmute(
      phase, year_mid = (year_start + year_end) / 2,
      series = sprintf("H3: FP_within (%s), per IQR of within-deviation", h3_slope_model),
      pct_gap = pct_gap_change_iqr,
      pct_gap_lo = pct_gap_change_iqr_lo,
      pct_gap_hi = pct_gap_change_iqr_hi
    )
)
h2_series_lab <- sprintf("H2: FP_between (%s), per doubling", h2_slope_model)
h3_series_lab <- sprintf("H3: FP_within (%s), per IQR of within-deviation", h3_slope_model)
plot_df <- plot_df %>%
  mutate(
    series = factor(
      series,
      levels = c(h2_series_lab, h3_series_lab)
    )
  )

col_map <- setNames(
  c("#d73027", "#4575b4"),
  c(h2_series_lab, h3_series_lab)
)

p_combined <- ggplot(plot_df, aes(x = year_mid, y = pct_gap, colour = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_vline(
    xintercept = PHASE_V2_BREAKS, linetype = "dotted",
    colour = "grey40", linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = pct_gap_lo, ymax = pct_gap_hi),
    position = position_dodge(width = 1.8), width = 1.0, linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 1.8), size = 2.6) +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_x_continuous(breaks = seq(1985, 2015, by = 5), limits = c(1984, 2016)) +
  labs(
    x = "Year (phase midpoint)",
    y = "Percent of remaining overprediction gap\nclosed (+) or widened (−)",
    title = "phase_v2 within-between H2/H3: proportional gap change by phase",
    subtitle = paste0(
      "Policy-anchored phases (1992 / 2002 / 2008). ",
      "H2 = ", h2_slope_model, " between slopes (doubling); ",
      "H3 = ", h3_slope_model, " within slopes (IQR). Dotted lines = policy breakpoints."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.subtitle = element_text(size = 7.8, colour = "grey30"),
    panel.grid.minor = element_blank()
  )
ggsave(path_out_fig_combined, p_combined, width = 10, height = 6, dpi = 150)
logmsg("Saved: ", path_out_fig_combined)

# ---------------------------------------------------------------------------
# Presentation figures (standalone H2 / H3; IQR convention; pooled H2 overlay)
# ---------------------------------------------------------------------------
Y_LABEL <- paste0(
  "Change in overprediction gap (% of remaining gap)\n",
  "(+ above zero = gap closed; \u2212 below zero = gap widened)"
)

pct_gap_change_from_slope <- function(slope, dx, r0) {
  gap0 <- 1 - r0
  gap1 <- 1 - r0 * exp(slope * dx)
  ((gap0 - gap1) / gap0) * 100
}

prepare_phase_plot_data <- function(df) {
  df %>%
    mutate(
      phase = factor(phase, levels = PHASE_V2),
      pct_gap = pct_gap_change_iqr,
      pct_gap_lo = pct_gap_change_iqr_lo,
      pct_gap_hi = pct_gap_change_iqr_hi,
      significant = (pct_gap_lo > 0) | (pct_gap_hi < 0)
    )
}

plot_gap_by_phase <- function(plot_df, title, subtitle = NULL, pooled_overlay = NULL) {
  p <- ggplot(plot_df, aes(x = phase, y = pct_gap)) +
    geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey45") +
    geom_errorbar(
      aes(ymin = pct_gap_lo, ymax = pct_gap_hi, colour = significant),
      width = 0.18, linewidth = 0.75
    ) +
    geom_point(
      aes(fill = significant, colour = significant),
      size = 4, shape = 21, stroke = 0.9
    ) +
    scale_colour_manual(
      values = c("TRUE" = "#2166ac", "FALSE" = "#2166ac"),
      guide = "none"
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("TRUE" = "#2166ac", "FALSE" = "white"),
      labels = c("TRUE" = "95% CI excludes zero", "FALSE" = "95% CI includes zero")
    ) +
    labs(x = "Phase (policy-anchored)", y = Y_LABEL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 10, angle = 15, hjust = 1),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5, colour = "grey30")
    )

  if (!is.null(pooled_overlay)) {
    p <- p +
      geom_hline(
        data = pooled_overlay, aes(yintercept = pct_gap),
        linetype = "dashed", colour = "#b2182b", linewidth = 0.8
      ) +
      geom_hline(
        data = pooled_overlay, aes(yintercept = pct_gap_lo),
        linetype = "dotted", colour = "#b2182b", linewidth = 0.45, alpha = 0.7
      ) +
      geom_hline(
        data = pooled_overlay, aes(yintercept = pct_gap_hi),
        linetype = "dotted", colour = "#b2182b", linewidth = 0.45, alpha = 0.7
      ) +
      annotate(
        "text",
        x = 4.35, y = pooled_overlay$pct_gap,
        label = "Pooled estimate\n(all years)",
        hjust = 0, vjust = -0.3, size = 3.2, colour = "#b2182b", fontface = "italic"
      )
  }
  p
}

h2_plot <- prepare_phase_plot_data(h2_out)
h3_plot <- prepare_phase_plot_data(h3_out)

dx_h2 <- h2_out$delta_x_iqr_FP_between[1]
r0_pooled <- mean(h2_out$baseline_ratio)
pooled_overlay <- NULL
if (pooled_fitted) {
  pooled_overlay <- tibble(
    pct_gap = pct_gap_change_from_slope(pooled_out$estimate, dx_h2, r0_pooled),
    pct_gap_lo = pct_gap_change_from_slope(pooled_out$ci_lo, dx_h2, r0_pooled),
    pct_gap_hi = pct_gap_change_from_slope(pooled_out$ci_hi, dx_h2, r0_pooled)
  )
  logmsg(sprintf(
    "H2 pooled overlay (IQR scale): %+.3f%% [%.3f%%, %+.3f%%]",
    pooled_overlay$pct_gap, pooled_overlay$pct_gap_lo, pooled_overlay$pct_gap_hi
  ))
}

p_h2 <- plot_gap_by_phase(
  h2_plot,
  title = "H2: fishing pressure\u2019s spatial effect by policy-anchored phase",
  subtitle = paste0(
    "Between-rectangle effect (", h2_slope_model, ", phase_v2); IQR of FP_between (",
    sprintf("%.2f", dx_h2),
    " log-hours).",
    if (pooled_fitted) " Dashed red = pooled FP_between (no phase interaction)." else ""
  ),
  pooled_overlay = pooled_overlay
)
ggsave(path_out_fig_h2, p_h2, width = 8.5, height = 5.5, dpi = 150)
logmsg("Saved: ", path_out_fig_h2)

p_h3 <- plot_gap_by_phase(
  h3_plot,
  title = "H3: fishing pressure\u2019s temporal effect by policy-anchored phase",
  subtitle = paste0(
    "Within-rectangle effect (", h3_slope_model, ", phase_v2); phase-specific IQR of FP_within."
  )
)
ggsave(path_out_fig_h3, p_h3, width = 8.5, height = 5.5, dpi = 150)
logmsg("Saved: ", path_out_fig_h3)

for (i in seq_len(nrow(h2_plot))) {
  r <- h2_plot[i, ]
  logmsg(sprintf(
    "  H2 fig %s: gap %+.2f%% [%.2f, %+.2f]; %s",
    r$phase, r$pct_gap, r$pct_gap_lo, r$pct_gap_hi,
    if (r$significant) "significant" else "not significant"
  ))
}
for (i in seq_len(nrow(h3_plot))) {
  r <- h3_plot[i, ]
  logmsg(sprintf(
    "  H3 fig %s: gap %+.2f%% [%.2f, %+.2f]; %s",
    r$phase, r$pct_gap, r$pct_gap_lo, r$pct_gap_hi,
    if (r$significant) "significant" else "not significant"
  ))
}

# ---------------------------------------------------------------------------
# Persist model objects
# ---------------------------------------------------------------------------
saveRDS(
  list(
    primary_model_v2 = fit_wb,
    fit_wb_car_v2 = fit_car,
    fit_pooled_between_car_v2 = fit_pooled,
    data = dat,
    adjMatrix = if (car_fitted) adjMatrix else NULL,
    formula_v2 = residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec),
    formula_car_v2 = if (car_fitted) {
      residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
    } else {
      NULL
    },
    formula_pooled_v2 = if (pooled_fitted) {
      residual ~ FP_between + FP_within * phase_v2 + adjacency(1 | stat_rec)
    } else {
      NULL
    },
    phase_v2_levels = PHASE_V2,
    car_fitted = car_fitted,
    pooled_fitted = pooled_fitted,
    slopes = slopes_all,
    proportional_H2 = h2_out,
    proportional_H3 = h3_out,
    pooled_coef = if (pooled_fitted) pooled_out else NULL
  ),
  path_out_models
)
logmsg("Saved: ", path_out_models)

# Also enrich primary_model_v2.rds with CAR if present (non-destructive add)
v2$fit_wb_car_v2 <- fit_car
v2$fit_pooled_between_car_v2 <- fit_pooled
v2$car_fitted <- car_fitted
v2$reporting_source <- path_out_models
saveRDS(v2, path_v2)
logmsg("Updated ", path_v2, " with CAR / pooled fits.")

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_slopes)
logmsg("- ", path_out_fe)
if (car_fitted) logmsg("- ", path_out_car_fe)
logmsg("- ", path_out_h2)
logmsg("- ", path_out_h3)
if (pooled_fitted) logmsg("- ", path_out_pooled)
logmsg("- ", path_out_models)
logmsg("- ", path_out_fig_combined)
logmsg("- ", path_out_fig_h2)
logmsg("- ", path_out_fig_h3)
logmsg("- ", path_out_session)
logmsg("- ", path_out_run_log, " (this file)")
logmsg("")
logmsg(
  "display_discussion drafts not rewritten by this script — use these phase_v2_* ",
  "artifacts when updating One page / H2 / H3 drafts."
)
logmsg("")
logmsg(
  "Supplementary CAR diagnostics (leave-one-rectangle-out, CAR R², signed vs ",
  "unsigned residual audit) are a separate step that reads this fit_car and ",
  "does not refit the signed primary model: ",
  "Rscript --vanilla pipeline/run_h2h3_car_reporting_diagnostics.R"
)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== phase_v2 reporting complete. ===\n")
