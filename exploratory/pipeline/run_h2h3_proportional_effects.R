# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H2/H3 proportional effect-size reporting (supplementary)
#
# PURPOSE: re-express already-fitted primary-model phase slopes on
# proportional / gap-closure scales for all four phases. No model refitting.
#
# Source slopes: outputs/h2h3_results_fp_slopes_by_phase.csv (primary_plain_re)
# Source data:   outputs/h2h3_results_model_objects.rds (analysis frame)
#
# For each phase and each predictor contrast Δx in {log(2), phase IQR of
# log(hours+1)}:
#   1. % change in B_obs/B_pred ratio: (exp(slope * Δx) - 1) * 100
#   2. percentage-point shift in the ratio, from phase baseline residual
#   3. % of remaining gap closed: ((1-r0)-(1-r1))/(1-r0)*100
# CIs: same transforms applied to slope 95% CI bounds (Δx > 0, monotonic).
#
# Run: Rscript --vanilla pipeline/run_h2h3_proportional_effects.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_proportional_effects.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_models <- file.path(project_root, "outputs", "h2h3_results_model_objects.rds")
path_slopes <- file.path(project_root, "outputs", "h2h3_results_fp_slopes_by_phase.csv")
stopifnot(file.exists(path_models), file.exists(path_slopes))

path_out_table <- file.path(project_root, "outputs", "h2h3_results_proportional_effects.csv")
path_out_fig <- file.path(fig_dir, "h2h3_results_gap_closed_by_phase.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_results_proportional_effects_run_log.md")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 proportional effect-size reporting — run log")
logmsg("")
logmsg(
  "Reporting/formatting only. No models refit. Re-expresses primary-model phase slopes ",
  "from the biomass-free results run on proportional and gap-closure scales."
)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
mod <- readRDS(path_models)
dat <- mod$data
slopes <- read_csv(path_slopes, show_col_types = FALSE) %>%
  filter(model_id == "primary_plain_re")

if (nrow(slopes) != 4L) {
  stop("Expected 4 primary phase-slope rows; found ", nrow(slopes))
}
if (isTRUE(mod$biomass_removed) || !"mean_ln_B_obs" %in% all.vars(stats::formula(mod$fit_primary))) {
  logmsg("Confirmed source model is biomass-free primary specification.")
} else {
  logmsg("FLAG: source model appears to include biomass — check upstream results run.")
}

# ---------------------------------------------------------------------------
# Phase baseline choice (item 2/3)
# ---------------------------------------------------------------------------
# Use each phase's MEDIAN residual as the baseline for ratio / gap calculations.
# Reason: residual is left-skewed within phases (long negative tail); the median
# matches the worked example already in the interpretation note (overall median)
# and is a more typical "starting haul" than the mean, which is pulled toward
# less-negative values by the upper tail. Mean residual is logged alongside for
# transparency but is NOT used in the reported metrics.
phase_base <- dat %>%
  group_by(phase) %>%
  summarise(
    n_hauls = dplyr::n(),
    median_residual = median(residual, na.rm = TRUE),
    mean_residual = mean(residual, na.rm = TRUE),
    baseline_ratio = exp(median_residual),
    mean_ratio_for_reference = exp(mean_residual),
    log_hours_q25 = as.numeric(stats::quantile(log_hours_total, 0.25, na.rm = TRUE)),
    log_hours_q75 = as.numeric(stats::quantile(log_hours_total, 0.75, na.rm = TRUE)),
    log_hours_iqr = log_hours_q75 - log_hours_q25,
    hours_q25 = as.numeric(stats::quantile(hours_total, 0.25, na.rm = TRUE)),
    hours_q75 = as.numeric(stats::quantile(hours_total, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

logmsg("")
logmsg("## Phase-specific baseline choice")
logmsg(
  "Baseline for percentage-point ratio shift and percent-gap-closed: each phase's ",
  "MEDIAN residual (ratio = exp(median residual)). Mean residual is recorded but not used."
)
logmsg("Why median: residual is left-skewed within phases; median is the typical haul and ",
       "matches the interpretation note's worked-example convention.")
for (i in seq_len(nrow(phase_base))) {
  r <- phase_base[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: n=%d; median residual=%.4f (ratio=%.4f); mean residual=%.4f (ratio=%.4f); ",
      "log(hours+1) IQR=%.3f (hours Q25–Q75: %.0f–%.0f)"
    ),
    as.character(r$phase), r$n_hauls, r$median_residual, r$baseline_ratio,
    r$mean_residual, r$mean_ratio_for_reference, r$log_hours_iqr, r$hours_q25, r$hours_q75
  ))
}

# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------
#' Apply slope (and CI bounds) through proportional metrics for a given Δx.
prop_from_slope <- function(slope, slope_lo, slope_hi, dx, r0) {
  # % change in ratio
  pct_ratio <- function(b) (exp(b * dx) - 1) * 100
  # ending ratio and percentage-point shift
  r1 <- function(b) r0 * exp(b * dx)
  pp_shift <- function(b) (r1(b) - r0) * 100
  # % of remaining gap closed (gap = 1 - ratio)
  gap0 <- 1 - r0
  pct_gap <- function(b) {
    g1 <- 1 - r1(b)
    ((gap0 - g1) / gap0) * 100
  }
  list(
    pct_ratio_change = pct_ratio(slope),
    pct_ratio_change_lo = pct_ratio(slope_lo),
    pct_ratio_change_hi = pct_ratio(slope_hi),
    ratio_start = r0,
    ratio_end = r1(slope),
    ratio_end_lo = r1(slope_lo),
    ratio_end_hi = r1(slope_hi),
    pp_ratio_shift = pp_shift(slope),
    pp_ratio_shift_lo = pp_shift(slope_lo),
    pp_ratio_shift_hi = pp_shift(slope_hi),
    pct_gap_closed = pct_gap(slope),
    pct_gap_closed_lo = pct_gap(slope_lo),
    pct_gap_closed_hi = pct_gap(slope_hi)
  )
}

joined <- slopes %>%
  mutate(phase = as.character(phase)) %>%
  left_join(phase_base %>% mutate(phase = as.character(phase)), by = "phase")

rows <- lapply(seq_len(nrow(joined)), function(i) {
  r <- joined[i, ]
  r0 <- r$baseline_ratio
  dx_double <- log(2)
  dx_iqr <- r$log_hours_iqr

  d_dbl <- prop_from_slope(r$fp_slope, r$fp_slope_lo, r$fp_slope_hi, dx_double, r0)
  d_iqr <- prop_from_slope(r$fp_slope, r$fp_slope_lo, r$fp_slope_hi, dx_iqr, r0)

  data.frame(
    phase = r$phase,
    year_start = r$year_start,
    year_end = r$year_end,
    n_hauls = r$n_hauls,
    fp_slope = r$fp_slope,
    fp_slope_lo = r$fp_slope_lo,
    fp_slope_hi = r$fp_slope_hi,
    baseline = "phase_median_residual",
    baseline_residual = r$median_residual,
    baseline_ratio = r0,
    log_hours_iqr = dx_iqr,
    hours_q25 = r$hours_q25,
    hours_q75 = r$hours_q75,
    # --- doubling ---
    delta_x_doubling = dx_double,
    pct_ratio_change_doubling = d_dbl$pct_ratio_change,
    pct_ratio_change_doubling_lo = d_dbl$pct_ratio_change_lo,
    pct_ratio_change_doubling_hi = d_dbl$pct_ratio_change_hi,
    ratio_start_pct = r0 * 100,
    ratio_end_doubling_pct = d_dbl$ratio_end * 100,
    ratio_end_doubling_pct_lo = d_dbl$ratio_end_lo * 100,
    ratio_end_doubling_pct_hi = d_dbl$ratio_end_hi * 100,
    pp_ratio_shift_doubling = d_dbl$pp_ratio_shift,
    pp_ratio_shift_doubling_lo = d_dbl$pp_ratio_shift_lo,
    pp_ratio_shift_doubling_hi = d_dbl$pp_ratio_shift_hi,
    pct_gap_closed_doubling = d_dbl$pct_gap_closed,
    pct_gap_closed_doubling_lo = d_dbl$pct_gap_closed_lo,
    pct_gap_closed_doubling_hi = d_dbl$pct_gap_closed_hi,
    # --- phase IQR ---
    delta_x_iqr = dx_iqr,
    pct_ratio_change_iqr = d_iqr$pct_ratio_change,
    pct_ratio_change_iqr_lo = d_iqr$pct_ratio_change_lo,
    pct_ratio_change_iqr_hi = d_iqr$pct_ratio_change_hi,
    ratio_end_iqr_pct = d_iqr$ratio_end * 100,
    ratio_end_iqr_pct_lo = d_iqr$ratio_end_lo * 100,
    ratio_end_iqr_pct_hi = d_iqr$ratio_end_hi * 100,
    pp_ratio_shift_iqr = d_iqr$pp_ratio_shift,
    pp_ratio_shift_iqr_lo = d_iqr$pp_ratio_shift_lo,
    pp_ratio_shift_iqr_hi = d_iqr$pp_ratio_shift_hi,
    pct_gap_closed_iqr = d_iqr$pct_gap_closed,
    pct_gap_closed_iqr_lo = d_iqr$pct_gap_closed_lo,
    pct_gap_closed_iqr_hi = d_iqr$pct_gap_closed_hi,
    stringsAsFactors = FALSE
  )
})

out <- bind_rows(rows) %>%
  mutate(
    phase = factor(phase, levels = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"))
  ) %>%
  arrange(phase)

write_csv(out, path_out_table)
logmsg("")
logmsg("## Table")
logmsg("Saved: ", path_out_table)
logmsg("")
logmsg("Column groups:")
logmsg("  - pct_ratio_change_* : percent change in B_obs/B_pred (= (exp(slope*Δx)-1)*100)")
logmsg("  - pp_ratio_shift_*   : percentage-point change in the ratio (e.g. 25.0% → 26.4% = +1.4 pp)")
logmsg("  - pct_gap_closed_*   : percent of remaining gap (1 - ratio) closed; negative = gap widens")
logmsg("  - *_doubling / *_iqr : Δx = log(2) vs phase-specific IQR of log(hours+1)")
logmsg("  - *_lo / *_hi        : same transforms of the slope 95% CI bounds")

logmsg("")
logmsg("### Summary (doubling of fishing hours)")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: slope=%+.4f; %% ratio change=%+.2f%% [%+.2f, %+.2f]; ",
      "ratio %.1f%% → %.1f%% (%+.2f pp); gap closed=%+.2f%% [%+.2f, %+.2f]"
    ),
    as.character(r$phase), r$fp_slope,
    r$pct_ratio_change_doubling, r$pct_ratio_change_doubling_lo, r$pct_ratio_change_doubling_hi,
    r$ratio_start_pct, r$ratio_end_doubling_pct, r$pp_ratio_shift_doubling,
    r$pct_gap_closed_doubling, r$pct_gap_closed_doubling_lo, r$pct_gap_closed_doubling_hi
  ))
}
logmsg("")
logmsg("### Summary (phase log-hours IQR)")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: Δx(IQR)=%.3f; %% ratio change=%+.2f%% [%+.2f, %+.2f]; ",
      "ratio %.1f%% → %.1f%% (%+.2f pp); gap closed=%+.2f%% [%+.2f, %+.2f]"
    ),
    as.character(r$phase), r$delta_x_iqr,
    r$pct_ratio_change_iqr, r$pct_ratio_change_iqr_lo, r$pct_ratio_change_iqr_hi,
    r$ratio_start_pct, r$ratio_end_iqr_pct, r$pp_ratio_shift_iqr,
    r$pct_gap_closed_iqr, r$pct_gap_closed_iqr_lo, r$pct_gap_closed_iqr_hi
  ))
}

# ---------------------------------------------------------------------------
# Figure: % gap closed by phase (doubling), CI bars
# ---------------------------------------------------------------------------
# Place phases on a year axis (phase midpoints) so break-year markers align
# visually with h2h3_results_fp_effect_by_phase.png.
plot_df <- out %>%
  mutate(
    year_mid = (year_start + year_end) / 2,
    phase_label = as.character(phase)
  )

p <- ggplot(plot_df, aes(x = year_mid, y = pct_gap_closed_doubling)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = c(1989, 2001, 2008), linetype = "dotted", colour = "grey40", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = pct_gap_closed_doubling_lo, ymax = pct_gap_closed_doubling_hi),
    width = 1.2, colour = "#d73027", linewidth = 0.7
  ) +
  geom_point(size = 2.8, colour = "#d73027") +
  geom_text(
    aes(label = phase_label, y = pct_gap_closed_doubling_hi),
    vjust = -0.6, size = 3, colour = "grey25"
  ) +
  scale_x_continuous(breaks = seq(1985, 2015, by = 5), limits = c(1984, 2016)) +
  labs(
    x = "Year (phase midpoint)",
    y = "Percent of remaining overprediction gap closed\nper doubling of fishing hours",
    title = "H2/H3: proportional gap closure by phase (doubling convention)",
    subtitle = paste0(
      "Companion to the absolute-scale FP-effect figure. Positive = overprediction gap shrinks ",
      "when fishing hours double; negative = gap widens. Error bars = 95% CI from slope CI. ",
      "Dotted lines = structural breaks 1989, 2001, 2008. Baseline = phase median residual."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.subtitle = element_text(size = 8.2, colour = "grey30"),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = {
    yr <- range(c(plot_df$pct_gap_closed_doubling_lo, plot_df$pct_gap_closed_doubling_hi), finite = TRUE)
    pad <- diff(yr) * 0.18
    c(yr[1] - pad * 0.3, yr[2] + pad)
  })

ggsave(path_out_fig, p, width = 10, height = 5.5, dpi = 150)
logmsg("")
logmsg("## Figure")
logmsg("Saved: ", path_out_fig)
logmsg(
  "Percent-of-gap-closed (doubling) by phase with CI bars; year-midpoint x-axis and ",
  "break-year markers match the absolute-scale FP-effect figure layout."
)

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_table)
logmsg("- ", path_out_fig)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 proportional effect-size reporting complete. ===\n")
