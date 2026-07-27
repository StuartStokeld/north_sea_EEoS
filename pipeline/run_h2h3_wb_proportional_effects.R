# H2/H3 proportional effect sizes for the WITHIN-BETWEEN decomposed slopes
#
# PURPOSE: re-express already-fitted FP_between (H2) and FP_within (H3) phase
# slopes on proportional / gap-change scales. No model refitting.
#
# SUPERSEDES: outputs/h2h3_results_proportional_effects.csv (blended-term table).
#
# H2 source: CAR (wb_car) FP_between phase slopes
# H3 source: primary (wb_primary) FP_within phase slopes (CAR within slopes
#            agree closely; primary used as default — noted in run log)
#
# Benchmarks:
#   H2: doubling (Δx = log(2)) AND IQR of FP_between across rectangles
#   H3: IQR of FP_within only — doubling NOT applicable (mean-zero deviations)
#
# Baseline ratio: exp(phase-specific median residual)
#
# Run: Rscript --vanilla pipeline/run_h2h3_wb_proportional_effects.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_wb_proportional_effects.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_slopes <- file.path(project_root, "outputs", "h2h3_wb_fp_slopes_by_phase.csv")
stopifnot(file.exists(path_models), file.exists(path_slopes))

path_out_h2 <- file.path(project_root, "outputs", "h2h3_wb_proportional_effects_H2.csv")
path_out_h3 <- file.path(project_root, "outputs", "h2h3_wb_proportional_effects_H3.csv")
path_out_fig <- file.path(fig_dir, "h2h3_wb_gap_change_by_phase.png")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_wb_proportional_effects_run_log.md")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 within-between proportional effect sizes — run log")
logmsg("")
logmsg(
  "Reporting only — no models refit. Re-expresses decomposed FP_between (H2) and ",
  "FP_within (H3) phase slopes on proportional / gap-change scales."
)
logmsg("")
logmsg(
  "SUPERSEDES the earlier blended-term proportional table ",
  "(`outputs/h2h3_results_proportional_effects.csv` / ",
  "`h2h3_results_gap_closed_by_phase.png`). That table used the undecomposed ",
  "`log_hours_total` slopes and is no longer the primary proportional reporting ",
  "artefact. Use the H2/H3 tables produced here instead."
)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
mod <- readRDS(path_models)
dat <- mod$data
slopes <- read_csv(path_slopes, show_col_types = FALSE)

h2_slopes <- slopes %>%
  filter(model_id == "wb_car", component == "FP_between")
h3_slopes <- slopes %>%
  filter(model_id == "wb_primary", component == "FP_within")
h3_car <- slopes %>%
  filter(model_id == "wb_car", component == "FP_within")

if (nrow(h2_slopes) != 4L || nrow(h3_slopes) != 4L) {
  stop("Expected 4 H2 (CAR between) and 4 H3 (primary within) phase slopes.")
}

logmsg("")
logmsg("## Coefficient sources")
logmsg("H2: CAR model (`wb_car`) FP_between phase slopes (per design).")
logmsg(
  "H3: primary model (`wb_primary`) FP_within phase slopes. CAR within slopes agree ",
  "closely (max |primary − CAR| phase-slope difference logged below); primary used as ",
  "the reported H3 source."
)
cmp_h3 <- h3_slopes %>%
  select(phase, primary = fp_slope) %>%
  left_join(h3_car %>% select(phase, car = fp_slope), by = "phase") %>%
  mutate(abs_diff = abs(primary - car))
logmsg(sprintf(
  "Max |primary − CAR| FP_within phase slope = %.5f (phases agree for reporting).",
  max(cmp_h3$abs_diff)
))

# ---------------------------------------------------------------------------
# Baselines and Δx benchmarks
# ---------------------------------------------------------------------------
phase_base <- dat %>%
  group_by(phase) %>%
  summarise(
    n_hauls = dplyr::n(),
    median_residual = median(residual, na.rm = TRUE),
    baseline_ratio = exp(median_residual),
    iqr_FP_within = IQR(FP_within, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(phase = as.character(phase))

# FP_between is time-invariant: one value per rectangle
rect_between <- dat %>%
  distinct(stat_rec, FP_between)
iqr_FP_between <- IQR(rect_between$FP_between, na.rm = TRUE)

logmsg("")
logmsg("## Baseline and Δx benchmarks")
logmsg(
  "Baseline for percentage-point ratio shift and gap-change: each phase's MEDIAN ",
  "residual; starting ratio = exp(median residual). Same convention as the (now ",
  "superseded) blended-term proportional table."
)
for (i in seq_len(nrow(phase_base))) {
  r <- phase_base[i, ]
  logmsg(sprintf(
    "  - %s: median residual=%.4f; baseline ratio=%.4f (%.1f%%); phase IQR(FP_within)=%.3f",
    r$phase, r$median_residual, r$baseline_ratio, 100 * r$baseline_ratio, r$iqr_FP_within
  ))
}
logmsg(sprintf(
  paste0(
    "H2 IQR benchmark: IQR(FP_between) across %d rectangles = %.4f ",
    "(applied to all phases; FP_between is time-invariant)."
  ),
  nrow(rect_between), iqr_FP_between
))
logmsg(
  "H3 IQR benchmark: phase-specific IQR(FP_within) among hauls in that phase ",
  "(realistic within-rectangle year-to-year spread in that period)."
)
logmsg(
  "H3 doubling convention: NOT APPLICABLE. FP_within is a mean-zero deviation ",
  "(can be negative); 'doubling fishing hours' is not a meaningful contrast on this ",
  "scale. Only the IQR-of-deviation benchmark is reported for H3."
)

# ---------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------
prop_metrics <- function(slope, slope_lo, slope_hi, dx, r0) {
  pct_ratio <- function(b) (exp(b * dx) - 1) * 100
  r1 <- function(b) r0 * exp(b * dx)
  pp_shift <- function(b) (r1(b) - r0) * 100
  gap0 <- 1 - r0
  # signed: positive = gap shrinks (closed); negative = gap widens
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
    ratio_end_lo = r1(slope_lo),
    ratio_end_hi = r1(slope_hi),
    pp_ratio_shift = pp_shift(slope),
    pp_ratio_shift_lo = pp_shift(slope_lo),
    pp_ratio_shift_hi = pp_shift(slope_hi),
    pct_gap_change = pct_gap_change(slope),
    pct_gap_change_lo = pct_gap_change(slope_lo),
    pct_gap_change_hi = pct_gap_change(slope_hi),
    gap_direction = gap_dir(slope)
  )
}

# ---------------------------------------------------------------------------
# H2 table (CAR between): doubling + IQR
# ---------------------------------------------------------------------------
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
    model_id = "wb_car",
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
    # doubling
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
    # IQR of FP_between across rectangles
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
  mutate(phase = factor(phase, levels = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"))) %>%
  arrange(phase)
write_csv(h2_out, path_out_h2)

logmsg("")
logmsg("## H2 table (FP_between, CAR)")
logmsg("Saved: ", path_out_h2)
logmsg("Gap labelling: positive pct_gap_change = gap closed/narrowed; negative = gap widened.")
for (i in seq_len(nrow(h2_out))) {
  r <- h2_out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: slope=%+.4f; doubling: ratio change=%+.2f%%, gap %s by %.2f%% ",
      "[%.2f, %.2f]; IQR(Δx=%.3f): ratio change=%+.2f%%, gap %s by %.2f%%"
    ),
    as.character(r$phase), r$fp_slope,
    r$pct_ratio_change_doubling, r$gap_direction_doubling, abs(r$pct_gap_change_doubling),
    r$pct_gap_change_doubling_lo, r$pct_gap_change_doubling_hi,
    r$delta_x_iqr_FP_between, r$pct_ratio_change_iqr, r$gap_direction_iqr,
    abs(r$pct_gap_change_iqr)
  ))
}

# ---------------------------------------------------------------------------
# H3 table (FP_within): IQR only — no doubling
# ---------------------------------------------------------------------------
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
    model_id = "wb_primary",
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
  mutate(phase = factor(phase, levels = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"))) %>%
  arrange(phase)
write_csv(h3_out, path_out_h3)

logmsg("")
logmsg("## H3 table (FP_within, primary) — IQR-of-deviation only")
logmsg("Saved: ", path_out_h3)
logmsg(
  "NOTE: No doubling-convention columns for H3. FP_within is mean-zero by construction; ",
  "a doubling benchmark does not apply. See doubling_note column in the CSV."
)
for (i in seq_len(nrow(h3_out))) {
  r <- h3_out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: slope=%+.4f; IQR(Δx=%.3f): ratio change=%+.2f%% [%+.2f, %+.2f]; ",
      "ratio %.1f%% → %.1f%% (%+.2f pp); gap %s by %.2f%% [%.2f, %.2f]"
    ),
    as.character(r$phase), r$fp_slope, r$delta_x_iqr_FP_within,
    r$pct_ratio_change_iqr, r$pct_ratio_change_iqr_lo, r$pct_ratio_change_iqr_hi,
    r$ratio_start_pct, r$ratio_end_iqr_pct, r$pp_ratio_shift_iqr,
    r$gap_direction_iqr, abs(r$pct_gap_change_iqr),
    r$pct_gap_change_iqr_lo, r$pct_gap_change_iqr_hi
  ))
}

# ---------------------------------------------------------------------------
# Figure: gap change by phase — H2 (doubling) and H3 (IQR) as two series
# ---------------------------------------------------------------------------
# Use comparable "headline" benchmarks: H2 doubling (natural for between means);
# H3 IQR-of-within (only applicable benchmark). Label axes/legend accordingly.
plot_df <- bind_rows(
  h2_out %>%
    transmute(
      phase, year_mid = (year_start + year_end) / 2,
      series = "H2: FP_between (CAR), per doubling",
      pct_gap = pct_gap_change_doubling,
      pct_gap_lo = pct_gap_change_doubling_lo,
      pct_gap_hi = pct_gap_change_doubling_hi,
      gap_direction = gap_direction_doubling
    ),
  h3_out %>%
    transmute(
      phase, year_mid = (year_start + year_end) / 2,
      series = "H3: FP_within (primary), per IQR of within-deviation",
      pct_gap = pct_gap_change_iqr,
      pct_gap_lo = pct_gap_change_iqr_lo,
      pct_gap_hi = pct_gap_change_iqr_hi,
      gap_direction = gap_direction_iqr
    )
) %>%
  mutate(
    series = factor(
      series,
      levels = c(
        "H2: FP_between (CAR), per doubling",
        "H3: FP_within (primary), per IQR of within-deviation"
      )
    )
  )

col_map <- c(
  "H2: FP_between (CAR), per doubling" = "#d73027",
  "H3: FP_within (primary), per IQR of within-deviation" = "#4575b4"
)

p <- ggplot(plot_df, aes(x = year_mid, y = pct_gap, colour = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = c(1989, 2001, 2008), linetype = "dotted", colour = "grey40", linewidth = 0.4) +
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
    title = "Within-between H2/H3: proportional gap change by phase",
    subtitle = paste0(
      "Positive = gap closed/narrowed; negative = gap widened. ",
      "H2 benchmark = doubling of rectangle-mean fishing hours (CAR between slopes). ",
      "H3 benchmark = phase IQR of FP_within (doubling not applicable). ",
      "Baseline = phase median residual. Dotted lines = breaks 1989, 2001, 2008."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.subtitle = element_text(size = 7.8, colour = "grey30"),
    panel.grid.minor = element_blank()
  )

ggsave(path_out_fig, p, width = 10, height = 6, dpi = 150)
logmsg("")
logmsg("## Figure")
logmsg("Saved: ", path_out_fig)
logmsg(
  "Two series: H2 gap-change under doubling; H3 gap-change under phase IQR(FP_within). ",
  "Sign convention labelled on the y-axis (closed vs widened)."
)

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_h2)
logmsg("- ", path_out_h3)
logmsg("- ", path_out_fig)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 within-between proportional effect-size reporting complete. ===\n")
