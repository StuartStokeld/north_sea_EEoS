# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Final presentation figures — one standalone figure per hypothesis (H2, H3).
#
# Reads already-computed proportional effect tables only; no model refitting.
# H2 pooled overlay: slope + CI from h2h3_wb_pooled_between_coef.csv, converted
# to the same IQR gap-change scale using delta_x and baseline conventions from
# h2h3_wb_proportional_effects_H2.csv (same transform as run_h2h3_wb_proportional_effects.R).
#
# Run: Rscript --vanilla pipeline/run_h2h3_presentation_figures.R

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
  stop("Run from pipeline/ or Rscript --vanilla pipeline/run_h2h3_presentation_figures.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_h2 <- file.path(project_root, "outputs", "h2h3_wb_proportional_effects_H2.csv")
path_h3 <- file.path(project_root, "outputs", "h2h3_wb_proportional_effects_H3.csv")
path_pooled <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_coef.csv")

path_fig_h2 <- file.path(fig_dir, "h2h3_presentation_H2_gap_change_by_phase.png")
path_fig_h3 <- file.path(fig_dir, "h2h3_presentation_H3_gap_change_by_phase.png")
path_run_log <- file.path(project_root, "outputs", "h2h3_presentation_figures_run_log.md")

stopifnot(file.exists(path_h2), file.exists(path_h3), file.exists(path_pooled))

PHASE_LEVELS <- c("1985-1988", "1989-2000", "2001-2007", "2008-2015")
Y_LABEL <- paste0(
  "Change in overprediction gap (% of remaining gap)\n",
  "(+ above zero = gap closed; \u2212 below zero = gap widened)"
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

pct_gap_change_from_slope <- function(slope, dx, r0) {
  gap0 <- 1 - r0
  gap1 <- 1 - r0 * exp(slope * dx)
  ((gap0 - gap1) / gap0) * 100
}

prepare_phase_plot_data <- function(df) {
  df %>%
    mutate(
      phase = factor(phase, levels = PHASE_LEVELS),
      pct_gap = pct_gap_change_iqr,
      pct_gap_lo = pct_gap_change_iqr_lo,
      pct_gap_hi = pct_gap_change_iqr_hi,
      significant = (pct_gap_lo > 0) | (pct_gap_hi < 0)
    )
}

plot_gap_by_phase <- function(
    plot_df,
    title,
    subtitle = NULL,
    pooled_overlay = NULL
) {
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
    labs(x = "Phase", y = Y_LABEL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5, colour = "grey30")
    )

  if (!is.null(pooled_overlay)) {
    p <- p +
      geom_hline(
        data = pooled_overlay,
        aes(yintercept = pct_gap),
        linetype = "dashed", colour = "#b2182b", linewidth = 0.8
      ) +
      geom_hline(
        data = pooled_overlay,
        aes(yintercept = pct_gap_lo),
        linetype = "dotted", colour = "#b2182b", linewidth = 0.45, alpha = 0.7
      ) +
      geom_hline(
        data = pooled_overlay,
        aes(yintercept = pct_gap_hi),
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

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
h2 <- read_csv(path_h2, show_col_types = FALSE) %>% prepare_phase_plot_data()
h3 <- read_csv(path_h3, show_col_types = FALSE) %>% prepare_phase_plot_data()
pooled_coef <- read_csv(path_pooled, show_col_types = FALSE) %>%
  filter(term == "FP_between")

if (nrow(h2) != 4L || nrow(h3) != 4L || nrow(pooled_coef) != 1L) {
  stop("Expected 4 H2 rows, 4 H3 rows, and 1 pooled FP_between row.")
}

logmsg("# Final presentation figures (H2 / H3) — run log")
logmsg("")
logmsg("Standalone supervisor figures. No models refit. All phase bar values read directly from proportional-effects CSVs (IQR convention).")

# ---------------------------------------------------------------------------
# H2 pooled overlay — convert pooled slope using H2 CSV conventions
# ---------------------------------------------------------------------------
dx_h2 <- h2$delta_x_iqr_FP_between[1]
r0_pooled <- mean(h2$baseline_ratio)
pooled_overlay <- tibble(
  pct_gap = pct_gap_change_from_slope(pooled_coef$estimate, dx_h2, r0_pooled),
  pct_gap_lo = pct_gap_change_from_slope(pooled_coef$ci_lo, dx_h2, r0_pooled),
  pct_gap_hi = pct_gap_change_from_slope(pooled_coef$ci_hi, dx_h2, r0_pooled)
)

logmsg("")
logmsg("## Source files")
logmsg("- H2 phase bars: ", path_h2)
logmsg("- H3 phase bars: ", path_h3)
logmsg("- H2 pooled overlay: ", path_pooled)

logmsg("")
logmsg("## H2 pooled overlay conversion (same IQR scale as phase bars)")
logmsg(sprintf(
  "Pooled slope = %+.6f (95%% CI [%+.6f, %+.6f]) from pooled coef CSV.",
  pooled_coef$estimate, pooled_coef$ci_lo, pooled_coef$ci_hi
))
logmsg(sprintf(
  "Delta-x = IQR(FP_between) = %.4f from H2 proportional CSV; baseline ratio = mean of 4 phase baselines = %.4f from H2 proportional CSV.",
  dx_h2, r0_pooled
))
logmsg(sprintf(
  "Converted gap change (IQR): %+.3f%% [%.3f%%, %+.3f%%] (transform from run_h2h3_wb_proportional_effects.R; not a refit).",
  pooled_overlay$pct_gap, pooled_overlay$pct_gap_lo, pooled_overlay$pct_gap_hi
))

# ---------------------------------------------------------------------------
# Figure 1: H2
# ---------------------------------------------------------------------------
p_h2 <- plot_gap_by_phase(
  h2,
  title = "H2: fishing pressure\u2019s spatial effect reverses direction over time",
  subtitle = paste0(
    "Between-rectangle effect (CAR model); IQR of FP_between across rectangles (",
    sprintf("%.2f", dx_h2), " log-hours units). Dashed red = pooled FP_between (no phase interaction)."
  ),
  pooled_overlay = pooled_overlay
)
ggsave(path_fig_h2, p_h2, width = 8.5, height = 5.5, dpi = 150)

logmsg("")
logmsg("## Figure 1 (H2)")
logmsg("Saved: ", path_fig_h2)
for (i in seq_len(nrow(h2))) {
  r <- h2[i, ]
  logmsg(sprintf(
    "  - %s: gap change %+.2f%% [%.2f%%, %+.2f%%]; %s",
    r$phase, r$pct_gap, r$pct_gap_lo, r$pct_gap_hi,
    if (r$significant) "significant" else "not significant"
  ))
}

# ---------------------------------------------------------------------------
# Figure 2: H3
# ---------------------------------------------------------------------------
p_h3 <- plot_gap_by_phase(
  h3,
  title = "H3: fishing pressure\u2019s temporal effect is intermittent",
  subtitle = "Within-rectangle effect (primary model); phase-specific IQR of FP_within (log-hours deviation)."
)
ggsave(path_fig_h3, p_h3, width = 8.5, height = 5.5, dpi = 150)

logmsg("")
logmsg("## Figure 2 (H3)")
logmsg("Saved: ", path_fig_h3)
logmsg("No pooled overlay (FP_within pooled contrast not computed).")
for (i in seq_len(nrow(h3))) {
  r <- h3[i, ]
  logmsg(sprintf(
    "  - %s: gap change %+.2f%% [%.2f%%, %+.2f%%]; IQR(FP_within)=%.3f; %s",
    r$phase, r$pct_gap, r$pct_gap_lo, r$pct_gap_hi, r$delta_x_iqr_FP_within,
    if (r$significant) "significant" else "not significant"
  ))
}

logmsg("")
logmsg("## Outputs")
logmsg("- ", path_fig_h2)
logmsg("- ", path_fig_h3)
logmsg("- ", path_run_log, " (this file)")

writeLines(run_log, path_run_log)
cat("\nSaved run log:", path_run_log, "\n")
cat("=== Presentation figures complete. ===\n")
