# Final presentation figures — H1.
#
# Figure 1 (Test 1): three-panel unfitted EEoS vs productivity 1:1 (E × m_min).
# Figure 2: fitted ln(E) OLS ceiling + Harte Fig. 2 productivity-ratio structure.
#
# Reads haul benchmarks / comparison / ln(E) coef / Fig 2 leverage tables;
# no model refitting.
#
# Run: Rscript --vanilla pipeline/run_h1_presentation_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript --vanilla pipeline/run_h1_presentation_figures.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

fig_dir <- file.path(project_root, "outputs", "figures")
display_fig_dir <- file.path(project_root, "display_discussion", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(display_fig_dir, recursive = TRUE, showWarnings = FALSE)

path_haul <- file.path(project_root, "outputs", "haul_h1_benchmarks.rds")
path_cmp <- file.path(project_root, "outputs", "h1_model_comparison.rds")
path_lne <- file.path(project_root, "outputs", "h1_lne_coefficients.csv")
path_fig2_lev <- file.path(project_root, "outputs", "h1_fig2_leverage_hauls.csv")
path_fig <- file.path(fig_dir, "h1_presentation_test1_eeos_vs_productivity.png")
path_fig_display <- file.path(
  display_fig_dir, "h1_presentation_test1_eeos_vs_productivity.png"
)
path_fig2 <- file.path(
  fig_dir, "h1_presentation_relative_structure_vs_absolute.png"
)
path_fig2_display <- file.path(
  display_fig_dir, "h1_presentation_relative_structure_vs_absolute.png"
)
path_run_log <- file.path(project_root, "outputs", "h1_presentation_figures_run_log.md")

stopifnot(
  file.exists(path_haul),
  file.exists(path_cmp),
  file.exists(path_lne),
  file.exists(path_fig2_lev)
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

haul <- readRDS(path_haul)
compare <- readRDS(path_cmp)
cmp <- compare$comparison

row_eeos <- cmp[cmp$model_id == "eeos_biomass", , drop = FALSE]
row_prod1 <- cmp[cmp$model_id == "productivity_1to1", , drop = FALSE]
stopifnot(nrow(row_eeos) == 1L, nrow(row_prod1) == 1L)
stopifnot(
  "ln_B_obs" %in% names(haul),
  "ln_B_pred" %in% names(haul),
  "ln_E_calibrated" %in% names(haul),
  "S" %in% names(haul)
)

r2_eeos <- as.numeric(row_eeos$log_r2)
cor2_eeos <- as.numeric(row_eeos$cor2)
r2_prod1 <- as.numeric(row_prod1$log_r2)
cor2_prod1 <- as.numeric(row_prod1$cor2)
med_ratio <- as.numeric(row_eeos$median_ratio)
ss_ratio <- as.numeric(compare$ss_ratio_eeos_over_productivity_1to1)
n_hauls <- nrow(haul)

theme_h1 <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    legend.position = "bottom",
    legend.key.width = unit(1.1, "cm"),
    legend.key.height = unit(0.35, "cm")
  )

# Shared axis range so panels (a) and (b) are visually comparable
axis_vals <- c(haul$ln_B_obs, haul$ln_B_pred, haul$ln_E_calibrated)
axis_vals <- axis_vals[is.finite(axis_vals)]
pad <- 0.04 * (max(axis_vals) - min(axis_vals))
lims <- c(min(axis_vals) - pad, max(axis_vals) + pad)

annot_lab <- function(r2, c2) {
  sprintf("log_r\u00b2 = %.3f\ncor\u00b2 = %.3f", r2, c2)
}

p_eeos <- ggplot(haul, aes(x = ln_B_obs, y = ln_B_pred, colour = log(S))) +
  geom_point(alpha = 0.08, size = 0.45, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey25", linewidth = 0.85) +
  scale_colour_viridis_c(option = "viridis", name = "ln(S)") +
  coord_fixed(xlim = lims, ylim = lims, expand = FALSE) +
  annotate(
    "label",
    x = lims[1] + 0.04 * diff(lims),
    y = lims[2] - 0.04 * diff(lims),
    hjust = 0, vjust = 1, size = 3.1,
    label = annot_lab(r2_eeos, cor2_eeos),
    fill = "white", alpha = 0.92
  ) +
  labs(
    x = expression(log(B[obs]) ~ "[g]"),
    y = expression(log(B[pred]) ~ "[g]"),
    title = "(a) EEoS — Harte Fig 1",
    subtitle = "Unfitted 1:1; median B_pred/B_obs \u2248 4.1\u00d7"
  ) +
  theme_h1 +
  theme(legend.position = "none")

p_prod1 <- ggplot(haul, aes(x = ln_B_obs, y = ln_E_calibrated, colour = log(S))) +
  geom_point(alpha = 0.08, size = 0.45, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey25", linewidth = 0.85) +
  scale_colour_viridis_c(option = "viridis", name = "ln(S)") +
  coord_fixed(xlim = lims, ylim = lims, expand = FALSE) +
  annotate(
    "label",
    x = lims[1] + 0.04 * diff(lims),
    y = lims[2] - 0.04 * diff(lims),
    hjust = 0, vjust = 1, size = 3.1,
    label = annot_lab(r2_prod1, cor2_prod1),
    fill = "white", alpha = 0.92
  ) +
  labs(
    x = expression(log(B[obs]) ~ "[g]"),
    y = expression(log(E %*% m[min]) ~ "[g-equiv]"),
    title = "(b) Productivity 1:1 baseline",
    subtitle = "Unfitted E \u00d7 m_min (same E, m_min as EEoS)"
  ) +
  theme_h1

r2_df <- data.frame(
  model = factor(
    c("EEoS (S, N, E)", "Productivity 1:1"),
    levels = c("Productivity 1:1", "EEoS (S, N, E)")
  ),
  r2 = c(r2_eeos, r2_prod1),
  colour = c("below", "above")
)

p_bar <- ggplot(r2_df, aes(x = model, y = r2, fill = colour)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.8, colour = "grey30") +
  geom_text(
    aes(
      label = sprintf("%.3f", r2),
      vjust = ifelse(r2 < 0, 1.6, -0.45)
    ),
    size = 3.6, fontface = "bold"
  ) +
  scale_fill_manual(values = c("below" = "#d6604d", "above" = "#4393c3")) +
  annotate(
    "text",
    x = 1.5,
    y = max(r2_df$r2) * 0.78,
    size = 3.0,
    colour = "grey35",
    label = sprintf(
      "SS_res ratio = %.1f\u00d7\nHarte criterion: NOT met",
      ss_ratio
    )
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.18))) +
  labs(
    x = NULL,
    y = expression(log_r^2),
    title = "(c) Predictive R\u00b2",
    subtitle = "Zero = no better than the mean"
  ) +
  theme_h1 +
  theme(plot.margin = margin(5.5, 12, 5.5, 5.5))

# Extract shared colour legend from productivity panel
legend_lnS <- get_legend(
  p_prod1 +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.margin = margin(0, 0, 0, 0)
    )
)
p_prod1_noleg <- p_prod1 + theme(legend.position = "none")

panels <- plot_grid(
  p_eeos, p_prod1_noleg, p_bar,
  nrow = 1,
  rel_widths = c(1, 1, 0.95),
  align = "h",
  axis = "tb"
)

p_body <- plot_grid(
  panels,
  legend_lnS,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

title_block <- ggdraw() +
  draw_label(
    "H1: Does EEoS predict haul-level biomass?",
    fontface = "bold",
    size = 14,
    x = 0, hjust = 0, y = 0.72
  ) +
  draw_label(
    sprintf(
      "NS-IBTS Q1 1985\u20132015 · n = %s hauls · unfitted 1:1 maps (no slope/intercept fit)",
      format(n_hauls, big.mark = ",")
    ),
    size = 10,
    colour = "grey30",
    x = 0, hjust = 0, y = 0.28
  ) +
  theme(plot.margin = margin(4, 10, 2, 10))

caption_block <- ggdraw() +
  draw_label(
    paste0(
      "Primary metric: log_r\u00b2 = 1 \u2212 SS_res/SS_tot on log scale. ",
      "cor\u00b2 is diagnostic only (rank agreement). ",
      sprintf("Median B_pred/B_obs = %.1f\u00d7.", med_ratio)
    ),
    size = 8.5,
    colour = "grey40",
    x = 0, hjust = 0, y = 0.5
  ) +
  theme(plot.margin = margin(2, 10, 4, 10))

p_main <- plot_grid(
  title_block,
  p_body,
  caption_block,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.07)
)

ggsave(path_fig, p_main, width = 13.2, height = 5.6, dpi = 180, bg = "white")
ok_copy <- file.copy(path_fig, path_fig_display, overwrite = TRUE)
if (!isTRUE(ok_copy)) {
  warning("Failed to copy figure to display_discussion/figures/")
}

logmsg("Wrote ", path_fig)
logmsg("Copied to ", path_fig_display)
logmsg(sprintf(
  "Fig1 metrics: EEoS log_r2=%.3f cor2=%.3f | prod1 log_r2=%.3f cor2=%.3f | SS ratio=%.2f | median ratio=%.2f",
  r2_eeos, cor2_eeos, r2_prod1, cor2_prod1, ss_ratio, med_ratio
))

# ---------------------------------------------------------------------------
# Figure 2 — relative structure holds; absolute prediction fails
# ---------------------------------------------------------------------------

row_lne <- cmp[cmp$model_id == "ln_e_ols", , drop = FALSE]
row_ratio <- cmp[cmp$model_id == "productivity_ratio", , drop = FALSE]
stopifnot(nrow(row_lne) == 1L, nrow(row_ratio) == 1L)

lne_coef <- read.csv(path_lne, stringsAsFactors = FALSE)
fig2_excl <- read.csv(path_fig2_lev, stringsAsFactors = FALSE)

lne_intercept <- as.numeric(lne_coef$Estimate[lne_coef$term == "(Intercept)"])
lne_slope <- as.numeric(lne_coef$Estimate[lne_coef$term == "log_E"])
r2_lne <- as.numeric(row_lne$log_r2)
ss_ratio_lne <- as.numeric(compare$ss_ratio_eeos_over_lne)

r2_ratio_all <- as.numeric(row_ratio$fig2_r2_pearson_all)
r2_ratio_trim <- as.numeric(row_ratio$fig2_r2_pearson_trimmed)
n_fig2_excl <- as.integer(row_ratio$fig2_r2_pearson_n_excluded)
if (is.na(n_fig2_excl) || n_fig2_excl < 1L) n_fig2_excl <- nrow(fig2_excl)
n_trimmed <- n_hauls - n_fig2_excl

stopifnot(
  "E" %in% names(haul),
  "fig2_predicted_ratio" %in% names(haul),
  "fig2_observed_ratio" %in% names(haul),
  "haul_id" %in% names(haul),
  is.finite(lne_intercept),
  is.finite(lne_slope)
)

haul_f2 <- haul %>%
  mutate(
    ln_E = log(E),
    ln_S = log(S),
    is_leverage = haul_id %in% fig2_excl$haul_id
  )

# (A) Fitted ln(E) OLS ceiling ------------------------------------------------
p_lne <- ggplot(haul_f2, aes(x = ln_E, y = ln_B_obs, colour = ln_S)) +
  geom_point(alpha = 0.08, size = 0.45, stroke = 0) +
  geom_abline(
    intercept = lne_intercept,
    slope = lne_slope,
    colour = "#b2182b",
    linewidth = 0.95
  ) +
  scale_colour_viridis_c(option = "viridis", name = "ln(S)") +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.35,
    size = 3.1,
    label = sprintf(
      "log_r\u00b2 = %.3f\nlog(B_obs) = %.2f + %.3f \u00d7 log(E)\nSS_res ratio EEoS/OLS = %.1f\u00d7",
      r2_lne, lne_intercept, lne_slope, ss_ratio_lne
    ),
    fill = "white",
    alpha = 0.92
  ) +
  labs(
    x = expression(log(E)),
    y = expression(log(B[obs]) ~ "[g]"),
    title = "(A) Fitted ln(E) OLS ceiling",
    subtitle = "Absorbs the scale offset that unfitted Test 1 penalises"
  ) +
  theme_h1 +
  theme(legend.position = "none")

# (B) Harte Fig. 2 productivity ratio -----------------------------------------
fig2_xlim <- c(0, 60)
fig2_ylim <- c(0, 60)
fig2_full_max <- max(
  c(haul_f2$fig2_predicted_ratio, haul_f2$fig2_observed_ratio),
  na.rm = TRUE
) * 1.05
n_fig2_offframe <- sum(
  haul_f2$fig2_predicted_ratio > fig2_xlim[2] |
    haul_f2$fig2_observed_ratio > fig2_ylim[2],
  na.rm = TRUE
)
haul_lev <- haul_f2 %>% filter(is_leverage)

fig2_base <- ggplot(
  haul_f2,
  aes(x = fig2_predicted_ratio, y = fig2_observed_ratio, colour = ln_S)
) +
  geom_point(alpha = 0.10, size = 0.45, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, colour = "grey25", linewidth = 0.85) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_colour_viridis_c(option = "viridis", name = "ln(S)") +
  theme_h1

p_fig2_main <- fig2_base +
  coord_cartesian(xlim = fig2_xlim, ylim = fig2_ylim, expand = FALSE) +
  annotate(
    "label",
    x = fig2_xlim[1] + diff(fig2_xlim) * 0.04,
    y = fig2_ylim[2] - diff(fig2_ylim) * 0.04,
    hjust = 0,
    vjust = 1,
    size = 3.1,
    label = sprintf(
      "Pearson R\u00b2 = %.3f (trimmed, N=%d)\nPearson R\u00b2 = %.3f (all hauls)\nn trimmed = %s",
      r2_ratio_trim, n_fig2_excl, r2_ratio_all,
      format(n_trimmed, big.mark = ",")
    ),
    fill = "white",
    alpha = 0.92
  ) +
  annotate(
    "label",
    x = fig2_xlim[2] - diff(fig2_xlim) * 0.02,
    y = fig2_ylim[1] + diff(fig2_ylim) * 0.04,
    hjust = 1,
    vjust = 0,
    size = 2.7,
    label = sprintf(
      "%d haul(s) outside panel\n(red = leverage; see inset)",
      n_fig2_offframe
    ),
    fill = "white",
    alpha = 0.88
  ) +
  labs(
    x = expression(Predicted ~ ratio ~ E/B^{3/4} ~ (sqrt ~ scale)),
    y = expression(Observed ~ ratio ~ E/B^{3/4} ~ (sqrt ~ scale)),
    title = "(B) Harte Fig. 2 — productivity ratio",
    subtitle = "Relative structure holds; pattern associated with disturbance"
  ) +
  theme(legend.position = "none")

p_fig2_inset <- fig2_base +
  geom_point(
    data = haul_lev,
    aes(x = fig2_predicted_ratio, y = fig2_observed_ratio),
    colour = "#b2182b",
    size = 2.0,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0, fig2_full_max), ylim = c(0, fig2_full_max)) +
  labs(x = NULL, y = NULL, title = "Full range", subtitle = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6),
    plot.background = element_rect(fill = "white", colour = "grey55", linewidth = 0.4),
    panel.background = element_rect(fill = "grey98", colour = NA),
    plot.margin = margin(2, 4, 2, 2)
  )

p_fig2_panel <- ggdraw(p_fig2_main) +
  draw_plot(p_fig2_inset, x = 0.52, y = 0.08, width = 0.45, height = 0.40)

legend_lnS_f2 <- get_legend(
  fig2_base +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

panels_f2 <- plot_grid(
  p_lne,
  p_fig2_panel,
  nrow = 1,
  rel_widths = c(1, 1.05),
  align = "h",
  axis = "tb"
)

p_body_f2 <- plot_grid(
  panels_f2,
  legend_lnS_f2,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

title_f2 <- ggdraw() +
  draw_label(
    "Relative productivity structure holds while absolute biomass prediction fails",
    fontface = "bold",
    size = 13,
    x = 0,
    hjust = 0,
    y = 0.72
  ) +
  draw_label(
    sprintf(
      "NS-IBTS Q1 1985\u20132015 · n = %s (trimmed n = %s) · offset is not wholesale theory failure",
      format(n_hauls, big.mark = ","),
      format(n_trimmed, big.mark = ",")
    ),
    size = 9.5,
    colour = "grey30",
    x = 0,
    hjust = 0,
    y = 0.28
  ) +
  theme(plot.margin = margin(4, 10, 2, 10))

caption_f2 <- ggdraw() +
  draw_label(
    paste0(
      "(A) Fitted log(B_obs) vs log(E) OLS with fitted line. ",
      "(B) Harte Fig. 2 productivity-ratio predicted vs observed (E/B^0.75), ",
      "highlighting two leverage hauls. ",
      "Tests: OLS; Pearson R\u00b2 all vs trimmed. ",
      "Sources: h1_lne_coefficients.csv, h1_harte_baseline_metrics.csv, ",
      "h1_fig2_leverage_hauls.csv."
    ),
    size = 7.8,
    colour = "grey40",
    x = 0,
    hjust = 0,
    y = 0.5,
    lineheight = 1.15
  ) +
  theme(plot.margin = margin(2, 10, 4, 10))

p_fig2_out <- plot_grid(
  title_f2,
  p_body_f2,
  caption_f2,
  ncol = 1,
  rel_heights = c(0.11, 1, 0.10)
)

ggsave(path_fig2, p_fig2_out, width = 12.5, height = 6.2, dpi = 180, bg = "white")
ok_copy2 <- file.copy(path_fig2, path_fig2_display, overwrite = TRUE)
if (!isTRUE(ok_copy2)) {
  warning("Failed to copy Figure 2 to display_discussion/figures/")
}

logmsg("Wrote ", path_fig2)
logmsg("Copied to ", path_fig2_display)
logmsg(sprintf(
  "Fig2 metrics: ln(E) OLS log_r2=%.3f (%.2f + %.3f*log E) | SS EEoS/OLS=%.2f | Pearson all=%.3f trimmed=%.3f (N=%d)",
  r2_lne, lne_intercept, lne_slope, ss_ratio_lne,
  r2_ratio_all, r2_ratio_trim, n_fig2_excl
))

writeLines(
  c(
    "# H1 presentation figures — run log",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Figures",
    "",
    paste0("- Fig 1 (Test 1): `", path_fig, "`"),
    paste0("  - Display copy: `", path_fig_display, "`"),
    paste0("- Fig 2 (relative structure): `", path_fig2, "`"),
    paste0("  - Display copy: `", path_fig2_display, "`"),
    "",
    "## Fig 1 metrics",
    "",
    sprintf("- n hauls: %s", format(n_hauls, big.mark = ",")),
    sprintf("- EEoS log_r2: %.6f", r2_eeos),
    sprintf("- EEoS cor2: %.6f", cor2_eeos),
    sprintf("- Productivity 1:1 (E × m_min) log_r2: %.6f", r2_prod1),
    sprintf("- Productivity 1:1 cor2: %.6f", cor2_prod1),
    sprintf("- SS_res ratio (EEoS / prod1): %.4f", ss_ratio),
    sprintf("- Median B_pred/B_obs: %.4f", med_ratio),
    "",
    "## Fig 2 metrics",
    "",
    sprintf("- ln(E) OLS log_r2: %.6f", r2_lne),
    sprintf("- Equation: log(B_obs) = %.4f + %.6f × log(E)", lne_intercept, lne_slope),
    sprintf("- SS_res ratio (EEoS / ln(E) OLS): %.4f", ss_ratio_lne),
    sprintf("- Pearson R² (all): %.6f", r2_ratio_all),
    sprintf("- Pearson R² (trimmed, N=%d): %.6f", n_fig2_excl, r2_ratio_trim),
    sprintf("- Trimmed n: %s", format(n_trimmed, big.mark = ",")),
    "",
    "## Console",
    "",
    paste0("    ", run_log)
  ),
  path_run_log
)
logmsg("Wrote ", path_run_log)
