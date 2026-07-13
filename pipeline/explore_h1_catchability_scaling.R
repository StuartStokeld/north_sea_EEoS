# H1 catchability scaling exploration (no correction applied)
# Documents systematic B_pred/B_obs offset and whether a simple scaling would be stable.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/explore_h1_catchability_scaling.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
project_root <- get_project_root_from(script_dir)

path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_out_csv <- file.path(project_root, "outputs", "h1_catchability_scaling_summary.csv")
path_out_quartile <- file.path(project_root, "outputs", "h1_catchability_by_quartile.csv")
path_out_fig <- file.path(project_root, "outputs", "figures", "h1_catchability_ratio_vs_bobs.png")
path_out_md <- file.path(project_root, "display_discussion", "H1_catchability_scaling_exploration.md")

stopifnot(file.exists(path_preds))
dir.create(dirname(path_out_fig), recursive = TRUE, showWarnings = FALSE)

haul <- readRDS(path_preds) %>%
  mutate(
    ln_B_obs = log(B_obs),
    ln_B_pred = log(B_pred),
    ratio = B_pred / B_obs,
    ln_ratio = ln_B_pred - ln_B_obs,
    b_quartile = ntile(ln_B_obs, 4L)
  )

# --- Constant multiplier (median ratio) ---
median_ratio <- median(haul$ratio, na.rm = TRUE)
mean_ln_ratio <- mean(haul$ln_ratio, na.rm = TRUE)

# --- Power-law: log(B_obs) ~ a + b * log(B_pred) ---
fit_power <- lm(ln_B_obs ~ ln_B_pred, data = haul)
power_intercept <- coef(fit_power)[["(Intercept)"]]
power_slope <- coef(fit_power)[["ln_B_pred"]]

# --- Would a single divisor fix log_r2? (exploratory only) ---
haul_scaled <- haul %>%
  mutate(
    B_pred_div_median = B_pred / median_ratio,
    ln_B_pred_div = log(B_pred_div_median)
  )
r2_raw <- log_r2(haul$ln_B_obs, haul$ln_B_pred)
r2_div <- log_r2(haul_scaled$ln_B_obs, haul_scaled$ln_B_pred_div)

quartile_tbl <- haul %>%
  group_by(b_quartile) %>%
  summarise(
    n = n(),
    ln_B_obs_min = min(ln_B_obs, na.rm = TRUE),
    ln_B_obs_max = max(ln_B_obs, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    mean_ln_ratio = mean(ln_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# Correlation of ratio with state variables (confound check for H2)
cor_tbl <- tibble(
  predictor = c("ln_B_obs", "ln_E_raw", "log_N", "log_S"),
  cor_with_ln_ratio = c(
    cor(haul$ln_B_obs, haul$ln_ratio, use = "complete.obs"),
    cor(log(haul$E_raw), haul$ln_ratio, use = "complete.obs"),
    cor(log(haul$N), haul$ln_ratio, use = "complete.obs"),
    cor(log(haul$S), haul$ln_ratio, use = "complete.obs")
  )
)

summary_tbl <- tibble(
  metric = c(
    "n_hauls",
    "median_B_pred_over_B_obs",
    "mean_log_ratio",
    "power_law_intercept",
    "power_law_slope",
    "log_r2_unadjusted",
    "log_r2_after_median_division",
    "interpretation"
  ),
  value = c(
    nrow(haul),
    median_ratio,
    mean_ln_ratio,
    power_intercept,
    power_slope,
    r2_raw,
    r2_div,
    "exploratory_only_no_correction_applied"
  )
)

write.csv(summary_tbl, path_out_csv, row.names = FALSE)
write.csv(quartile_tbl, path_out_quartile, row.names = FALSE)

p <- ggplot(haul, aes(x = ln_B_obs, y = ln_ratio)) +
  geom_point(alpha = 0.05, size = 0.4, colour = "grey35") +
  geom_hline(yintercept = log(median_ratio), linetype = 2, colour = "#4393c3") +
  geom_smooth(method = "lm", se = FALSE, colour = "#d6604d", linewidth = 0.9) +
  labs(
    x = "log(B_obs) [g]",
    y = "log(B_pred / B_obs)",
    title = "Catchability scaling exploration (not applied)",
    subtitle = sprintf(
      "Median ratio = %.2f×; power-law slope = %.3f; log_r2 if divided = %.3f",
      median_ratio, power_slope, r2_div
    )
  ) +
  theme_minimal(base_size = 11)

ggsave(path_out_fig, p, width = 8, height = 5.5, dpi = 120)

md <- c(
  "# H1 catchability scaling exploration",
  "",
  paste("Generated:", Sys.Date(), "— **exploratory only; no correction applied to the pipeline.**"),
  "",
  "## Question",
  "",
  "EEoS systematically overpredicts catch biomass (median `B_pred/B_obs` ≈ **",
  sprintf("%.2f", median_ratio),
  "×**). Is this a stable multiplicative catchability offset that could be removed with a single scaler?",
  "",
  "## Findings",
  "",
  sprintf("- **Median ratio:** %.2f× (constant multiplicative overprediction if catchability were uniform)", median_ratio),
  sprintf("- **Mean log-ratio:** %.3f (≈ log(%.2f))", mean_ln_ratio, exp(mean_ln_ratio)),
  sprintf("- **Power-law** `log(B_obs) ~ a + b × log(B_pred)`: intercept = %.3f, slope = **%.3f** (1.0 = proportional scaling)", power_intercept, power_slope),
  sprintf("- **log_r2** without correction: **%.3f**", r2_raw),
  sprintf("- **log_r2** if every `B_pred` were divided by the median ratio (exploratory): **%.3f**", r2_div),
  "",
  "### Quartile instability (why a naive fix is risky)",
  "",
  "| Quartile | log(B_obs) range | Median B_pred/B_obs |",
  "|----------|------------------|---------------------|",
  sprintf("| Q1 | %.1f–%.1f | %.2f× |", quartile_tbl$ln_B_obs_min[1], quartile_tbl$ln_B_obs_max[1], quartile_tbl$median_ratio[1]),
  sprintf("| Q2 | %.1f–%.1f | %.2f× |", quartile_tbl$ln_B_obs_min[2], quartile_tbl$ln_B_obs_max[2], quartile_tbl$median_ratio[2]),
  sprintf("| Q3 | %.1f–%.1f | %.2f× |", quartile_tbl$ln_B_obs_min[3], quartile_tbl$ln_B_obs_max[3], quartile_tbl$median_ratio[3]),
  sprintf("| Q4 | %.1f–%.1f | %.2f× |", quartile_tbl$ln_B_obs_min[4], quartile_tbl$ln_B_obs_max[4], quartile_tbl$median_ratio[4]),
  "",
  "Overprediction **increases with observed biomass** (3.1× → 5.5×), so a single catchability multiplier is **not** adequate.",
  "",
  "### Correlation of log-ratio with predictors (H2 confound check)",
  "",
  "| Predictor | cor with log(B_pred/B_obs) |",
  "|-----------|---------------------------:|"
)

for (i in seq_len(nrow(cor_tbl))) {
  md <- c(md, sprintf("| %s | %.3f |", cor_tbl$predictor[i], cor_tbl$cor_with_ln_ratio[i]))
}

md <- c(
  md,
  "",
  "## Recommendation (no implementation)",
  "",
  "1. **Do not** apply a global `B_pred / median_ratio` correction in H1 — it would artificially inflate log_r2 without validating the catchability mechanism.",
  "2. **Do not** use biomass-dependent scaling in H2 without explicit fishing-pressure adjustment — ratio correlates with `log(B_obs)`.",
  "3. **Carry forward** the scale-dependent offset as a design constraint: H2/H3 use **mean absolute residual** magnitude, not calibrated absolute biomass.",
  "",
  sprintf("![Ratio vs log B_obs](../outputs/figures/h1_catchability_ratio_vs_bobs.png)"),
  "",
  "*Outputs: `outputs/h1_catchability_scaling_summary.csv`, `outputs/h1_catchability_by_quartile.csv`.*"
)

writeLines(md, path_out_md)

cat("=== Catchability scaling exploration ===\n\n")
print(summary_tbl)
cat("\nQuartile ratios:\n")
print(quartile_tbl)
cat("\nSaved:\n")
cat(" ", path_out_csv, "\n")
cat(" ", path_out_quartile, "\n")
cat(" ", path_out_fig, "\n")
cat(" ", path_out_md, "\n")
