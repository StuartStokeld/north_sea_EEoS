# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Step 0 ROBUSTNESS CHECK — is the D-vs-fishing-pressure cross-sectional
# correlation (r = 0.183, R2 = 0.033, p = 0.020, n = 161) a stable signal or
# an artefact of a few high-leverage rectangles / biased missingness?
#
# ROBUSTNESS CHECK, NOT A RE-RUN OF STEP 0: does NOT change the Step 0
# pipeline, thresholds, or outputs (h2_dominance_diagnostic_helpers.R,
# run_h2_dominance_fishing_pressure_diagnostic.R, step0_*.csv/.rds). Reads
# only the already-built step0_rectangle_panel.csv.
#
# Reports findings uninterpreted, as with Step 0 — does NOT decide whether D
# should become a required H2 control. That is Stuart's call.
#
# Primary focus is D per the briefing (the significant Step 0 result); D-only
# checks match the spec exactly. size_CV rows are also computed throughout
# as a cheap, clearly-labelled supplementary extension for consistency (not
# requested, but not requiring extra data / assumptions beyond D).

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2_dominance_robustness_check.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_dominance_diagnostic_helpers.R"))
source(file.path(script_dir, "R", "h2_dominance_robustness_helpers.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_rect_panel <- file.path(project_root, "outputs", "step0_rectangle_panel.csv")
path_exclusions <- file.path(project_root, "outputs", "step0_exclusion_counts.csv")

path_out_md <- file.path(project_root, "outputs", "step0_robustness_check.md")
path_out_pearson_spearman <- file.path(project_root, "outputs", "step0_robustness_pearson_spearman.csv")
path_out_cooks <- file.path(project_root, "outputs", "step0_robustness_cooks_distance.csv")
path_out_refit <- file.path(project_root, "outputs", "step0_robustness_refit_after_cooks_removal.csv")
path_out_loess <- file.path(project_root, "outputs", "step0_robustness_linear_vs_loess.csv")
path_out_missingness <- file.path(project_root, "outputs", "step0_robustness_missingness.csv")
path_out_fig_cooks <- file.path(project_root, "outputs", "figures", "step0_robustness_cooks_distance.png")
path_out_fig_overlay <- file.path(project_root, "outputs", "figures", "step0_robustness_linear_vs_loess_overlay.png")

stopifnot(file.exists(path_rect_panel), file.exists(path_exclusions))
dir.create(dirname(path_out_fig_cooks), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load Step 0 rectangle panel (187 rows: 161 with paired fishing hours, 26
# dropped for missing Couce coverage — verified against step0_exclusion_counts.csv)
# ---------------------------------------------------------------------------
rect_panel <- read.csv(path_rect_panel, stringsAsFactors = FALSE)
retained <- rect_panel %>% filter(!is.na(mean_annual_hours_total))
dropped <- rect_panel %>% filter(is.na(mean_annual_hours_total))

cat("=== Step 0 robustness check: D vs fishing-pressure correlation ===\n")
cat("Rectangle panel loaded:", nrow(rect_panel), "total\n")
cat("  Retained (paired D + fishing hours, used for the original Step 0 correlation):", nrow(retained), "\n")
cat("  Dropped (no Couce fishing-hours coverage):", nrow(dropped), "\n\n")

# ===========================================================================
# Check 1: Spearman alongside Pearson
# ===========================================================================
cat("--- Check 1: Pearson vs Spearman ---\n")

corr_D <- pearson_spearman_comparison(retained, "mean_D")
corr_cv <- pearson_spearman_comparison(retained, "mean_size_CV") # supplementary, not requested

print(corr_D)
cat("(supplementary, not requested — size_CV for consistency)\n")
print(corr_cv)

write_csv(bind_rows(corr_D, corr_cv), path_out_pearson_spearman)
cat("Saved", path_out_pearson_spearman, "\n\n")

# ===========================================================================
# Check 2: leverage / influence diagnostics (Cook's distance)
# ===========================================================================
cat("--- Check 2: Cook's distance / leverage ---\n")

cooks_D <- cooks_distance_table(retained, "mean_D")
top5_D <- head(cooks_D, 5L)
cat("Top 5 rectangles by Cook's distance (mean_D ~ fishing_hours):\n")
print(top5_D)

refit1_D <- refit_excluding_top_cooks(retained, "mean_D", n_remove = 1L)
refit3_D <- refit_excluding_top_cooks(retained, "mean_D", n_remove = 3L)
cat(sprintf(
  "\nRefit after removing top 1 by Cook's D (%s): r=%.3f, R2=%.3f, p=%.3g, n=%d\n",
  refit1_D$removed_stat_rec, refit1_D$correlation, refit1_D$r_squared, refit1_D$p_value, refit1_D$n
))
cat(sprintf(
  "Refit after removing top 3 by Cook's D (%s): r=%.3f, R2=%.3f, p=%.3g, n=%d\n",
  refit3_D$removed_stat_rec, refit3_D$correlation, refit3_D$r_squared, refit3_D$p_value, refit3_D$n
))

write_csv(cooks_D, path_out_cooks)
refit_table <- bind_rows(
  cross_sectional_association(retained, "mean_D") %>% mutate(n_removed = 0L, removed_stat_rec = NA_character_),
  refit1_D,
  refit3_D
) %>%
  mutate(refit_label = c("original (n_removed=0)", "top1_removed", "top3_removed")) %>%
  relocate(refit_label)
write_csv(refit_table, path_out_refit)
cat("Saved", path_out_cooks, "and", path_out_refit, "\n\n")

p_cooks <- ggplot(
  cooks_D %>% mutate(rank = row_number()),
  aes(x = rank, y = cooks_distance)
) +
  geom_point(colour = "grey30", size = 1.6) +
  geom_segment(aes(x = rank, xend = rank, y = 0, yend = cooks_distance), colour = "grey60", linewidth = 0.3) +
  geom_hline(
    aes(yintercept = cooks_threshold_4_over_n[1]),
    colour = "#d6604d", linetype = "dashed", linewidth = 0.7
  ) +
  labs(
    x = "Rectangle rank (by Cook's distance, descending)",
    y = "Cook's distance",
    title = "Step 0 robustness check: Cook's distance — mean_D ~ fishing_hours",
    subtitle = sprintf(
      "n = %d rectangles; dashed line = conventional threshold 4/n = %.4f",
      nrow(cooks_D), cooks_D$cooks_threshold_4_over_n[1]
    )
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_cooks, p_cooks, width = 8, height = 5.5, dpi = 120)
cat("Saved", path_out_fig_cooks, "\n\n")

# ===========================================================================
# Check 3: linear vs loess fit comparison
# ===========================================================================
cat("--- Check 3: linear vs loess fit comparison ---\n")

loess_dev_D <- linear_vs_loess_deviation(retained, "mean_D")
cat(sprintf(
  "mean_D ~ fishing_hours: max |linear - loess| = %.4f at fishing_hours = %.0f (linear=%.4f, loess=%.4f)\n",
  loess_dev_D$max_abs_diff, loess_dev_D$at_fishing_hours,
  loess_dev_D$linear_fitted_at_max, loess_dev_D$loess_fitted_at_max
))

loess_dev_cv <- linear_vs_loess_deviation(retained, "mean_size_CV") # supplementary
write_csv(bind_rows(loess_dev_D, loess_dev_cv), path_out_loess)
cat("Saved", path_out_loess, "(size_CV row is supplementary, not requested)\n\n")

fit_lm_D <- lm(mean_D ~ mean_annual_hours_total, data = retained)
fit_loess_D <- loess(mean_D ~ mean_annual_hours_total, data = retained, span = 0.75)
grid_x <- seq(min(retained$mean_annual_hours_total), max(retained$mean_annual_hours_total), length.out = 200L)
overlay_df <- bind_rows(
  tibble(mean_annual_hours_total = grid_x, fitted = predict(fit_lm_D, newdata = data.frame(mean_annual_hours_total = grid_x)), fit = "linear"),
  tibble(mean_annual_hours_total = grid_x, fitted = predict(fit_loess_D, newdata = data.frame(mean_annual_hours_total = grid_x)), fit = "loess")
)

p_overlay <- ggplot(retained, aes(x = mean_annual_hours_total, y = mean_D)) +
  geom_point(alpha = 0.5, size = 1.6, colour = "grey40") +
  geom_line(
    data = overlay_df, aes(x = mean_annual_hours_total, y = fitted, colour = fit, linetype = fit),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = c(linear = "#d6604d", loess = "#1a1a1a")) +
  scale_linetype_manual(values = c(linear = "solid", loess = "dashed")) +
  labs(
    x = "Mean annual fishing hours (Couce et al. 2020)",
    y = "Rectangle-mean D (Berger-Parker dominance)",
    colour = "Fit", linetype = "Fit",
    title = "Step 0 robustness check: linear vs loess fit — mean_D ~ fishing_hours",
    subtitle = sprintf(
      "n = %d rectangles; max |linear - loess| = %.4f at fishing_hours = %.0f",
      nrow(retained), loess_dev_D$max_abs_diff, loess_dev_D$at_fishing_hours
    )
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_overlay, p_overlay, width = 8, height = 5.5, dpi = 120)
cat("Saved", path_out_fig_overlay, "\n\n")

# ===========================================================================
# Check 4: missingness check on the 26 dropped rectangles
# ===========================================================================
cat("--- Check 4: missingness check (26 dropped vs", nrow(retained), "retained) ---\n")

miss_D <- missingness_comparison(dropped, retained, "mean_D")
miss_cv <- missingness_comparison(dropped, retained, "mean_size_CV")
print(bind_rows(miss_D, miss_cv))

write_csv(bind_rows(miss_D, miss_cv), path_out_missingness)
cat("Saved", path_out_missingness, "\n\n")

# ===========================================================================
# Single summary document
# ===========================================================================
md_table <- function(df, digits = 4L) {
  df <- as.data.frame(df)
  num_cols <- vapply(df, is.numeric, logical(1))
  for (col in names(df)[num_cols]) {
    df[[col]] <- format(round(df[[col]], digits), nsmall = 0, trim = TRUE)
  }
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, body)
}

md <- c(
  "# Step 0 robustness check — D vs fishing-pressure correlation",
  "",
  paste(
    "Generated:", Sys.Date(),
    "— robustness check on the Step 0 cross-sectional D-fishing-pressure",
    "correlation (r = 0.183, R2 = 0.033, p = 0.020, n = 161). Does not modify",
    "the Step 0 pipeline/thresholds/outputs. Reported uninterpreted, as with",
    "Step 0 — interpretation is Stuart's call."
  ),
  "",
  sprintf(
    "Rectangle counts: %d total in step0_rectangle_panel.csv -> %d retained (paired D + fishing hours) / %d dropped (no Couce coverage).",
    nrow(rect_panel), nrow(retained), nrow(dropped)
  ),
  "",
  "## Check 1 — Pearson vs Spearman (D vs fishing_hours)",
  "",
  md_table(corr_D %>% select(-variable)),
  "",
  "Supplementary (not requested) — size_CV vs fishing_hours, for consistency:",
  "",
  md_table(corr_cv %>% select(-variable)),
  "",
  "## Check 2 — Cook's distance / leverage (mean_D ~ fishing_hours)",
  "",
  sprintf("Top 5 rectangles by Cook's distance (n = %d total):", nrow(cooks_D)),
  "",
  md_table(top5_D),
  "",
  "Re-fit correlation after removing the highest-Cook's-D rectangle(s):",
  "",
  md_table(
    refit_table %>%
      select(refit_label, n, correlation, r_squared, slope, p_value, n_removed, removed_stat_rec)
  ),
  "",
  "## Check 3 — linear vs loess fit comparison (mean_D ~ fishing_hours)",
  "",
  md_table(loess_dev_D %>% select(-variable)),
  "",
  sprintf(
    "![Cook's distance](figures/%s)",
    basename(path_out_fig_cooks)
  ),
  "",
  sprintf(
    "![Linear vs loess overlay](figures/%s)",
    basename(path_out_fig_overlay)
  ),
  "",
  "## Check 4 — missingness check (dropped vs retained rectangles)",
  "",
  sprintf(
    "%d rectangles dropped for no Couce fishing-hours coverage vs %d retained rectangles, compared on rectangle-mean D and rectangle-mean size_CV (same haul-level source / aggregation as Step 0):",
    nrow(dropped), nrow(retained)
  ),
  "",
  md_table(bind_rows(miss_D, miss_cv)),
  "",
  "## Notes on scope",
  "",
  "- Checks 1-3 are D-only per the briefing (the significant Step 0 result under review); size_CV rows in Checks 1 and 3 are a supplementary, clearly-labelled extension using the same already-built panel, not requested in the brief.",
  "- Check 4 reports both D and size_CV per the brief.",
  "- No conclusion is drawn here about whether D should become a required H2 control.",
  "",
  "*Outputs: `outputs/step0_robustness_pearson_spearman.csv`, `outputs/step0_robustness_cooks_distance.csv`,",
  "`outputs/step0_robustness_refit_after_cooks_removal.csv`, `outputs/step0_robustness_linear_vs_loess.csv`,",
  "`outputs/step0_robustness_missingness.csv`, `outputs/figures/step0_robustness_cooks_distance.png`,",
  "`outputs/figures/step0_robustness_linear_vs_loess_overlay.png`.*"
)

writeLines(md, path_out_md)
cat("Saved", path_out_md, "\n")

cat("\n=== Step 0 robustness check complete — findings uninterpreted; no conclusion drawn. ===\n")
