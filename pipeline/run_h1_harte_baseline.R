# Harte et al. (2022) unfitted baseline (Tier 1) — run after build_eeos_predictions.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h1_harte_baseline.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
source(file.path(r_dir, "h1_harte_baseline_helpers.R"))
project_root <- get_project_root_from(script_dir)

path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_metrics_csv <- file.path(project_root, "outputs", "h1_harte_baseline_metrics.csv")
path_out_rds <- file.path(project_root, "outputs", "h1_harte_baseline.rds")
path_fig_dir <- file.path(project_root, "outputs", "figures")

stopifnot(file.exists(path_preds))

haul <- readRDS(path_preds)
haul <- augment_haul_harte_baseline(haul, project_root)
metrics <- harte_baseline_metrics(haul)
fig2_lev <- fig2_leverage_diagnostic(haul, n_exclude = 2L)
path_fig2_lev <- file.path(project_root, "outputs", "h1_fig2_leverage_hauls.csv")

dir.create(dirname(path_metrics_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, path_metrics_csv, row.names = FALSE)
if (nrow(fig2_lev$excluded_hauls) > 0L) {
  write.csv(fig2_lev$excluded_hauls, path_fig2_lev, row.names = FALSE)
}
saveRDS(
  list(
    metrics = metrics,
    n_hauls = nrow(haul),
    fig2_leverage = fig2_lev[c(
      "fig2_r2_pearson_all", "fig2_r2_pearson_trimmed",
      "fig2_r2_pearson_n_excluded"
    )],
    productivity_1to1_note = paste(
      "Productivity 1:1 (headline) uses E_calibrated = E × m_min — the same",
      "normalised E passed to biomass() and the same independently-derived m_min",
      "used for B_pred = B_pred_norm × m_min. Unfitted log-log identity vs B_obs.",
      "E_raw-based metrics are retained as model_id productivity_1to1_uncalibrated",
      "(diagnostic only; not the Test 1 headline baseline)."
    )
  ),
  path_out_rds
)
save_harte_baseline_figures(haul, path_fig_dir)

cat("=== Harte baseline (Tier 1, unfitted) ===\n")
cat("Hauls:", nrow(haul), "\n\n")
print(metrics[, c("model", "log_r2", "cor2", "log_rmse", "median_ratio")],
      row.names = FALSE)
ratio_row <- metrics[metrics$model_id == "productivity_ratio", , drop = FALSE]
cat("\nFig 2 R² all / trimmed (N=",
    ratio_row$fig2_r2_pearson_n_excluded, "): ",
    round(ratio_row$fig2_r2_pearson_all, 3), " / ",
    round(ratio_row$fig2_r2_pearson_trimmed, 3), "\n", sep = "")
cat("\nSaved:\n")
cat(" ", path_metrics_csv, "\n")
cat(" ", path_out_rds, "\n")
cat(" ", path_fig2_lev, "\n")
cat(" ", file.path(path_fig_dir, "harte_fig1_logB_pred_vs_obs.png"), "\n")
cat(" ", file.path(path_fig_dir, "harte_fig2_productivity_ratio.png"), "\n")
cat(" ", file.path(path_fig_dir, "productivity_1to1_logB_vs_logE_calibrated.png"), "\n")
cat(" ", file.path(path_fig_dir, "productivity_1to1_uncalibrated_logB_vs_logEraw.png"), "\n")
