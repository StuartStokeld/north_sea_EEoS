# H1 model comparison: Harte baseline (Tier 1) + ln(E) OLS (Tier 2)
# Run after build_eeos_predictions.R (and optionally run_h1_harte_baseline.R for figures).

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h1_lne_reference.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
source(file.path(r_dir, "h1_harte_baseline_helpers.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "h1_lne_helpers.R"))

path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_haul_out <- file.path(project_root, "outputs", "haul_h1_benchmarks.rds")
path_coef <- file.path(project_root, "outputs", "h1_lne_coefficients.csv")
path_compare_csv <- file.path(project_root, "outputs", "h1_model_comparison.csv")
path_compare_rds <- file.path(project_root, "outputs", "h1_model_comparison.rds")
path_fig_dir <- file.path(project_root, "outputs", "figures")

dir.create(dirname(path_haul_out), recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(path_preds))

haul <- readRDS(path_preds)
haul <- augment_haul_harte_baseline(haul, project_root)
aug <- augment_haul_ln_e(haul)
haul <- aug$haul
fit_info <- aug$fit

comparison <- compare_all_h1_models(haul)

# Fig 2 leverage diagnostic (Test 2 only — does not alter haul table / other tests)
fig2_lev <- fig2_leverage_diagnostic(haul, n_exclude = 2L)
path_fig2_lev <- file.path(project_root, "outputs", "h1_fig2_leverage_hauls.csv")
if (nrow(fig2_lev$excluded_hauls) > 0L) {
  write.csv(fig2_lev$excluded_hauls, path_fig2_lev, row.names = FALSE)
}

ss_eeos <- comparison$ss_res_log[comparison$model_id == "eeos_biomass"]
ss_prod1 <- comparison$ss_res_log[comparison$model_id == "productivity_1to1"]
ss_lne <- comparison$ss_res_log[comparison$model_id == "ln_e_ols"]
ss_ratio <- comparison$ss_res_log[comparison$model_id == "productivity_ratio"]

comparison$ss_ratio_vs_eeos <- comparison$ss_res_log / ss_eeos
comparison$harte_criterion_vs_prod1 <- comparison$ss_res_log < 0.5 * ss_prod1
comparison$harte_criterion_vs_lne <- comparison$ss_res_log < 0.5 * ss_lne
comparison$harte_criterion_vs_prod1[comparison$model_id == "eeos_biomass"] <- NA
comparison$harte_criterion_vs_lne[comparison$model_id == "eeos_biomass"] <- NA

coef_df <- as.data.frame(fit_info$coef)
coef_df$term <- rownames(coef_df)
rownames(coef_df) <- NULL
coef_df <- coef_df[, c("term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]

saveRDS(
  list(
    fit = fit_info,
    coefficients = coef_df,
    comparison = comparison,
    fig2_leverage = fig2_lev[c(
      "fig2_r2_pearson_all", "fig2_r2_pearson_trimmed",
      "fig2_r2_pearson_n_excluded"
    )],
    ss_res_eeos = ss_eeos,
    ss_res_productivity_1to1 = ss_prod1,
    ss_res_lne = ss_lne,
    ss_ratio_eeos_over_productivity_1to1 = ss_eeos / ss_prod1,
    ss_ratio_eeos_over_lne = ss_eeos / ss_lne,
    harte_criterion_met_vs_prod1 = ss_eeos < 0.5 * ss_prod1,
    harte_criterion_met_vs_lne = ss_eeos < 0.5 * ss_lne,
    n_hauls = nrow(haul)
  ),
  path_compare_rds
)

write.csv(coef_df, path_coef, row.names = FALSE)
write.csv(comparison, path_compare_csv, row.names = FALSE)
saveRDS(haul, path_haul_out)

fig_path <- file.path(path_fig_dir, "lnE_reference_logB_vs_logE.png")
png(fig_path, width = 900, height = 600, res = 120)
plot(
  log(haul$E),
  log(haul$B_obs),
  pch = 16,
  col = adjustcolor("black", 0.15),
  cex = 0.5,
  xlab = expression(log(E) ~ "(normalised)"),
  ylab = expression(log(B[obs]) ~ "(g)"),
  main = "ln(E) OLS (fitted correlative): ln(B_obs) ~ ln(E)"
)
abline(fit_info$fit, col = "red", lwd = 2)
legend(
  "topleft",
  legend = sprintf("R² = %.3f", fit_info$r_squared),
  bty = "n"
)
dev.off()

cat("=== ln(E) OLS (Tier 2, fitted correlative) ===\n")
cat("Hauls:", nrow(haul), "\n")
cat("R² (log scale):", round(fit_info$r_squared, 4), "\n")
cat("Slope (log E):", round(coef_df$Estimate[coef_df$term == "log_E"], 4), "\n")
cat("Intercept:", round(coef_df$Estimate[coef_df$term == "(Intercept)"], 4), "\n\n")

cat("=== Unified H1 model comparison ===\n")
print(
  comparison[, c("model", "fitted", "tier", "log_r2", "cor2", "log_rmse", "median_ratio")],
  row.names = FALSE
)

ratio_row <- comparison[comparison$model_id == "productivity_ratio", , drop = FALSE]
cat("\n=== Fig 2 leverage diagnostic (Test 2 only) ===\n")
cat("R² Pearson all:    ", round(ratio_row$fig2_r2_pearson_all, 3), "\n")
cat("R² Pearson trimmed:", round(ratio_row$fig2_r2_pearson_trimmed, 3),
    "(N =", ratio_row$fig2_r2_pearson_n_excluded, "excluded)\n")
if (nrow(fig2_lev$excluded_hauls) > 0L) {
  cat("Excluded hauls:\n")
  print(fig2_lev$excluded_hauls, row.names = FALSE)
}

cat("\nHarte criterion vs productivity 1:1 (EEoS SS < 0.5 × prod1 SS):",
    if (ss_eeos < 0.5 * ss_prod1) "MET" else "NOT MET", "\n")
cat("SS_res ratio (EEoS / productivity 1:1):", round(ss_eeos / ss_prod1, 3), "\n")
cat("Harte criterion vs ln(E) OLS (EEoS SS < 0.5 × ln(E) SS):",
    if (ss_eeos < 0.5 * ss_lne) "MET" else "NOT MET", "\n")
cat("SS_res ratio (EEoS / ln(E) OLS):", round(ss_eeos / ss_lne, 3), "\n\n")

cat("Saved:\n")
cat(" ", path_haul_out, "\n")
cat(" ", path_coef, "\n")
cat(" ", path_compare_csv, "\n")
cat(" ", path_compare_rds, "\n")
cat(" ", path_fig2_lev, "\n")
cat(" ", fig_path, "\n")
