# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Spec A (HISTORICAL — RE baseline): lagged FP_between under (1 | stat_rec).
#
# SUPERSEDED as the presented Spec A claim. Current Spec A is CAR-scoped:
#   pipeline/run_h2h3_spec_a_car_identifiability.R
#   pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R
# See outputs/spec_a_car_routing_note.md.
#
# Formula (historical):
#   residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2
#              + FP_within * phase_v2 + (1 | stat_rec)  [REML]
#
# Run from repo root:
#   Rscript --vanilla pipeline/build_knn_spatial_weights.R   # if weights absent
#   Rscript --vanilla pipeline/run_h2h3_spec_a_lag_refit.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_spec_a_lag_refit.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_spatial_autocorr_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_spec_a_lag_refit.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Package 'foreign' required for queen listw BLUP Moran comparison.")
}

path_primary_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_knn <- file.path(project_root, "outputs", "knn_listw_k4.rds")
path_lag_rect <- file.path(project_root, "outputs", "fp_between_lag_rectangle.rds")
path_out_model <- file.path(project_root, "outputs", "primary_model_v2_spec_a.rds")
path_out_vif <- file.path(project_root, "outputs", "spec_a_vif_fp_between_lag.csv")
path_out_vif_md <- file.path(project_root, "outputs", "spec_a_vif_fp_between_lag.md")
path_out_blup <- file.path(project_root, "outputs", "spec_a_blup_spatial_diagnostic.csv")
path_out_run_log <- file.path(project_root, "outputs", "spec_a_lag_refit_run_log.md")
path_out_session <- file.path(project_root, "outputs", "spec_a_lag_refit_sessionInfo.txt")
paths <- h2_output_paths(project_root)

stopifnot(file.exists(path_primary_v2))

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# Spec A — lagged FP_between refit — run log")
logmsg("")
logmsg(
  "Adds FP_between_lag * phase_v2 to the current primary model (single-change principle). ",
  "No year term; (1 | stat_rec) unchanged."
)

logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)

# ---------------------------------------------------------------------------
# Load primary v2 + k-NN weights
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

primary_v2 <- readRDS(path_primary_v2)
if (is.null(primary_v2$primary_model_v2) || is.null(primary_v2$data)) {
  stop("primary_model_v2.rds missing primary_model_v2 or data.")
}
dat <- primary_v2$data
fit_v2 <- primary_v2$primary_model_v2
phases <- if (!is.null(primary_v2$phase_v2_labels)) {
  as.character(primary_v2$phase_v2_labels)
} else {
  c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
}

dat$stat_rec <- factor(dat$stat_rec)
dat$phase_v2 <- factor(dat$phase_v2, levels = phases)

logmsg("Loaded primary v2: ", path_primary_v2)
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  nrow(dat), dplyr::n_distinct(dat$stat_rec), min(dat$year), max(dat$year)
))

if (!file.exists(path_knn)) {
  stop(
    "knn_listw_k4.rds not found. Run: Rscript pipeline/build_knn_spatial_weights.R"
  )
}
knn_pack <- readRDS(path_knn)
nb_knn <- knn_pack$nb
logmsg("Loaded k-NN weights: ", path_knn)
logmsg(sprintf("k = %d; panel N = %d", knn_pack$k, knn_pack$panel_n))

# ---------------------------------------------------------------------------
# FP_between_lag (rectangle level → haul broadcast)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## FP_between_lag construction")

fp_rect <- unique(dat[, c("stat_rec", "FP_between"), drop = FALSE])
fp_rect$stat_rec <- normalize_stat_rec(fp_rect$stat_rec)
stopifnot(nrow(fp_rect) == length(unique(dat$stat_rec)))

lag_rect <- build_fp_between_lag_rectangle(fp_rect, nb_knn)
saveRDS(lag_rect, path_lag_rect)
write_csv(lag_rect, sub("\\.rds$", ".csv", path_lag_rect))

logmsg(sprintf(
  "FP_between_lag range [%.4f, %.4f]; cor(FP_between, FP_between_lag) = %.4f",
  min(lag_rect$FP_between_lag), max(lag_rect$FP_between_lag),
  stats::cor(lag_rect$FP_between, lag_rect$FP_between_lag)
))
logmsg("Saved: ", path_lag_rect)

dat <- join_fp_between_lag_hauls(dat, lag_rect)

# ---------------------------------------------------------------------------
# VIF diagnostic
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## VIF (FP_between vs FP_between_lag, rectangle level)")

vif_tab <- vif_fp_between_lag(lag_rect)
write_csv(vif_tab, path_out_vif)

vif_lines <- c(
  "# Spec A — VIF diagnostic (FP_between vs FP_between_lag)",
  "",
  "Rectangle-level check (n = 158). High VIF expected by construction ",
  "(spatially clustered fishing pressure); informative, not disqualifying.",
  "",
  sprintf("- Pearson r = **%.4f**", vif_tab$cor[1]),
  sprintf("- r² = **%.4f**", vif_tab$r_squared[1]),
  sprintf("- VIF (either term) = **%.3f**", vif_tab$vif[1]),
  "",
  "| Term | VIF |",
  "|------|-----|",
  sprintf("| FP_between | %.3f |", vif_tab$vif[1]),
  sprintf("| FP_between_lag | %.3f |", vif_tab$vif[2]),
  ""
)
writeLines(vif_lines, path_out_vif_md)
logmsg(sprintf("VIF = %.3f (r = %.4f)", vif_tab$vif[1], vif_tab$cor[1]))
logmsg("Saved: ", path_out_vif)
logmsg("Saved: ", path_out_vif_md)

# ---------------------------------------------------------------------------
# Refit Spec A
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model refit (Spec A)")

formula_spec_a <- residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +
  FP_within * phase_v2 + (1 | stat_rec)
logmsg("Formula: ", paste(deparse(formula_spec_a), collapse = " "), " [REML]")

time_fit <- system.time({
  fit_spec_a <- glmmTMB(formula_spec_a, data = dat, REML = TRUE)
})
logmsg(sprintf("Refit elapsed: %.2f sec.", time_fit[["elapsed"]]))
logmsg(sprintf("Fitted formula: %s", paste(deparse(formula(fit_spec_a)), collapse = " ")))

# ---------------------------------------------------------------------------
# BLUP Moran / Geary (queen listw — same baseline as primary v2)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## BLUP spatial diagnostic (queen listw; baseline Moran I ≈ 0.552)")

panel <- readRDS(paths$panel)
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
nb_queen <- nb_to_listw_W(
  build_queen_nb_from_ices_dbf(panel, project_root),
  zero.policy = TRUE
)
listw_queen <- nb_queen$listw

blups_v2 <- extract_primary_blups_from_fit(fit_v2, panel$stat_rec)
blups_spec_a <- extract_primary_blups_from_fit(fit_spec_a, panel$stat_rec)

test_v2 <- run_spatial_tests(blups_v2$blup, listw_queen, "primary_v2_BLUPs")
test_spec_a <- run_spatial_tests(blups_spec_a$blup, listw_queen, "spec_a_BLUPs")

blup_diag <- data.frame(
  model = c("primary_v2", "spec_a"),
  variable = "(1|stat_rec) intercepts",
  morans_i = c(test_v2$row$morans_i, test_spec_a$row$morans_i),
  geary_c = c(test_v2$row$geary_c, test_spec_a$row$geary_c),
  p_moran = c(test_v2$row$p_moran, test_spec_a$row$p_moran),
  p_geary = c(test_v2$row$p_geary, test_spec_a$row$p_geary),
  n = c(test_v2$row$n, test_spec_a$row$n),
  delta_morans_i_vs_v2 = c(0, test_spec_a$row$morans_i - test_v2$row$morans_i),
  note = c(
    "baseline primary_model_v2 (archived Moran I ≈ 0.552)",
    "Spec A with FP_between_lag * phase_v2"
  ),
  stringsAsFactors = FALSE
)
write_csv(blup_diag, path_out_blup)

logmsg(sprintf(
  "primary_v2 BLUP Moran I = %.6f; Spec A BLUP Moran I = %.6f (Δ = %+.6f)",
  test_v2$row$morans_i, test_spec_a$row$morans_i,
  test_spec_a$row$morans_i - test_v2$row$morans_i
))
logmsg("Saved: ", path_out_blup)

# ---------------------------------------------------------------------------
# Save artifact
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")

saveRDS(
  list(
    primary_model_spec_a = fit_spec_a,
    data = dat,
    formula_spec_a = formula_spec_a,
    phase_v2_labels = phases,
    fp_between_lag_rectangle = lag_rect,
    knn_weights_path = path_knn,
    vif = vif_tab,
    blup_diagnostic = blup_diag,
    primary_v2_source = path_primary_v2,
    comparison_note = paste0(
      "Single change vs primary_model_v2: added FP_between_lag * phase_v2 only."
    )
  ),
  path_out_model
)

logmsg("Saved: ", path_out_model)
logmsg("- ", path_lag_rect)
logmsg("- ", path_out_vif)
logmsg("- ", path_out_vif_md)
logmsg("- ", path_out_blup)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Spec A lag refit complete. ===\n")
