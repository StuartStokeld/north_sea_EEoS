# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Spec B: lagged-outcome (neighbour mean ln_B_obs) + Spec A+B joint fit.
#
# Spec B (vs primary):
#   residual ~ FP_between * phase_v2 + FP_within * phase_v2 + B_lag_neighbour
#              + (1 | stat_rec)  [REML]
#
# Spec A+B joint:
#   residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2
#              + FP_within * phase_v2 + B_lag_neighbour + (1 | stat_rec)  [REML]
#
# B_lag_neighbour: contemporaneous k-NN (k=4) mean of neighbours' rectangle-year
# mean ln_B_obs; missing neighbours dropped and re-normalised.
#
# Run from repo root:
#   Rscript --vanilla pipeline/run_h2h3_spec_b_lag_refit.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_spec_b_lag_refit.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_spatial_autocorr_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("Package 'glmmTMB' required. Run with: Rscript --vanilla pipeline/run_h2h3_spec_b_lag_refit.R")
}
suppressPackageStartupMessages(library(glmmTMB))

if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Package 'foreign' required for queen listw Moran comparison.")
}

paths <- h2_output_paths(project_root)
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_primary_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_spec_a <- file.path(project_root, "outputs", "primary_model_v2_spec_a.rds")
path_knn <- file.path(project_root, "outputs", "knn_listw_k4.rds")
path_fp_lag <- file.path(project_root, "outputs", "fp_between_lag_rectangle.rds")

path_out_lag_ry <- file.path(project_root, "outputs", "b_lag_neighbour_rectangle_year.rds")
path_out_model_b <- file.path(project_root, "outputs", "primary_model_v2_spec_b.rds")
path_out_model_ab <- file.path(project_root, "outputs", "primary_model_v2_spec_ab.rds")
path_out_vif_b <- file.path(project_root, "outputs", "spec_b_vif.csv")
path_out_vif_b_md <- file.path(project_root, "outputs", "spec_b_vif.md")
path_out_vif_ab <- file.path(project_root, "outputs", "spec_ab_vif.csv")
path_out_vif_ab_md <- file.path(project_root, "outputs", "spec_ab_vif.md")
path_out_spatial <- file.path(project_root, "outputs", "spec_b_spatial_diagnostic.csv")
path_out_coefs <- file.path(project_root, "outputs", "spec_b_coefficient_comparison.csv")
path_out_run_log <- file.path(project_root, "outputs", "spec_b_lag_refit_run_log.md")
path_out_session <- file.path(project_root, "outputs", "spec_b_lag_refit_sessionInfo.txt")
path_out_summary <- file.path(project_root, "outputs", "spec_b_lag_refit_summary.md")

stopifnot(
  file.exists(path_haul),
  file.exists(path_primary_v2),
  file.exists(path_knn),
  file.exists(paths$panel)
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# Spec B / Spec A+B — lagged-outcome refit — run log")
logmsg("")
logmsg(
  "Spec B adds pooled B_lag_neighbour (neighbour mean ln_B_obs) to the primary. ",
  "Spec A+B adds both FP_between_lag * phase_v2 and B_lag_neighbour. ",
  "No year term; (1 | stat_rec) unchanged."
)

logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

panel <- readRDS(paths$panel)
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
knn <- readRDS(path_knn)
nb <- knn$nb
if (!identical(normalize_stat_rec(knn$stat_rec_order), panel$stat_rec)) {
  stop("knn_listw_k4.rds stat_rec_order does not match panel$stat_rec.")
}

primary_v2 <- readRDS(path_primary_v2)
fit_v2 <- primary_v2$primary_model_v2
dat <- primary_v2$data
phases <- if (!is.null(primary_v2$phase_v2_labels)) {
  as.character(primary_v2$phase_v2_labels)
} else {
  c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
}
dat$stat_rec <- factor(dat$stat_rec)
dat$phase_v2 <- factor(dat$phase_v2, levels = phases)

logmsg("Primary v2: ", path_primary_v2)
logmsg("Haul predictions: ", path_haul)
logmsg("k-NN weights: ", path_knn, sprintf(" (k=%d)", knn$k))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  nrow(dat), dplyr::n_distinct(dat$stat_rec), min(dat$year), max(dat$year)
))

# Spec A fit (for comparison); rebuild FP_between_lag join if needed
fit_spec_a <- NULL
if (file.exists(path_spec_a)) {
  spec_a_obj <- readRDS(path_spec_a)
  fit_spec_a <- spec_a_obj$primary_model_spec_a
  logmsg("Loaded Spec A fit: ", path_spec_a)
} else {
  logmsg("WARNING: Spec A artifact missing — Spec A columns in comparison will be NA.")
}

# ---------------------------------------------------------------------------
# Build B_lag_neighbour (rectangle-year)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## B_lag_neighbour construction")

haul <- readRDS(path_haul)
lag_ry <- build_b_lag_neighbour_rectangle_year(
  haul = haul,
  panel_stat_rec = knn$stat_rec_order,
  nb = nb,
  year_min = H2_YEAR_MIN,
  year_max = H2_YEAR_MAX
)

# Audit: analysis rect-years vs neighbour coverage
analysis_keys <- unique(paste(normalize_stat_rec(dat$stat_rec), dat$year, sep = "\r"))
lag_keys <- paste(lag_ry$stat_rec, lag_ry$year, sep = "\r")
analysis_lag <- lag_ry[lag_keys %in% analysis_keys, , drop = FALSE]
n_partial <- sum(analysis_lag$n_neighbours_used < knn$k & analysis_lag$n_neighbours_used > 0L)
n_full <- sum(analysis_lag$n_neighbours_used == knn$k)
n_none <- sum(analysis_lag$n_neighbours_used == 0L)
n_couce_gap_source <- sum(
  is.finite(lag_ry$mean_ln_B_obs_ry) &
    !(lag_keys %in% analysis_keys)
)

logmsg(sprintf(
  "Rectangle-year table: %d rows (158 × %d years); with biomass: %d",
  nrow(lag_ry), H2_YEAR_MAX - H2_YEAR_MIN + 1L,
  sum(is.finite(lag_ry$mean_ln_B_obs_ry))
))
logmsg(sprintf(
  "Analysis rect-years: %d; k-NN full=%d; partial=%d; none=%d",
  nrow(analysis_lag), n_full, n_partial, n_none
))
logmsg(sprintf(
  "Neighbour-source-only rect-years (biomass present, not in analysis data): %d",
  n_couce_gap_source
))
if (n_none > 0L) {
  stop("Unexpected: analysis rect-years with zero available neighbours.")
}

saveRDS(lag_ry, path_out_lag_ry)
write_csv(lag_ry, sub("\\.rds$", ".csv", path_out_lag_ry))
logmsg("Saved: ", path_out_lag_ry)

dat_b <- join_b_lag_neighbour_hauls(dat, lag_ry)
logmsg(sprintf(
  "Joined B_lag_neighbour to hauls: range [%.4f, %.4f]; mean n_neighbours_used=%.3f",
  min(dat_b$B_lag_neighbour), max(dat_b$B_lag_neighbour),
  mean(dat_b$n_neighbours_used)
))

# FP_between_lag for A+B
if (file.exists(path_fp_lag)) {
  fp_lag <- readRDS(path_fp_lag)
} else if (!is.null(fit_spec_a) && !is.null(spec_a_obj$fp_between_lag_rectangle)) {
  fp_lag <- spec_a_obj$fp_between_lag_rectangle
} else {
  fp_rect <- unique(dat_b[, c("stat_rec", "FP_between"), drop = FALSE])
  fp_rect$stat_rec <- normalize_stat_rec(fp_rect$stat_rec)
  fp_lag <- build_fp_between_lag_rectangle(fp_rect, nb)
}
dat_ab <- join_fp_between_lag_hauls(dat_b, fp_lag)

# ---------------------------------------------------------------------------
# VIF
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## VIF diagnostics")

# Use unique rect-year rows for B_lag (time-varying) + rectangle-constant FP terms
ry_vif <- dat_ab %>%
  distinct(stat_rec, year, .keep_all = TRUE) %>%
  as.data.frame()

vif_b <- vif_continuous_terms(
  ry_vif,
  c("FP_between", "FP_within", "B_lag_neighbour")
)
vif_ab <- vif_continuous_terms(
  ry_vif,
  c("FP_between", "FP_within", "FP_between_lag", "B_lag_neighbour")
)
write_csv(vif_b, path_out_vif_b)
write_csv(vif_ab, path_out_vif_ab)

write_vif_md <- function(path, title, vif_tab) {
  lines <- c(
    paste0("# ", title),
    "",
    "Rectangle-year grain (unique `stat_rec` × `year` in analysis data).",
    "VIF = 1/(1−R²) from each term ~ remaining continuous terms.",
    "",
    "| Term | VIF | R² (others) | n |",
    "|------|-----|-------------|---|"
  )
  for (i in seq_len(nrow(vif_tab))) {
    r <- vif_tab[i, ]
    lines <- c(lines, sprintf(
      "| %s | %.3f | %.4f | %d |",
      r$term, r$vif, r$r_squared_others, r$n
    ))
  }
  lines <- c(lines, "")
  writeLines(lines, path)
}
write_vif_md(path_out_vif_b_md, "Spec B — VIF diagnostic", vif_b)
write_vif_md(path_out_vif_ab_md, "Spec A+B — VIF diagnostic", vif_ab)
logmsg(sprintf(
  "Spec B max VIF = %.3f (%s)",
  max(vif_b$vif), vif_b$term[which.max(vif_b$vif)]
))
logmsg(sprintf(
  "Spec A+B max VIF = %.3f (%s)",
  max(vif_ab$vif), vif_ab$term[which.max(vif_ab$vif)]
))
logmsg("Saved: ", path_out_vif_b, " / ", path_out_vif_ab)

# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fits")

formula_b <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 +
  B_lag_neighbour + (1 | stat_rec)
formula_ab <- residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +
  FP_within * phase_v2 + B_lag_neighbour + (1 | stat_rec)

logmsg("Spec B formula: ", paste(deparse(formula_b), collapse = " "))
t_b <- system.time({
  fit_b <- glmmTMB(formula_b, data = dat_b, REML = TRUE)
})
logmsg(sprintf("Spec B elapsed: %.2f sec.", t_b[["elapsed"]]))

logmsg("Spec A+B formula: ", paste(deparse(formula_ab), collapse = " "))
t_ab <- system.time({
  fit_ab <- glmmTMB(formula_ab, data = dat_ab, REML = TRUE)
})
logmsg(sprintf("Spec A+B elapsed: %.2f sec.", t_ab[["elapsed"]]))

b_lag_b <- unname(glmmTMB::fixef(fit_b)$cond[["B_lag_neighbour"]])
b_lag_ab <- unname(glmmTMB::fixef(fit_ab)$cond[["B_lag_neighbour"]])
logmsg(sprintf("B_lag_neighbour coef — Spec B: %+.6f; Spec A+B: %+.6f", b_lag_b, b_lag_ab))

# ---------------------------------------------------------------------------
# Spatial diagnostics (BLUP + residual-level Moran)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Spatial diagnostics (queen listw; BLUP + rectangle-mean residual)")

nb_queen <- nb_to_listw_W(
  build_queen_nb_from_ices_dbf(panel, project_root),
  zero.policy = TRUE
)
listw_queen <- nb_queen$listw

spatial_one <- function(fit, dat_fit, model_id) {
  blups <- extract_primary_blups_from_fit(fit, panel$stat_rec)
  resid_pack <- extract_resid_by_rect_glmmTMB(fit, dat_fit, panel$stat_rec)
  t_blup <- run_spatial_tests(blups$blup, listw_queen, paste0(model_id, "_BLUP"))
  t_resid <- run_spatial_tests(
    resid_pack$resid_by_rect$resid, listw_queen, paste0(model_id, "_resid")
  )
  rbind(
    data.frame(
      model = model_id, level = "BLUP",
      morans_i = t_blup$row$morans_i, geary_c = t_blup$row$geary_c,
      p_moran = t_blup$row$p_moran, p_geary = t_blup$row$p_geary,
      n = t_blup$row$n, stringsAsFactors = FALSE
    ),
    data.frame(
      model = model_id, level = "rectangle_mean_residual",
      morans_i = t_resid$row$morans_i, geary_c = t_resid$row$geary_c,
      p_moran = t_resid$row$p_moran, p_geary = t_resid$row$p_geary,
      n = t_resid$row$n, stringsAsFactors = FALSE
    )
  )
}

spatial_rows <- list(
  spatial_one(fit_v2, dat, "primary_v2")
)
if (!is.null(fit_spec_a)) {
  # Spec A data may already have FP_between_lag; rebuild from primary + lag join
  dat_a <- join_fp_between_lag_hauls(dat, fp_lag)
  dat_a$stat_rec <- factor(dat_a$stat_rec, levels = levels(dat$stat_rec))
  dat_a$phase_v2 <- factor(dat_a$phase_v2, levels = phases)
  spatial_rows <- c(spatial_rows, list(spatial_one(fit_spec_a, dat_a, "spec_a")))
}
spatial_rows <- c(
  spatial_rows,
  list(
    spatial_one(fit_b, dat_b, "spec_b"),
    spatial_one(fit_ab, dat_ab, "spec_ab")
  )
)
spatial_diag <- dplyr::bind_rows(spatial_rows)
write_csv(spatial_diag, path_out_spatial)

logmsg("| Model | Level | Moran's I | Geary's C | p (Moran) |")
logmsg("|---|---|---|---|---|")
for (i in seq_len(nrow(spatial_diag))) {
  r <- spatial_diag[i, ]
  logmsg(sprintf(
    "| %s | %s | %.4f | %.4f | %.3g |",
    r$model, r$level, r$morans_i, r$geary_c, r$p_moran
  ))
}
logmsg("Saved: ", path_out_spatial)

# ---------------------------------------------------------------------------
# Coefficient comparison (H2 / H3 phase slopes)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Coefficient comparison (FP_between / FP_within phase slopes)")

slope_block <- function(fit, model_id) {
  if (is.null(fit)) return(NULL)
  dplyr::bind_rows(
    extract_wb_phase_slopes(
      fit, "FP_between", model_id, "H2_spatial_between", phases = phases
    ),
    extract_wb_phase_slopes(
      fit, "FP_within", model_id, "H3_temporal_within", phases = phases
    )
  )
}

coef_cmp <- dplyr::bind_rows(
  slope_block(fit_v2, "primary_v2"),
  slope_block(fit_spec_a, "spec_a"),
  slope_block(fit_b, "spec_b"),
  slope_block(fit_ab, "spec_ab")
)
write_csv(coef_cmp, path_out_coefs)
logmsg("Saved: ", path_out_coefs)

# Log H2 slope deltas vs primary for Spec B / A+B
for (mid in c("spec_b", "spec_ab")) {
  for (ph in phases) {
    base <- coef_cmp$fp_slope[
      coef_cmp$model_id == "primary_v2" &
        coef_cmp$component == "FP_between" & coef_cmp$phase == ph
    ]
    new <- coef_cmp$fp_slope[
      coef_cmp$model_id == mid &
        coef_cmp$component == "FP_between" & coef_cmp$phase == ph
    ]
    if (length(base) == 1L && length(new) == 1L) {
      logmsg(sprintf(
        "  %s H2 %s: primary=%+.4f → %s=%+.4f (Δ=%+.4f)",
        mid, ph, base, mid, new, new - base
      ))
    }
  }
}

# ---------------------------------------------------------------------------
# Save model artifacts
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")

saveRDS(
  list(
    primary_model_spec_b = fit_b,
    data = dat_b,
    formula_spec_b = formula_b,
    phase_v2_labels = phases,
    b_lag_neighbour_rectangle_year = lag_ry,
    knn_weights_path = path_knn,
    vif = vif_b,
    comparison_note = "Single change vs primary_model_v2: added pooled B_lag_neighbour only."
  ),
  path_out_model_b
)
saveRDS(
  list(
    primary_model_spec_ab = fit_ab,
    data = dat_ab,
    formula_spec_ab = formula_ab,
    phase_v2_labels = phases,
    b_lag_neighbour_rectangle_year = lag_ry,
    fp_between_lag_rectangle = fp_lag,
    knn_weights_path = path_knn,
    vif = vif_ab,
    comparison_note = paste0(
      "Adds FP_between_lag * phase_v2 (Spec A) and pooled B_lag_neighbour (Spec B) ",
      "to primary_model_v2."
    )
  ),
  path_out_model_ab
)
logmsg("Saved: ", path_out_model_b)
logmsg("Saved: ", path_out_model_ab)

# ---------------------------------------------------------------------------
# Summary markdown
# ---------------------------------------------------------------------------
blup_i <- function(mid) {
  spatial_diag$morans_i[spatial_diag$model == mid & spatial_diag$level == "BLUP"]
}
resid_i <- function(mid) {
  spatial_diag$morans_i[
    spatial_diag$model == mid & spatial_diag$level == "rectangle_mean_residual"
  ]
}

summary_lines <- c(
  "# Spec B / Spec A+B — lagged-outcome summary",
  "",
  "## Formulas",
  "",
  sprintf("- Spec B: `%s`", paste(deparse(formula_b), collapse = " ")),
  sprintf("- Spec A+B: `%s`", paste(deparse(formula_ab), collapse = " ")),
  "",
  "## B_lag_neighbour",
  "",
  sprintf("- Analysis rect-years: %d (full k=%d: %d; partial: %d)", nrow(analysis_lag), knn$k, n_full, n_partial),
  sprintf("- Neighbour-source-only rect-years (incl. Couce gaps): %d", n_couce_gap_source),
  sprintf("- Spec B `B_lag_neighbour` coef: **%+.6f**", b_lag_b),
  sprintf("- Spec A+B `B_lag_neighbour` coef: **%+.6f**", b_lag_ab),
  "",
  "## Moran's I (queen listw)",
  "",
  "| Model | BLUP I | Residual I |",
  "|-------|--------|------------|",
  sprintf("| primary_v2 | %.4f | %.4f |", blup_i("primary_v2"), resid_i("primary_v2")),
  if (!is.null(fit_spec_a)) {
    sprintf("| spec_a | %.4f | %.4f |", blup_i("spec_a"), resid_i("spec_a"))
  } else {
    "| spec_a | NA | NA |"
  },
  sprintf("| spec_b | %.4f | %.4f |", blup_i("spec_b"), resid_i("spec_b")),
  sprintf("| spec_ab | %.4f | %.4f |", blup_i("spec_ab"), resid_i("spec_ab")),
  "",
  "Baselines: Finding 1 BLUP I ≈ 0.552; Spec A BLUP I ≈ 0.532.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", path_out_lag_ry),
  sprintf("- `%s`", path_out_model_b),
  sprintf("- `%s`", path_out_model_ab),
  sprintf("- `%s`", path_out_spatial),
  sprintf("- `%s`", path_out_coefs),
  ""
)
writeLines(unlist(summary_lines), path_out_summary)
logmsg("Saved: ", path_out_summary)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Spec B / Spec A+B refit complete. ===\n")
