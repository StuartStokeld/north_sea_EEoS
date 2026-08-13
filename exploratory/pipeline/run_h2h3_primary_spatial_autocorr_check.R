# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Residual spatial autocorrelation check on the primary H2/H3 model
# (`residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)`).
#
# Compares Moran/Geary on primary-model BLUPs and rectangle-collapsed residuals
# against the archived H2c diagnostics in outputs/h2_spatial_diagnostics.csv
# (verbatim; not recomputed), using the same 158-rectangle panel and queen
# contiguity structure as pipeline/run_h2_models.R.
#
# Run: Rscript --vanilla pipeline/run_h2h3_primary_spatial_autocorr_check.R
#
# Dependencies: base R + foreign (for ICES DBF). Does not require a live
# glmmTMB/spdep/sf install — recovers BLUPs/residuals from the saved model
# artifact and rebuilds queen weights from the same ICES SOUTH/WEST grid
# (validated to reproduce archived Moran/Geary statistics to machine precision).

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_primary_spatial_autocorr_check.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_spatial_autocorr_helpers.R"))

if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Package 'foreign' required to read the ICES rectangle DBF.")
}

paths <- h2_output_paths(project_root)
path_archived_diag <- paths$spatial_diagnostics
path_panel <- paths$panel
path_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")

path_out_comparison <- file.path(project_root, "outputs", "spatial_diagnostics_comparison.csv")
path_out_blups <- file.path(project_root, "outputs", "blups_by_rectangle.csv")
path_out_resids <- file.path(project_root, "outputs", "residuals_by_rectangle.csv")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_primary_spatial_autocorr_run_log.md")
path_out_session <- file.path(project_root, "outputs", "h2h3_primary_spatial_autocorr_sessionInfo.txt")

stopifnot(
  file.exists(path_archived_diag),
  file.exists(path_panel),
  file.exists(path_models)
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 primary-model residual spatial autocorrelation check — run log")
logmsg("")
logmsg(
  "Tests whether the primary hierarchical model `(1 | stat_rec)` has absorbed the spatial ",
  "clustering documented in archived H2c (`outputs/h2_spatial_diagnostics.csv`), or whether ",
  "clustering remains in BLUPs and/or rectangle-collapsed model residuals."
)

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(
  "Runtime path: base R + foreign. BLUPs/residuals recovered from saved glmmTMB artifact ",
  "without loading glmmTMB; queen weights from ICES DBF SOUTH/WEST (equivalent to ",
  "poly2nb queen=TRUE on this regular lattice); Moran/Geary via Cliff–Ord randomization ",
  "moments matching h2_global_spatial_tests() / spdep defaults (zero.policy, alternative greater)."
)

# ---------------------------------------------------------------------------
# Step 1b finding — archived residual definition (before any new numbers)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1b — Archived residual definition (checked before producing numbers)")
logmsg(archived_ols_residual_definition())
logmsg("")
logmsg(
  "Primary-model residual extraction: conditional response residual = ",
  "y − (Xβ + b_i), recovered from fit_wb$obj$env$last.par.best. This matches ",
  "`residuals(glmmTMB_fit, type = \"response\")` (includes the random intercept). ",
  "Pearson residuals are not used."
)

# ---------------------------------------------------------------------------
# Load archived diagnostics, panel, primary model
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

archived_rows <- load_archived_spatial_diagnostics(path_archived_diag)
logmsg("Loaded archived diagnostics verbatim from: ", path_archived_diag)
logmsg(sprintf(
  "Archived rows: %s",
  paste(sprintf("%s (N=%d)", archived_rows$variable, archived_rows$n), collapse = "; ")
))

panel <- readRDS(path_panel)
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
if (nrow(panel) != 158L) {
  stop("Expected archived H2 panel N = 158; got ", nrow(panel))
}
logmsg(sprintf("Archived H2 panel: %s (%d rectangles).", path_panel, nrow(panel)))

models <- readRDS(path_models)
if (is.null(models$fit_wb)) {
  stop("h2h3_wb_model_objects.rds missing fit_wb — primary model artifact not found.")
}
fit_wb <- models$fit_wb
dat <- models$data
if (is.null(dat)) stop("h2h3_wb_model_objects.rds missing data.")
logmsg(sprintf(
  "Loaded primary model artifact: %s (fit_wb; formula: %s).",
  path_models, paste(deparse(models$formula_wb), collapse = " ")
))
logmsg("Model was NOT refit — using saved glmmTMB object and its associated analysis data.")
logmsg(sprintf(
  "Analysis data: %d hauls, %d unique rectangles, years %d–%d.",
  nrow(dat), length(unique(dat$stat_rec)), min(dat$year), max(dat$year)
))

# ---------------------------------------------------------------------------
# Spatial weights
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Spatial weights (matched to H2c)")
logmsg(
  "No persisted nb/listw RDS was found in outputs/ from run_h2_models.R. ",
  "Rebuilding queen contiguity from the same ICES shapefile DBF ",
  "(gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.dbf) via SOUTH/WEST grid ",
  "indices (0.5° lat × 1° lon), row-standardised style W, zero.policy=TRUE — ",
  "the regular-lattice equivalent of build_h2_spatial_weights()/poly2nb(queen=TRUE)."
)
nb_pack <- nb_to_listw_W(build_queen_nb_from_ices_dbf(panel, project_root), zero.policy = TRUE)
listw <- nb_pack$listw
logmsg(sprintf(
  "Queen contiguity weights: N=%d rectangles; n_isolated=%d; mean degree=%.3f.",
  length(listw$neighbours), nb_pack$n_isolated, mean(lengths(nb_pack$nb))
))
arch_iso <- utils::read.csv(path_archived_diag, stringsAsFactors = FALSE)$n_isolated[1]
logmsg(sprintf("Archived diagnostics n_isolated=%s (must match rebuild: %d).", arch_iso, nb_pack$n_isolated))
if (as.integer(arch_iso) != nb_pack$n_isolated) {
  stop("n_isolated mismatch vs archived diagnostics — weights not comparable.")
}

# Validation: recomputed archived variables must match CSV statistics
logmsg("")
logmsg("### Weights validation against archived statistics")
ols_fit <- stats::lm(mean_abs_residual ~ mean_annual_hours_total, data = panel)
val_vars <- list(
  mean_abs_residual = panel$mean_abs_residual,
  log_mean_annual_hours_total = panel$log_mean_annual_hours_total,
  ols_primary_abs_residuals = stats::residuals(ols_fit)
)
for (vn in names(val_vars)) {
  m <- moran_test_listw(val_vars[[vn]], listw)
  arch_i <- archived_rows$morans_i[archived_rows$variable == vn]
  diff_i <- abs(m$statistic - arch_i)
  logmsg(sprintf(
    "  %s: recomputed Moran I=%.8f vs archived %.8f (|diff|=%.3g) %s",
    vn, m$statistic, arch_i, diff_i,
    if (diff_i < 1e-10) "OK" else "FLAG"
  ))
  if (diff_i >= 1e-10) {
    stop("Weights validation failed for ", vn, " — aborting before new tests.")
  }
}
logmsg("Weights validated: recomputed Moran I matches archived CSV to machine precision.")

# ---------------------------------------------------------------------------
# Step 1 — Extract BLUPs and residuals
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1 — Extract BLUPs and rectangle-collapsed residuals")

blups <- extract_primary_blups_from_fit(fit_wb, panel$stat_rec)
# attach n_years / n_hauls
counts <- as.data.frame(dat, stringsAsFactors = FALSE)
counts$stat_rec <- normalize_stat_rec(counts$stat_rec)
n_years <- tapply(counts$year, counts$stat_rec, function(y) length(unique(y)))
n_hauls <- tapply(counts$year, counts$stat_rec, length)
blups$n_years <- as.integer(n_years[blups$stat_rec])
blups$n_hauls <- as.integer(n_hauls[blups$stat_rec])
if (anyNA(blups$n_years) || anyNA(blups$n_hauls)) {
  stop("Failed to attach n_years/n_hauls to BLUPs — ID mismatch.")
}
logmsg(sprintf(
  "BLUPs: nrow=%d; setequal(panel IDs)=%s; aligned to panel row order for listw.",
  nrow(blups), setequal(blups$stat_rec, panel$stat_rec)
))

resid_pack <- extract_primary_resid_by_rect_from_fit(fit_wb, dat, panel$stat_rec)
resid_by_rect <- resid_pack$resid_by_rect
logmsg(sprintf(
  "Residuals: nrow(resid_by_rect)=%d; columns resid (signed mean) and abs_resid (mean |residual|).",
  nrow(resid_by_rect)
))
logmsg(
  "FLAG — archived-comparable residual row: ",
  "**Primary model residuals — signed mean** (`rectangle-mean residual`). ",
  "Basis: script/CSV check shows ols_primary_abs_residuals = signed OLS residuals of ",
  "mean_abs_residual ~ hours, not an absolute residual; role-matched leftover after the ",
  "primary model is the signed response residual. Mean-absolute residual is retained as ",
  "secondary (magnitude clustering)."
)

# ---------------------------------------------------------------------------
# Step 2 — Moran / Geary
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 2 — Moran / Geary on primary-model objects")

test_blup <- run_spatial_tests(blups$blup, listw, "(1|stat_rec) intercepts")
test_sgn <- run_spatial_tests(resid_by_rect$resid, listw, "rectangle-mean residual")
test_abs <- run_spatial_tests(resid_by_rect$abs_resid, listw, "rectangle-mean |residual|")

print_test <- function(t) {
  logmsg(sprintf(
    "  %s: Moran I=%.6f (E=%.6f, Var=%.6g, z=%.3f, p=%.6g); Geary C=%.6f (E=%.6f, p=%.6g); N=%d",
    t$label,
    t$row$morans_i, t$row$expected_moran, t$row$variance_moran, t$row$z_moran, t$row$p_moran,
    t$row$geary_c, t$row$expected_geary, t$row$p_geary,
    t$row$n
  ))
}
print_test(test_blup)
print_test(test_sgn)
print_test(test_abs)

# ---------------------------------------------------------------------------
# Step 3 — Comparison table
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 3 — Comparison table")

new_rows <- rbind(
  cbind(stage = "Primary model BLUPs", test_blup$row, archived_comparable_residual = FALSE),
  cbind(stage = "Primary model residuals — signed mean", test_sgn$row, archived_comparable_residual = TRUE),
  cbind(stage = "Primary model residuals — mean absolute", test_abs$row, archived_comparable_residual = FALSE)
)

comparison <- rbind(archived_rows, new_rows)
comp_flag <- as.logical(comparison$archived_comparable_residual)
comparison$comparable_to_ols_primary_abs_residuals <- ifelse(
  comparison$variable == "ols_primary_abs_residuals",
  "Archived reference row (signed OLS residuals of mean_abs_residual ~ hours)",
  ifelse(
    comp_flag %in% TRUE,
    "YES — role-matched leftover after primary model (signed mean response residual); see Step 1b finding",
    ifelse(
      comparison$stage == "Primary model residuals — mean absolute",
      "NO — secondary magnitude check only (not the archived OLS residual definition)",
      ifelse(
        comparison$stage == "Primary model BLUPs",
        "NO — BLUP map, not residual map",
        "Archived raw-variable row (not a residual)"
      )
    )
  )
)
comparison$footnote_archived_comparable <- ifelse(
  comp_flag %in% TRUE,
  "FLAG: this new residual row is the one used as primary input to Decision Rule (b)",
  NA_character_
)

# Column order
comparison <- comparison[, c(
  "stage", "variable", "morans_i", "geary_c", "p_moran", "p_geary", "n",
  "expected_moran", "expected_geary", "variance_moran", "z_moran", "source",
  "comparable_to_ols_primary_abs_residuals", "footnote_archived_comparable"
)]

utils::write.csv(comparison, path_out_comparison, row.names = FALSE)
logmsg("Saved: ", path_out_comparison)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
logmsg("")
logmsg("| Stage | Variable | Moran's I | Geary's C | p (Moran) | p (Geary) | N |")
logmsg("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  logmsg(sprintf(
    "| %s | `%s` | %.3f | %.3f | %s | %s | %d |",
    r$stage, r$variable, r$morans_i, r$geary_c, fmt_p(r$p_moran), fmt_p(r$p_geary), r$n
  ))
}
logmsg("")
logmsg(paste0(
  "Footnote: archived-comparable new residual row = ",
  "`Primary model residuals — signed mean` / `rectangle-mean residual`. ",
  archived_ols_residual_definition()
))

utils::write.csv(blups, path_out_blups, row.names = FALSE)
utils::write.csv(resid_by_rect, path_out_resids, row.names = FALSE)
logmsg("Saved: ", path_out_blups)
logmsg("Saved: ", path_out_resids)
logmsg(
  "Caveat: BLUPs are shrinkage estimates (attenuated vs true rectangle effects), ",
  "especially for sparse rectangles — see n_years / n_hauls in blups_by_rectangle.csv. ",
  "Non-significant Moran on BLUPs = no strong evidence of remaining clustering, ",
  "not proof that clustering is absent."
)

# ---------------------------------------------------------------------------
# Step 4 — Decision rule
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 4 — Decision-rule classification")

decision <- classify_spatial_decision(
  p_blup_moran = test_blup$row$p_moran,
  i_blup = test_blup$row$morans_i,
  p_resid_primary_moran = test_sgn$row$p_moran,
  i_resid_primary = test_sgn$row$morans_i,
  p_resid_secondary_moran = test_abs$row$p_moran,
  i_resid_secondary = test_abs$row$morans_i,
  primary_resid_label = "signed mean residual",
  secondary_resid_label = "mean absolute residual",
  alpha = ALPHA_SPATIAL
)

logmsg(sprintf("Alpha for significance: %.2f (Moran's I one-sided greater, matching archived).", ALPHA_SPATIAL))
logmsg(sprintf(
  "BLUPs: Moran I=%.4f, p=%.6g → %s",
  decision$blup_moran_i, test_blup$row$p_moran,
  if (decision$blup_significant) "SIGNIFICANT" else "non-significant"
))
logmsg(sprintf(
  "Primary residual input (signed mean): Moran I=%.4f, p=%.6g → %s",
  decision$primary_resid_moran_i, test_sgn$row$p_moran,
  if (decision$primary_residual_significant) "SIGNIFICANT" else "non-significant"
))
logmsg(sprintf(
  "Secondary residual (mean absolute): Moran I=%.4f, p=%.6g → %s",
  test_abs$row$morans_i, test_abs$row$p_moran,
  if (decision$secondary_residual_significant) "SIGNIFICANT" else "non-significant"
))
logmsg("")
logmsg(sprintf("**Decision rule %d:**", decision$rule_id))
logmsg(paste0('"', decision$classification, '"'))
logmsg("")
logmsg(decision$disagreement_note)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("Decision rule", decision$rule_id, "\n")
