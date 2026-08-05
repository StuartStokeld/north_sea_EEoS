# H2/H3 primary model refit with policy-anchored phase breakpoints
#
# PURPOSE: replace data-driven phase breaks (1989 / 2001 / 2008) with
# policy-anchored breaks (1992 CFP / 2002 CFP / 2008 LTMP–MSFD) and refit
# the primary within-between model. Original `phase` column left intact.
#
# phase_v2:
#   1985–1991  Pre-reform baseline
#   1992–2001  1992 CFP reform
#   2002–2007  2002 CFP reform (in force Jan 2003)
#   2008–2015  2008 long-term management plan / MSFD
#
# PRIMARY (revised): residual ~ FP_between * phase_v2 + FP_within * phase_v2
#                    + (1 | stat_rec)   [REML]
#
# Run: Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R
#
# Independent of bootstrap / bathymetry tracks.

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_phase_v2_refit.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_spatial_autocorr_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

if (!requireNamespace("foreign", quietly = TRUE)) {
  stop("Package 'foreign' required to rebuild the established queen listw from the ICES DBF.")
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
paths <- h2_output_paths(project_root)
path_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_panel <- paths$panel
path_archived_blup <- file.path(project_root, "outputs", "spatial_diagnostics_comparison.csv")

path_out_model <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_out_compare <- file.path(project_root, "outputs", "phase_v2_vs_original_comparison.csv")
path_out_blup <- file.path(project_root, "outputs", "phase_v2_blup_diagnostic.csv")
path_out_run_log <- file.path(project_root, "outputs", "phase_v2_refit_run_log.md")
path_out_session <- file.path(project_root, "outputs", "phase_v2_refit_sessionInfo.txt")

stopifnot(file.exists(path_models), file.exists(path_panel))

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 primary model — policy-anchored phase_v2 refit — run log")
logmsg("")
logmsg(
  "Refits the primary within-between model with policy-anchored phase breakpoints ",
  "(1992 / 2002 / 2008) instead of the data-driven structural breaks (1989 / 2001 / 2008). ",
  "Original `phase` column retained for comparison. No bootstrap or bathymetry logic."
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
logmsg(sprintf("glmmTMB %s", as.character(utils::packageVersion("glmmTMB"))))

# ---------------------------------------------------------------------------
# Load original primary model + data
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")

models <- readRDS(path_models)
if (is.null(models$fit_wb)) {
  stop("h2h3_wb_model_objects.rds missing fit_wb — original primary model not found.")
}
primary_model <- models$fit_wb
dat <- models$data
if (is.null(dat)) stop("h2h3_wb_model_objects.rds missing data.")

stopifnot(
  all(c("year", "phase", "residual", "FP_between", "FP_within", "stat_rec") %in% names(dat))
)

logmsg(sprintf(
  "Loaded original primary model: %s (formula: %s).",
  path_models, paste(deparse(formula(primary_model)), collapse = " ")
))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d.",
  nrow(dat), dplyr::n_distinct(dat$stat_rec), min(dat$year), max(dat$year)
))
logmsg(sprintf(
  "Original phase levels: %s.",
  paste(levels(dat$phase), collapse = " | ")
))

# ---------------------------------------------------------------------------
# Define phase_v2 (keep original phase intact)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## phase_v2 definition (policy-anchored)")
logmsg(
  "breaks = c(-Inf, 1991, 2001, 2007, Inf); right = TRUE → ",
  "1985–1991 | 1992–2001 | 2002–2007 | 2008–2015."
)
logmsg(
  "Anchors: pre-reform baseline; 1992 CFP reform; 2002 CFP reform (in force Jan 2003); ",
  "2008 LTMP / MSFD."
)

dat$phase_v2 <- cut(
  dat$year,
  breaks = c(-Inf, 1991, 2001, 2007, Inf),
  labels = c("1985-1991", "1992-2001", "2002-2007", "2008-2015"),
  right = TRUE
)
if (anyNA(dat$phase_v2)) {
  stop("phase_v2 introduced NAs — check year coverage.")
}

# Boundary-year verification
logmsg("")
logmsg("### Boundary-year check: table(year, phase_v2)")
year_phase_tab <- table(dat$year, dat$phase_v2)
print(year_phase_tab)
# Capture table into run log as markdown
tab_df <- as.data.frame.matrix(year_phase_tab)
tab_df$year <- as.integer(rownames(tab_df))
tab_df <- tab_df[, c("year", setdiff(names(tab_df), "year"))]
# Spot-check boundary years land in the intended phase
boundary_expect <- c(
  "1991" = "1985-1991",
  "1992" = "1992-2001",
  "2001" = "1992-2001",
  "2002" = "2002-2007",
  "2007" = "2002-2007",
  "2008" = "2008-2015"
)
for (yr in names(boundary_expect)) {
  ph <- as.character(unique(dat$phase_v2[dat$year == as.integer(yr)]))
  ok <- identical(ph, boundary_expect[[yr]])
  logmsg(sprintf(
    "  year %s → phase_v2 %s (expected %s) %s",
    yr, paste(ph, collapse = ","), boundary_expect[[yr]],
    if (ok) "OK" else "FAIL"
  ))
  if (!ok) stop("Boundary year misassigned: ", yr)
}

# Phase sizes
phase_sizes <- dat %>%
  group_by(phase_v2) %>%
  summarise(
    n_years = dplyr::n_distinct(year),
    n_obs = dplyr::n(),
    year_min = min(year),
    year_max = max(year),
    .groups = "drop"
  )
logmsg("")
logmsg("### Phase sizes (expect roughly 7 / 10 / 6 / 8 years)")
for (i in seq_len(nrow(phase_sizes))) {
  r <- phase_sizes[i, ]
  logmsg(sprintf(
    "  %s: n_years=%d, n_obs=%d (years %d–%d)",
    as.character(r$phase_v2), r$n_years, r$n_obs, r$year_min, r$year_max
  ))
}

# ---------------------------------------------------------------------------
# Step 1 — Refit
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1 — Refit primary model with phase_v2")

formula_v2 <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)
logmsg("Target formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec) [REML]")
logmsg(
  "Note: stats::update(. ~ . - phase + phase_v2) leaves interaction terms on the old ",
  "`phase` factor; refitting explicitly with the full formula."
)

# update(. ~ . - phase + phase_v2) does not cleanly rewrite interaction terms
# (leaves FP_between:phase / phase:FP_within). Always refit explicitly.
time_v2 <- system.time({
  primary_model_v2 <- glmmTMB(formula_v2, data = dat, REML = TRUE)
})
fitted_fml <- paste(deparse(formula(primary_model_v2)), collapse = " ")
logmsg(sprintf("Refit elapsed: %.2f sec.", time_v2[["elapsed"]]))
logmsg(sprintf("Fitted formula: %s", fitted_fml))
stripped <- gsub("phase_v2", "", fitted_fml)
if (!grepl("phase_v2", fitted_fml, fixed = TRUE) || grepl("\\bphase\\b", stripped)) {
  stop("Refit formula incorrect (missing phase_v2 or still contains phase): ", fitted_fml)
}

# ---------------------------------------------------------------------------
# Step 2 — Compare against original (from saved fit object)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 2 — Compare against original model (values from saved fit)")

extract_ref_coef <- function(fit, term) {
  fe <- summary(fit)$coefficients$cond
  if (!term %in% rownames(fe)) {
    stop("Term '", term, "' not found in fixed effects of ", deparse(formula(fit)))
  }
  list(coef = unname(fe[term, "Estimate"]), se = unname(fe[term, "Std. Error"]))
}

rect_re_var <- function(fit) {
  vc <- VarCorr(fit)
  as.numeric(attr(vc$cond$stat_rec, "stddev")[1])^2
}

fe_between_orig <- extract_ref_coef(primary_model, "FP_between")
fe_within_orig <- extract_ref_coef(primary_model, "FP_within")
fe_between_v2 <- extract_ref_coef(primary_model_v2, "FP_between")
fe_within_v2 <- extract_ref_coef(primary_model_v2, "FP_within")

aic_orig <- as.numeric(stats::AIC(primary_model))
bic_orig <- as.numeric(stats::BIC(primary_model))
aic_v2 <- as.numeric(stats::AIC(primary_model_v2))
bic_v2 <- as.numeric(stats::BIC(primary_model_v2))
re_var_orig <- rect_re_var(primary_model)
re_var_v2 <- rect_re_var(primary_model_v2)

comparison <- data.frame(
  quantity = c(
    "FP_between_coef",
    "FP_between_SE",
    "FP_within_coef",
    "FP_within_SE",
    "AIC",
    "BIC",
    "rectangle_RE_variance"
  ),
  note = c(
    "reference-phase main effect (original phase 1985-1988)",
    "SE of FP_between reference-phase main effect",
    "reference-phase main effect (original phase 1985-1988)",
    "SE of FP_within reference-phase main effect",
    "Akaike information criterion (REML)",
    "Bayesian information criterion (REML)",
    "VarCorr (1|stat_rec) intercept variance"
  ),
  original_phase = c(
    fe_between_orig$coef,
    fe_between_orig$se,
    fe_within_orig$coef,
    fe_within_orig$se,
    aic_orig,
    bic_orig,
    re_var_orig
  ),
  revised_phase_v2 = c(
    fe_between_v2$coef,
    fe_between_v2$se,
    fe_within_v2$coef,
    fe_within_v2$se,
    aic_v2,
    bic_v2,
    re_var_v2
  ),
  stringsAsFactors = FALSE
)
# Clarify that v2 reference phase is 1985-1991
comparison$note[comparison$quantity == "FP_between_coef"] <-
  "reference-phase main effect (orig 1985-1988; v2 1985-1991)"
comparison$note[comparison$quantity == "FP_within_coef"] <-
  "reference-phase main effect (orig 1985-1988; v2 1985-1991)"

write_csv(comparison, path_out_compare)
logmsg("Saved: ", path_out_compare)
logmsg("")
logmsg("| Quantity | Original (phase) | Revised (phase_v2) |")
logmsg("|---|---|---|")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  logmsg(sprintf(
    "| %s | %.6g | %.6g |",
    r$quantity, r$original_phase, r$revised_phase_v2
  ))
}
logmsg(sprintf(
  "Delta AIC (v2 − original) = %.3f; Delta BIC = %.3f.",
  aic_v2 - aic_orig, bic_v2 - bic_orig
))

# ---------------------------------------------------------------------------
# Step 3 — BLUP spatial diagnostic (same listw construction as original check)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 3 — BLUP Moran / Geary diagnostic")
logmsg(
  "No persisted listw RDS exists in outputs/. Reusing the established 158-rectangle ",
  "queen-contiguity construction from run_h2h3_primary_spatial_autocorr_check.R ",
  "(ICES DBF SOUTH/WEST grid; style W; zero.policy=TRUE) — identical to the weights ",
  "that produced the archived original BLUP Moran I ≈ 0.552 / Geary C ≈ 0.441. ",
  "Not a new weights scheme."
)

panel <- readRDS(path_panel)
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
if (nrow(panel) != 158L) {
  stop("Expected archived H2 panel N = 158; got ", nrow(panel))
}

nb_pack <- nb_to_listw_W(build_queen_nb_from_ices_dbf(panel, project_root), zero.policy = TRUE)
listw <- nb_pack$listw
logmsg(sprintf(
  "listw: N=%d; n_isolated=%d; mean degree=%.3f.",
  length(listw$neighbours), nb_pack$n_isolated, mean(lengths(nb_pack$nb))
))

# Align BLUPs to panel order (same as original diagnostic)
blups_orig <- extract_primary_blups_from_fit(primary_model, panel$stat_rec)
blups_v2 <- extract_primary_blups_from_fit(primary_model_v2, panel$stat_rec)

test_orig <- run_spatial_tests(blups_orig$blup, listw, "original_phase_BLUPs")
test_v2 <- run_spatial_tests(blups_v2$blup, listw, "phase_v2_BLUPs")

# Cross-check original against archived comparison CSV if present
if (file.exists(path_archived_blup)) {
  arch <- utils::read.csv(path_archived_blup, stringsAsFactors = FALSE)
  arch_blup <- arch[arch$stage == "Primary model BLUPs", , drop = FALSE]
  if (nrow(arch_blup) == 1L) {
    logmsg(sprintf(
      "Archived original BLUP: Moran I=%.6f, Geary C=%.6f (source: spatial_diagnostics_comparison.csv).",
      arch_blup$morans_i, arch_blup$geary_c
    ))
    logmsg(sprintf(
      "Recomputed original BLUP on same listw: Moran I=%.6f, Geary C=%.6f (|ΔI|=%.3g, |ΔC|=%.3g).",
      test_orig$row$morans_i, test_orig$row$geary_c,
      abs(test_orig$row$morans_i - arch_blup$morans_i),
      abs(test_orig$row$geary_c - arch_blup$geary_c)
    ))
  }
}

blup_diag <- data.frame(
  model = c("original_phase", "phase_v2"),
  variable = "(1|stat_rec) intercepts",
  morans_i = c(test_orig$row$morans_i, test_v2$row$morans_i),
  geary_c = c(test_orig$row$geary_c, test_v2$row$geary_c),
  p_moran = c(test_orig$row$p_moran, test_v2$row$p_moran),
  p_geary = c(test_orig$row$p_geary, test_v2$row$p_geary),
  n = c(test_orig$row$n, test_v2$row$n),
  expected_moran = c(test_orig$row$expected_moran, test_v2$row$expected_moran),
  expected_geary = c(test_orig$row$expected_geary, test_v2$row$expected_geary),
  variance_moran = c(test_orig$row$variance_moran, test_v2$row$variance_moran),
  z_moran = c(test_orig$row$z_moran, test_v2$row$z_moran),
  note = c(
    "data-driven phase 1989/2001/2008; matches archived Moran≈0.552 Geary≈0.441",
    "policy-anchored phase_v2 1992/2002/2008"
  ),
  stringsAsFactors = FALSE
)
write_csv(blup_diag, path_out_blup)
logmsg("Saved: ", path_out_blup)
logmsg("")
logmsg("| Model | Moran's I | Geary's C | p (Moran) | p (Geary) | N |")
logmsg("|---|---|---|---|---|---|")
fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.4g", p)))
}
for (i in seq_len(nrow(blup_diag))) {
  r <- blup_diag[i, ]
  logmsg(sprintf(
    "| %s | %.4f | %.4f | %s | %s | %d |",
    r$model, r$morans_i, r$geary_c, fmt_p(r$p_moran), fmt_p(r$p_geary), r$n
  ))
}

# ---------------------------------------------------------------------------
# Save model artifact
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")

saveRDS(
  list(
    primary_model_v2 = primary_model_v2,
    data = dat,
    formula_v2 = formula_v2,
    phase_v2_breaks = c(-Inf, 1991, 2001, 2007, Inf),
    phase_v2_labels = c("1985-1991", "1992-2001", "2002-2007", "2008-2015"),
    phase_v2_anchors = c(
      "Pre-reform baseline",
      "1992 CFP reform",
      "2002 CFP reform (in force Jan 2003)",
      "2008 long-term management plan / MSFD"
    ),
    original_model_source = path_models,
    phase_sizes = as.data.frame(phase_sizes),
    comparison = comparison,
    blup_diagnostic = blup_diag
  ),
  path_out_model
)
logmsg("Saved: ", path_out_model)
logmsg("- ", path_out_compare)
logmsg("- ", path_out_blup)
logmsg("- ", path_out_session)
logmsg("- ", path_out_run_log, " (this file)")
logmsg("")
logmsg(
  "Note: original `phase` and `h2h3_wb_model_objects.rds` are unchanged. ",
  "`primary_model_v2.rds` is the policy-anchored primary refit."
)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== phase_v2 primary refit complete. ===\n")
