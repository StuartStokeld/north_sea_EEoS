# Spec A on CAR — Task 1 identifiability diagnostic (single fit)
#
# Fit once:
#   residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2
#              + FP_within * phase_v2 + adjacency(1 | stat_rec)
#
# Gate: rho stability vs primary CAR, SE inflation on FP_between_lag,
#       |cor(FP_between_lag, CAR BLUPs)|, convergence.
#
# Run: Rscript --vanilla pipeline/run_h2h3_spec_a_car_identifiability.R

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_spec_a_car_identifiability.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop("Package 'spaMM' required (check .R_libs/spaMM).")
}
suppressPackageStartupMessages(library(spaMM))

PHASE_V2_LEVELS <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

path_primary <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_round2 <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
path_fp_lag <- file.path(project_root, "outputs", "fp_between_lag_rectangle.rds")
path_spec_a_re <- file.path(project_root, "outputs", "primary_model_v2_spec_a.rds")
path_knn <- file.path(project_root, "outputs", "knn_listw_k4.rds")

path_out_fit <- file.path(project_root, "outputs", "primary_model_v2_spec_a_car.rds")
path_out_diag <- file.path(project_root, "outputs", "spec_a_car_identifiability.csv")
path_out_md <- file.path(project_root, "outputs", "spec_a_car_identifiability.md")
path_out_run_log <- file.path(project_root, "outputs", "spec_a_car_identifiability_run_log.md")
path_out_session <- file.path(project_root, "outputs", "spec_a_car_identifiability_sessionInfo.txt")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

extract_car_rho <- function(fit) {
  ran <- spaMM::get_ranPars(fit)
  rho <- ran$corrPars[["1"]]$rho
  rhorange <- attr(ran, "moreargs")[["1"]]$rhorange
  list(
    rho = unname(rho),
    rho_lo = unname(rhorange[1]),
    rho_hi = unname(rhorange[2]),
    lambda = unname(ran$lambda[[1]]),
    phi = unname(ran$phi)
  )
}

extract_car_blups <- function(fit, panel_ids) {
  re <- spaMM::ranef(fit)[[1]]
  ids <- normalize_stat_rec(names(re))
  if (is.null(ids) || anyNA(ids) || !length(ids)) {
    # spaMM may return unnamed vector ordered by factor levels
    levs <- levels(fit$data$stat_rec)
    if (is.null(levs)) levs <- levels(model.frame(fit)$stat_rec)
    ids <- normalize_stat_rec(levs)
  }
  blups <- data.frame(
    stat_rec = ids,
    blup = as.numeric(re),
    stringsAsFactors = FALSE
  )
  panel_ids <- normalize_stat_rec(panel_ids)
  blups <- blups[match(panel_ids, blups$stat_rec), , drop = FALSE]
  if (anyNA(blups$stat_rec) || anyNA(blups$blup)) {
    stop("CAR BLUP alignment to panel failed.")
  }
  blups
}

phase_slope_table <- function(fit, term_name, phases) {
  b <- spaMM::fixef(fit)
  V <- as.matrix(stats::vcov(fit))
  rows <- lapply(seq_along(phases), function(i) {
    ph <- phases[[i]]
    if (i == 1L) {
      est <- unname(b[[term_name]])
      se <- sqrt(V[[term_name, term_name]])
    } else {
      int_term <- find_fp_phase_interaction_name(names(b), term_name, ph)
      cvec <- stats::setNames(rep(0, length(b)), names(b))
      cvec[[term_name]] <- 1
      cvec[[int_term]] <- 1
      est <- as.numeric(sum(cvec * b))
      se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    }
    data.frame(
      term = term_name, phase = ph, estimate = est, se = se,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

logmsg("# Spec A on CAR — Task 1 identifiability diagnostic")
logmsg("")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg(sprintf("spaMM %s", as.character(utils::packageVersion("spaMM"))))

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
stopifnot(file.exists(path_primary), file.exists(path_round2), file.exists(path_fp_lag))

primary <- readRDS(path_primary)
fit_car <- primary$fit_wb_car_v2
if (is.null(fit_car)) stop("primary_model_v2.rds missing fit_wb_car_v2")
dat <- primary$data
if (is.null(dat)) stop("primary_model_v2.rds missing data")

round2 <- readRDS(path_round2)
adjMatrix <- round2$adjMatrix
if (is.null(adjMatrix)) stop("round2 RDS missing adjMatrix")

fp_lag <- readRDS(path_fp_lag)
dat <- join_fp_between_lag_hauls(dat, fp_lag)
dat$stat_rec <- factor(as.character(normalize_stat_rec(dat$stat_rec)),
                       levels = rownames(adjMatrix))
if (anyNA(dat$stat_rec)) stop("stat_rec levels do not match adjMatrix")
dat$phase_v2 <- factor(dat$phase_v2, levels = PHASE_V2_LEVELS)

logmsg("Primary CAR: ", path_primary)
logmsg("adjMatrix: ", path_round2)
logmsg("FP_between_lag: ", path_fp_lag)
logmsg(sprintf("n hauls = %d; n rectangles = %d", nrow(dat), nlevels(dat$stat_rec)))

# Primary CAR rho
rho_primary <- extract_car_rho(fit_car)
logmsg(sprintf(
  "Primary CAR rho = %.6f (admissible [%.6f, %.6f]); lambda = %.6g; phi = %.6g",
  rho_primary$rho, rho_primary$rho_lo, rho_primary$rho_hi,
  rho_primary$lambda, rho_primary$phi
))

# ---------------------------------------------------------------------------
# Fit CAR + lag
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Fit CAR + FP_between_lag")

formula_car_lag <- residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +
  FP_within * phase_v2 + adjacency(1 | stat_rec)
logmsg("Formula: ", paste(deparse(formula_car_lag), collapse = " "))

warns <- character(0)
fit_car_lag <- NULL
t_fit <- system.time({
  withCallingHandlers(
    {
      fit_car_lag <- spaMM::fitme(
        formula_car_lag,
        data = dat,
        adjMatrix = adjMatrix,
        method = "REML"
      )
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
})
logmsg(sprintf("Fit elapsed: %.2f sec", t_fit[["elapsed"]]))
if (length(warns)) {
  logmsg("Warnings:")
  for (w in warns) logmsg("  - ", w)
} else {
  logmsg("Warnings: none")
}

rho_lag <- extract_car_rho(fit_car_lag)
logmsg(sprintf(
  "CAR+lag rho = %.6f (admissible [%.6f, %.6f]); lambda = %.6g; phi = %.6g",
  rho_lag$rho, rho_lag$rho_lo, rho_lag$rho_hi, rho_lag$lambda, rho_lag$phi
))
delta_rho <- rho_lag$rho - rho_primary$rho
frac_primary <- (rho_primary$rho - rho_primary$rho_lo) /
  (rho_primary$rho_hi - rho_primary$rho_lo)
frac_lag <- (rho_lag$rho - rho_lag$rho_lo) / (rho_lag$rho_hi - rho_lag$rho_lo)
logmsg(sprintf(
  "Δ rho (lag − primary) = %+.6f; frac of upper bound: primary=%.4f, lag=%.4f",
  delta_rho, frac_primary, frac_lag
))

# ---------------------------------------------------------------------------
# Coefficients / SEs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Phase coefficients and SEs")

slopes_fp <- phase_slope_table(fit_car_lag, "FP_between", PHASE_V2_LEVELS)
slopes_lag <- phase_slope_table(fit_car_lag, "FP_between_lag", PHASE_V2_LEVELS)

# RE Spec A SEs for comparison if available
se_re_lag <- rep(NA_real_, length(PHASE_V2_LEVELS))
names(se_re_lag) <- PHASE_V2_LEVELS
if (file.exists(path_spec_a_re)) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    logmsg("glmmTMB not available — skipping RE Spec A SE comparison")
  } else {
    suppressPackageStartupMessages(library(glmmTMB))
    re_obj <- readRDS(path_spec_a_re)
    re_fit <- re_obj$primary_model_spec_a
    re_tab <- extract_wb_phase_slopes(
      re_fit, "FP_between_lag", "spec_a_re", "lag", phases = PHASE_V2_LEVELS
    )
    se_re_lag[re_tab$phase] <- re_tab$fp_slope_se
    logmsg("Loaded RE Spec A for lag SE comparison")
  }
}

coef_tab <- dplyr::bind_rows(
  slopes_fp %>% mutate(model = "CAR+lag", role = "FP_between"),
  slopes_lag %>% mutate(model = "CAR+lag", role = "FP_between_lag")
)
coef_tab$se_re_spec_a_lag <- ifelse(
  coef_tab$role == "FP_between_lag",
  unname(se_re_lag[coef_tab$phase]),
  NA_real_
)
# SE ratio: lag SE / FP_between SE in same CAR+lag model
fp_se_map <- stats::setNames(slopes_fp$se, slopes_fp$phase)
coef_tab$se_ratio_vs_fp_between <- ifelse(
  coef_tab$role == "FP_between_lag",
  coef_tab$se / fp_se_map[coef_tab$phase],
  NA_real_
)
coef_tab$se_ratio_vs_re_lag <- ifelse(
  coef_tab$role == "FP_between_lag" & is.finite(coef_tab$se_re_spec_a_lag),
  coef_tab$se / coef_tab$se_re_spec_a_lag,
  NA_real_
)

for (i in seq_len(nrow(slopes_fp))) {
  logmsg(sprintf(
    "  FP_between %s: est=%+.6f SE=%.6f",
    slopes_fp$phase[i], slopes_fp$estimate[i], slopes_fp$se[i]
  ))
}
for (i in seq_len(nrow(slopes_lag))) {
  logmsg(sprintf(
    "  FP_between_lag %s: est=%+.6f SE=%.6f; SE/FP_between=%.2f; SE/RE_lag=%s",
    slopes_lag$phase[i], slopes_lag$estimate[i], slopes_lag$se[i],
    slopes_lag$se[i] / slopes_fp$se[i],
    if (is.finite(se_re_lag[slopes_lag$phase[i]])) {
      sprintf("%.2f", slopes_lag$se[i] / se_re_lag[slopes_lag$phase[i]])
    } else {
      "NA"
    }
  ))
}

# ---------------------------------------------------------------------------
# cor(FP_between_lag, CAR BLUPs)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## cor(FP_between_lag, CAR spatial BLUPs)")

panel_ids <- rownames(adjMatrix)
blups_primary <- extract_car_blups(fit_car, panel_ids)
blups_lag <- extract_car_blups(fit_car_lag, panel_ids)
lag_map <- fp_lag[, c("stat_rec", "FP_between_lag")]
lag_map$stat_rec <- normalize_stat_rec(lag_map$stat_rec)
lag_aligned <- lag_map$FP_between_lag[match(normalize_stat_rec(panel_ids), lag_map$stat_rec)]

r_primary <- stats::cor(lag_aligned, blups_primary$blup)
r_lagfit <- stats::cor(lag_aligned, blups_lag$blup)
logmsg(sprintf(
  "cor(FP_between_lag, primary-CAR BLUPs) = %+.4f",
  r_primary
))
logmsg(sprintf(
  "cor(FP_between_lag, CAR+lag BLUPs) = %+.4f",
  r_lagfit
))

# ---------------------------------------------------------------------------
# Gate decision
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Gate decision")

# Flags (conservative; brief guide |r|>0.7; large rho shift; SE inflation)
flag_rho <- abs(delta_rho) > 0.02 || abs(frac_lag - frac_primary) > 0.15
flag_corr <- abs(r_primary) > 0.7 || abs(r_lagfit) > 0.7
se_ratios_fp <- coef_tab$se_ratio_vs_fp_between[
  coef_tab$role == "FP_between_lag" & is.finite(coef_tab$se_ratio_vs_fp_between)
]
se_ratios_re <- coef_tab$se_ratio_vs_re_lag[
  coef_tab$role == "FP_between_lag" & is.finite(coef_tab$se_ratio_vs_re_lag)
]
flag_se <- any(se_ratios_fp > 2.5) || any(se_ratios_re > 2.5)
flag_conv <- length(warns) > 0L
# Also flag if rho pinned at bound in both (already near-bound primary) AND
# lag absorbs spatial structure via high BLUP correlation
near_bound_both <- frac_primary > 0.95 && frac_lag > 0.95

gate_clear <- !flag_rho && !flag_corr && !flag_se && !flag_conv

logmsg(sprintf("Flag rho shift: %s (Δrho=%+.6f)", flag_rho, delta_rho))
logmsg(sprintf("Flag |cor(lag, BLUP)| > 0.7: %s (primary=%+.4f, lagfit=%+.4f)",
               flag_corr, r_primary, r_lagfit))
logmsg(sprintf(
  "Flag SE inflation: %s (max SE_lag/SE_FP=%.2f; max SE_lag/SE_RE=%.2f)",
  flag_se,
  if (length(se_ratios_fp)) max(se_ratios_fp) else NA_real_,
  if (length(se_ratios_re)) max(se_ratios_re) else NA_real_
))
logmsg(sprintf("Flag convergence warnings: %s", flag_conv))
logmsg(sprintf("Both fits near CAR rho upper bound: %s", near_bound_both))

if (gate_clear) {
  gate <- "CLEAR"
  gate_note <- paste0(
    "Task 1 clears: proceed to Task 2 (CAR+lag permutation bootstrap)."
  )
} else {
  gate <- "GATED"
  reasons <- c(
    if (flag_rho) "rho shifted materially vs primary CAR",
    if (flag_corr) "|cor(FP_between_lag, CAR BLUPs)| > 0.7",
    if (flag_se) "FP_between_lag SEs inflated vs FP_between and/or RE Spec A",
    if (flag_conv) "fit produced warnings"
  )
  gate_note <- paste0(
    "Task 1 does not clear (", paste(reasons, collapse = "; "),
    "). Skip Task 2; run Task 2b fallback."
  )
}
logmsg("")
logmsg("**GATE: ", gate, "**")
logmsg(gate_note)

# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
diag_row <- data.frame(
  gate = gate,
  rho_primary = rho_primary$rho,
  rho_car_lag = rho_lag$rho,
  delta_rho = delta_rho,
  frac_upper_primary = frac_primary,
  frac_upper_car_lag = frac_lag,
  cor_lag_blup_primary = r_primary,
  cor_lag_blup_car_lag = r_lagfit,
  max_se_ratio_vs_fp_between = if (length(se_ratios_fp)) max(se_ratios_fp) else NA_real_,
  max_se_ratio_vs_re_lag = if (length(se_ratios_re)) max(se_ratios_re) else NA_real_,
  n_warnings = length(warns),
  fit_sec = unname(t_fit[["elapsed"]]),
  stringsAsFactors = FALSE
)
write_csv(diag_row, path_out_diag)
write_csv(coef_tab, sub("identifiability.csv", "identifiability_coefs.csv", path_out_diag))

saveRDS(
  list(
    fit_car_lag = fit_car_lag,
    data = dat,
    formula = formula_car_lag,
    adjMatrix = adjMatrix,
    fp_between_lag_rectangle = fp_lag,
    rho_primary = rho_primary,
    rho_car_lag = rho_lag,
    blups_primary = blups_primary,
    blups_car_lag = blups_lag,
    coef_table = coef_tab,
    diagnostic = diag_row,
    warnings = warns,
    gate = gate,
    gate_note = gate_note
  ),
  path_out_fit
)

md <- c(
  "# Spec A on CAR — Task 1 identifiability diagnostic",
  "",
  sprintf("**GATE: %s**", gate),
  "",
  gate_note,
  "",
  "## Rho",
  "",
  sprintf("- Primary CAR rho: **%.6f** (%.1f%% of admissible upper bound)",
          rho_primary$rho, 100 * frac_primary),
  sprintf("- CAR+lag rho: **%.6f** (%.1f%% of admissible upper bound)",
          rho_lag$rho, 100 * frac_lag),
  sprintf("- Δ rho: **%+.6f**", delta_rho),
  "",
  "## cor(FP_between_lag, CAR BLUPs)",
  "",
  sprintf("- vs primary-CAR BLUPs: **%+.4f**", r_primary),
  sprintf("- vs CAR+lag BLUPs: **%+.4f**", r_lagfit),
  "",
  "## FP_between_lag phase SEs (CAR+lag)",
  "",
  "| Phase | est | SE | SE / FP_between SE | SE / RE Spec A lag SE |",
  "|-------|-----|----|--------------------|------------------------|"
)
for (i in seq_len(nrow(slopes_lag))) {
  ph <- slopes_lag$phase[i]
  md <- c(md, sprintf(
    "| %s | %+.4f | %.4f | %.2f | %s |",
    ph, slopes_lag$estimate[i], slopes_lag$se[i],
    slopes_lag$se[i] / slopes_fp$se[i],
    if (is.finite(se_re_lag[ph])) sprintf("%.2f", slopes_lag$se[i] / se_re_lag[ph]) else "NA"
  ))
}
md <- c(
  md, "",
  sprintf("- Warnings: %s", if (length(warns)) paste(warns, collapse = "; ") else "none"),
  sprintf("- Fit time: %.2f sec", t_fit[["elapsed"]]),
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", path_out_fit),
  sprintf("- `%s`", path_out_diag),
  ""
)
writeLines(md, path_out_md)
writeLines(run_log, path_out_run_log)
cat("\nSaved:", path_out_md, "\n")
cat("GATE:", gate, "\n")
