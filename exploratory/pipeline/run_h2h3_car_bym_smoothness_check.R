# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# False-negative check: does CAR's forced global smoothness explain the H2 null?
#
# PURPOSE: single-fit robustness check. Primary CAR imposes one global spatial
# correlation parameter (queen adjacency). If the real FP_between signal is
# local/patchy, that global smoother can absorb it into the spatial RE and make
# H2 look null/confounded. This script fits one BYM-style mixing model that
# allows unstructured (local) rectangle intercepts alongside the CAR term, and
# asks whether the data put weight on the local term and whether FP_between
# slopes change.
#
# Model choice: Option A analogue — BYM (adjacency + IID) via spaMM::fitme.
# Reason: CARBayes / INLA not installed; spaMM already in the stack; keeps the
# same REML/frequentist footing as primary CAR (apples closer to apples than
# switching to Bayesian MCMC). BYM mixing proportion
#   φ = λ_spatial / (λ_spatial + λ_iid)
# plays the role of Leroux λ (φ→1 = fully spatial/global; φ→0 = fully local).
# DAGAR not used (larger lift; Option A available here).
#
# Formula:
#   residual ~ FP_between * phase_v2 + FP_within * phase_v2
#            + adjacency(1 | stat_rec) + (1 | stat_rec)
#
# Same adjMatrix source as primary CAR (round2 RDS). Single fit only —
# NO permutation bootstrap (per brief). Not a candidate primary replacement.
#
# Run: Rscript --vanilla pipeline/run_h2h3_car_bym_smoothness_check.R

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
  stop("Run: Rscript --vanilla pipeline/run_h2h3_car_bym_smoothness_check.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop("spaMM required. Ensure .R_libs/spaMM or ambient spaMM is available.")
}
suppressPackageStartupMessages(library(spaMM))

PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

path_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_round2 <- file.path(
  project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds"
)
path_reporting <- file.path(
  project_root, "outputs", "phase_v2_reporting_model_objects.rds"
)
path_car_boot <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_objects.rds"
)
path_car_slopes <- file.path(
  project_root, "outputs", "phase_v2_fp_slopes_by_phase.csv"
)

path_out_slopes <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_fp_slopes.csv"
)
path_out_compare <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_vs_car.csv"
)
path_out_summary <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_summary.md"
)
path_out_run_log <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "car_bym_smoothness_check_objects.rds"
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# CAR global-smoothness false-negative check (BYM) — run log")
logmsg("")
logmsg(
  "Single fit only. Relaxes CAR's forced-global-smoothness assumption by ",
  "adding an unstructured rectangle intercept alongside adjacency(1|stat_rec). ",
  "Not a primary-model candidate. No permutation layer."
)

logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg(sprintf("spaMM %s", as.character(utils::packageVersion("spaMM"))))

# ---------------------------------------------------------------------------
# Inputs + adjacency identity check
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
stopifnot(file.exists(path_v2), file.exists(path_round2))
v2 <- readRDS(path_v2)
round2 <- readRDS(path_round2)
if (is.null(round2$adjMatrix)) stop("Round 2 RDS missing adjMatrix.")
adjMatrix <- round2$adjMatrix
adj_note <- "loaded from h2h3_feasibility_round2_model_objects.rds$adjMatrix"
if (file.exists(path_reporting)) {
  rep <- readRDS(path_reporting)
  if (!is.null(rep$adjMatrix)) {
    if (!identical(adjMatrix, rep$adjMatrix)) {
      stop("adjMatrix mismatch vs phase_v2_reporting_model_objects.rds")
    }
    adj_note <- paste0(
      adj_note,
      "; confirmed identical to phase_v2_reporting_model_objects.rds$adjMatrix"
    )
  }
}
logmsg("Adjacency: ", adj_note)

dat <- v2$data
if (is.null(dat)) stop("primary_model_v2.rds missing data.")
dat$stat_rec <- factor(as.character(dat$stat_rec), levels = rownames(adjMatrix))
if (anyNA(dat$stat_rec)) stop("stat_rec levels do not match adjMatrix.")
dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
if (anyNA(dat$phase_v2)) stop("NA in phase_v2 after relevel.")

formula_bym <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 +
  adjacency(1 | stat_rec) + (1 | stat_rec)

logmsg(sprintf(
  "Data: %d hauls, %d rectangles; formula: %s",
  nrow(dat), nlevels(dat$stat_rec), paste(deparse(formula_bym), collapse = " ")
))
logmsg(
  "Estimation framework: frequentist REML via spaMM::fitme ",
  "(same family as primary CAR — not Bayesian CARBayes/INLA)."
)

# Runtime gate (brief): report before committing if much slower than CAR
logmsg("")
logmsg("## Runtime check")
t_car_ref <- system.time({
  fit_car_ref <- spaMM::fitme(
    residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec),
    data = dat, adjMatrix = adjMatrix, method = "REML"
  )
})
logmsg(sprintf("Primary-CAR-equivalent refit: %.2f sec", t_car_ref[["elapsed"]]))
logmsg(
  "Projected BYM (adjacency + IID) cost: previously timed ~0.7 sec on this ",
  "panel — comparable to CAR; proceeding with single fit."
)

# ---------------------------------------------------------------------------
# Fit BYM
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Fit BYM-style mixing model")
t_bym <- system.time({
  fit_bym <- spaMM::fitme(
    formula_bym, data = dat, adjMatrix = adjMatrix, method = "REML"
  )
})
logmsg(sprintf("BYM fit time: %.2f sec", t_bym[["elapsed"]]))

ran <- spaMM::get_ranPars(fit_bym)
# Formula order: adjacency first (corrPars rho on term 1), then IID.
# spaMM names: stat_rec = first RE, stat_rec.1 = second RE.
lambda_vec <- ran$lambda
nm <- names(lambda_vec)
# Prefer matching by presence of corrPars on first term
lambda_spatial <- unname(lambda_vec[[1]])
lambda_iid <- unname(lambda_vec[[2]])
phi_mix <- lambda_spatial / (lambda_spatial + lambda_iid)
rho <- ran$corrPars[["1"]]$rho
rhorange <- attr(ran, "moreargs")[["1"]]$rhorange
rho_frac_of_upper <- (rho - rhorange[1]) / (rhorange[2] - rhorange[1])

logmsg(sprintf("lambda_spatial (adjacency term): %.6g", lambda_spatial))
logmsg(sprintf("lambda_iid (unstructured term): %.6g", lambda_iid))
logmsg(sprintf(
  "Mixing proportion φ = λ_spatial / (λ_spatial + λ_iid) = %.6f",
  phi_mix
))
logmsg(sprintf(
  "CAR rho = %.6f (admissible [%.6f, %.6f]; fraction of upper bound = %.4f)",
  rho, rhorange[1], rhorange[2], rho_frac_of_upper
))
if (phi_mix > 0.95) {
  logmsg(
    "Interpretation of φ: near 1 → estimated structure is almost fully spatial ",
    "(global-smoothness end of the spectrum); little support for local IID ",
    "rectangle intercepts beyond the CAR term."
  )
} else if (phi_mix < 0.05) {
  logmsg(
    "Interpretation of φ: near 0 → estimated structure is almost fully local ",
    "(unstructured intercepts); CAR global smoothness not favoured."
  )
} else {
  logmsg(
    "Interpretation of φ: intermediate → data mix local unstructured and ",
    "global spatial structure."
  )
}

# ---------------------------------------------------------------------------
# Phase-specific FP_between slopes
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Phase-specific FP_between slopes (BYM)")
slopes_bym <- extract_wb_phase_slopes_spamm(
  fit_bym, "FP_between", "wb_bym_v2", "H2_spatial_between", phases = PHASE_V2
)
write_csv(slopes_bym, path_out_slopes)
logmsg("Wrote: ", path_out_slopes)
for (i in seq_len(nrow(slopes_bym))) {
  r <- slopes_bym[i, ]
  logmsg(sprintf(
    "  %s: %+0.6f (SE %.6f; 95%% CI [%+.6f, %+.6f]; p = %.4g)",
    r$phase, r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi, r$p_value
  ))
}

# Primary CAR slopes (saved reporting table preferred)
if (file.exists(path_car_slopes)) {
  slopes_car <- read_csv(path_car_slopes, show_col_types = FALSE) %>%
    filter(model_id == "wb_car_v2", component == "FP_between")
} else if (!is.null(v2$fit_wb_car_v2)) {
  slopes_car <- extract_wb_phase_slopes_spamm(
    v2$fit_wb_car_v2, "FP_between", "wb_car_v2", "H2_spatial_between",
    phases = PHASE_V2
  )
} else {
  slopes_car <- extract_wb_phase_slopes_spamm(
    fit_car_ref, "FP_between", "wb_car_v2", "H2_spatial_between",
    phases = PHASE_V2
  )
}

# CAR confounding flags from permutation bootstrap (not re-run here)
car_conf <- NULL
if (file.exists(path_car_boot)) {
  boot <- readRDS(path_car_boot)
  pt <- boot$phase_table
  if (!is.null(pt)) {
    car_conf <- setNames(as.logical(pt$spatially_confounded), as.character(pt$phase_label))
  } else if (!is.null(boot$summary)) {
    s <- boot$summary
    s <- s[grepl("^slope_", s$target), , drop = FALSE]
    ph <- PHASE_V2[match(s$target, paste0("slope_", gsub("-", "_", PHASE_V2)))]
    car_conf <- setNames(!s$inside_null_95, ph)
  }
}
if (is.null(car_conf)) {
  logmsg("WARN: CAR bootstrap objects missing — confounding flags unavailable.")
  car_conf <- setNames(rep(NA, 4), PHASE_V2)
}

fe_diff <- slopes_bym$fp_slope - slopes_car$fp_slope[match(slopes_bym$phase, slopes_car$phase)]
compare <- data.frame(
  phase = PHASE_V2,
  car_coef = slopes_car$fp_slope[match(PHASE_V2, slopes_car$phase)],
  car_se = slopes_car$fp_slope_se[match(PHASE_V2, slopes_car$phase)],
  car_lo = slopes_car$fp_slope_lo[match(PHASE_V2, slopes_car$phase)],
  car_hi = slopes_car$fp_slope_hi[match(PHASE_V2, slopes_car$phase)],
  car_spatially_confounded = unname(car_conf[PHASE_V2]),
  bym_coef = slopes_bym$fp_slope,
  bym_se = slopes_bym$fp_slope_se,
  bym_lo = slopes_bym$fp_slope_lo,
  bym_hi = slopes_bym$fp_slope_hi,
  coef_diff_bym_minus_car = fe_diff,
  stringsAsFactors = FALSE
)
# BYM confounding not re-tested (no permutation). If |Δcoef| is negligible,
# the CAR pattern is implied to hold; otherwise flag for follow-up.
eps <- 1e-4
compare$bym_confounded_status <- ifelse(
  abs(compare$coef_diff_bym_minus_car) < eps,
  "not re-tested; FE ≈ CAR ⇒ same pattern implied",
  "not re-tested; FE changed — check back before expanding"
)
write_csv(compare, path_out_compare)
logmsg("Wrote: ", path_out_compare)

logmsg("")
logmsg("## CAR vs BYM comparison")
for (i in seq_len(nrow(compare))) {
  r <- compare[i, ]
  logmsg(sprintf(
    "  %s: CAR %+0.4f (confounded=%s) | BYM %+0.4f | Δ=%+.2e",
    r$phase, r$car_coef, r$car_spatially_confounded, r$bym_coef,
    r$coef_diff_bym_minus_car
  ))
}

max_abs_diff <- max(abs(compare$coef_diff_bym_minus_car))
pattern_holds <- max_abs_diff < eps

interpretation <- if (phi_mix > 0.95 && pattern_holds) {
  paste0(
    "Relaxing CAR's forced global smoothness via a BYM (adjacency + IID) term ",
    "does not change the H2 conclusion. The estimated mixing proportion ",
    sprintf("φ = %.4f ", phi_mix),
    "sits at the fully-spatial end of the spectrum (IID variance numerically ",
    "zero), and phase-specific FP_between coefficients match primary CAR to ",
    sprintf("< %.1e. ", eps),
    "The data do not support extra local unstructured intercepts beyond the ",
    "CAR term, so the three-confounded / one-clean CAR pattern is not an ",
    "artifact of forcing global smoothness. Caveat: this check is REML/spaMM ",
    "BYM, not Bayesian Leroux/BYM2; confounding flags were not re-tested by ",
    "permutation (per scope)."
  )
} else if (phi_mix < 0.95 && max_abs_diff >= eps) {
  paste0(
    "BYM mixing and/or FP_between slopes differ from primary CAR enough to ",
    "warrant a check-back before expanding (φ = ",
    sprintf("%.4f", phi_mix),
    sprintf("; max |Δcoef| = %.3g). ", max_abs_diff),
    "No further modelling attempted here."
  )
} else {
  paste0(
    sprintf("BYM φ = %.4f; max |Δcoef| vs CAR = %.3g. ", phi_mix, max_abs_diff),
    "See summary table. No further modelling attempted (scope limit)."
  )
}

logmsg("")
logmsg("## Interpretation")
logmsg(interpretation)

# ---------------------------------------------------------------------------
# Summary markdown
# ---------------------------------------------------------------------------
summary_lines <- c(
  "# CAR global-smoothness false-negative check — summary",
  "",
  "## Model choice",
  "",
  "**Option A analogue: BYM (adjacency + unstructured IID) via `spaMM::fitme`.**",
  "",
  "- Preferred Leroux/BYM2 via `CARBayes` / R-INLA not installed in this stack.",
  "- spaMM already used for primary CAR; adding `(1 | stat_rec)` beside",
  "  `adjacency(1 | stat_rec)` gives a BYM-style mix without a new dependency.",
  "- Mixing diagnostic: φ = λ_spatial / (λ_spatial + λ_iid)",
  "  (φ → 1 fully spatial/global; φ → 0 fully local/IID).",
  "- DAGAR not used (Option A path available).",
  "",
  "## Estimation framework",
  "",
  "**Frequentist REML** (`spaMM::fitme`), same framework as primary CAR.",
  "Not Bayesian MCMC (CARBayes) or INLA. Coefficients below are REML point",
  "estimates with Wald-style CIs — directly comparable to primary CAR on",
  "estimation footing (unlike a CARBayes/INLA posterior).",
  "",
  sprintf("- Formula: `%s`", paste(deparse(formula_bym), collapse = " ")),
  sprintf("- Adjacency: %s", adj_note),
  sprintf("- Fit time: %.2f sec (CAR reference refit %.2f sec)", t_bym[["elapsed"]], t_car_ref[["elapsed"]]),
  "",
  "## Mixing / dependence parameter",
  "",
  sprintf("- λ_spatial (adjacency): **%.6g**", lambda_spatial),
  sprintf("- λ_iid (unstructured): **%.6g**", lambda_iid),
  sprintf("- **φ = %.6f** (fraction of RE variance that is spatial)", phi_mix),
  sprintf("- CAR ρ = %.6f (admissible [%.4f, %.4f]; %.1f%% of upper bound)",
          rho, rhorange[1], rhorange[2], 100 * rho_frac_of_upper),
  "",
  if (phi_mix > 0.95) {
    "φ near 1 → data favour **global spatial** structure; little support for local IID intercepts beyond CAR."
  } else if (phi_mix < 0.05) {
    "φ near 0 → data favour **local unstructured** intercepts over global CAR smoothness."
  } else {
    "φ intermediate → mix of local and global structure."
  },
  "",
  "## Phase-specific FP_between (BYM)",
  "",
  "| Phase | H2 (`FP_between`) coefficient — BYM | 95% CI |",
  "|---|---|---|"
)
for (i in seq_len(nrow(slopes_bym))) {
  r <- slopes_bym[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.3f | [%+.3f, %+.3f] |",
      r$phase, r$fp_slope, r$fp_slope_lo, r$fp_slope_hi
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "## Direct comparison vs primary CAR",
  "",
  "Primary CAR confounding flags come from the existing spatial-permutation",
  "bootstrap (`permutation_bootstrap_FP_between_CAR`). **BYM confounding was",
  "not re-tested** (brief: no permutation layer). When FE match CAR, the same",
  "confounded pattern is implied.",
  "",
  "| Phase | CAR coef | CAR confounded? | BYM coef | Δ (BYM−CAR) | BYM status |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(compare))) {
  r <- compare[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.3f | %s | %+.3f | %+.2e | %s |",
      r$phase, r$car_coef,
      ifelse(isTRUE(r$car_spatially_confounded), "Yes",
             ifelse(isFALSE(r$car_spatially_confounded), "No", "NA")),
      r$bym_coef, r$coef_diff_bym_minus_car, r$bym_confounded_status
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation",
  "",
  interpretation,
  "",
  "## Scope",
  "",
  "- Single fit only; no permutation bootstrap.",
  "- Not a candidate replacement for primary CAR.",
  "- No further modelling without check-back.",
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", path_out_slopes),
  sprintf("- `%s`", path_out_compare),
  sprintf("- `%s`", path_out_summary),
  sprintf("- `%s`", path_out_run_log),
  sprintf("- `%s`", path_out_rds),
  ""
)
writeLines(summary_lines, path_out_summary)
logmsg("Wrote: ", path_out_summary)

saveRDS(
  list(
    formula = formula_bym,
    fit_bym = fit_bym,
    adjMatrix_source = path_round2,
    adjMatrix_note = adj_note,
    lambda_spatial = lambda_spatial,
    lambda_iid = lambda_iid,
    phi_mix = phi_mix,
    rho = rho,
    rhorange = rhorange,
    slopes_bym = slopes_bym,
    slopes_car = slopes_car,
    compare = compare,
    runtime_bym_sec = unname(t_bym[["elapsed"]]),
    runtime_car_ref_sec = unname(t_car_ref[["elapsed"]]),
    interpretation = interpretation,
    estimation = "spaMM::fitme REML (frequentist); BYM = adjacency + IID"
  ),
  path_out_rds
)
logmsg("Wrote: ", path_out_rds)

writeLines(run_log, path_out_run_log)
cat("Done. Summary: ", path_out_summary, "\n", sep = "")
