# Audit helpers: H1 metrics and the primary CAR model refit for a modified haul table.
#
# Used by the audit diagnosis scripts (audit/scripts/*.R) to compare data treatments
# without touching the live pipeline. Reuses the live pipeline functions for the panel,
# Couce join, FP decomposition and slope extraction, so the only thing that changes is
# the haul table passed in. Run from the repo root.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB) # numFactor(), used inside build_feasibility_data()
})

for (f in c(
  "h1_common.R", "h2_common.R", "h2_panel_helpers.R",
  "h2h3_feasibility_helpers.R", "h2h3_within_between_helpers.R"
)) {
  source(file.path("pipeline", "R", f))
}

AUDIT_PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

#' Haul key (Survey_Year_Quarter_Country_Ship_HaulNo) from a FishGlob haul_id.
audit_haul_key <- function(haul_id) {
  vapply(strsplit(haul_id, " "), function(p) paste(p[c(1:5, 8)], collapse = "_"), "")
}

#' Recompute ln B and residual columns after B_obs or B_pred changed.
audit_recompute_residuals <- function(haul) {
  haul %>%
    mutate(
      ln_B_obs = log(B_obs),
      ln_B_pred = log(B_pred),
      residual = ln_B_obs - ln_B_pred,
      abs_residual = abs(residual)
    )
}

#' H1 metrics on the natural-log scale.
#' bias = mean(ln B_obs - ln B_pred) (negative = over-prediction);
#' bias-removed R2 = 1 - sum((r - mean r)^2) / SS_tot (1:1 line shifted by the mean bias).
audit_h1_metrics <- function(haul) {
  r <- haul$ln_B_obs - haul$ln_B_pred
  ss_tot <- sum((haul$ln_B_obs - mean(haul$ln_B_obs))^2)
  tibble(
    n_hauls = nrow(haul),
    log_r2_1to1 = log_r2(haul$ln_B_obs, haul$ln_B_pred),
    bias_mean_ln = mean(r),
    log_r2_bias_removed = 1 - sum((r - mean(r))^2) / ss_tot,
    median_ratio_pred_obs = median(haul$B_pred / haul$B_obs),
    cor2_ln = stats::cor(haul$ln_B_obs, haul$ln_B_pred)^2
  )
}

#' Build the H2/H3 haul-level analysis data exactly as the live pipeline does
#' (build_h2_rectangle_residuals -> build_h2_analysis_panel -> build_feasibility_data ->
#' add_fp_within_between -> phase_v2), from a haul table.
audit_build_h2h3_data <- function(haul,
                                  couce_rectangle = readRDS("outputs/h2_couce_rectangle_effort.rds"),
                                  couce_year = readRDS("outputs/h2_couce_year_effort.rds")) {
  residual_panel <- build_h2_rectangle_residuals(haul)
  panel <- build_h2_analysis_panel(
    residual_panel, couce_rectangle,
    min_hauls = H2_MIN_HAULS_DEFAULT, require_fishing = TRUE
  )
  built <- build_feasibility_data(haul, panel, couce_year, H2_YEAR_MIN, H2_YEAR_MAX)
  dat <- built$data
  dat$stat_rec <- factor(as.character(dat$stat_rec), levels = as.character(panel$stat_rec))
  dat <- add_fp_within_between(dat)
  dat$phase_v2 <- cut(
    dat$year,
    breaks = c(-Inf, 1991, 2001, 2007, Inf),
    labels = AUDIT_PHASE_V2, right = TRUE
  )
  dat
}

#' Fit the primary CAR model with response `response` and return the 8 phase slopes
#' (4 H2 FP_between, 4 H3 FP_within) with raw and Bonferroni (x4 per family) p-values.
#' The adjacency is the live one, restricted to the rectangles present (reported).
audit_fit_primary_car <- function(dat, response = c("residual", "abs_residual"),
                                  adj_full = readRDS(
                                    "outputs/h2h3_feasibility_round2_model_objects.rds"
                                  )$adjMatrix) {
  response <- match.arg(response)
  rects <- intersect(rownames(adj_full), unique(as.character(dat$stat_rec)))
  adj <- adj_full[rects, rects, drop = FALSE]
  d <- dat
  d$stat_rec <- factor(as.character(d$stat_rec), levels = rects)
  d$y <- d[[response]]
  fit <- spaMM::fitme(
    y ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec),
    data = d, adjMatrix = adj, method = "REML"
  )
  slopes <- bind_rows(
    extract_wb_phase_slopes_spamm(fit, "FP_between", "audit", "H2", phases = AUDIT_PHASE_V2),
    extract_wb_phase_slopes_spamm(fit, "FP_within", "audit", "H3", phases = AUDIT_PHASE_V2)
  ) %>%
    transmute(
      response = response,
      hypothesis = hypothesis_group, component, phase,
      slope = fp_slope, se = fp_slope_se, p_raw = p_value,
      p_bonf = pmin(1, 4 * p_value)
    )
  attr(slopes, "n_hauls") <- nrow(d)
  attr(slopes, "n_rect") <- length(rects)
  attr(slopes, "rho") <- tryCatch(spaMM::get_ranPars(fit, which = "corrPars")[[1]]$rho,
                                  error = function(e) NA_real_)
  attr(slopes, "fit") <- fit
  slopes
}

#' One treatment: H1 metrics + CAR slopes for signed and absolute residual.
audit_run_treatment <- function(haul, label) {
  haul <- audit_recompute_residuals(haul)
  h1 <- audit_h1_metrics(haul) %>% mutate(treatment = label, .before = 1)
  dat <- audit_build_h2h3_data(haul)
  sl <- bind_rows(lapply(c("residual", "abs_residual"), function(rsp) {
    s <- audit_fit_primary_car(dat, rsp)
    s %>% mutate(n_hauls_h2h3 = attr(s, "n_hauls"), n_rect = attr(s, "n_rect"),
                 rho = attr(s, "rho"))
  })) %>% mutate(treatment = label, .before = 1)
  list(h1 = h1, slopes = sl, data = dat)
}

#' Compact markdown table of slopes: one row per treatment x response x component x phase.
audit_format_slopes <- function(sl) {
  sl %>%
    mutate(
      cell = sprintf(
        "%+.3f (%.3f)%s", slope, se,
        ifelse(p_bonf < 0.05, " **", ifelse(p_raw < 0.05, " *", ""))
      )
    )
}

#' Write a data frame as a GitHub markdown table.
audit_md_table <- function(df, digits = 3) {
  df <- as.data.frame(df)
  fmt <- function(x) {
    if (is.numeric(x)) {
      ifelse(is.na(x), "NA", ifelse(abs(x - round(x)) < 1e-12 & abs(x) < 1e7,
                                     format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE),
                                     formatC(x, digits = digits, format = "fg", flag = "#")))
    } else {
      as.character(x)
    }
  }
  body <- vapply(seq_len(nrow(df)), function(i) {
    paste0("| ", paste(vapply(df[i, , drop = FALSE], fmt, ""), collapse = " | "), " |")
  }, "")
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body
  )
}

# ---------------------------------------------------------------------------
# EEoS solver (same reticulate path as pipeline/build_eeos_predictions.R) + convergence check
# ---------------------------------------------------------------------------

#' Initialise reticulate with the project venv and source biomass.py once.
audit_eeos_init <- function() {
  if (!exists(".audit_eeos_ready", envir = .GlobalEnv)) {
    reticulate::use_python(".venv/bin/python", required = TRUE)
    reticulate::source_python("equation_of_state/biomass.py", envir = .GlobalEnv)
    reticulate::py_run_string("
def audit_predict(S, N, E):
    import pandas as pd, numpy as np
    s = pd.Series({'S': float(S), 'N': float(N), 'E': float(E)})
    b = biomass(s)
    ls = mete_lambdas(s, thresh=False)
    con = constraints(ls, s)
    return [float(b), float(np.max(np.abs(con)))]
")
    assign(".audit_eeos_ready", TRUE, envir = .GlobalEnv)
  }
  invisible(TRUE)
}

#' Dimensionless EEoS biomass and the max |relative constraint deviation| for each (S, N, E).
audit_eeos_predict <- function(S, N, E) {
  audit_eeos_init()
  out <- t(vapply(seq_along(S), function(i) {
    r <- tryCatch(unlist(reticulate::py$audit_predict(S[i], N[i], E[i])),
                  error = function(e) c(NA_real_, NA_real_))
    r
  }, c(0, 0)))
  tibble(B_pred_norm = out[, 1], max_abs_constraint_dev = out[, 2])
}
