# Helpers for the H2/H3 shared hierarchical model RESULTS run
# (run_h2h3_shared_model_results.R). Unlike the feasibility helpers, these
# support interpreting fishing-pressure and fishing-pressure x phase
# coefficients as the substantive H2/H3 findings.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Label fixed-effect terms by which hypothesis they answer.
#' Reference-phase log_hours_total = H2 (FP effect in 1985-1988).
#' log_hours_total:phase* interactions = H3 (change in FP effect relative
#' to the reference phase). Phase main effects are controls, not H2/H3
#' answer terms. (Biomass is not in the corrected primary specification.)
label_h2_h3_terms <- function(term) {
  dplyr::case_when(
    term == "log_hours_total" ~ "H2 (FP effect, reference phase 1985-1988)",
    grepl("^log_hours_total:phase", term) ~ "H3 (FP x phase interaction vs reference)",
    grepl("^phase", term) ~ "control (phase main effect)",
    term == "(Intercept)" ~ "control (intercept)",
    TRUE ~ "other"
  )
}

#' Attach hypothesis labels and a short plain-language term description to a
#' tidy fixed-effects table (columns: term, estimate, ...).
annotate_fixed_effects <- function(fe_table) {
  fe_table %>%
    mutate(
      hypothesis = label_h2_h3_terms(term),
      term_plain = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "log_hours_total" ~ "Fishing pressure (log(hours+1)), phase 1985-1988",
        term == "phase1989-2000" ~ "Phase 1989-2000 (vs 1985-1988)",
        term == "phase2001-2007" ~ "Phase 2001-2007 (vs 1985-1988)",
        term == "phase2008-2015" ~ "Phase 2008-2015 (vs 1985-1988)",
        term == "log_hours_total:phase1989-2000" ~
          "FP x phase 1989-2000 (change in FP slope vs 1985-1988)",
        term == "log_hours_total:phase2001-2007" ~
          "FP x phase 2001-2007 (change in FP slope vs 1985-1988)",
        term == "log_hours_total:phase2008-2015" ~
          "FP x phase 2008-2015 (change in FP slope vs 1985-1988)",
        TRUE ~ term
      )
    )
}

#' Nakagawa & Schielzeth (2013) marginal / conditional R^2 for a gaussian
#' glmmTMB mixed model with a single random intercept. Implemented here
#' without the performance/MuMIn packages (neither is installed in this
#' project's ambient library).
nakagawa_r2_glmmtmb <- function(fit) {
  fe_pred <- as.numeric(stats::predict(fit, re.form = NA))
  var_f <- stats::var(fe_pred)
  vc <- glmmTMB::VarCorr(fit)
  var_r <- as.numeric(attr(vc$cond$stat_rec, "stddev")[1])^2
  var_e <- as.numeric(stats::sigma(fit))^2
  denom <- var_f + var_r + var_e
  list(
    r2_marginal = var_f / denom,
    r2_conditional = (var_f + var_r) / denom,
    var_fixed = var_f,
    var_random = var_r,
    var_residual = var_e
  )
}

#' Per-phase fishing-pressure slopes from a spaMM CAR fit (same delta-method
#' construction as extract_phase_fp_slopes() for glmmTMB in the temporal
#' robustness helpers).
extract_phase_fp_slopes_spamm <- function(fit_car) {
  b <- spaMM::fixef(fit_car)
  V <- as.matrix(stats::vcov(fit_car))
  # Align names: spaMM may drop "(Intercept)" naming inconsistently; use names(b)
  bounds <- data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    year_start = c(1985L, 1989L, 2001L, 2008L),
    year_end = c(1988L, 2000L, 2007L, 2015L),
    stringsAsFactors = FALSE
  )
  ref_phase <- bounds$phase[1]
  rows <- lapply(seq_len(nrow(bounds)), function(i) {
    ph <- bounds$phase[i]
    if (identical(ph, ref_phase)) {
      est <- unname(b[["log_hours_total"]])
      se <- sqrt(V[["log_hours_total", "log_hours_total"]])
    } else {
      term <- paste0("log_hours_total:phase", ph)
      if (!term %in% names(b)) {
        stop("Expected interaction term not found in CAR model: ", term)
      }
      cvec <- stats::setNames(rep(0, length(b)), names(b))
      cvec[["log_hours_total"]] <- 1
      cvec[[term]] <- 1
      est <- as.numeric(sum(cvec * b))
      se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    }
    data.frame(
      model_id = "car_adjacency",
      phase = ph,
      year_start = bounds$year_start[i],
      year_end = bounds$year_end[i],
      fp_slope = est,
      fp_slope_se = se,
      fp_slope_lo = est - 1.96 * se,
      fp_slope_hi = est + 1.96 * se,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Summarise GAM year-by-year FP slopes within each categorical phase
#' (mean slope and mean CI bounds over years in the phase) for side-by-side
#' comparison with the phase-model slopes. Not a parametric phase effect —
#' a descriptive aggregation of the smooth, flagged as such.
summarise_gam_slopes_by_phase <- function(gam_year_slopes) {
  bounds <- data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    year_start = c(1985L, 1989L, 2001L, 2008L),
    year_end = c(1988L, 2000L, 2007L, 2015L),
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(lapply(seq_len(nrow(bounds)), function(i) {
    ph <- bounds$phase[i]
    ys <- gam_year_slopes %>%
      filter(year >= bounds$year_start[i], year <= bounds$year_end[i])
    data.frame(
      model_id = "smooth_year_gam",
      phase = ph,
      year_start = bounds$year_start[i],
      year_end = bounds$year_end[i],
      fp_slope = mean(ys$fp_slope),
      fp_slope_se = mean(ys$fp_slope_se),
      fp_slope_lo = mean(ys$fp_slope_lo),
      fp_slope_hi = mean(ys$fp_slope_hi),
      note = "mean of year-by-year GAM slopes within phase (not a parametric phase coefficient)",
      stringsAsFactors = FALSE
    )
  }))
}

#' Approximate observation-level Cook's distance for a gaussian glmmTMB fit,
#' using Pearson residuals and the fixed-effects hat matrix (random effects
#' ignored in the hat calculation — a standard first-pass influence screen,
#' not a full conditional Cook's D). Builds X from the known fixed-effects
#' formula so lme4 is not required.
approximate_cooks_glmmtmb_noforte <- function(fit, data) {
  X <- stats::model.matrix(~ log_hours_total * phase, data = data)
  hat <- tryCatch({
    Q <- qr(X)
    rowSums(qr.Q(Q)^2)
  }, error = function(e) rep(NA_real_, nrow(data)))
  mu <- as.numeric(stats::predict(fit, re.form = NULL))
  resid_pearson <- (data$residual - mu) / stats::sigma(fit)
  p <- ncol(X)
  cooks <- (resid_pearson^2) * hat / (p * pmax(1 - hat, 1e-8)^2)
  data.frame(
    row_id = seq_len(nrow(data)),
    stat_rec = as.character(data$stat_rec),
    year = data$year,
    hat = hat,
    pearson_resid = resid_pearson,
    cooks_approx = cooks,
    stringsAsFactors = FALSE
  )
}

#' Joint Wald test that a set of named coefficients are jointly zero.
#' Uses the coefficient vector and vcov; returns chi-square, df, p.
wald_joint_zero <- function(beta, V, term_names) {
  missing <- setdiff(term_names, names(beta))
  if (length(missing) > 0L) {
    stop("Terms not found in coefficient vector: ", paste(missing, collapse = ", "))
  }
  b <- beta[term_names]
  v <- V[term_names, term_names, drop = FALSE]
  # chi^2 = b' V^{-1} b
  w <- as.numeric(t(b) %*% solve(v, b))
  df <- length(term_names)
  p <- stats::pchisq(w, df = df, lower.tail = FALSE)
  list(statistic = w, df = df, p_value = p, terms = term_names)
}
