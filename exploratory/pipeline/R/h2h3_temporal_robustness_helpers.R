# Helpers for the H2/H3 TEMPORAL ROBUSTNESS CHECK
# (run_h2h3_temporal_robustness.R): categorical phase vs continuous linear year
# vs smooth year. This is a feasibility/robustness check, not a results run —
# nothing here interprets fishing-pressure or fishing-pressure x time
# coefficients as answering H2/H3.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Centre year at the analysis-data mean for numerical stability. Returns a
#' list with the centred vector and the centre used (so it can be recovered
#' for prediction / slope-vs-year extraction). Same mean-centring for every
#' model that uses year_c within a single run.
build_year_centred <- function(year) {
  centre <- mean(year, na.rm = TRUE)
  list(year_c = year - centre, year_centre = centre)
}

#' Phase calendar bounds matching build_phase_factor() in
#' h2h3_feasibility_helpers.R (break year = first year of the new phase).
phase_calendar_bounds <- function() {
  data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    year_start = c(1985L, 1989L, 2001L, 2008L),
    year_end = c(1988L, 2000L, 2007L, 2015L),
    stringsAsFactors = FALSE
  )
}

#' Per-phase fishing-pressure slope (partial derivative of predicted residual
#' w.r.t. log_hours_total) from the categorical-phase glmmTMB model, with
#' delta-method SE. Reference phase slope = coef(log_hours_total); other
#' phases = that plus the corresponding interaction. Returns one row per
#' phase, with year_start/year_end for step-function plotting.
extract_phase_fp_slopes <- function(fit_phase) {
  b <- glmmTMB::fixef(fit_phase)$cond
  V <- stats::vcov(fit_phase)$cond
  bounds <- phase_calendar_bounds()
  ref_phase <- bounds$phase[1]

  rows <- lapply(seq_len(nrow(bounds)), function(i) {
    ph <- bounds$phase[i]
    if (identical(ph, ref_phase)) {
      est <- unname(b[["log_hours_total"]])
      se <- sqrt(V[["log_hours_total", "log_hours_total"]])
    } else {
      term <- paste0("log_hours_total:phase", ph)
      if (!term %in% names(b)) {
        stop("Expected interaction term not found in phase model: ", term)
      }
      cvec <- stats::setNames(rep(0, length(b)), names(b))
      cvec[["log_hours_total"]] <- 1
      cvec[[term]] <- 1
      est <- as.numeric(sum(cvec * b))
      se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    }
    data.frame(
      model_id = "categorical_phase",
      year = NA_integer_,
      year_start = bounds$year_start[i],
      year_end = bounds$year_end[i],
      phase = ph,
      fp_slope = est,
      fp_slope_se = se,
      fp_slope_lo = est - 1.96 * se,
      fp_slope_hi = est + 1.96 * se,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Year-by-year fishing-pressure slope from the continuous-linear-year
#' glmmTMB model: slope(year) = b_log + b_interact * (year - year_centre),
#' with delta-method SE. Returns one row per year in `years`.
extract_linear_fp_slopes <- function(fit_linear, years, year_centre) {
  b <- glmmTMB::fixef(fit_linear)$cond
  V <- stats::vcov(fit_linear)$cond
  year_c <- years - year_centre
  b_log <- unname(b[["log_hours_total"]])
  b_int <- unname(b[["log_hours_total:year_c"]])
  v_ll <- V[["log_hours_total", "log_hours_total"]]
  v_ii <- V[["log_hours_total:year_c", "log_hours_total:year_c"]]
  v_li <- V[["log_hours_total", "log_hours_total:year_c"]]

  est <- b_log + b_int * year_c
  se <- sqrt(v_ll + year_c^2 * v_ii + 2 * year_c * v_li)
  data.frame(
    model_id = "continuous_linear_year",
    year = as.integer(years),
    year_start = NA_integer_,
    year_end = NA_integer_,
    phase = NA_character_,
    fp_slope = est,
    fp_slope_se = se,
    fp_slope_lo = est - 1.96 * se,
    fp_slope_hi = est + 1.96 * se,
    stringsAsFactors = FALSE
  )
}

#' Year-by-year fishing-pressure slope from the varying-coefficient GAM:
#' formula uses s(year, by = log_hours_total), which contributes
#' log_hours_total * f(year) to the linear predictor. Evaluating
#' predict(..., type = "terms") at log_hours_total = 1 therefore recovers
#' f(year) itself — the partial derivative of predicted residual w.r.t.
#' log_hours_total as a function of year. Uses a dummy rectangle for the
#' random intercept (its contribution is zero under type = "terms" for the
#' by-smooth column, but mgcv still requires a valid factor level).
extract_gam_fp_slopes <- function(fit_gam, years, data) {
  # Include mean_ln_B_obs only when present in the fitted formula (older
  # with-biomass fits); omit for the corrected no-biomass GAM.
  nd <- data.frame(
    year = as.integer(years),
    log_hours_total = 1,
    stat_rec = factor(levels(data$stat_rec)[1], levels = levels(data$stat_rec))
  )
  if ("mean_ln_B_obs" %in% all.vars(stats::formula(fit_gam)) &&
        "mean_ln_B_obs" %in% names(data)) {
    nd$mean_ln_B_obs <- mean(data$mean_ln_B_obs, na.rm = TRUE)
  }
  pr <- stats::predict(fit_gam, newdata = nd, type = "terms", se.fit = TRUE)
  by_col <- grep("log_hours_total", colnames(pr$fit), value = TRUE)
  if (length(by_col) != 1L) {
    stop(
      "Expected exactly one by-smooth column matching log_hours_total in ",
      "predict(type='terms'); found: ", paste(colnames(pr$fit), collapse = ", ")
    )
  }
  est <- as.numeric(pr$fit[, by_col])
  se <- as.numeric(pr$se.fit[, by_col])
  data.frame(
    model_id = "smooth_year_gam",
    year = as.integer(years),
    year_start = NA_integer_,
    year_end = NA_integer_,
    phase = NA_character_,
    fp_slope = est,
    fp_slope_se = se,
    fp_slope_lo = est - 1.96 * se,
    fp_slope_hi = est + 1.96 * se,
    stringsAsFactors = FALSE
  )
}

#' Expand the phase-model's per-phase slopes onto a year grid as a step
#' function (same value for every year within a phase), for overlay plotting
#' against the continuous/smooth curves.
expand_phase_slopes_to_years <- function(phase_slopes, years) {
  bounds <- phase_calendar_bounds()
  dplyr::bind_rows(lapply(years, function(y) {
    row <- bounds[y >= bounds$year_start & y <= bounds$year_end, , drop = FALSE]
    if (nrow(row) != 1L) stop("Year ", y, " does not fall in exactly one phase.")
    ph <- phase_slopes[phase_slopes$phase == row$phase, , drop = FALSE]
    data.frame(
      model_id = "categorical_phase",
      year = as.integer(y),
      year_start = ph$year_start,
      year_end = ph$year_end,
      phase = ph$phase,
      fp_slope = ph$fp_slope,
      fp_slope_se = ph$fp_slope_se,
      fp_slope_lo = ph$fp_slope_lo,
      fp_slope_hi = ph$fp_slope_hi,
      stringsAsFactors = FALSE
    )
  }))
}

#' Unified effects summary row(s) for a glmmTMB parametric model (phase or
#' linear). Same columns as the GAM summary below where possible; terms that
#' don't apply (edf, Ref.df) are NA.
tidy_parametric_effects <- function(fit, model_id) {
  s <- summary(fit)
  coef_mat <- s$coefficients$cond
  data.frame(
    model_id = model_id,
    term_type = "parametric",
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, "z value"],
    p_value = coef_mat[, "Pr(>|z|)"],
    edf = NA_real_,
    ref_df = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Unified effects summary for an mgcv GAM: parametric coefficients PLUS
#' smooth-term rows (edf / Ref.df / F / approx p). Smooth rows leave
#' estimate/std_error as NA — they are not a single slope.
tidy_gam_effects <- function(fit, model_id) {
  s <- summary(fit)
  p_rows <- if (!is.null(s$p.table) && nrow(s$p.table) > 0L) {
    data.frame(
      model_id = model_id,
      term_type = "parametric",
      term = rownames(s$p.table),
      estimate = s$p.table[, "Estimate"],
      std_error = s$p.table[, "Std. Error"],
      statistic = s$p.table[, "t value"],
      p_value = s$p.table[, "Pr(>|t|)"],
      edf = NA_real_,
      ref_df = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
  s_rows <- if (!is.null(s$s.table) && nrow(s$s.table) > 0L) {
    data.frame(
      model_id = model_id,
      term_type = "smooth",
      term = rownames(s$s.table),
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = s$s.table[, "F"],
      p_value = s$s.table[, "p-value"],
      edf = s$s.table[, "edf"],
      ref_df = s$s.table[, "Ref.df"],
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
  dplyr::bind_rows(p_rows, s_rows)
}

#' Convergence / numerical diagnostics for a glmmTMB fit (same fields as the
#' Round 1 feasibility helper, kept local so this script does not depend on
#' Round 1 helpers beyond data prep).
glmmtmb_convergence_brief <- function(fit, model_id) {
  conv_code <- fit$fit$convergence
  conv_msg <- fit$fit$message
  pdHess <- isTRUE(fit$sdr$pdHess)
  grad <- fit$sdr$gradient.fixed
  max_abs_grad <- if (!is.null(grad)) max(abs(grad), na.rm = TRUE) else NA_real_
  data.frame(
    model_id = model_id,
    framework = "glmmTMB",
    converged = isTRUE(conv_code == 0),
    optimizer_convergence_code = conv_code,
    optimizer_message = ifelse(is.null(conv_msg) || conv_msg == "", "none", conv_msg),
    hessian_positive_definite = pdHess,
    max_abs_gradient = max_abs_grad,
    gam_iterations = NA_integer_,
    gam_optimizer = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Convergence / numerical diagnostics for an mgcv::gam fit (outer-newton
#' REML): converged flag, iteration count, gradient range, Hessian positive-
#' definiteness, and the basis-dimension (k) check summary as a single
#' text field for the run log.
mgcv_convergence_brief <- function(fit, model_id) {
  conv <- isTRUE(fit$converged)
  # mgcv stores outer-iteration info in fit$outer.info when available
  n_iter <- if (!is.null(fit$iter)) as.integer(fit$iter) else NA_integer_
  opt_name <- if (!is.null(fit$method)) as.character(fit$method) else "REML"
  grad_range <- tryCatch({
    oi <- fit$outer.info
    if (!is.null(oi) && !is.null(oi$grad)) max(abs(oi$grad), na.rm = TRUE) else NA_real_
  }, error = function(e) NA_real_)
  pdHess <- tryCatch({
    oi <- fit$outer.info
    if (!is.null(oi) && !is.null(oi$hess)) {
      ev <- eigen(oi$hess, only.values = TRUE, symmetric = TRUE)$values
      all(ev > 0)
    } else {
      NA
    }
  }, error = function(e) NA)

  k_check_text <- tryCatch({
    # Capture the printed k-index table from gam.check without plotting
    out <- utils::capture.output({
      old_dev <- grDevices::dev.cur()
      grDevices::pdf(NULL)
      on.exit({
        grDevices::dev.off()
        if (old_dev > 1) grDevices::dev.set(old_dev)
      }, add = TRUE)
      mgcv::gam.check(fit, rep = 0)
    })
    paste(out, collapse = " | ")
  }, error = function(e) paste("gam.check() failed:", conditionMessage(e)))

  data.frame(
    model_id = model_id,
    framework = "mgcv::gam",
    converged = conv,
    optimizer_convergence_code = NA_integer_,
    optimizer_message = if (conv) "full convergence" else "did not report full convergence",
    hessian_positive_definite = pdHess,
    max_abs_gradient = grad_range,
    gam_iterations = n_iter,
    gam_optimizer = opt_name,
    gam_check_text = k_check_text,
    stringsAsFactors = FALSE
  )
}

#' Information criteria for a fitted model. For glmmTMB, AIC()/BIC() report
#' whatever criterion the fit was estimated under (REML or ML); for mgcv,
#' AIC() uses the usual marginal AIC approximation. Callers should refit
#' under ML when they want a fair cross-model AIC/BIC comparison across
#' different fixed-effect structures (REML AIC is not strictly comparable
#' across models with different fixed effects).
extract_ic <- function(fit, model_id, estimation_method) {
  data.frame(
    model_id = model_id,
    estimation_method = estimation_method,
    aic = as.numeric(stats::AIC(fit)),
    bic = as.numeric(stats::BIC(fit)),
    logLik = as.numeric(stats::logLik(fit)),
    n_params_reported = attr(stats::logLik(fit), "df"),
    stringsAsFactors = FALSE
  )
}
