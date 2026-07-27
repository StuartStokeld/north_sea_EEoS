# Helpers for the H2/H3 shared hierarchical model FEASIBILITY CHECK
# (run_h2h3_shared_model_feasibility.R). This is a feasibility check, not a
# results run — nothing here interprets the fishing-pressure or
# fishing-pressure x phase coefficients as answering H2/H3.

suppressPackageStartupMessages({
  library(dplyr)
})

#' 4-level phase factor from the 3 statistically robust structural breaks
#' (1989, 2001, 2008; see outputs/h2h3_designA4_structbreak_years.csv). The
#' marginal 1997 break is excluded from this variable per the brief.
#' CONVENTION (a modelling choice, not verified against strucchange's exact
#' segment-boundary indexing): each break year is treated as the FIRST year
#' of the new phase, e.g. phase 2 runs 1989-2000 inclusive. Not re-derived
#' from strucchange's internal breakpoint indexing here.
build_phase_factor <- function(year) {
  breaks <- c(-Inf, 1989, 2001, 2008, Inf)
  labels <- c("1985-1988", "1989-2000", "2001-2007", "2008-2015")
  cut(year, breaks = breaks, labels = labels, right = FALSE)
}

#' Build the haul-level analysis dataset for the feasibility model.
#' Universe = the ESTABLISHED H2 analysis panel's 158 rectangles
#' (outputs/h2_rectangle_panel.rds: n_hauls >= H2_MIN_HAULS_DEFAULT AND
#' has Couce fishing-pressure coverage) — the same universe the original H2
#' SEM/SAR analysis used. No new rectangle-inclusion rule is introduced here.
build_feasibility_data <- function(haul, panel, couce_year, year_min, year_max) {
  required_haul <- c("stat_rec", "year", "residual")
  missing_haul <- setdiff(required_haul, names(haul))
  if (length(missing_haul) > 0L) {
    stop("haul predictions missing columns: ", paste(missing_haul, collapse = ", "))
  }
  required_panel <- c("stat_rec", "mean_ln_B_obs", "rect_lon", "rect_lat")
  missing_panel <- setdiff(required_panel, names(panel))
  if (length(missing_panel) > 0L) {
    stop("h2_rectangle_panel missing columns: ", paste(missing_panel, collapse = ", "))
  }

  h <- haul %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(
      stat_rec %in% panel$stat_rec,
      is.finite(residual),
      year >= year_min,
      year <= year_max
    )
  n_after_universe <- nrow(h)

  h <- h %>%
    inner_join(
      couce_year %>%
        mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
        select(stat_rec, year, hours_total),
      by = c("stat_rec", "year")
    )
  n_after_couce_year_join <- nrow(h)

  h <- h %>%
    inner_join(
      panel %>% select(stat_rec, mean_ln_B_obs, rect_lon, rect_lat),
      by = "stat_rec"
    ) %>%
    mutate(
      log_hours_total = log(hours_total + 1),
      phase = build_phase_factor(year),
      pos = numFactor(rect_lon, rect_lat),
      dummy = factor(1L)
    )

  list(
    data = h,
    n_hauls_in_h2_universe = n_after_universe,
    n_hauls_dropped_no_rect_year_couce = n_after_universe - n_after_couce_year_join,
    n_hauls_final = nrow(h)
  )
}

#' Tidy fixed-effects table (estimate, SE, statistic, p-value) from a
#' glmmTMB fit, tagged with model_id.
tidy_fixed_effects <- function(fit, model_id) {
  s <- summary(fit)
  coef_mat <- s$coefficients$cond
  data.frame(
    model_id = model_id,
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, "z value"],
    p_value = coef_mat[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

#' Convergence diagnostics for a glmmTMB fit: optimizer convergence code,
#' Hessian positive-definiteness, max abs gradient, and any NA/singular
#' variance-component flags. Does not force convergence by simplifying
#' the model — reports whatever the fit produced.
glmmtmb_convergence_report <- function(fit, model_id) {
  conv_code <- fit$fit$convergence
  conv_msg <- fit$fit$message
  pdHess <- isTRUE(fit$sdr$pdHess)
  grad <- fit$sdr$gradient.fixed
  max_abs_grad <- if (!is.null(grad)) max(abs(grad), na.rm = TRUE) else NA_real_
  nas_in_se <- any(is.na(sqrt(diag(fit$sdr$cov.fixed))))

  diagnose_out <- tryCatch(
    {
      out <- capture.output(glmmTMB::diagnose(fit))
      paste(out, collapse = " | ")
    },
    error = function(e) paste("diagnose() failed:", conditionMessage(e)),
    warning = function(w) paste("diagnose() warning:", conditionMessage(w))
  )

  data.frame(
    model_id = model_id,
    optimizer_convergence_code = conv_code,
    optimizer_message = ifelse(is.null(conv_msg) || conv_msg == "", "none", conv_msg),
    hessian_positive_definite = pdHess,
    max_abs_gradient = max_abs_grad,
    any_na_std_error = nas_in_se,
    diagnose_output = diagnose_out,
    stringsAsFactors = FALSE
  )
}

#' Variance-components table. For the spatial model, the exp() structure's
#' internal parameters are theta1 = log(spatial SD), theta2 = log(range);
#' recovered via fit$sdr$par.fixed / cov.fixed and reported on the natural
#' scale via a delta-method (log-normal) 95% CI. For the plain model,
#' reports the (1|stat_rec) SD via VarCorr with its glmmTMB-reported CI.
extract_variance_components <- function(fit, model_id, spatial = TRUE) {
  vc <- VarCorr(fit)
  resid_sd <- attr(vc$cond, "sc")
  if (is.null(resid_sd)) resid_sd <- sigma(fit)

  if (spatial) {
    pf <- fit$sdr$par.fixed
    theta_idx <- which(names(pf) == "theta")
    if (length(theta_idx) < 2L) {
      stop("Expected 2 theta parameters (log sd, log range) for exp() structure; found ", length(theta_idx))
    }
    theta_idx <- theta_idx[1:2]
    se <- sqrt(diag(fit$sdr$cov.fixed))[theta_idx]
    log_sd <- pf[theta_idx[1]]
    log_range <- pf[theta_idx[2]]
    se_log_sd <- se[1]
    se_log_range <- se[2]

    data.frame(
      model_id = model_id,
      component = c("rectangle_spatial_sd", "spatial_range"),
      estimate = c(exp(log_sd), exp(log_range)),
      ci_low = c(exp(log_sd - 1.96 * se_log_sd), exp(log_range - 1.96 * se_log_range)),
      ci_high = c(exp(log_sd + 1.96 * se_log_sd), exp(log_range + 1.96 * se_log_range)),
      se_on_log_scale = c(se_log_sd, se_log_range),
      residual_sd = resid_sd,
      stringsAsFactors = FALSE
    )
  } else {
    sd_df <- as.data.frame(vc$cond$stat_rec)
    rect_sd <- attr(vc$cond$stat_rec, "stddev")[1]
    ci <- tryCatch(
      confint(fit, parm = "theta_"),
      error = function(e) matrix(NA_real_, nrow = 1L, ncol = 2L, dimnames = list("Std.Dev.(Intercept)|stat_rec", c("2.5 %", "97.5 %")))
    )
    row_match <- grep("stat_rec", rownames(ci), fixed = TRUE)
    ci_low <- if (length(row_match) > 0L) ci[row_match[1], 1] else NA_real_
    ci_high <- if (length(row_match) > 0L) ci[row_match[1], 2] else NA_real_

    data.frame(
      model_id = model_id,
      component = "rectangle_intercept_sd",
      estimate = rect_sd,
      ci_low = ci_low,
      ci_high = ci_high,
      se_on_log_scale = NA_real_,
      residual_sd = resid_sd,
      stringsAsFactors = FALSE
    )
  }
}

#' Partial-pooling table: per-rectangle random-intercept estimate joined to
#' haul count, PLUS an UNPOOLED reference point (the per-rectangle mean of
#' the "partial residual" = observed residual minus the model's own
#' population-level [fixed-effects-only, re.form=NA] prediction — i.e. what
#' each rectangle's intercept would be if fit with NO pooling at all, using
#' the SAME fixed effects as the mixed model). shrinkage_ratio =
#' random_intercept / unpooled_intercept is therefore on a comparable scale
#' (both are deviations from the fixed-effects prediction). Expected
#' pattern: shrinkage_ratio closer to 0 (heavy shrinkage toward the
#' population mean) for low-haul-count rectangles, closer to 1 (little
#' shrinkage) for high-haul-count rectangles. Works for both the spatial
#' (`exp(pos+0|dummy)`) and plain (`(1|stat_rec)`) random-effect structures.
extract_partial_pooling <- function(fit, data, model_id, spatial = TRUE) {
  partial_resid <- data$residual - stats::predict(fit, newdata = data, re.form = NA)
  raw_means <- data %>%
    mutate(.partial_resid = partial_resid) %>%
    group_by(stat_rec) %>%
    summarise(n_hauls = dplyr::n(), unpooled_intercept = mean(.partial_resid, na.rm = TRUE), .groups = "drop")

  if (spatial) {
    re <- ranef(fit)$cond$dummy
    lookup <- data %>%
      distinct(stat_rec, pos) %>%
      mutate(pos_label = paste0("pos", as.character(pos)))
    values <- as.numeric(re[1, ])
    re_df <- data.frame(
      pos_label = colnames(re),
      random_intercept = values,
      stringsAsFactors = FALSE
    )
    out <- lookup %>%
      inner_join(re_df, by = "pos_label") %>%
      select(stat_rec, random_intercept)
  } else {
    re <- ranef(fit)$cond$stat_rec
    out <- data.frame(
      stat_rec = rownames(re),
      random_intercept = re[, "(Intercept)"],
      stringsAsFactors = FALSE
    )
  }

  out %>%
    inner_join(raw_means, by = "stat_rec") %>%
    mutate(
      model_id = model_id,
      shrinkage_ratio = random_intercept / unpooled_intercept
    ) %>%
    arrange(n_hauls)
}

# ---------------------------------------------------------------------------
# Round 2: adjacency-based (CAR) spatial structure via spaMM::fitme()
# ---------------------------------------------------------------------------

#' Build the binary (0/1) queen-adjacency matrix for spaMM's adjacency() CAR
#' term FROM THE REUSED spdep neighbour list (`weights$nb`, as returned by
#' build_h2_spatial_weights() — the exact same function/parameters used by
#' the original H2 SEM/SAR analysis). This does NOT rebuild the adjacency
#' RULE: `nb` is the identical neighbour-list object; only the matrix
#' *representation* differs from the row-standardised "W"-style `listw` used
#' by errorsarlm()/lagsarlm() (style = "W"), because a proper CAR model
#' (spaMM's adjacency() family, and CAR models generally) is conventionally
#' parameterised on the raw symmetric binary adjacency matrix style = "B",
#' with the correlation strength captured by a single free `rho`, not by
#' row-standardising the matrix itself. Row order is fixed to `panel$stat_rec`
#' so it aligns with the analysis data's `stat_rec` factor levels.
build_car_adjacency_matrix <- function(nb, panel) {
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  rownames(W) <- colnames(W) <- panel$stat_rec
  W
}

#' Admissible range for the CAR rho parameter: (1/lambda_min, 1/lambda_max)
#' of the adjacency matrix's eigenvalues (standard CAR identifiability
#' bound; confirmed against spaMM's own reported range on a toy example).
car_rho_admissible_range <- function(W) {
  eig <- eigen(W, only.values = TRUE, symmetric = TRUE)$values
  c(lower = 1 / min(eig), upper = 1 / max(eig))
}

#' Tidy fixed-effects table from a spaMM fitme() fit, tagged with model_id.
#' spaMM's default summary() does not report a p-value column (only
#' Estimate, Cond. SE, t-value) for a gaussian-family HLfit; a Wald normal-
#' approximation p-value is added here (2 * pnorm(-|t|)) for direct
#' comparability with the glmmTMB z-value p-values used elsewhere in this
#' feasibility check (both are large-sample Wald approximations, no formal
#' finite-sample t-distribution df adjustment in either case).
tidy_fixed_effects_spamm <- function(fit, model_id) {
  bt <- summary(fit, verbose = FALSE)$beta_table
  data.frame(
    model_id = model_id,
    term = rownames(bt),
    estimate = bt[, "Estimate"],
    std_error = bt[, "Cond. SE"],
    statistic = bt[, "t-value"],
    p_value = 2 * stats::pnorm(-abs(bt[, "t-value"])),
    stringsAsFactors = FALSE
  )
}

#' Convergence / robustness report for the spaMM CAR fit. spaMM's fitme()
#' does not expose an iterative-optimizer status code the way glmmTMB does;
#' instead this reports (a) any warnings/errors raised during the fit
#' (fit$warnings; empty list = none raised), (b) whether refitting from
#' several different starting values for rho (near the lower admissible
#' bound, zero, near the upper admissible bound) converges to the SAME
#' fitted rho (a standard robustness check for local-optimum issues when no
#' single "convergence code" is available), and (c) whether all fixed-effect
#' standard errors are finite.
spamm_car_convergence_report <- function(fit, model_id, formula_car, data, adjMatrix, rhorange) {
  warn_list <- fit$warnings
  n_warnings <- length(warn_list)

  init_rhos <- c(
    near_lower = unname(rhorange["lower"]) * 0.9,
    zero = 0,
    near_upper = unname(rhorange["upper"]) * 0.9
  )
  refit_rhos <- vapply(init_rhos, function(r0) {
    f <- tryCatch(
      spaMM::fitme(
        formula_car, data = data, adjMatrix = adjMatrix, method = "REML",
        init = list(corrPars = list("1" = list(rho = r0)))
      ),
      error = function(e) NULL
    )
    if (is.null(f)) return(NA_real_)
    spaMM::get_ranPars(f, which = "corrPars")[["1"]]$rho
  }, numeric(1))
  fitted_rho <- spaMM::get_ranPars(fit, which = "corrPars")[["1"]]$rho
  robust_to_start <- all(is.finite(refit_rhos)) && all(abs(refit_rhos - fitted_rho) < 1e-4)

  bt <- summary(fit, verbose = FALSE)$beta_table
  any_na_se <- any(!is.finite(bt[, "Cond. SE"]))

  data.frame(
    model_id = model_id,
    n_warnings_during_fit = n_warnings,
    warnings_text = if (n_warnings > 0L) paste(unlist(warn_list), collapse = " | ") else "none",
    refit_near_lower_bound_rho = unname(refit_rhos["near_lower"]),
    refit_zero_start_rho = unname(refit_rhos["zero"]),
    refit_near_upper_bound_rho = unname(refit_rhos["near_upper"]),
    fitted_rho = fitted_rho,
    robust_to_starting_value = robust_to_start,
    any_na_or_infinite_fixed_effect_se = any_na_se,
    stringsAsFactors = FALSE
  )
}

#' CAR spatial-parameter table: rho estimate, admissible range, parametric
#' bootstrap 95% CI (via spaMM::confint with parm = a function extracting
#' rho, since profile-likelihood confint() in spaMM only covers fixed-effect
#' coefficients directly), and a boundary-pinning flag (rho within 1% of
#' either admissible bound). Also reports the CAR variance component
#' ('lambda').
extract_car_spatial_param <- function(fit, adjMatrix, rhorange, model_id, nsim = 99, seed = 123) {
  fitted_rho <- spaMM::get_ranPars(fit, which = "corrPars")[["1"]]$rho
  lambda_est <- fit$lambda[[1]]

  # NOTE: must call the bare S3 generic confint() (dispatching to spaMM's
  # confint.HLfit method), NOT spaMM::confint() — spaMM does not export a
  # `confint` object itself (only the `confint.HLfit` method), so
  # `spaMM::confint(...)` fails with "not an exported object" and would
  # silently be swallowed by tryCatch below without this fix.
  boot_out <- tryCatch(
    stats::confint(
      fit,
      parm = function(f) spaMM::get_ranPars(f, which = "corrPars")[["1"]]$rho,
      boot_args = list(nsim = nsim, seed = seed, nb_cores = 1),
      verbose = FALSE
    ),
    error = function(e) {
      warning("CAR rho bootstrap CI failed: ", conditionMessage(e))
      NULL
    }
  )
  # attr(boot_out, "table") is a NAMED LIST of 1x2 matrices, one per
  # boot.ci() interval type ("normal", "percent", "basic"); percentile is
  # used as the primary (most standard, no normality assumption) interval.
  ci_tab <- if (!is.null(boot_out)) attr(boot_out, "table") else NULL
  if (!is.null(ci_tab) && "percent" %in% names(ci_tab)) {
    ci_low <- ci_tab$percent[1, 1]
    ci_high <- ci_tab$percent[1, 2]
  } else if (!is.null(ci_tab) && length(ci_tab) > 0L) {
    ci_low <- ci_tab[[1]][1, 1]
    ci_high <- ci_tab[[1]][1, 2]
  } else {
    ci_low <- NA_real_
    ci_high <- NA_real_
  }

  range_width <- rhorange["upper"] - rhorange["lower"]
  dist_to_lower <- (fitted_rho - rhorange["lower"]) / range_width
  dist_to_upper <- (rhorange["upper"] - fitted_rho) / range_width
  pinned <- dist_to_lower < 0.01 || dist_to_upper < 0.01

  data.frame(
    model_id = model_id,
    rho_estimate = fitted_rho,
    rho_ci_low_boot = ci_low,
    rho_ci_high_boot = ci_high,
    rho_admissible_lower = unname(rhorange["lower"]),
    rho_admissible_upper = unname(rhorange["upper"]),
    fraction_of_range_from_lower = dist_to_lower,
    fraction_of_range_from_upper = dist_to_upper,
    boundary_pinned_flag = pinned,
    car_lambda_variance = lambda_est,
    nsim_bootstrap = nsim,
    stringsAsFactors = FALSE
  )
}

#' Partial-pooling table for the spaMM CAR fit, on the SAME
#' unpooled_intercept / shrinkage_ratio definition as extract_partial_pooling()
#' (Round 1), for direct comparability: unpooled_intercept = per-rectangle
#' mean of (observed residual - population-level [re.form = NA] prediction).
extract_partial_pooling_spamm <- function(fit, data, model_id) {
  partial_resid <- data$residual - as.numeric(stats::predict(fit, newdata = data, re.form = NA))
  raw_means <- data %>%
    mutate(.partial_resid = partial_resid) %>%
    group_by(stat_rec) %>%
    summarise(n_hauls = dplyr::n(), unpooled_intercept = mean(.partial_resid, na.rm = TRUE), .groups = "drop")

  re <- spaMM::ranef(fit)[[1]]
  re_df <- data.frame(
    stat_rec = names(re),
    random_intercept = as.numeric(re),
    stringsAsFactors = FALSE
  )

  re_df %>%
    inner_join(raw_means, by = "stat_rec") %>%
    mutate(
      model_id = model_id,
      shrinkage_ratio = random_intercept / unpooled_intercept
    ) %>%
    arrange(n_hauls)
}
