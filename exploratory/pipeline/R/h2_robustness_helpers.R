# H2 robustness diagnostics (leverage, LM tests, spatial sensitivity)
# Requires h2_spatial_helpers.R, h2_model_helpers.R

suppressPackageStartupMessages({
  library(dplyr)
})

PRIMARY_OLS_FORMULA <- mean_abs_residual ~ mean_annual_hours_total
SEM_PRIMARY_FORMULA <- mean_abs_residual ~ mean_annual_hours_total
SEM_BIOMASS_FORMULA <- mean_abs_residual ~ mean_annual_hours_total + mean_ln_B_obs

#' Run Moran's I / Geary's C + SEM on one analysis panel with fresh weights.
h2_spatial_panel_summary <- function(panel,
                                     rectangles_sf,
                                     project_root,
                                     panel_label,
                                     ols_formula = PRIMARY_OLS_FORMULA,
                                     sem_formula = SEM_PRIMARY_FORMULA,
                                     sem_model_id = "sem_primary_abs") {
  weights <- build_h2_spatial_weights(panel, rectangles_sf, project_root)
  ols_fit <- lm(ols_formula, data = panel)

  spatial_resid <- h2_global_spatial_tests(residuals(ols_fit), weights$listw) %>%
    mutate(
      panel = panel_label,
      variable = "ols_primary_abs_residuals",
      n_rectangles = nrow(panel),
      n_isolated = weights$n_isolated
    )

  sem_coef <- fit_h2_sem(
    panel,
    weights$listw,
    formula = sem_formula,
    model_id = sem_model_id
  ) %>%
    mutate(panel = panel_label)

  list(
    panel = panel_label,
    n_rectangles = nrow(panel),
    n_isolated = weights$n_isolated,
    spatial_diagnostics = spatial_resid,
    sem_results = sem_coef,
    ols_fit = ols_fit,
    weights = weights
  )
}

#' Flatten spatial sensitivity results to one summary table.
h2_spatial_sensitivity_table <- function(spatial_runs) {
  bind_rows(lapply(spatial_runs, function(run) {
    moran <- run$spatial_diagnostics %>%
      filter(test == "morans_i") %>%
      slice(1)
    geary <- run$spatial_diagnostics %>%
      filter(test == "geary_c") %>%
      slice(1)
    sem_hours <- run$sem_results %>%
      filter(term == "mean_annual_hours_total") %>%
      slice(1)

    data.frame(
      panel = run$panel,
      n_rectangles = run$n_rectangles,
      n_isolated = run$n_isolated,
      morans_i = moran$statistic,
      morans_i_p = moran$p_value,
      geary_c = geary$statistic,
      geary_c_p = geary$p_value,
      sem_beta_hours = sem_hours$estimate,
      sem_p_hours = sem_hours$p_value,
      sem_lambda = sem_hours$lambda[1L],
      stringsAsFactors = FALSE
    )
  }))
}

#' Cook's distance and leverage for primary OLS; refit excluding flagged points.
h2_leverage_diagnostics <- function(panel,
                                    ols_formula = PRIMARY_OLS_FORMULA,
                                    cook_threshold_multiplier = 4) {
  fit <- lm(ols_formula, data = panel)
  n <- nrow(panel)
  cooks <- stats::cooks.distance(fit)
  leverage <- stats::hatvalues(fit)
  cook_threshold <- cook_threshold_multiplier / n

  diag_df <- panel %>%
    mutate(
      stat_rec = stat_rec,
      cooks_distance = as.numeric(cooks),
      leverage = as.numeric(leverage),
      flagged_cooks = cooks_distance > cook_threshold
    ) %>%
    select(
      stat_rec,
      n_hauls,
      mean_abs_residual,
      mean_annual_hours_total,
      cooks_distance,
      leverage,
      flagged_cooks
    ) %>%
    arrange(desc(cooks_distance))

  flagged <- diag_df %>% filter(flagged_cooks)
  panel_trim <- panel %>% filter(!stat_rec %in% flagged$stat_rec)

  refit <- if (nrow(panel_trim) >= 3L && nrow(flagged) > 0L) {
    fit_h2_ols(
      panel_trim,
      ols_formula,
      "primary_abs_min10_no_high_leverage",
      term = "mean_annual_hours_total"
    )
  } else {
    NULL
  }

  full <- fit_h2_ols(
    panel,
    ols_formula,
    "primary_abs_min10_full",
    term = "mean_annual_hours_total"
  )

  list(
    diagnostics = diag_df,
    flagged_rectangles = flagged,
    cook_threshold = cook_threshold,
    full_panel = full,
    trimmed_panel = refit,
    n_flagged = nrow(flagged),
    n_trimmed = nrow(panel_trim)
  )
}

#' Lagrange Multiplier tests for spatial model choice on OLS residuals.
h2_lm_spatial_tests <- function(ols_fit, listw) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for lm.LMtests().")
  }

  tests <- spdep::lm.LMtests(
    ols_fit,
    listw,
    test = c("LMerr", "LMlag", "RLMerr", "RLMlag")
  )

  name_map <- c(
    LMerr = "LM_error",
    LMlag = "LM_lag",
    RLMerr = "RLM_error",
    RLMlag = "RLM_lag",
    RSerr = "LM_error",
    RSlag = "LM_lag",
    adjRSerr = "RLM_error",
    adjRSlag = "RLM_lag"
  )

  raw_names <- names(tests)
  pretty_names <- unname(name_map[raw_names])
  pretty_names[is.na(pretty_names)] <- raw_names[is.na(pretty_names)]

  data.frame(
    test = pretty_names,
    test_raw = raw_names,
    statistic = vapply(tests, function(x) unname(x$statistic), numeric(1)),
    p_value = vapply(tests, function(x) unname(x$p.value), numeric(1)),
    df = vapply(tests, function(x) unname(x$parameter), numeric(1)),
    stringsAsFactors = FALSE
  )
}
