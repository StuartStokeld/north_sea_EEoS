# OLS model helpers for H2
# Requires h2_common.R

#' Fit a simple OLS model and return a one-row coefficient summary.
fit_h2_ols <- function(data,
                       formula,
                       model_id,
                       term = "mean_annual_hours_total") {
  fit <- lm(formula, data = data)
  sm <- summary(fit)
  ci <- confint(fit)
  coef_name <- term
  if (!coef_name %in% rownames(sm$coefficients)) {
    stop("Term not found in model: ", coef_name)
  }

  data.frame(
    model_id = model_id,
    term = coef_name,
    estimate = sm$coefficients[coef_name, "Estimate"],
    std_error = sm$coefficients[coef_name, "Std. Error"],
    statistic = sm$coefficients[coef_name, "t value"],
    p_value = sm$coefficients[coef_name, "Pr(>|t|)"],
    ci_low = ci[coef_name, 1L],
    ci_high = ci[coef_name, 2L],
    r_squared = sm$r.squared,
    adj_r_squared = sm$adj.r.squared,
    n = nrow(model.frame(fit)),
    stringsAsFactors = FALSE
  )
}

#' Fit the pre-specified H2 model set for one analysis panel.
fit_h2_model_set <- function(panel, min_hauls_label = "default") {
  models <- list(
    fit_h2_ols(
      panel,
      mean_abs_residual ~ mean_annual_hours_total,
      paste0("primary_abs_", min_hauls_label)
    ),
    fit_h2_ols(
      panel,
      mean_residual ~ mean_annual_hours_total,
      paste0("secondary_signed_", min_hauls_label)
    ),
    fit_h2_ols(
      panel,
      mean_abs_residual ~ mean_annual_hours_otter,
      paste0("sensitivity_otter_", min_hauls_label),
      term = "mean_annual_hours_otter"
    ),
    fit_h2_ols(
      panel,
      mean_abs_residual ~ mean_annual_hours_beam,
      paste0("sensitivity_beam_", min_hauls_label),
      term = "mean_annual_hours_beam"
    ),
    fit_h2_ols(
      panel,
      mean_abs_residual ~ log_mean_annual_hours_total,
      paste0("sensitivity_log_hours_", min_hauls_label),
      term = "log_mean_annual_hours_total"
    ),
    fit_h2_ols(
      panel,
      mean_abs_residual ~ mean_annual_hours_total + mean_ln_B_obs,
      paste0("sensitivity_biomass_covariate_", min_hauls_label)
    )
  )

  bind_rows(models)
}
