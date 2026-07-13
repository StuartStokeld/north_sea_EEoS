# ln(E) reference model helpers (Harte et al. 2022 benchmark)
# Requires h1_common.R to be sourced first.

EEOS_MODEL <- "EEoS (S, N, E)"
LNE_MODEL <- "ln(E) OLS (fitted correlative)"

#' Fit haul-level ln(B_obs) ~ ln(E) reference model.
fit_ln_e_reference <- function(B_obs_g, E_norm) {
  df <- data.frame(
    log_B_obs = log(B_obs_g),
    log_E = log(E_norm)
  )
  df <- df[is.finite(df$log_B_obs) & is.finite(df$log_E), , drop = FALSE]
  fit <- lm(log_B_obs ~ log_E, data = df)
  list(
    fit = fit,
    coef = coef(summary(fit)),
    r_squared = summary(fit)$r.squared,
    adj_r_squared = summary(fit)$adj.r.squared,
    sigma = summary(fit)$sigma,
    n = nrow(df)
  )
}

#' Build haul-level table with ln(E) predictions and residuals.
augment_haul_ln_e <- function(haul) {
  if (!"ln_B_obs" %in% names(haul)) {
    haul$ln_B_obs <- log(haul$B_obs)
  }
  fit_info <- fit_ln_e_reference(haul$B_obs, haul$E)
  ln_B_pred_E <- predict(fit_info$fit, newdata = data.frame(log_E = log(haul$E)))
  haul$ln_B_pred_E <- ln_B_pred_E
  haul$B_pred_lnE <- exp(ln_B_pred_E)
  haul$residual_lnE <- haul$ln_B_obs - ln_B_pred_E
  haul$abs_residual_lnE <- abs(haul$residual_lnE)
  list(haul = haul, fit = fit_info)
}

#' One model row for the unified comparison table.
model_metrics_row <- function(
    model, model_id, fitted, tier, comparison_type,
    obs, pred, log_scale = TRUE) {
  m <- evaluate_prediction(obs, pred, log_scale = log_scale)
  m$fig2_r2_pearson <- NA_real_
  m$fig2_r2_pearson_all <- NA_real_
  m$fig2_r2_pearson_trimmed <- NA_real_
  m$fig2_r2_pearson_n_excluded <- NA_integer_
  m$fig2_r2_cod_extended <- NA_real_
  cbind(
    data.frame(
      model = model,
      model_id = model_id,
      fitted = fitted,
      tier = tier,
      comparison_type = comparison_type,
      stringsAsFactors = FALSE
    ),
    m
  )
}

#' Fig 2 model row (Harte Pearson R² + extended cod diagnostic).
fig2_metrics_row <- function(haul) {
  m <- harte_fig2_metrics(haul$fig2_predicted_ratio, haul$fig2_observed_ratio)
  cbind(
    data.frame(
      model = "Productivity ratio (Harte Fig 2)",
      model_id = "productivity_ratio",
      fitted = FALSE,
      tier = 1L,
      comparison_type = "E/B^(3/4) predicted vs observed (normalised E)",
      stringsAsFactors = FALSE
    ),
    m
  )
}

#' Unified H1 model comparison (Tier 1 unfitted + Tier 2 fitted ln(E)).
#'
#' Headline productivity_1to1 uses E × m_min (same E as biomass(), same m_min
#' as B_pred grams conversion). productivity_1to1_uncalibrated (E_raw) is
#' retained as a catchability / unit diagnostic only — not the Test 1 baseline.
compare_all_h1_models <- function(haul) {
  rows <- list(
    model_metrics_row(
      EEOS_MODEL, "eeos_biomass", FALSE, 1L,
      "log B_pred vs log B_obs",
      haul$B_obs, haul$B_pred
    ),
    model_metrics_row(
      "Productivity 1:1 (log(E × m_min) vs log B_obs)",
      "productivity_1to1", FALSE, 1L,
      "unfitted 1:1 productivity map (E_norm × m_min)",
      haul$B_obs, haul$E_calibrated
    ),
    fig2_metrics_row(haul),
    model_metrics_row(
      "Productivity 1:1 uncalibrated (log E_raw vs log B_obs)",
      "productivity_1to1_uncalibrated", FALSE, 1L,
      "diagnostic: uncalibrated E_raw 1:1 (not headline)",
      haul$B_obs, haul$E_raw
    )
  )

  if ("ln_B_pred_E" %in% names(haul)) {
    rows[[length(rows) + 1L]] <- model_metrics_row(
      LNE_MODEL, "ln_e_ols", TRUE, 2L,
      "fitted lm(log B_obs ~ log E)",
      haul$B_obs, haul$B_pred_lnE
    )
  }

  do.call(rbind, rows)
}

#' Compare EEoS vs ln(E) reference on log-scale predictive metrics.
compare_eeos_vs_lne <- function(haul) {
  cmp <- compare_all_h1_models(haul)
  cmp$ss_ratio_vs_eeos <- cmp$ss_res_log / cmp$ss_res_log[cmp$model_id == "eeos_biomass"]
  cmp$harte_criterion_met <- cmp$ss_res_log < 0.5 * cmp$ss_res_log[cmp$model_id == "eeos_biomass"]
  cmp$harte_criterion_met[cmp$model_id == "eeos_biomass"] <- NA
  cmp
}
