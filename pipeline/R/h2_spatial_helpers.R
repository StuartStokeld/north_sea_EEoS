# Spatial autocorrelation and spatial error models for H2
# Requires h2_common.R, sf

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

#' Load ICES rectangles and attach stat_rec from ICESNAME.
load_ices_rectangles_sf <- function(project_root) {
  path <- h2_ices_shapefile_path(project_root)
  stopifnot(file.exists(path))
  st_read(path, quiet = TRUE) %>%
    mutate(stat_rec = normalize_stat_rec(ICESNAME))
}

#' Build queen contiguity weights for analysis rectangles.
build_h2_spatial_weights <- function(panel, rectangles_sf, project_root) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for spatial analysis. Install via renv::install('spdep').")
  }

  panel_sf <- rectangles_sf %>%
    filter(stat_rec %in% panel$stat_rec) %>%
    select(stat_rec) %>%
    inner_join(panel, by = "stat_rec")

  if (nrow(panel_sf) != nrow(panel)) {
    stop("Could not match all analysis rectangles to ICES shapefile.")
  }

  # Preserve panel row order for aligning model residuals.
  panel_sf <- panel_sf[match(panel$stat_rec, panel_sf$stat_rec), , drop = FALSE]

  nb <- spdep::poly2nb(panel_sf, queen = TRUE)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  list(
    panel_sf = panel_sf,
    nb = nb,
    listw = listw,
    n_isolated = sum(spdep::card(nb) == 0L)
  )
}

#' Global Moran's I and Geary's C on a numeric vector aligned to weights.
h2_global_spatial_tests <- function(x, listw) {
  x <- as.numeric(x)
  moran <- spdep::moran.test(x, listw, zero.policy = TRUE)
  geary <- spdep::geary.test(x, listw, zero.policy = TRUE)

  data.frame(
    test = c("morans_i", "geary_c"),
    statistic = c(unname(moran$estimate[["Moran I statistic"]]), geary$estimate[["Geary C statistic"]]),
    expected = c(unname(moran$estimate[["Expectation"]]), 1),
    p_value = c(moran$p.value, geary$p.value),
    stringsAsFactors = FALSE
  )
}

#' Fit spatial error model; returns one row per coefficient term.
fit_h2_sem <- function(panel,
                       listw,
                       formula = mean_abs_residual ~ mean_annual_hours_total,
                       model_id = "sem_primary_abs") {
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("Package 'spatialreg' is required. Install via renv::install('spatialreg').")
  }

  fit <- spatialreg::errorsarlm(
    formula,
    data = panel,
    listw = listw,
    zero.policy = TRUE
  )

  coef_mat <- coef(summary(fit))
  ci <- confint(fit)
  terms <- setdiff(rownames(coef_mat), "(Intercept)")
  lambda_val <- unname(fit$lambda)

  bind_rows(lapply(terms, function(term) {
    data.frame(
      model_id = model_id,
      term = term,
      estimate = coef_mat[term, "Estimate"],
      std_error = coef_mat[term, "Std. Error"],
      statistic = coef_mat[term, "z value"],
      p_value = coef_mat[term, "Pr(>|z|)"],
      ci_low = ci[term, 1L],
      ci_high = ci[term, 2L],
      lambda = lambda_val,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      n = nrow(panel),
      stringsAsFactors = FALSE
    )
  }))
}

#' Fit spatial lag (SAR) model; returns coefficient rows plus model fit summary.
fit_h2_sar_lag <- function(panel,
                           listw,
                           formula = mean_abs_residual ~ mean_annual_hours_total,
                           model_id = "sar_lag_primary_abs") {
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("Package 'spatialreg' is required. Install via renv::install('spatialreg').")
  }

  fit <- spatialreg::lagsarlm(
    formula,
    data = panel,
    listw = listw,
    zero.policy = TRUE
  )

  sm <- summary(fit)
  coef_mat <- coef(sm)
  ci <- confint(fit)
  terms <- setdiff(rownames(coef_mat), "(Intercept)")
  rho_val <- unname(fit$rho)
  rho_se <- unname(fit$rho.se)
  rho_z <- if (is.finite(rho_se) && rho_se > 0) rho_val / rho_se else NA_real_
  rho_p <- if (is.finite(rho_z)) 2 * stats::pnorm(-abs(rho_z)) else NA_real_

  coef_rows <- bind_rows(lapply(terms, function(term) {
    data.frame(
      model_id = model_id,
      term = term,
      estimate = coef_mat[term, "Estimate"],
      std_error = coef_mat[term, "Std. Error"],
      statistic = coef_mat[term, "z value"],
      p_value = coef_mat[term, "Pr(>|z|)"],
      ci_low = ci[term, 1L],
      ci_high = ci[term, 2L],
      rho = rho_val,
      rho_se = rho_se,
      rho_z = rho_z,
      rho_p_value = rho_p,
      log_likelihood = as.numeric(logLik(fit)),
      aic = AIC(fit),
      n = nrow(panel),
      stringsAsFactors = FALSE
    )
  }))

  coef_rows
}

#' Extract SEM log-likelihood and AIC by refitting (read-only comparison helper).
fit_h2_sem_fit_stats <- function(panel,
                                 listw,
                                 formula = mean_abs_residual ~ mean_annual_hours_total,
                                 model_id = "sem_primary_abs") {
  fit <- spatialreg::errorsarlm(
    formula,
    data = panel,
    listw = listw,
    zero.policy = TRUE
  )
  data.frame(
    model_id = model_id,
    spatial_spec = "sem",
    log_likelihood = as.numeric(logLik(fit)),
    aic = AIC(fit),
    spatial_param = unname(fit$lambda),
    spatial_param_name = "lambda",
    n = nrow(panel),
    stringsAsFactors = FALSE
  )
}

#' Save H2 diagnostic and map figures.
save_h2_figures <- function(panel, panel_sf, weights, ols_fit, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping H2 figures.")
    return(invisible(NULL))
  }

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  p_scatter <- ggplot2::ggplot(
    panel,
    ggplot2::aes(x = mean_annual_hours_total, y = mean_abs_residual)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#b2182b") +
    ggplot2::labs(
      x = "Mean annual fishing hours (1985–2015)",
      y = "Mean |EEoS residual| per rectangle",
      title = "H2: Fishing pressure vs EEoS residual magnitude"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    file.path(fig_dir, "h2_scatter_pressure_vs_abs_residual.png"),
    p_scatter,
    width = 8,
    height = 5.5,
    dpi = 150
  )

  p_map_dv <- ggplot2::ggplot(panel_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = mean_abs_residual), colour = NA) +
    ggplot2::scale_fill_viridis_c(option = "C") +
    ggplot2::labs(
      fill = "Mean |residual|",
      title = "Mean |EEoS residual| by ICES rectangle"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    file.path(fig_dir, "h2_map_abs_residual.png"),
    p_map_dv,
    width = 8,
    height = 6,
    dpi = 150
  )

  p_map_iv <- ggplot2::ggplot(panel_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = mean_annual_hours_total), colour = NA) +
    ggplot2::scale_fill_viridis_c(option = "A", trans = "log1p") +
    ggplot2::labs(
      fill = "Mean annual hours",
      title = "Mean annual fishing hours by ICES rectangle"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    file.path(fig_dir, "h2_map_fishing_pressure.png"),
    p_map_iv,
    width = 8,
    height = 6,
    dpi = 150
  )

  if (!is.null(ols_fit)) {
    grDevices::png(
      file.path(fig_dir, "h2_moran_scatter.png"),
      width = 800,
      height = 600,
      res = 150
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    invisible(spdep::moran.plot(
      residuals(ols_fit),
      weights$listw,
      zero.policy = TRUE,
      main = "Moran plot: OLS residuals"
    ))
  }

  invisible(fig_dir)
}

#' Format p-values for figure annotations.
h2_format_p_value <- function(p) {
  if (!is.finite(p)) {
    return("p = NA")
  }
  if (p < 0.001) {
    return("p < 0.001")
  }
  sprintf("p = %.3f", p)
}

#' Save headline two-panel figure: OLS scatter + OLS vs SEM coefficient comparison.
save_h2_topline_figure <- function(panel, ols_row, sem_row, spatial_diagnostics, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; skipping H2 topline figure.")
    return(invisible(NULL))
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork not available; skipping H2 topline figure.")
    return(invisible(NULL))
  }

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  hours_scale <- 10000
  moran_row <- spatial_diagnostics[
    spatial_diagnostics$variable == "ols_primary_abs_residuals" &
      spatial_diagnostics$test == "morans_i",
    ,
    drop = FALSE
  ]
  moran_i <- if (nrow(moran_row) > 0L) moran_row$statistic[1L] else NA_real_

  coef_df <- data.frame(
    model = factor(
      c("OLS", "Spatial error model"),
      levels = c("OLS", "Spatial error model")
    ),
    estimate = c(ols_row$estimate, sem_row$estimate) * hours_scale,
    ci_low = c(ols_row$ci_low, sem_row$ci_low) * hours_scale,
    ci_high = c(ols_row$ci_high, sem_row$ci_high) * hours_scale,
    p_value = c(ols_row$p_value, sem_row$p_value),
    stringsAsFactors = FALSE
  )
  coef_df$label <- vapply(
    coef_df$p_value,
    function(p) {
      if (!is.finite(p) || p >= 0.05) {
        h2_format_p_value(p)
      } else {
        paste0(h2_format_p_value(p), "*")
      }
    },
    character(1)
  )
  coef_df$fill <- ifelse(coef_df$p_value < 0.05, "#b2182b", "#4575b4")

  theme_h2 <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(colour = "grey30", size = 9.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )

  p_scatter <- ggplot2::ggplot(
    panel,
    ggplot2::aes(x = mean_annual_hours_total, y = mean_abs_residual)
  ) +
    ggplot2::geom_point(colour = "grey35", alpha = 0.75, size = 2) +
    ggplot2::geom_smooth(
      method = "lm",
      se = TRUE,
      colour = "#b2182b",
      fill = "#ef8a62",
      linewidth = 0.9
    ) +
    ggplot2::labs(
      x = "Mean annual fishing hours (1985–2015)",
      y = "Mean |EEoS residual| per rectangle",
      title = "Rectangle-level relationship",
      subtitle = sprintf(
        "N = %d rectangles · OLS slope negative (p = %.3f) · Moran's I = %.2f",
        nrow(panel),
        ols_row$p_value,
        moran_i
      )
    ) +
    theme_h2

  p_coef <- ggplot2::ggplot(
    coef_df,
    ggplot2::aes(x = model, y = estimate, colour = fill)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      width = 0.08,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, y = ci_high),
      vjust = -0.7,
      size = 3.3,
      colour = "grey20",
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 2.35,
      y = max(coef_df$ci_high, na.rm = TRUE) * 0.85,
      label = "Pre-registered\nexpectation: β > 0",
      hjust = 1,
      size = 3,
      colour = "grey40",
      lineheight = 0.95
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(
      x = NULL,
      y = sprintf("Change in mean |residual| per %s annual hours", format(hours_scale, big.mark = ",")),
      title = "Fishing-pressure effect",
      subtitle = "OLS association is not robust to spatial error correction (SEM)"
    ) +
    theme_h2 +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))

  p_topline <- p_scatter + p_coef +
    patchwork::plot_layout(widths = c(1.15, 1)) +
    patchwork::plot_annotation(
      title = "H2: Higher fishing pressure is not associated with larger EEoS residuals",
      subtitle = paste0(
        "ICES rectangles · NS-IBTS Q1 · 1985–2015 · ",
        "Opposite to pre-registered direction in OLS; null after spatial correction"
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(colour = "grey35", size = 10.5)
      )
    )

  out_path <- file.path(fig_dir, "h2_topline_result.png")
  ggplot2::ggsave(out_path, p_topline, width = 11, height = 5.2, dpi = 150)
  invisible(out_path)
}
