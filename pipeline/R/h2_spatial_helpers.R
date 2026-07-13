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

#' Fit spatial error model for primary H2 relationship.
fit_h2_sem <- function(panel, listw) {
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("Package 'spatialreg' is required. Install via renv::install('spatialreg').")
  }

  fit <- spatialreg::errorsarlm(
    mean_abs_residual ~ mean_annual_hours_total,
    data = panel,
    listw = listw,
    zero.policy = TRUE
  )

  coef_mat <- coef(summary(fit))
  ci <- confint(fit)
  term <- "mean_annual_hours_total"

  data.frame(
    model_id = "sem_primary_abs",
    term = term,
    estimate = coef_mat[term, "Estimate"],
    std_error = coef_mat[term, "Std. Error"],
    statistic = coef_mat[term, "z value"],
    p_value = coef_mat[term, "Pr(>|z|)"],
    ci_low = ci[term, 1L],
    ci_high = ci[term, 2L],
    lambda = unname(fit$lambda),
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
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
