# Harte et al. (2022) unfitted baseline comparisons (Tier 1)
# Requires h1_common.R to be sourced first.

#' Join E_raw from haul state if missing from EEoS predictions table.
ensure_haul_E_raw <- function(haul, project_root) {
  if ("E_raw" %in% names(haul) && !all(is.na(haul$E_raw))) {
    return(haul)
  }
  path_state <- file.path(project_root, "outputs", "haul_state_variables.rds")
  stopifnot(file.exists(path_state))
  state <- readRDS(path_state)
  key <- if ("haul_id" %in% names(haul)) "haul_id" else "haul_key"
  state_key <- if ("haul_id" %in% names(state)) "haul_id" else "haul_key"
  merge(
    haul,
    state[, c(state_key, "E_raw"), drop = FALSE],
    by = setNames(state_key, key),
    all.x = TRUE
  )
}

#' Add Harte baseline columns: calibrated productivity proxy, Fig 2 ratios.
#'
#' Productivity 1:1 uses E_calibrated = E × m_min (same normalised E as
#' biomass() / Test 2, same m_min as B_pred = B_pred_norm × m_min). This is
#' unfitted — m_min is from the smallest-individual convention, not fit to
#' B_obs. E_raw is retained only as an uncalibrated diagnostic.
augment_haul_harte_baseline <- function(haul, project_root = NULL) {
  if (!is.null(project_root)) {
    haul <- ensure_haul_E_raw(haul, project_root)
  }
  if (!"ln_B_obs" %in% names(haul)) {
    haul$ln_B_obs <- log(haul$B_obs)
  }
  stopifnot("E" %in% names(haul), "m_min" %in% names(haul))
  # Same independently-derived m_min used to put B_pred into grams.
  haul$E_calibrated <- haul$E * haul$m_min
  haul$ln_E_calibrated <- log(haul$E_calibrated)
  haul$ln_E_raw <- log(haul$E_raw)
  # Harte Fig 2: E / B^(3/4) using normalised E (same E as biomass() input), not E_raw.
  haul$fig2_predicted_ratio <- productivity_ratio(haul$E, haul$B_pred)
  haul$fig2_observed_ratio <- productivity_ratio(haul$E, haul$B_obs)
  haul$ratio_pred <- haul$fig2_predicted_ratio
  haul$ratio_obs <- haul$fig2_observed_ratio
  haul
}

#' Tier-1 unfitted model metrics (one row per model).
harte_baseline_metrics <- function(haul) {
  eeos <- evaluate_prediction(haul$B_obs, haul$B_pred)
  eeos$fig2_r2_pearson <- NA_real_
  eeos$fig2_r2_pearson_all <- NA_real_
  eeos$fig2_r2_pearson_trimmed <- NA_real_
  eeos$fig2_r2_pearson_n_excluded <- NA_integer_
  eeos$fig2_r2_cod_extended <- NA_real_
  eeos$model <- "EEoS biomass (Harte Fig 1)"
  eeos$model_id <- "eeos_biomass"
  eeos$fitted <- FALSE
  eeos$tier <- 1L
  eeos$comparison_type <- "log B_pred vs log B_obs"

  # Headline productivity baseline: E_norm × m_min (grams-equivalent, unfitted).
  # Not E_raw — that quantity is uncalibrated and kept only as a diagnostic below.
  prod1 <- evaluate_prediction(haul$B_obs, haul$E_calibrated)
  prod1$fig2_r2_pearson <- NA_real_
  prod1$fig2_r2_pearson_all <- NA_real_
  prod1$fig2_r2_pearson_trimmed <- NA_real_
  prod1$fig2_r2_pearson_n_excluded <- NA_integer_
  prod1$fig2_r2_cod_extended <- NA_real_
  prod1$model <- "Productivity 1:1 (log(E × m_min) vs log B_obs)"
  prod1$model_id <- "productivity_1to1"
  prod1$fitted <- FALSE
  prod1$tier <- 1L
  prod1$comparison_type <- "unfitted 1:1 productivity map (E_norm × m_min)"

  # Diagnostic only — pre-calibration E_raw baseline; not used in headline Test 1.
  prod1_uncal <- evaluate_prediction(haul$B_obs, haul$E_raw)
  prod1_uncal$fig2_r2_pearson <- NA_real_
  prod1_uncal$fig2_r2_pearson_all <- NA_real_
  prod1_uncal$fig2_r2_pearson_trimmed <- NA_real_
  prod1_uncal$fig2_r2_pearson_n_excluded <- NA_integer_
  prod1_uncal$fig2_r2_cod_extended <- NA_real_
  prod1_uncal$model <- "Productivity 1:1 uncalibrated (log E_raw vs log B_obs)"
  prod1_uncal$model_id <- "productivity_1to1_uncalibrated"
  prod1_uncal$fitted <- FALSE
  prod1_uncal$tier <- 1L
  prod1_uncal$comparison_type <- "diagnostic: uncalibrated E_raw 1:1 (not headline)"

  ratio <- harte_fig2_metrics(haul$fig2_predicted_ratio, haul$fig2_observed_ratio)
  ratio$model <- "Productivity ratio (Harte Fig 2)"
  ratio$model_id <- "productivity_ratio"
  ratio$fitted <- FALSE
  ratio$tier <- 1L
  ratio$comparison_type <- "E/B^(3/4) predicted vs observed (normalised E)"

  rbind(eeos, prod1, ratio, prod1_uncal)
}

#' Base R scatter with 1:1 line and metric annotation (Harte-style).
plot_harte_1to1 <- function(
    x, y,
    xlab, ylab, main,
    annotate_cor2 = TRUE) {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  rmin <- min(c(x, y), na.rm = TRUE)
  rmax <- max(c(x, y), na.rm = TRUE)
  pad <- 0.05 * (rmax - rmin)
  rmin <- rmin - pad
  rmax <- rmax + pad

  plot(
    x, y,
    pch = 16,
    col = adjustcolor("black", 0.12),
    cex = 0.45,
    xlab = xlab,
    ylab = ylab,
    main = main,
    xlim = c(rmin, rmax),
    ylim = c(rmin, rmax),
    asp = 1
  )
  abline(0, 1, lty = 2, lwd = 2, col = "grey40")
  if (annotate_cor2 && length(x) >= 2L) {
    r2_cod <- if (min(x) > 0 && min(y) > 0) {
      log_r2(log(x), log(y))
    } else {
      r2(x, y)
    }
    c2 <- cor(x, y)^2
    legend(
      "topleft",
      legend = c(
        sprintf("log R^2 = %.3f", r2_cod),
        sprintf("cor^2 = %.3f", c2)
      ),
      bty = "n",
      cex = 0.9
    )
  }
  invisible(list(x = x, y = y))
}

save_harte_baseline_figures <- function(haul, fig_dir) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  png(
    file.path(fig_dir, "harte_fig1_logB_pred_vs_obs.png"),
    width = 900, height = 900, res = 120
  )
  plot_harte_1to1(
    haul$ln_B_pred, haul$ln_B_obs,
    xlab = expression(log(B[pred]) ~ "[g]"),
    ylab = expression(log(B[obs]) ~ "[g]"),
    main = "Harte Fig 1: EEoS predicted vs observed biomass"
  )
  dev.off()

  png(
    file.path(fig_dir, "productivity_1to1_logB_vs_logE_calibrated.png"),
    width = 900, height = 900, res = 120
  )
  plot_harte_1to1(
    haul$ln_E_calibrated, haul$ln_B_obs,
    xlab = expression(log(E %*% m[min]) ~ "[g-equiv]"),
    ylab = expression(log(B[obs]) ~ "[g]"),
    main = "Unfitted productivity 1:1: log(E × m_min) vs log(B_obs)"
  )
  dev.off()

  # Diagnostic: pre-calibration E_raw baseline (not headline Test 1).
  png(
    file.path(fig_dir, "productivity_1to1_uncalibrated_logB_vs_logEraw.png"),
    width = 900, height = 900, res = 120
  )
  plot_harte_1to1(
    haul$ln_E_raw, haul$ln_B_obs,
    xlab = expression(log(E[raw])),
    ylab = expression(log(B[obs]) ~ "[g]"),
    main = "Diagnostic (uncalibrated): log(E_raw) vs log(B_obs)"
  )
  dev.off()

  png(
    file.path(fig_dir, "harte_fig2_productivity_ratio.png"),
    width = 900, height = 900, res = 120
  )
  plot_harte_fig2(haul$fig2_predicted_ratio, haul$fig2_observed_ratio)
  dev.off()
}

#' Harte Fig 2 scatter: sqrt axes (visual only), tight zoom + full-range inset.
#' Stats still use raw ratios; annotates trimmed Pearson R² as primary.
plot_harte_fig2 <- function(
    predicted_ratio, observed_ratio,
    xlim = c(0, 60), ylim = c(0, 60),
    n_exclude = 2L) {
  valid <- is.finite(predicted_ratio) & is.finite(observed_ratio) &
    predicted_ratio > 0 & observed_ratio > 0
  x <- predicted_ratio[valid]
  y <- observed_ratio[valid]
  m <- harte_fig2_metrics(predicted_ratio, observed_ratio, n_exclude = n_exclude)
  n_off <- sum(x > xlim[2] | y > ylim[2])

  sx <- sqrt(x)
  sy <- sqrt(y)
  sxlim <- sqrt(xlim)
  sylim <- sqrt(ylim)

  plot(
    sx, sy,
    pch = 16,
    col = adjustcolor("black", 0.12),
    cex = 0.45,
    xlab = expression(Predicted ~ ratio ~ E/B^{3/4} ~ (sqrt ~ scale)),
    ylab = expression(Observed ~ ratio ~ E/B^{3/4} ~ (sqrt ~ scale)),
    main = "Harte Fig 2: productivity ratio (sqrt display)",
    xlim = sxlim,
    ylim = sylim,
    asp = 1,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE
  )
  # Tick labels in original ratio units
  xt <- pretty(xlim, n = 5)
  yt <- pretty(ylim, n = 5)
  axis(1, at = sqrt(xt), labels = xt)
  axis(2, at = sqrt(yt), labels = yt)
  box()
  abline(0, 1, lty = 1, lwd = 2, col = "grey40")
  legend(
    "topleft",
    legend = c(
      sprintf(
        "R² (Pearson, trimmed N=%d) = %.3f",
        m$fig2_r2_pearson_n_excluded,
        m$fig2_r2_pearson_trimmed
      ),
      sprintf("R² (Pearson, all) = %.3f", m$fig2_r2_pearson_all),
      if (n_off > 0L) {
        sprintf("%d haul(s) outside panel (see inset)", n_off)
      } else {
        NULL
      }
    ),
    bty = "n",
    cex = 0.8
  )

  # Full-range inset (same sqrt transform)
  full_max <- max(c(x, y), na.rm = TRUE) * 1.05
  sfull <- c(0, sqrt(full_max))
  op <- par(fig = c(0.55, 0.95, 0.08, 0.42), new = TRUE, mar = c(2, 2, 1.2, 0.4))
  on.exit(par(op), add = TRUE)
  plot(
    sx, sy,
    pch = 16,
    col = adjustcolor("black", 0.10),
    cex = 0.25,
    xlim = sfull,
    ylim = sfull,
    asp = 1,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    main = "Full range",
    cex.main = 0.75
  )
  abline(0, 1, col = "grey40", lwd = 1)
  # Highlight top-N leverage points
  if (n_exclude > 0L && length(x) >= n_exclude + 2L) {
    mag_ord <- order(x * y, decreasing = TRUE)
    hi <- mag_ord[seq_len(n_exclude)]
    points(sx[hi], sy[hi], pch = 16, col = "#b2182b", cex = 0.9)
  }
  ft <- pretty(c(0, full_max), n = 3)
  axis(1, at = sqrt(ft), labels = ft, cex.axis = 0.55)
  axis(2, at = sqrt(ft), labels = ft, cex.axis = 0.55)
  box(col = "grey40")
  invisible(m)
}
