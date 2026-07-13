# Helpers for H1 null-model analysis (B randomisation)
# Source R/h1_common.R before this file (provides log_r2, etc.).

distribution_skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  if (s == 0) {
    return(0)
  }
  mean((x - m)^3) / s^3
}

distribution_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  if (s == 0) {
    return(0)
  }
  mean((x - m)^4) / s^4 - 3
}

#' Assess shape of log(B_obs) and choose primary null sampling scheme.
#'
#' Returns list with decision, metrics, and plain-language rationale.
assess_b_null_sampling <- function(B_obs_g, n_hist_bins = 30L) {
  logB <- log(B_obs_g[is.finite(B_obs_g) & B_obs_g > 0])
  if (length(logB) < 10L) {
    stop("Too few valid B_obs values for distribution assessment.")
  }

  rng <- diff(range(logB))
  iqr_range_ratio <- if (rng > 0) IQR(logB) / rng else NA_real_
  skew <- distribution_skewness(logB)
  kurt <- distribution_excess_kurtosis(logB)

  h <- hist(logB, breaks = n_hist_bins, plot = FALSE)
  bin_props <- h$counts / sum(h$counts)
  peak_ratio <- max(bin_props) / (1 / n_hist_bins)

  # Peaked: central mass tight vs full range, histogram peak >> uniform,
  # or positive excess kurtosis (supervisor concern: resampling ~ circular).
  peaked <- isTRUE(iqr_range_ratio < 0.25) ||
    isTRUE(peak_ratio > 2.5) ||
    isTRUE(kurt > 0.75)

  primary_method <- if (peaked) "uniform_log_b_95" else "b_shuffle"
  robustness_method <- if (peaked) "b_shuffle" else "uniform_log_b_95"

  rationale <- if (peaked) {
    paste0(
      "log(B_obs) is concentrated (IQR/range = ", round(iqr_range_ratio, 3),
      ", peak ratio = ", round(peak_ratio, 2),
      ", excess kurtosis = ", round(kurt, 2),
      "). Empirical resampling would largely recreate the observed biomass scale; ",
      "primary null uses a uniform draw on the central 95% of log(B_obs)."
    )
  } else {
    paste0(
      "log(B_obs) is not strongly peaked (IQR/range = ", round(iqr_range_ratio, 3),
      ", peak ratio = ", round(peak_ratio, 2),
      "). Primary null permutes observed B_obs across hauls (B shuffle)."
    )
  }

  q025 <- as.numeric(quantile(logB, 0.025))
  q975 <- as.numeric(quantile(logB, 0.975))

  list(
    primary_method = primary_method,
    robustness_method = robustness_method,
    peaked = peaked,
    rationale = rationale,
    metrics = data.frame(
      n = length(logB),
      log_median = median(logB),
      log_iqr = IQR(logB),
      log_range = rng,
      iqr_range_ratio = iqr_range_ratio,
      skewness = skew,
      excess_kurtosis = kurt,
      hist_peak_ratio = peak_ratio,
      log_q025 = q025,
      log_q975 = q975,
      stringsAsFactors = FALSE
    ),
    log_bounds_95 = c(q025, q975)
  )
}

#' One null replicate: randomise B_obs, keep B_pred fixed.
null_r2_replicate <- function(
  B_obs_g,
  B_pred_g,
  method = c("uniform_log_b_95", "b_shuffle"),
  log_bounds_95 = NULL
) {
  method <- match.arg(method)
  log_pred <- log(B_pred_g)
  n <- length(B_obs_g)

  B_null <- if (method == "b_shuffle") {
    sample(B_obs_g, size = n, replace = FALSE)
  } else {
    if (is.null(log_bounds_95)) {
      logB <- log(B_obs_g)
      log_bounds_95 <- quantile(logB, c(0.025, 0.975), names = FALSE)
    }
    exp(runif(n, min = log_bounds_95[1L], max = log_bounds_95[2L]))
  }

  log_r2(log(B_null), log_pred)
}

#' Run n_perm null replicates; return vector of null R² values.
run_null_simulation <- function(
  B_obs_g,
  B_pred_g,
  method,
  n_perm = 999L,
  log_bounds_95 = NULL,
  seed = 42L
) {
  set.seed(seed)
  vapply(
    seq_len(n_perm),
    function(i) {
      null_r2_replicate(B_obs_g, B_pred_g, method, log_bounds_95)
    },
    numeric(1)
  )
}

#' Summarise observed vs null R² distribution.
summarise_null_test <- function(T_obs, T_null, method_label) {
  T_null <- T_null[is.finite(T_null)]
  p_one_sided <- mean(T_null >= T_obs)
  p_two_sided <- mean(abs(T_null - mean(T_null)) >= abs(T_obs - mean(T_null)))

  data.frame(
    method = method_label,
    T_obs = T_obs,
    null_median = median(T_null),
    null_mean = mean(T_null),
    null_q025 = quantile(T_null, 0.025),
    null_q975 = quantile(T_null, 0.975),
    p_value_one_sided = p_one_sided,
    p_value_two_sided = p_two_sided,
    stringsAsFactors = FALSE
  )
}

plot_null_r2_distribution <- function(
  T_obs,
  T_null,
  method_label,
  path_out
) {
  png(path_out, width = 900, height = 600, res = 120)
  on.exit(dev.off(), add = TRUE)

  hist(
    T_null,
    breaks = 40,
    col = adjustcolor("steelblue", 0.6),
    border = "white",
    main = paste0("Null distribution of R² (log scale): ", method_label),
    xlab = expression(R^2 ~ "(log B_obs vs log B_pred)"),
    ylab = "Count"
  )
  abline(v = T_obs, col = "red", lwd = 2, lty = 1)
  abline(v = median(T_null), col = "darkblue", lwd = 2, lty = 2)
  legend(
    "topleft",
    legend = c(
      paste0("Observed R² = ", round(T_obs, 3)),
      paste0("Null median = ", round(median(T_null), 3))
    ),
    col = c("red", "darkblue"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
}
