# Shared helpers for H1 analysis scripts (north_sea_eeos)

#' Walk up from start_dir until FishGlob_data/ is found (repo workspace root).
get_project_root_from <- function(start_dir) {
  dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(dir, "FishGlob_data"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      break
    }
    dir <- parent
  }
  stop(
    "Could not find project root (directory containing FishGlob_data/). ",
    "Started from: ", start_dir
  )
}

#' Resolve project root when run via Rscript or interactively.
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  start <- if (length(file_arg) > 0L) {
    dirname(sub("^--file=", "", file_arg))
  } else {
    getwd()
  }
  get_project_root_from(start)
}

#' Source pipeline R helpers and return project_root + script_dir.
load_pipeline_context <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg) > 0L) {
    normalizePath(
      dirname(sub("^--file=", "", file_arg)),
      winslash = "/",
      mustWork = TRUE
    )
  } else if (dir.exists("R")) {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  } else {
    stop("Run scripts with Rscript pipeline/script.R or set wd to pipeline/")
  }
  source(file.path(script_dir, "R", "h1_common.R"))
  list(
    project_root = get_project_root_from(script_dir),
    script_dir = script_dir
  )
}

#' R-squared on log scale (same metric as build_eeos_predictions.R).
log_r2 <- function(log_obs, log_pred) {
  valid <- is.finite(log_obs) & is.finite(log_pred)
  log_obs <- log_obs[valid]
  log_pred <- log_pred[valid]
  if (length(log_obs) < 2L) {
    return(NA_real_)
  }
  ss_res <- sum((log_obs - log_pred)^2)
  ss_tot <- sum((log_obs - mean(log_obs))^2)
  if (ss_tot == 0) {
    return(NA_real_)
  }
  1 - ss_res / ss_tot
}

log_ss_res <- function(log_obs, log_pred) {
  sum((log_obs - log_pred)^2, na.rm = TRUE)
}

median_abs_log_residual <- function(log_obs, log_pred) {
  median(abs(log_obs - log_pred), na.rm = TRUE)
}

#' Squared Pearson correlation on log scale (diagnostic; not the primary H1 R²).
log_cor2 <- function(log_obs, log_pred) {
  valid <- is.finite(log_obs) & is.finite(log_pred)
  log_obs <- log_obs[valid]
  log_pred <- log_pred[valid]
  if (length(log_obs) < 2L) {
    return(NA_real_)
  }
  stats::cor(log_obs, log_pred)^2
}

#' Squared Pearson correlation on arbitrary scale.
cor2 <- function(obs, pred) {
  valid <- is.finite(obs) & is.finite(pred)
  obs <- obs[valid]
  pred <- pred[valid]
  if (length(obs) < 2L) {
    return(NA_real_)
  }
  stats::cor(obs, pred)^2
}

#' Coefficient of determination on arbitrary scale (can be negative).
r2 <- function(obs, pred) {
  valid <- is.finite(obs) & is.finite(pred)
  obs <- obs[valid]
  pred <- pred[valid]
  if (length(obs) < 2L) {
    return(NA_real_)
  }
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  if (ss_tot == 0) {
    return(NA_real_)
  }
  1 - ss_res / ss_tot
}

#' Root mean squared error on log scale.
log_rmse <- function(log_obs, log_pred) {
  valid <- is.finite(log_obs) & is.finite(log_pred)
  log_obs <- log_obs[valid]
  log_pred <- log_pred[valid]
  if (length(log_obs) == 0L) {
    return(NA_real_)
  }
  sqrt(mean((log_obs - log_pred)^2))
}

#' Root mean squared error on raw scale.
raw_rmse <- function(obs, pred) {
  valid <- is.finite(obs) & is.finite(pred)
  obs <- obs[valid]
  pred <- pred[valid]
  if (length(obs) == 0L) {
    return(NA_real_)
  }
  sqrt(mean((obs - pred)^2))
}

#' Harte-style productivity ratio E / B^(3/4).
productivity_ratio <- function(E_mass075, B_g) {
  E_mass075 / B_g^(3 / 4)
}

#' Rank Fig 2 hauls by leverage proxy (predicted_ratio × observed_ratio).
#' Returns indices into the original vectors (not the valid subset).
fig2_leverage_rank_idx <- function(predicted_ratio, observed_ratio) {
  valid <- which(
    is.finite(predicted_ratio) & is.finite(observed_ratio) &
      predicted_ratio > 0 & observed_ratio > 0
  )
  if (length(valid) == 0L) {
    return(integer())
  }
  mag <- predicted_ratio[valid] * observed_ratio[valid]
  valid[order(mag, decreasing = TRUE)]
}

#' Fig 2 leverage diagnostic: Pearson R² all vs trimmed (top-N magnitude points).
#' Does not drop outliers from the haul table — diagnostic only for the ratio regression.
fig2_leverage_diagnostic <- function(
    haul,
    n_exclude = 2L,
    pred_col = "fig2_predicted_ratio",
    obs_col = "fig2_observed_ratio") {
  stopifnot(pred_col %in% names(haul), obs_col %in% names(haul))
  n_exclude <- as.integer(n_exclude)
  pred <- haul[[pred_col]]
  obs <- haul[[obs_col]]
  ranked <- fig2_leverage_rank_idx(pred, obs)
  n_valid <- length(ranked)
  n_drop <- min(n_exclude, max(0L, n_valid - 2L))
  drop_idx <- if (n_drop > 0L) ranked[seq_len(n_drop)] else integer()

  r2_all <- if (n_valid >= 2L) {
    cor2(pred[ranked], obs[ranked])
  } else {
    NA_real_
  }
  keep <- setdiff(ranked, drop_idx)
  r2_trim <- if (length(keep) >= 2L) {
    cor2(pred[keep], obs[keep])
  } else {
    NA_real_
  }

  id_col <- if ("haul_id" %in% names(haul)) {
    "haul_id"
  } else if ("haul_key" %in% names(haul)) {
    "haul_key"
  } else {
    NA_character_
  }
  flag_cols <- c(id_col, "S", "N", "E", "B_obs", "B_pred", pred_col, obs_col)
  flag_cols <- flag_cols[!is.na(flag_cols) & flag_cols %in% names(haul)]
  excluded <- if (length(drop_idx) > 0L) {
    out <- haul[drop_idx, flag_cols, drop = FALSE]
    out$leverage_rank <- seq_along(drop_idx)
    out$leverage_product <- pred[drop_idx] * obs[drop_idx]
    out
  } else {
    data.frame()
  }

  list(
    fig2_r2_pearson_all = r2_all,
    fig2_r2_pearson_trimmed = r2_trim,
    fig2_r2_pearson_n_excluded = n_drop,
    excluded_hauls = excluded,
    drop_idx = drop_idx
  )
}

#' Harte Fig 2 metrics on raw ratio scale (make_plots.py Fig 2 block).
#' x = predicted_ratio, y = observed_ratio; headline = Pearson r² (linregress rvalue²).
#' Also returns leverage-trimmed Pearson R² (top n_exclude by |pred×obs| magnitude).
harte_fig2_metrics <- function(predicted_ratio, observed_ratio, n_exclude = 2L) {
  valid <- is.finite(predicted_ratio) & is.finite(observed_ratio) &
    predicted_ratio > 0 & observed_ratio > 0
  x <- predicted_ratio[valid]
  y <- observed_ratio[valid]
  n <- length(x)
  empty <- data.frame(
    fig2_r2_pearson = NA_real_,
    fig2_r2_pearson_all = NA_real_,
    fig2_r2_pearson_trimmed = NA_real_,
    fig2_r2_pearson_n_excluded = NA_integer_,
    fig2_r2_cod_extended = NA_real_,
    cor2 = NA_real_,
    log_r2 = NA_real_,
    log_rmse = NA_real_,
    raw_rmse = NA_real_,
    median_ratio = NA_real_,
    median_abs_log_resid = NA_real_,
    ss_res_log = NA_real_,
    n = n,
    stringsAsFactors = FALSE
  )
  if (n < 2L) {
    return(empty)
  }
  pearson_all <- cor2(x, y)
  n_drop <- min(as.integer(n_exclude), n - 2L)
  if (n_drop > 0L) {
    mag_ord <- order(x * y, decreasing = TRUE)
    keep <- -mag_ord[seq_len(n_drop)]
    pearson_trim <- cor2(x[keep], y[keep])
  } else {
    pearson_trim <- pearson_all
    n_drop <- 0L
  }
  cod_ext <- r2(y, x)
  data.frame(
    fig2_r2_pearson = pearson_all,
    fig2_r2_pearson_all = pearson_all,
    fig2_r2_pearson_trimmed = pearson_trim,
    fig2_r2_pearson_n_excluded = n_drop,
    fig2_r2_cod_extended = cod_ext,
    cor2 = pearson_all,
    log_r2 = cod_ext,
    log_rmse = NA_real_,
    raw_rmse = raw_rmse(y, x),
    median_ratio = median(x / y, na.rm = TRUE),
    median_abs_log_resid = median(abs(y - x), na.rm = TRUE),
    ss_res_log = sum((y - x)^2),
    n = n,
    stringsAsFactors = FALSE
  )
}

#' Full metric bundle for one obs/pred pairing.
evaluate_prediction <- function(obs, pred, log_scale = TRUE) {
  if (log_scale) {
    log_obs <- log(obs)
    log_pred <- log(pred)
    data.frame(
      log_r2 = log_r2(log_obs, log_pred),
      cor2 = log_cor2(log_obs, log_pred),
      log_rmse = log_rmse(log_obs, log_pred),
      raw_rmse = raw_rmse(obs, pred),
      median_ratio = median(pred / obs, na.rm = TRUE),
      median_abs_log_resid = median_abs_log_residual(log_obs, log_pred),
      ss_res_log = log_ss_res(log_obs, log_pred),
      n = sum(is.finite(obs) & is.finite(pred) & obs > 0 & pred > 0),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      log_r2 = r2(obs, pred),
      cor2 = cor2(obs, pred),
      log_rmse = NA_real_,
      raw_rmse = raw_rmse(obs, pred),
      median_ratio = median(pred / obs, na.rm = TRUE),
      median_abs_log_resid = median(abs(obs - pred), na.rm = TRUE),
      ss_res_log = sum((obs - pred)^2, na.rm = TRUE),
      n = sum(is.finite(obs) & is.finite(pred)),
      stringsAsFactors = FALSE
    )
  }
}
