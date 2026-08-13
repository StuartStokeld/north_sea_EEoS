# Helpers for bathymetry / directional residual anisotropy diagnostic.
# Design: display_discussion/Design_bathymetry_spatial_anisotropy.md

#' GEBCO paths for the Aug 2026 North Sea clip.
gebco_paths <- function(project_root) {
  base <- file.path(
    project_root,
    "data",
    "external",
    "GEBCO_02_Aug_2026_e390ca8d46b0"
  )
  list(
    dir = base,
    bathy = file.path(base, "gebco_2026_n65.0_s50.0_w-5.0_e10.0_geotiff.tif"),
    tid = file.path(base, "gebco_2026_tid_n65.0_s50.0_w-5.0_e10.0_geotiff.tif")
  )
}

#' Ensure .R_libs (local installs of circular/gstat) is on the search path.
ensure_bathymetry_libs <- function(project_root) {
  local_lib <- file.path(project_root, ".R_libs")
  if (dir.exists(local_lib)) {
    .libPaths(c(normalizePath(local_lib, winslash = "/"), .libPaths()))
  }
  need <- c("circular", "gstat", "sp")
  missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install into ", local_lib, " (see design §6)."
    )
  }
  invisible(NULL)
}

#' Fold a directed bearing (degrees) onto the axial half-circle [0, 180).
fold_axial_deg <- function(deg) {
  ((deg %% 180) + 180) %% 180
}

#' Pairwise geographic bearings (gstat plane convention: deg CCW from east)
#' and great-circle-ish Euclidean distance in degrees (lon compressed by cos lat).
pair_geometry <- function(lon, lat) {
  n <- length(lon)
  stopifnot(length(lat) == n)
  lat0 <- mean(lat, na.rm = TRUE)
  coslat <- cos(lat0 * pi / 180)
  i <- rep(seq_len(n), each = n)
  j <- rep(seq_len(n), times = n)
  keep <- i < j
  i <- i[keep]
  j <- j[keep]
  dx <- (lon[j] - lon[i]) * coslat
  dy <- lat[j] - lat[i]
  dist <- sqrt(dx * dx + dy * dy)
  # undirected pair bearing folded to [0, 180)
  bearing <- (atan2(dy, dx) * 180 / pi) %% 180
  data.frame(i = i, j = j, dist = dist, bearing = bearing)
}

#' Assign undirected bearing to directional bins (centres 0,45,90,135).
bearing_bin <- function(bearing_deg, centres = c(0, 45, 90, 135), tol = 22.5) {
  # circular distance on the 180° axis
  d <- vapply(centres, function(c) {
    raw <- abs(bearing_deg - c)
    pmin(raw, 180 - raw)
  }, numeric(length(bearing_deg)))
  if (is.null(dim(d))) {
    d <- matrix(d, nrow = 1L)
  }
  idx <- max.col(-d, ties.method = "first")
  ok <- d[cbind(seq_len(nrow(d)), idx)] <= tol + 1e-9
  out <- centres[idx]
  out[!ok] <- NA_real_
  out
}

#' Empirical directional variogram (semivariance by distance lag × direction).
#' @return data.frame(np, dist, gamma, dir.hor, id)
empirical_directional_variogram <- function(
    z,
    lon,
    lat,
    alphas = c(0, 45, 90, 135),
    tol = 22.5,
    width = 1.0,
    cutoff = 8.0
) {
  geom <- pair_geometry(lon, lat)
  geom$dir <- bearing_bin(geom$bearing, centres = alphas, tol = tol)
  geom <- geom[is.finite(geom$dir) & geom$dist > 0 & geom$dist <= cutoff, , drop = FALSE]
  geom$gamma_pair <- 0.5 * (z[geom$i] - z[geom$j])^2
  geom$lag <- width * (floor(geom$dist / width) + 0.5)
  # drop pairs beyond last full bin centre under cutoff
  agg <- aggregate(
    cbind(gamma = geom$gamma_pair, dist = geom$dist),
    by = list(dir.hor = geom$dir, lag = geom$lag),
    FUN = mean,
    na.rm = TRUE
  )
  np <- aggregate(
    geom$gamma_pair,
    by = list(dir.hor = geom$dir, lag = geom$lag),
    FUN = length
  )
  names(np)[3] <- "np"
  out <- merge(agg, np, by = c("dir.hor", "lag"))
  out <- out[order(out$dir.hor, out$lag), , drop = FALSE]
  rownames(out) <- NULL
  out$id <- paste0("dir_", out$dir.hor)
  out[, c("np", "dist", "gamma", "dir.hor", "id", "lag")]
}

#' Fit exponential variogram range/sill per direction (descriptive).
#' Uses gstat::fit.variogram when available; falls back to NA on failure.
fit_directional_exponentials <- function(vgram_df) {
  ensure_has <- function() {
    if (!requireNamespace("gstat", quietly = TRUE) || !requireNamespace("sp", quietly = TRUE)) {
      stop("gstat/sp required for fit_directional_exponentials")
    }
  }
  ensure_has()
  dirs <- sort(unique(vgram_df$dir.hor))
  rows <- lapply(dirs, function(d) {
    sub <- vgram_df[vgram_df$dir.hor == d & is.finite(vgram_df$gamma) & vgram_df$np >= 5, , drop = FALSE]
    if (nrow(sub) < 3L) {
      return(data.frame(
        dir.hor = d, nugget = NA_real_, sill = NA_real_, range = NA_real_,
        n_lags = nrow(sub), fit_ok = FALSE
      ))
    }
    # gstat expects class variogram with columns np, dist, gamma, dir.hor, id
    # and numeric (not integer) columns for the C backend
    vg <- data.frame(
      np = as.numeric(sub$np),
      dist = as.numeric(sub$dist),
      gamma = as.numeric(sub$gamma),
      dir.hor = as.numeric(sub$dir.hor),
      dir.ver = 0,
      id = as.character(sub$id),
      stringsAsFactors = FALSE
    )
    class(vg) <- c("gstatVariogram", "data.frame")
    psill0 <- stats::var(sub$gamma, na.rm = TRUE)
    if (!is.finite(psill0) || psill0 <= 0) psill0 <- max(sub$gamma, na.rm = TRUE)
    init <- gstat::vgm(
      psill = psill0,
      model = "Exp",
      range = max(stats::median(sub$dist, na.rm = TRUE), 1),
      nugget = max(min(sub$gamma, na.rm = TRUE), 0)
    )
    fit <- tryCatch(
      gstat::fit.variogram(vg, model = init, fit.method = 7),
      error = function(e) {
        tryCatch(
          gstat::fit.variogram(vg, model = init, fit.method = 1),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(fit) || nrow(fit) < 1L) {
      return(data.frame(
        dir.hor = d, nugget = NA_real_, sill = NA_real_, range = NA_real_,
        n_lags = nrow(sub), fit_ok = FALSE
      ))
    }
    nug <- if ("Nug" %in% fit$model) fit$psill[fit$model == "Nug"][1] else 0
    exp_row <- fit[fit$model == "Exp", , drop = FALSE]
    if (nrow(exp_row) < 1L) {
      return(data.frame(
        dir.hor = d, nugget = nug, sill = NA_real_, range = NA_real_,
        n_lags = nrow(sub), fit_ok = FALSE
      ))
    }
    data.frame(
      dir.hor = d,
      nugget = as.numeric(nug),
      sill = as.numeric(exp_row$psill[1]),
      range = as.numeric(exp_row$range[1]),
      n_lags = nrow(sub),
      fit_ok = TRUE
    )
  })
  do.call(rbind, rows)
}

#' Local residual-correlation bearing per rectangle.
#'
#' Among pairs involving rectangle k within [lag_min, lag_max], pick the
#' direction bin with lowest mean semivariance (strongest similarity).
local_correlation_bearing <- function(
    z,
    lon,
    lat,
    lag_min = 0.5,
    lag_max = 4.0,
    alphas = c(0, 45, 90, 135),
    tol = 22.5,
    min_pairs_per_dir = 2L
) {
  n <- length(z)
  geom <- pair_geometry(lon, lat)
  geom$dir <- bearing_bin(geom$bearing, centres = alphas, tol = tol)
  geom <- geom[
    is.finite(geom$dir) & geom$dist >= lag_min & geom$dist <= lag_max,
    ,
    drop = FALSE
  ]
  geom$gamma_pair <- 0.5 * (z[geom$i] - z[geom$j])^2

  bearings <- rep(NA_real_, n)
  n_pairs <- integer(n)
  for (k in seq_len(n)) {
    hit <- geom$i == k | geom$j == k
    sub <- geom[hit, , drop = FALSE]
    n_pairs[k] <- nrow(sub)
    if (nrow(sub) < length(alphas) * min_pairs_per_dir) {
      next
    }
    tab <- aggregate(
      sub$gamma_pair,
      by = list(dir = sub$dir),
      FUN = function(x) c(mean = mean(x), n = length(x))
    )
    # aggregate returns matrix column when FUN returns vector
    means <- tab$x[, "mean"]
    ns <- tab$x[, "n"]
    ok <- ns >= min_pairs_per_dir
    if (!any(ok)) {
      next
    }
    # lowest semivariance among admissible dirs
    pick <- which(ok)[which.min(means[ok])]
    bearings[k] <- tab$dir[pick]
  }
  data.frame(
    bearing_corr_deg = bearings,
    n_local_pairs = n_pairs
  )
}

#' Jammalamadaka–Sarma circular–circular correlation via circular::cor.circular.
#'
#' Bearings are axial (along-shelf lines): fold to [0,180), convert to radians,
#' double the angle (standard axial→circular map), then call cor.circular.
jammalamadaka_sarma_test <- function(bearing_a_deg, bearing_b_deg) {
  if (!requireNamespace("circular", quietly = TRUE)) {
    stop("Package 'circular' required")
  }
  ok <- is.finite(bearing_a_deg) & is.finite(bearing_b_deg)
  a <- fold_axial_deg(bearing_a_deg[ok])
  b <- fold_axial_deg(bearing_b_deg[ok])
  n <- length(a)
  if (n < 10L) {
    return(list(
      n = n, rho = NA_real_, statistic = NA_real_, p.value = NA_real_,
      note = "too few paired finite bearings"
    ))
  }
  # Axial → circular: double the angle
  a_circ <- circular::circular(2 * a * pi / 180, type = "angles", units = "radians",
                               template = "none", modulo = "2pi", zero = 0, rotation = "counter")
  b_circ <- circular::circular(2 * b * pi / 180, type = "angles", units = "radians",
                               template = "none", modulo = "2pi", zero = 0, rotation = "counter")
  fit <- circular::cor.circular(a_circ, b_circ, test = TRUE)
  # cor.circular with test=TRUE returns a list / htest-like object depending on version
  rho <- if (is.list(fit) && !is.null(fit$cor)) {
    as.numeric(fit$cor)
  } else if (is.numeric(fit)) {
    as.numeric(fit)
  } else {
    as.numeric(fit[[1]])
  }
  pval <- if (is.list(fit) && !is.null(fit$p.value)) {
    as.numeric(fit$p.value)
  } else {
    NA_real_
  }
  stat <- if (is.list(fit) && !is.null(fit$statistic)) {
    as.numeric(fit$statistic)
  } else {
    NA_real_
  }
  list(n = n, rho = rho, statistic = stat, p.value = pval, note = "axial-doubled")
}

#' Permutation p-value for Jammalamadaka–Sarma circular correlation.
#'
#' Shuffles the pairing of bearing_b relative to bearing_a (999+ draws).
#' Two-sided: fraction of |rho_perm| >= |rho_obs|.
jammalamadaka_sarma_permute <- function(
    bearing_a_deg,
    bearing_b_deg,
    n_perm = 999L,
    seed = 20260804L
) {
  ok <- is.finite(bearing_a_deg) & is.finite(bearing_b_deg)
  a <- bearing_a_deg[ok]
  b <- bearing_b_deg[ok]
  n <- length(a)
  obs <- jammalamadaka_sarma_test(a, b)
  if (!is.finite(obs$rho) || n < 10L) {
    return(list(
      n = n, rho = obs$rho, p_asymptotic = obs$p.value,
      p_permutation = NA_real_, n_perm = n_perm, note = obs$note
    ))
  }
  set.seed(seed)
  rho_perm <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    rho_perm[k] <- jammalamadaka_sarma_test(a, sample(b))$rho
  }
  p_perm <- mean(abs(rho_perm) >= abs(obs$rho) - 1e-15)
  # include observed among nulls (common conservative convention optional);
  # brief asks shuffle-based p — use strict count of extremes / n_perm
  list(
    n = n,
    rho = obs$rho,
    statistic = obs$statistic,
    p_asymptotic = obs$p.value,
    p_permutation = p_perm,
    n_perm = n_perm,
    note = paste0(obs$note, "; perm=", n_perm)
  )
}

#' Global compass along-shelf foil (depth-independent).
#'
#' North Sea piecewise compass model: N–S (90°) west of 4°E; E–W (0°) at/ east
#' of 4°E. Not confirmatory — contrast only.
global_compass_along_shelf_deg <- function(lon) {
  ifelse(lon < 4, 90, 0)
}

#' Plot residual-correlation bearing vs along-shelf depth axis.
save_bearing_comparison_plot <- function(bearing_tab, path) {
  grDevices::png(path, width = 900, height = 800, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)

  x <- fold_axial_deg(bearing_tab$bearing_depth)
  y <- fold_axial_deg(bearing_tab$bearing_resid)
  flagged <- as.logical(bearing_tab$tid_flag_lt50)
  flagged[is.na(flagged)] <- FALSE

  plot(
    NA, NA,
    xlim = c(0, 180), ylim = c(0, 180),
    xlab = "Local along-shelf depth axis (deg, axial; = depth-gradient + 90)",
    ylab = "Residual-correlation bearing (deg, axial)",
    main = "Residual vs depth-derived along-shelf bearing"
  )
  abline(0, 1, col = "grey70", lty = 2)
  # also show axial wrap of the 1:1 line at ±180 (same axis)
  points(
    x[!flagged], y[!flagged],
    pch = 16, col = grDevices::adjustcolor("#1b4f72", 0.75)
  )
  points(
    x[flagged], y[flagged],
    pch = 17, col = "#b03a2e", cex = 1.2
  )
  legend(
    "topleft",
    legend = c(
      sprintf("TID OK (n=%d)", sum(!flagged)),
      sprintf("TID flagged <50%% (n=%d; kept in plot)", sum(flagged)),
      "1:1 (aligned along-shelf axes)"
    ),
    pch = c(16, 17, NA),
    lty = c(NA, NA, 2),
    col = c("#1b4f72", "#b03a2e", "grey70"),
    bty = "n"
  )
}

save_directional_variogram_plot <- function(vgram_df, fits_df, path, title) {
  grDevices::png(path, width = 900, height = 700, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  dirs <- sort(unique(vgram_df$dir.hor))
  cols <- grDevices::hcl.colors(length(dirs), "Dark 3")
  ylim <- range(vgram_df$gamma, na.rm = TRUE)
  xlim <- range(vgram_df$dist, na.rm = TRUE)
  plot(
    NA, NA,
    xlim = xlim, ylim = ylim,
    xlab = "Lag distance (deg, lon-compressed)",
    ylab = "Semivariance",
    main = title
  )
  for (k in seq_along(dirs)) {
    d <- dirs[k]
    sub <- vgram_df[vgram_df$dir.hor == d, , drop = FALSE]
    points(sub$dist, sub$gamma, col = cols[k], pch = 16)
    lines(sub$dist[order(sub$dist)], sub$gamma[order(sub$dist)], col = cols[k])
    fr <- fits_df[fits_df$dir.hor == d & isTRUE(fits_df$fit_ok[1]), , drop = FALSE]
    if (nrow(fr) == 1L && is.finite(fr$range) && is.finite(fr$sill)) {
      xs <- seq(xlim[1], xlim[2], length.out = 100)
      ys <- fr$nugget + fr$sill * (1 - exp(-xs / fr$range))
      lines(xs, ys, col = cols[k], lty = 2)
    }
  }
  legend(
    "topleft",
    legend = paste0(dirs, "°"),
    col = cols, pch = 16, lty = 1, bty = "n", title = "dir.hor"
  )
}

#' Classify confirmatory verdict from primary + robustness p-values.
classify_anisotropy_verdict <- function(p_primary, p_robust, alpha = 0.05) {
  prim_sig <- is.finite(p_primary) && p_primary < alpha
  rob_sig <- is.finite(p_robust) && p_robust < alpha
  if (prim_sig && rob_sig) {
    list(
      verdict = "confirmed",
      next_step = paste0(
        "Replace raw lon/lat distance with an environment-informed distance ",
        "(cross-shelf weighted or depth-contour distance) in any future ",
        "spatial-covariance model. Do not re-attempt a plain isotropic model."
      )
    )
  } else if (!prim_sig) {
    list(
      verdict = "not_confirmed",
      next_step = paste0(
        "Deprioritize further distance-metric refinement; prioritize the ",
        "spatial-lag covariate test (Decision rule 2 from the Moran/BLUP run)."
      )
    )
  } else {
    # primary sig, robustness not
    list(
      verdict = "partial",
      next_step = paste0(
        "Do not treat as confirmed alignment; consider SST-front layer before ",
        "committing to a depth-based distance metric."
      )
    )
  }
}
