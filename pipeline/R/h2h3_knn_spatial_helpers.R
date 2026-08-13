# k-NN spatial weights and FP_between spatial-lag helpers (Spec A).
# Requires h2_common.R, h2_spatial_helpers.R, h2h3_spatial_autocorr_helpers.R

KNN_K_PRIMARY <- 4L
KNN_TIE_TOL_KM <- 1e-6

#' Great-circle distance (km) between two WGS84 points (degrees).
haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  lat1r <- lat1 * rad
  lat2r <- lat2 * rad
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

#' Pairwise great-circle distance matrix (km); diagonal set to Inf (exclude self).
great_circle_dist_matrix_km <- function(lon, lat) {
  n <- length(lon)
  d <- matrix(Inf, n, n)
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      dij <- haversine_km(lon[i], lat[i], lon[j], lat[j])
      d[i, j] <- dij
      d[j, i] <- dij
    }
  }
  d
}

#' k-NN neighbour list with explicit tie-breaking: rank by (distance, stat_rec).
#'
#' Returns list with nb (spdep-style), audit data.frame, and coordinates used.
build_knn_nb_explicit <- function(stat_rec,
                                  lon,
                                  lat,
                                  k = KNN_K_PRIMARY,
                                  tie_tol_km = KNN_TIE_TOL_KM) {
  stat_rec <- normalize_stat_rec(stat_rec)
  n <- length(stat_rec)
  if (length(lon) != n || length(lat) != n) {
    stop("stat_rec, lon, lat must have equal length.")
  }
  if (k >= n) {
    stop("k must be smaller than the number of rectangles (k=", k, ", n=", n, ").")
  }

  dist <- great_circle_dist_matrix_km(lon, lat)
  nb <- vector("list", n)
  audit <- vector("list", n)

  for (i in seq_len(n)) {
    candidates <- setdiff(seq_len(n), i)
    ranks <- data.frame(
      j = candidates,
      dist_km = dist[i, candidates],
      stat_rec_j = stat_rec[candidates],
      stringsAsFactors = FALSE
    )
    ranks <- ranks[order(ranks$dist_km, ranks$stat_rec_j), , drop = FALSE]
    nbr_idx <- ranks$j[seq_len(k)]
    nb[[i]] <- as.integer(nbr_idx)

    tie_at_k <- FALSE
    dist_k <- ranks$dist_km[k]
    if (nrow(ranks) > k) {
      dist_k1 <- ranks$dist_km[k + 1L]
      tie_at_k <- abs(dist_k - dist_k1) <= tie_tol_km
    }

    audit[[i]] <- data.frame(
      stat_rec = stat_rec[i],
      rect_lon = lon[i],
      rect_lat = lat[i],
      k = k,
      n_neighbours = length(nbr_idx),
      neighbour_stat_rec = paste(stat_rec[nbr_idx], collapse = "|"),
      neighbour_dist_km = paste(sprintf("%.6f", dist[i, nbr_idx]), collapse = "|"),
      dist_kth_km = dist_k,
      dist_kplus1_km = if (nrow(ranks) > k) ranks$dist_km[k + 1L] else NA_real_,
      tie_at_k = tie_at_k,
      tie_tol_km = tie_tol_km,
      stringsAsFactors = FALSE
    )
  }

  names(nb) <- stat_rec
  class(nb) <- "nb"
  attr(nb, "region.id") <- stat_rec
  attr(nb, "type") <- "knn"
  attr(nb, "k") <- k
  attr(nb, "call") <- match.call()

  list(
    nb = nb,
    audit = do.call(rbind, audit),
    stat_rec = stat_rec,
    lon = lon,
    lat = lat,
    k = k,
    tie_tol_km = tie_tol_km
  )
}

#' Rectangle centroids on the ICES 0.5°×1° lattice (WGS84 degrees).
#'
#' Uses the same ICES DBF (SOUTH/WEST) as build_queen_nb_from_ices_dbf().
#' Centroid = southwest corner + half cell dimensions (equivalent to
#' st_centroid() on axis-aligned ICES rectangle polygons).
panel_rectangle_centroids <- function(panel, project_root) {
  if (!requireNamespace("foreign", quietly = TRUE)) {
    stop("Package 'foreign' required to read ICES rectangle DBF.")
  }
  dbf_path <- file.path(
    project_root, "gis", "ICES_rectangles", "ICES_Statistical_Rectangles_Eco.dbf"
  )
  if (!file.exists(dbf_path)) {
    stop("ICES DBF not found: ", dbf_path)
  }
  dbf <- foreign::read.dbf(dbf_path)
  dbf$stat_rec <- normalize_stat_rec(dbf$ICESNAME)
  panel_ids <- normalize_stat_rec(panel$stat_rec)
  gi <- dbf[dbf$stat_rec %in% panel_ids, c("stat_rec", "SOUTH", "WEST"), drop = FALSE]
  gi <- gi[!duplicated(gi$stat_rec), , drop = FALSE]
  if (nrow(gi) != length(panel_ids)) {
    stop(
      "Could not match all panel rectangles to ICES DBF (matched ",
      nrow(gi), " of ", length(panel_ids), ")"
    )
  }
  gi <- gi[match(panel_ids, gi$stat_rec), , drop = FALSE]
  list(
    stat_rec = panel_ids,
    rect_lon = gi$WEST + 0.5,
    rect_lat = gi$SOUTH + 0.25,
    source = "ICES_DBF_SOUTH_WEST_centroid"
  )
}

#' Build row-standardised k-NN listw (style W) for the H2 panel.
build_knn_spatial_weights <- function(panel,
                                      project_root,
                                      k = KNN_K_PRIMARY,
                                      tie_tol_km = KNN_TIE_TOL_KM) {
  geo <- panel_rectangle_centroids(panel, project_root)
  knn_pack <- build_knn_nb_explicit(
    geo$stat_rec, geo$rect_lon, geo$rect_lat,
    k = k, tie_tol_km = tie_tol_km
  )
  w_pack <- nb_to_listw_W(knn_pack$nb, zero.policy = TRUE)
  list(
    listw = w_pack$listw,
    nb = w_pack$nb,
    n_isolated = w_pack$n_isolated,
    audit = knn_pack$audit,
    centroids = data.frame(
      stat_rec = geo$stat_rec,
      rect_lon = geo$rect_lon,
      rect_lat = geo$rect_lat,
      stringsAsFactors = FALSE
    ),
    k = k,
    tie_tol_km = tie_tol_km,
    crs = "GCS_WGS_1984 (geographic degrees; great-circle distances in km)"
  )
}

#' Row-standardised spatial lag: mean of neighbours' values (self excluded).
compute_spatial_lag <- function(values, nb) {
  values <- as.numeric(values)
  n <- length(values)
  if (length(nb) != n) {
    stop("values length (", n, ") != nb length (", length(nb), ").")
  }
  out <- numeric(n)
  for (i in seq_len(n)) {
    nbr <- nb[[i]]
    if (length(nbr) == 0L) {
      out[i] <- NA_real_
    } else {
      out[i] <- mean(values[nbr])
    }
  }
  out
}

#' Spatial lag averaging over available neighbours only (NA-aware).
#'
#' Returns list(lag, n_neighbours_used). When some neighbours are NA, the lag
#' is the mean of the finite ones (re-normalised); n_neighbours_used is their
#' count. Lag is NA only when no neighbours are available.
compute_spatial_lag_available <- function(values, nb) {
  values <- as.numeric(values)
  n <- length(values)
  if (length(nb) != n) {
    stop("values length (", n, ") != nb length (", length(nb), ").")
  }
  lag <- rep(NA_real_, n)
  n_used <- integer(n)
  for (i in seq_len(n)) {
    nbr <- nb[[i]]
    if (length(nbr) == 0L) {
      n_used[i] <- 0L
      next
    }
    v <- values[nbr]
    ok <- is.finite(v)
    n_used[i] <- sum(ok)
    if (n_used[i] > 0L) {
      lag[i] <- mean(v[ok])
    }
  }
  list(lag = lag, n_neighbours_used = n_used)
}

#' Rectangle-year mean ln_B_obs + B_lag_neighbour from panel-universe hauls.
#'
#' @param haul haul_eeos_predictions (or equivalent) with stat_rec, year, ln_B_obs
#' @param panel_stat_rec character vector of 158 panel IDs (region.id / panel order)
#' @param nb k-NN neighbour list aligned to panel_stat_rec
#' @param year_min, year_max inclusive year window
build_b_lag_neighbour_rectangle_year <- function(haul,
                                                panel_stat_rec,
                                                nb,
                                                year_min = 1985L,
                                                year_max = 2015L) {
  haul <- as.data.frame(haul)
  required <- c("stat_rec", "year", "ln_B_obs")
  missing <- setdiff(required, names(haul))
  if (length(missing) > 0L) {
    stop("haul missing columns: ", paste(missing, collapse = ", "))
  }
  panel_ids <- normalize_stat_rec(panel_stat_rec)
  region_id <- normalize_stat_rec(attr(nb, "region.id"))
  if (!identical(panel_ids, region_id)) {
    stop("panel_stat_rec must match nb region.id order exactly.")
  }

  h <- haul
  h$stat_rec <- normalize_stat_rec(h$stat_rec)
  h <- h[
    h$stat_rec %in% panel_ids &
      is.finite(h$ln_B_obs) &
      h$year >= year_min &
      h$year <= year_max,
    ,
    drop = FALSE
  ]

  split_keys <- paste(h$stat_rec, h$year, sep = "\r")
  split_idx <- split(seq_len(nrow(h)), split_keys)
  ry <- do.call(rbind, lapply(names(split_idx), function(k) {
    parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    ii <- split_idx[[k]]
    data.frame(
      stat_rec = parts[[1]],
      year = as.integer(parts[[2]]),
      mean_ln_B_obs_ry = mean(h$ln_B_obs[ii]),
      n_hauls = length(ii),
      stringsAsFactors = FALSE
    )
  }))

  years <- seq.int(year_min, year_max)
  out_rows <- vector("list", length(years))
  for (yi in seq_along(years)) {
    yr <- years[[yi]]
    vals <- rep(NA_real_, length(panel_ids))
    n_hauls <- rep(0L, length(panel_ids))
    sub <- ry[ry$year == yr, , drop = FALSE]
    if (nrow(sub) > 0L) {
      idx <- match(normalize_stat_rec(sub$stat_rec), panel_ids)
      ok <- !is.na(idx)
      vals[idx[ok]] <- sub$mean_ln_B_obs_ry[ok]
      n_hauls[idx[ok]] <- as.integer(sub$n_hauls[ok])
    }
    lag_pack <- compute_spatial_lag_available(vals, nb)
    out_rows[[yi]] <- data.frame(
      stat_rec = panel_ids,
      year = yr,
      mean_ln_B_obs_ry = vals,
      n_hauls_ry = n_hauls,
      B_lag_neighbour = lag_pack$lag,
      n_neighbours_used = lag_pack$n_neighbours_used,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out
}

#' Join B_lag_neighbour onto haul-level data by (stat_rec, year).
join_b_lag_neighbour_hauls <- function(dat, lag_ry) {
  dat <- as.data.frame(dat)
  lag_map <- lag_ry[, c("stat_rec", "year", "B_lag_neighbour", "n_neighbours_used"),
                    drop = FALSE]
  lag_map$stat_rec <- normalize_stat_rec(lag_map$stat_rec)
  dat$stat_rec_chr <- normalize_stat_rec(dat$stat_rec)
  key_dat <- paste(dat$stat_rec_chr, dat$year, sep = "\r")
  key_lag <- paste(lag_map$stat_rec, lag_map$year, sep = "\r")
  idx <- match(key_dat, key_lag)
  if (anyNA(idx)) {
    stop(
      "Haul (stat_rec, year) not found in B_lag_neighbour table (n missing = ",
      sum(is.na(idx)), ")."
    )
  }
  dat$B_lag_neighbour <- lag_map$B_lag_neighbour[idx]
  dat$n_neighbours_used <- lag_map$n_neighbours_used[idx]
  dat$stat_rec_chr <- NULL
  if (anyNA(dat$B_lag_neighbour)) {
    stop("B_lag_neighbour is NA for some analysis hauls — unexpected under locked NA policy.")
  }
  dat
}

#' Pairwise VIF among continuous predictors (1/(1-R²) from each ~ others).
vif_continuous_terms <- function(data, terms) {
  data <- as.data.frame(data)
  missing <- setdiff(terms, names(data))
  if (length(missing) > 0L) {
    stop("vif_continuous_terms: missing columns: ", paste(missing, collapse = ", "))
  }
  X <- data[, terms, drop = FALSE]
  keep <- stats::complete.cases(X)
  X <- X[keep, , drop = FALSE]
  if (nrow(X) < 3L) stop("Too few complete rows for VIF.")
  rows <- lapply(terms, function(t) {
    others <- setdiff(terms, t)
    if (length(others) == 0L) {
      return(data.frame(term = t, vif = 1, r_squared_others = 0, n = nrow(X),
                        stringsAsFactors = FALSE))
    }
    fml <- stats::as.formula(paste(t, "~", paste(others, collapse = " + ")))
    fit <- stats::lm(fml, data = X)
    r2 <- summary(fit)$r.squared
    data.frame(
      term = t,
      vif = 1 / (1 - r2),
      r_squared_others = r2,
      n = nrow(X),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Rectangle-mean response residuals from a fitted glmmTMB model (generic formula).
#'
#' Uses the model's own frame for `stat_rec` / `year` alignment so row order
#' matches `residuals(fit)` regardless of any separately passed data copy.
extract_resid_by_rect_glmmTMB <- function(fit, dat = NULL, panel_stat_rec) {
  panel_ids <- normalize_stat_rec(panel_stat_rec)
  resid_vec <- as.numeric(stats::residuals(fit, type = "response"))
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  if (is.null(mf) || !"stat_rec" %in% names(mf)) {
    if (is.null(dat)) stop("Need dat when model.frame(fit) lacks stat_rec.")
    dat <- as.data.frame(dat)
    if (length(resid_vec) != nrow(dat)) {
      stop("residuals length (", length(resid_vec), ") != nrow(dat) (", nrow(dat), ")")
    }
    mf <- dat
  }
  if (length(resid_vec) != nrow(mf)) {
    stop("residuals length (", length(resid_vec), ") != model frame n (", nrow(mf), ")")
  }
  resid_df <- data.frame(
    stat_rec = normalize_stat_rec(mf$stat_rec),
    year = if ("year" %in% names(mf)) mf$year else NA_integer_,
    resid = resid_vec,
    stringsAsFactors = FALSE
  )
  split_ids <- split(seq_len(nrow(resid_df)), resid_df$stat_rec)
  resid_by_rect <- do.call(rbind, lapply(panel_ids, function(id) {
    idx <- split_ids[[id]]
    if (is.null(idx) || length(idx) == 0L) {
      stop("No hauls for rectangle ", id)
    }
    r <- resid_df$resid[idx]
    data.frame(
      stat_rec = id,
      resid = mean(r),
      abs_resid = mean(abs(r)),
      n_years = length(unique(resid_df$year[idx])),
      n_hauls = length(idx),
      stringsAsFactors = FALSE
    )
  }))
  rownames(resid_by_rect) <- NULL
  list(resid_df = resid_df, resid_by_rect = resid_by_rect)
}

#' Build rectangle-level FP_between_lag table from FP_between and k-NN nb.
build_fp_between_lag_rectangle <- function(fp_rectangle, nb) {
  fp_rectangle <- as.data.frame(fp_rectangle)
  required <- c("stat_rec", "FP_between")
  missing <- setdiff(required, names(fp_rectangle))
  if (length(missing) > 0L) {
    stop("fp_rectangle missing columns: ", paste(missing, collapse = ", "))
  }
  fp_rectangle$stat_rec <- normalize_stat_rec(fp_rectangle$stat_rec)
  region_id <- normalize_stat_rec(attr(nb, "region.id"))
  if (!identical(fp_rectangle$stat_rec, region_id)) {
    fp_rectangle <- fp_rectangle[match(region_id, fp_rectangle$stat_rec), , drop = FALSE]
    if (anyNA(fp_rectangle$stat_rec)) {
      stop("FP_between rectangle table does not align to nb region.id order.")
    }
  }

  lag <- compute_spatial_lag(fp_rectangle$FP_between, nb)
  nbr_mat <- vapply(seq_along(nb), function(i) {
    paste(normalize_stat_rec(region_id[nb[[i]]]), collapse = "|")
  }, character(1))

  out <- data.frame(
    stat_rec = region_id,
    FP_between = fp_rectangle$FP_between,
    FP_between_lag = lag,
    neighbour_stat_rec = nbr_mat,
    stringsAsFactors = FALSE
  )
  out
}

#' Join rectangle-level FP_between_lag onto haul-level data by stat_rec.
join_fp_between_lag_hauls <- function(dat, lag_rectangle) {
  dat <- as.data.frame(dat)
  lag_map <- lag_rectangle[, c("stat_rec", "FP_between_lag"), drop = FALSE]
  lag_map$stat_rec <- normalize_stat_rec(lag_map$stat_rec)
  dat$stat_rec <- normalize_stat_rec(dat$stat_rec)
  idx <- match(dat$stat_rec, lag_map$stat_rec)
  if (anyNA(idx)) {
    stop("Haul rectangles not found in FP_between_lag table.")
  }
  dat$FP_between_lag <- lag_map$FP_between_lag[idx]
  dat
}

#' Recompute FP_between_lag after a rectangle-level FP_between shuffle (bootstrap).
recompute_fp_between_lag_on_data <- function(data, nb, fp_col = "FP_between",
                                             rectangle_col = "stat_rec",
                                             lag_col = "FP_between_lag") {
  data <- as.data.frame(data)
  fp_map <- unique(data[, c(rectangle_col, fp_col), drop = FALSE])
  region_id <- normalize_stat_rec(attr(nb, "region.id"))
  fp_map[[rectangle_col]] <- normalize_stat_rec(fp_map[[rectangle_col]])
  fp_map <- fp_map[match(region_id, fp_map[[rectangle_col]]), , drop = FALSE]
  if (anyNA(fp_map[[rectangle_col]])) {
    stop("Permuted data missing rectangles required by k-NN nb.")
  }
  lag_vec <- compute_spatial_lag(fp_map[[fp_col]], nb)
  names(lag_vec) <- region_id
  data[[lag_col]] <- lag_vec[normalize_stat_rec(data[[rectangle_col]])]
  data
}

#' VIF for FP_between vs FP_between_lag at rectangle level (158 rows).
vif_fp_between_lag <- function(lag_rectangle) {
  x <- lag_rectangle$FP_between
  y <- lag_rectangle$FP_between_lag
  r <- stats::cor(x, y)
  r2 <- r^2
  vif_x <- 1 / (1 - r2)
  vif_y <- 1 / (1 - r2)
  data.frame(
    term = c("FP_between", "FP_between_lag"),
    vif = c(vif_x, vif_y),
    cor = c(r, r),
    r_squared = c(r2, r2),
    n_rectangles = nrow(lag_rectangle),
    stringsAsFactors = FALSE
  )
}
