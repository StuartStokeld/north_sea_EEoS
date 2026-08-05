# Helpers for residual spatial-autocorrelation check on the primary
# H2/H3 (1 | stat_rec) model (run_h2h3_primary_spatial_autocorr_check.R).
#
# Designed to run without a live glmmTMB/spdep/sf install when needed:
# BLUPs/residuals are recovered from the saved glmmTMB artifact; queen
# weights are rebuilt from the same ICES rectangle DBF fields (SOUTH/WEST)
# used by the shapefile; Moran/Geary match the Cliff–Ord / spdep formulas
# used by h2_global_spatial_tests() (zero.policy / style W / greater).

ALPHA_SPATIAL <- 0.05

#' Load the three archived Moran/Geary rows from H2c verbatim.
load_archived_spatial_diagnostics <- function(path_csv) {
  raw <- utils::read.csv(path_csv, stringsAsFactors = FALSE)
  required <- c("test", "statistic", "expected", "p_value", "variable", "n_rectangles")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) {
    stop("Archived diagnostics missing columns: ", paste(missing, collapse = ", "))
  }

  vars_wanted <- c(
    "mean_abs_residual",
    "log_mean_annual_hours_total",
    "ols_primary_abs_residuals"
  )
  stage_map <- c(
    mean_abs_residual = "Raw DV",
    log_mean_annual_hours_total = "Raw IV",
    ols_primary_abs_residuals = "OLS residuals"
  )

  rows <- lapply(vars_wanted, function(v) {
    sub <- raw[raw$variable == v, , drop = FALSE]
    moran <- sub[sub$test == "morans_i", , drop = FALSE]
    geary <- sub[sub$test == "geary_c", , drop = FALSE]
    if (nrow(moran) != 1L || nrow(geary) != 1L) {
      stop("Expected one Moran's I and one Geary's C row for variable '", v, "'")
    }
    data.frame(
      stage = unname(stage_map[[v]]),
      variable = v,
      morans_i = moran$statistic[[1L]],
      geary_c = geary$statistic[[1L]],
      p_moran = moran$p_value[[1L]],
      p_geary = geary$p_value[[1L]],
      expected_moran = moran$expected[[1L]],
      expected_geary = geary$expected[[1L]],
      variance_moran = NA_real_,
      z_moran = NA_real_,
      n = as.integer(moran$n_rectangles[[1L]]),
      source = "archived_h2_spatial_diagnostics.csv",
      archived_comparable_residual = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Document how archived ols_primary_abs_residuals was constructed in H2c.
archived_ols_residual_definition <- function() {
  paste0(
    "FINDING (from pipeline/run_h2_models.R + outputs/h2_spatial_diagnostics.csv): ",
    "`ols_primary_abs_residuals` is residuals(lm(mean_abs_residual ~ mean_annual_hours_total)), ",
    "i.e. signed OLS residuals from a regression whose dependent variable was already ",
    "rectangle-mean absolute EEoS residual. It is NOT abs(signed residual). The 'abs' in the ",
    "variable name refers to the DV (mean_abs_residual), not to an absolute-value transform ",
    "of the OLS residuals themselves. Closest role-match among the new primary-model residual ",
    "rows is therefore the signed rectangle-mean response residual (leftover after the fitted ",
    "model), not the mean-absolute collapse — despite the archived variable name."
  )
}

#' Recover BLUPs from a saved glmmTMB fit without loading glmmTMB methods.
extract_primary_blups_from_fit <- function(fit, panel_stat_rec) {
  lp <- fit$obj$env$last.par.best
  b <- unname(lp[names(lp) == "b"])
  levs <- levels(fit$modelInfo$reTrms$cond$flist$stat_rec)
  if (length(b) != length(levs)) {
    stop("BLUP length (", length(b), ") != RE levels (", length(levs), ")")
  }
  blups <- data.frame(
    stat_rec = normalize_stat_rec(levs),
    blup = as.numeric(b),
    stringsAsFactors = FALSE
  )
  panel_ids <- normalize_stat_rec(panel_stat_rec)
  if (nrow(blups) != length(panel_ids)) {
    stop("BLUP count (", nrow(blups), ") != archived panel N (", length(panel_ids), ")")
  }
  if (!setequal(blups$stat_rec, panel_ids)) {
    stop("BLUP rectangle IDs do not match the archived H2 panel set exactly.")
  }
  blups <- blups[match(panel_ids, blups$stat_rec), , drop = FALSE]
  if (anyNA(blups$blup) || anyNA(blups$stat_rec)) {
    stop("NA introduced when aligning BLUPs to archived panel IDs — aborting.")
  }
  rownames(blups) <- NULL
  blups
}

#' Conditional response residuals: y - (Xβ + b_i), matching
#' residuals(glmmTMB_fit, type = "response").
extract_primary_resid_by_rect_from_fit <- function(fit, dat, panel_stat_rec) {
  dat <- as.data.frame(dat)
  dat$stat_rec <- normalize_stat_rec(dat$stat_rec)
  panel_ids <- normalize_stat_rec(panel_stat_rec)

  lp <- fit$obj$env$last.par.best
  beta <- unname(lp[names(lp) == "beta"])
  b <- unname(lp[names(lp) == "b"])
  levs <- normalize_stat_rec(levels(fit$modelInfo$reTrms$cond$flist$stat_rec))
  b_map <- stats::setNames(b, levs)

  X <- stats::model.matrix(~ FP_between * phase + FP_within * phase, data = dat)
  if (ncol(X) != length(beta)) {
    stop("model.matrix columns (", ncol(X), ") != beta length (", length(beta), ")")
  }
  # Align beta to X columns using saved FE table order verified at development time;
  # colnames(X) match glmmTMB beta order for this formula.
  eta <- as.numeric(X %*% beta) + as.numeric(b_map[dat$stat_rec])
  if (anyNA(eta)) {
    stop("NA in linear predictor — rectangle ID mismatch between data and RE levels.")
  }
  resid_vec <- as.numeric(dat$residual) - eta

  resid_df <- data.frame(
    stat_rec = dat$stat_rec,
    year = dat$year,
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

  if (nrow(resid_by_rect) != length(panel_ids)) {
    stop("resid_by_rect N (", nrow(resid_by_rect), ") != archived panel N (", length(panel_ids), ")")
  }
  if (anyNA(resid_by_rect$resid) || anyNA(resid_by_rect$abs_resid)) {
    stop("NA in rectangle-collapsed residuals — aborting.")
  }
  list(resid_df = resid_df, resid_by_rect = resid_by_rect)
}

#' Queen contiguity neighbour list from ICES DBF SOUTH/WEST grid indices.
#' For the regular 0.5°×1° ICES lattice this matches spdep::poly2nb(..., queen=TRUE).
build_queen_nb_from_ices_dbf <- function(panel, project_root) {
  if (!requireNamespace("foreign", quietly = TRUE)) {
    stop("Package 'foreign' required to read ICES rectangle DBF without sf.")
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
  gi$row_idx <- as.integer(round(gi$SOUTH / 0.5))
  gi$col_idx <- as.integer(round(gi$WEST / 1))

  key <- paste(gi$row_idx, gi$col_idx, sep = "_")
  lookup <- stats::setNames(seq_along(panel_ids), key)

  offsets <- expand.grid(dr = -1:1, dc = -1:1)
  offsets <- offsets[!(offsets$dr == 0 & offsets$dc == 0), , drop = FALSE]

  nb <- vector("list", length(panel_ids))
  for (i in seq_along(panel_ids)) {
    cand_keys <- paste(gi$row_idx[i] + offsets$dr, gi$col_idx[i] + offsets$dc, sep = "_")
    nbr <- unname(lookup[cand_keys])
    nbr <- sort(as.integer(nbr[!is.na(nbr)]))
    nb[[i]] <- nbr
  }
  names(nb) <- panel_ids
  class(nb) <- "nb"
  attr(nb, "region.id") <- panel_ids
  attr(nb, "type") <- "queen"
  attr(nb, "call") <- match.call()
  nb
}

#' Row-standardised (style W) spatial weights from an nb list.
nb_to_listw_W <- function(nb, zero.policy = TRUE) {
  n <- length(nb)
  neighbours <- nb
  weights <- vector("list", n)
  n_isolated <- 0L
  for (i in seq_len(n)) {
    k <- length(neighbours[[i]])
    if (k == 0L) {
      n_isolated <- n_isolated + 1L
      if (!zero.policy) stop("Region ", i, " has no neighbours")
      weights[[i]] <- numeric(0)
    } else {
      weights[[i]] <- rep(1 / k, k)
    }
  }
  listw <- list(style = "W", neighbours = neighbours, weights = weights)
  class(listw) <- c("listw", "nb")
  attr(listw, "region.id") <- attr(nb, "region.id")
  attr(listw, "zero.policy") <- zero.policy
  list(listw = listw, n_isolated = n_isolated, nb = nb)
}

#' Spatially lagged values for style-W listw.
lag_listw <- function(listw, x) {
  x <- as.numeric(x)
  n <- length(x)
  out <- numeric(n)
  for (i in seq_len(n)) {
    nbr <- listw$neighbours[[i]]
    if (length(nbr) == 0L) {
      out[i] <- 0
    } else {
      out[i] <- sum(listw$weights[[i]] * x[nbr])
    }
  }
  out
}

listw_S0 <- function(listw) {
  sum(vapply(listw$weights, sum, numeric(1)))
}

listw_S1 <- function(listw) {
  n <- length(listw$neighbours)
  # Build sparse W then S1 = 0.5 * sum_ij (w_ij + w_ji)^2
  # Efficient double loop over existing edges.
  s1 <- 0
  for (i in seq_len(n)) {
    nbr_i <- listw$neighbours[[i]]
    w_i <- listw$weights[[i]]
    if (length(nbr_i) == 0L) next
    for (k in seq_along(nbr_i)) {
      j <- nbr_i[k]
      wij <- w_i[k]
      # find w_ji
      nbr_j <- listw$neighbours[[j]]
      wji <- 0
      if (length(nbr_j) > 0L) {
        pos <- match(i, nbr_j)
        if (!is.na(pos)) wji <- listw$weights[[j]][pos]
      }
      s1 <- s1 + (wij + wji)^2
    }
  }
  0.5 * s1
}

listw_S2 <- function(listw) {
  n <- length(listw$neighbours)
  # For each i: (row sum_i + col sum_i)^2
  row_sum <- vapply(listw$weights, sum, numeric(1))
  col_sum <- numeric(n)
  for (i in seq_len(n)) {
    nbr <- listw$neighbours[[i]]
    w <- listw$weights[[i]]
    if (length(nbr) == 0L) next
    col_sum[nbr] <- col_sum[nbr] + w
  }
  sum((row_sum + col_sum)^2)
}

#' Moran's I — Cliff–Ord randomization moments; alternative = "greater".
#' Matches spdep::moran.test(..., zero.policy=TRUE) for finite x / style W.
moran_test_listw <- function(x, listw) {
  x <- as.numeric(x)
  n <- length(x)
  z <- x - mean(x)
  S0 <- listw_S0(listw)
  num <- sum(z * lag_listw(listw, z))
  den <- sum(z^2)
  I <- (n / S0) * (num / den)
  EI <- -1 / (n - 1)

  S1 <- listw_S1(listw)
  S2 <- listw_S2(listw)
  k <- (sum(z^4) / n) / (den / n)^2
  # Randomization variance (Cliff & Ord / spdep)
  num1 <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
  num2 <- k * ((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)
  den_v <- (n - 1) * (n - 2) * (n - 3) * S0^2
  VI <- (num1 - num2) / den_v - EI^2

  z_score <- (I - EI) / sqrt(VI)
  p_greater <- stats::pnorm(z_score, lower.tail = FALSE)
  list(
    statistic = I,
    expectation = EI,
    variance = VI,
    z = z_score,
    p.value = p_greater
  )
}

#' Geary's C — Cliff–Ord randomization moments; alternative = "greater"
#' (spdep default: greater means C smaller than expectation under positive SA).
geary_test_listw <- function(x, listw) {
  x <- as.numeric(x)
  n <- length(x)
  z <- x - mean(x)
  S0 <- listw_S0(listw)
  # C = ((n-1)/(2*S0)) * sum_ij w_ij (x_i - x_j)^2 / sum_i z_i^2
  num <- 0
  for (i in seq_len(n)) {
    nbr <- listw$neighbours[[i]]
    w <- listw$weights[[i]]
    if (length(nbr) == 0L) next
    num <- num + sum(w * (x[i] - x[nbr])^2)
  }
  den <- sum(z^2)
  C <- ((n - 1) / (2 * S0)) * (num / den)
  EC <- 1

  S1 <- listw_S1(listw)
  S2 <- listw_S2(listw)
  # spdep geary.test randomization variance
  k <- (sum(z^4) / n) / (den / n)^2
  n2 <- n^2
  # Formula from Cliff & Ord / spdep geary.test source
  VC <- (1 / (2 * (n + 1) * S0^2)) *
    ((n - 1) * (2 * (2 * n2 - 6 * n + 3) * S1 - (n - 1) * (n2 - n + 3) * S2 +
      3 * (n - 1)^2 * S0^2) -
      k * ((n^2 - n) * S1 - 2 * (n^2 - 2 * n + 3) * S0^2 + (n - 1) * (n2 - n + 3) * S2 -
        3 * (n - 1)^2 * S0^2) / ((n - 2) * (n - 3))) -
    EC^2
  # Fallback if VC unstable: use simpler form from spdep when needed
  if (!is.finite(VC) || VC <= 0) {
    VC <- ((2 * S1 + S2) * (n - 1) - 4 * S0^2) / (2 * (n + 1) * S0^2)
  }

  # For alternative="greater" (positive SA), C is expected smaller than 1,
  # so the test statistic is standardized toward lower tail of C... but
  # spdep reports p for "greater" as P(C <= observed) when testing positive SA?
  # Empirically: spdep geary.test alternative="greater" uses
  # pnorm((C - EC)/sqrt(VC), lower.tail = TRUE) for the usual positive-SA reading.
  # Check against archived values in the runner validation step.
  z_score <- (C - EC) / sqrt(VC)
  p_greater <- stats::pnorm(z_score, lower.tail = TRUE)
  list(
    statistic = C,
    expectation = EC,
    variance = VC,
    z = z_score,
    p.value = p_greater
  )
}

run_spatial_tests <- function(x, listw, label) {
  moran <- moran_test_listw(x, listw)
  geary <- geary_test_listw(x, listw)
  list(
    label = label,
    moran = moran,
    geary = geary,
    row = data.frame(
      variable = label,
      morans_i = moran$statistic,
      geary_c = geary$statistic,
      p_moran = moran$p.value,
      p_geary = geary$p.value,
      expected_moran = moran$expectation,
      expected_geary = geary$expectation,
      variance_moran = moran$variance,
      z_moran = moran$z,
      n = length(x),
      source = "primary_model_check",
      stringsAsFactors = FALSE
    )
  )
}

classify_spatial_decision <- function(p_blup_moran, i_blup,
                                     p_resid_primary_moran, i_resid_primary,
                                     p_resid_secondary_moran, i_resid_secondary,
                                     primary_resid_label, secondary_resid_label,
                                     alpha = ALPHA_SPATIAL) {
  blup_sig <- isTRUE(p_blup_moran < alpha)
  resid_pri_sig <- isTRUE(p_resid_primary_moran < alpha)
  resid_sec_sig <- isTRUE(p_resid_secondary_moran < alpha)

  if (blup_sig) {
    rule_id <- 2L
    text <- paste0(
      "Exchangeability directly falsified: rectangle intercepts are themselves spatially ",
      "clustered. Supports prioritizing a spatial-lag covariate over further error-covariance ",
      "modelling."
    )
  } else if (resid_pri_sig) {
    rule_id <- 3L
    text <- paste0(
      "Remaining structure is not coming from the rectangle-level mean. Investigate omitted ",
      "spatially-varying covariates (start with fishing-effort confounding) or ",
      "within-rectangle temporal-spatial interaction, rather than the RE structure itself."
    )
  } else {
    rule_id <- 1L
    text <- paste0(
      "Plain RE captures spatial pattern adequately in practice — exchangeability is a ",
      "practical simplification, not a mechanistic fit. Report as an explicit limitation."
    )
  }

  disagree <- xor(resid_pri_sig, resid_sec_sig)
  disagreement_note <- if (disagree) {
    paste0(
      "SIGNED vs ABSOLUTE residual tests DISAGREE at alpha=", alpha, ": ",
      primary_resid_label, " Moran p=", signif(p_resid_primary_moran, 4),
      " (I=", signif(i_resid_primary, 4), "); ",
      secondary_resid_label, " Moran p=", signif(p_resid_secondary_moran, 4),
      " (I=", signif(i_resid_secondary, 4), "). ",
      "Error magnitude and error direction are telling different spatial stories — ",
      "do not collapse to one residual conclusion."
    )
  } else {
    paste0(
      "Signed and absolute residual Moran tests agree on significance at alpha=", alpha,
      " (both ", if (resid_pri_sig) "significant" else "non-significant", ")."
    )
  }

  list(
    rule_id = rule_id,
    classification = text,
    blup_significant = blup_sig,
    primary_residual_significant = resid_pri_sig,
    secondary_residual_significant = resid_sec_sig,
    residual_tests_disagree = disagree,
    disagreement_note = disagreement_note,
    blup_moran_i = i_blup,
    primary_resid_moran_i = i_resid_primary
  )
}
