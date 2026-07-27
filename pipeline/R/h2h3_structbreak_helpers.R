# Structural-break check helpers — fishing-pressure time series.
# See CURSOR_BRIEFING "Structural Break Check — Fishing Pressure Time Series"
# (chat-supplied, not a repo file) for the full spec these helpers implement.
#
# REPORTING ONLY: no phase boundary (visual or statistical) is recommended
# here. The visual boundaries (1990/2000/2010, from Section A of the H2/H3
# design-support task) are NOT adjusted to fit the statistical result, and
# the statistical result is NOT forced to match the visual phases.
#
# Requires strucchange (not part of the project's renv-managed dependency
# set at the time this was written — installed ad hoc; see run log).

suppressPackageStartupMessages({
  library(dplyr)
})

#' Model-selection criterion table (RSS + BIC) for every candidate break
#' count strucchange's dynamic program evaluated (0 up to the h-determined
#' max), from a fitted `strucchange::breakpoints()` "breakpointsfull" object.
#' Returns one row per m, with `is_bic_optimal` flagging the BIC minimiser —
#' the ONLY basis for "optimal" used anywhere in this task, not chosen to
#' match any particular number of breaks in advance.
build_breakpoint_criterion_table <- function(bp_full) {
  sbp <- summary(bp_full)
  rss <- sbp$RSS
  m <- as.integer(colnames(rss))
  bic <- unname(rss["BIC", ])
  rss_val <- unname(rss["RSS", ])
  tibble::tibble(
    n_breaks = m,
    rss = rss_val,
    bic = bic,
    is_bic_optimal = bic == min(bic)
  )
}

#' Break-year table (with confidence interval, if available) for a
#' specific number of breaks extracted from a "breakpointsfull" object via
#' `strucchange::confint()`. `year_vector` maps observation index -> year
#' (same ordering as the data breakpoints() was fit on). CI columns are
#' NA if confint() cannot be computed (e.g. too few observations per
#' segment) — reported as NA, not silently dropped.
build_breakpoint_year_table <- function(bp_full, n_breaks_opt, year_vector) {
  bp_extracted <- strucchange::breakpoints(bp_full, breaks = n_breaks_opt)
  obs_idx <- as.integer(bp_extracted$breakpoints)

  ci_lower <- rep(NA_integer_, length(obs_idx))
  ci_upper <- rep(NA_integer_, length(obs_idx))
  ci_ok <- TRUE
  ci <- tryCatch(stats::confint(bp_full, breaks = n_breaks_opt), error = function(e) NULL)
  if (!is.null(ci)) {
    ci_mat <- ci$confint
    ci_lower <- as.integer(ci_mat[, "2.5 %"])
    ci_upper <- as.integer(ci_mat[, "97.5 %"])
  } else {
    ci_ok <- FALSE
  }

  tibble::tibble(
    break_number = seq_along(obs_idx),
    obs_index = obs_idx,
    break_year = year_vector[obs_idx],
    ci_lower_obs_index = ci_lower,
    ci_upper_obs_index = ci_upper,
    ci_lower_year = ifelse(is.na(ci_lower), NA_integer_, year_vector[ci_lower]),
    ci_upper_year = ifelse(is.na(ci_upper), NA_integer_, year_vector[ci_upper]),
    ci_computed = ci_ok
  )
}

#' Bidirectional nearest-neighbour distance between the statistically
#' selected break years and the visual phase boundaries — numbers only, no
#' judgement about which (if either) set is "right". Each row is one
#' directional comparison (visual->nearest stat break, and stat break->
#' nearest visual boundary), so neither list's perspective is privileged.
compare_breaks_to_visual <- function(break_years, visual_years) {
  visual_to_stat <- tibble::tibble(
    direction = "visual_boundary_to_nearest_statistical_break",
    reference_year = visual_years,
    nearest_other_year = vapply(visual_years, function(v) break_years[which.min(abs(break_years - v))], numeric(1)),
    abs_diff_years = vapply(visual_years, function(v) min(abs(break_years - v)), numeric(1))
  )
  stat_to_visual <- tibble::tibble(
    direction = "statistical_break_to_nearest_visual_boundary",
    reference_year = break_years,
    nearest_other_year = vapply(break_years, function(b) visual_years[which.min(abs(visual_years - b))], numeric(1)),
    abs_diff_years = vapply(break_years, function(b) min(abs(visual_years - b)), numeric(1))
  )
  dplyr::bind_rows(visual_to_stat, stat_to_visual)
}
