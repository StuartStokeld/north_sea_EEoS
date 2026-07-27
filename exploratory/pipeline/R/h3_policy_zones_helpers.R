# H3 strategy feasibility helpers — policy-period split x coarse spatial zones.
# See CURSOR_BRIEFING "H3 Strategy Feasibility — Policy-Period Split x Coarse
# Spatial Zones" (chat-supplied, not a repo file) for the full spec these
# helpers implement.
#
# FEASIBILITY/VISUALISATION ONLY: no H3 model is fit here, no spatial scheme
# or parameter value is selected as final. Descriptive pre/post comparisons
# (residual, fishing pressure) are plotting aids only — no significance
# testing is performed.
#
# Requires h1_common.R, h2_common.R (normalize_stat_rec, H2_YEAR_MIN/MAX).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Grid indexing — shared foundation for both Scheme A (block merge) and
# Scheme B (rook adjacency for contiguous zone construction). Derived from
# the ICES shapefile's own SOUTH/WEST corner fields (0.5 deg lat x 1 deg lon
# grid, confirmed against the shapefile directly), not re-parsed from the
# rectangle code string.
# ---------------------------------------------------------------------------

#' One row per stat_rec with integer grid row/col indices.
#' row_idx = round(SOUTH / 0.5); col_idx = round(WEST / 1).
#' Indices are anchored to absolute lat/lon (not to the subset's own min), so
#' block/adjacency assignment is stable regardless of which rectangles are
#' included in a given analysis universe.
build_rect_grid_index <- function(ices_sf) {
  df <- as.data.frame(ices_sf)
  required <- c("stat_rec", "SOUTH", "WEST")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop("ices_sf missing columns: ", paste(missing, collapse = ", "))
  }
  df %>%
    distinct(stat_rec, SOUTH, WEST) %>%
    mutate(
      row_idx = as.integer(round(SOUTH / 0.5)),
      col_idx = as.integer(round(WEST / 1))
    ) %>%
    select(stat_rec, row_idx, col_idx)
}

#' Rook (edge-sharing) adjacency list keyed by stat_rec, restricted to
#' `universe` (a character vector of stat_rec). Diagonal (corner-only)
#' touches are deliberately excluded — a rook-contiguous region is what
#' visually reads as a single connected block/zone on the map.
build_rook_adjacency <- function(grid_index, universe) {
  gi <- grid_index %>% filter(stat_rec %in% universe)
  key <- function(r, c) paste(r, c, sep = "_")
  lookup <- setNames(gi$stat_rec, key(gi$row_idx, gi$col_idx))

  adj <- lapply(seq_len(nrow(gi)), function(i) {
    r <- gi$row_idx[i]
    c <- gi$col_idx[i]
    cand <- unname(lookup[key(c(r + 1L, r - 1L, r, r), c(c, c, c + 1L, c - 1L))])
    cand[!is.na(cand)]
  })
  setNames(adj, gi$stat_rec)
}

# ---------------------------------------------------------------------------
# Scheme A — geographic block merge (blind to fishing pressure)
# ---------------------------------------------------------------------------

#' Assign each rectangle to a block_size x block_size grid block. Purely a
#' spatial-resolution reduction; block boundaries are anchored to absolute
#' grid indices (see build_rect_grid_index), not to the universe's own extent.
build_scheme_a_blocks <- function(grid_index, universe, block_size) {
  grid_index %>%
    filter(stat_rec %in% universe) %>%
    mutate(
      block_row = row_idx %/% block_size,
      block_col = col_idx %/% block_size,
      unit_id = sprintf("A%dx%d_r%02d_c%02d", block_size, block_size, block_row, block_col)
    ) %>%
    select(stat_rec, unit_id)
}

# ---------------------------------------------------------------------------
# Scheme B — contiguous fishing-pressure zones
# ---------------------------------------------------------------------------

#' Quantile-based pressure tier (1 = lowest) via dplyr::ntile — near-equal-
#' count groups, robust to ties, standard "quantile-based" split.
assign_pressure_tiers <- function(rect_pressure, n_tiers) {
  stopifnot(all(c("stat_rec", "mean_annual_hours_total") %in% names(rect_pressure)))
  rect_pressure %>%
    filter(!is.na(mean_annual_hours_total)) %>%
    mutate(tier = dplyr::ntile(mean_annual_hours_total, n_tiers)) %>%
    select(stat_rec, mean_annual_hours_total, tier)
}

#' Connected-component clustering of same-tier rectangles, restricted to
#' rook adjacency (see build_rook_adjacency). Method: breadth-first search
#' over the tier-filtered adjacency graph — plain BFS/union-find, no external
#' graph package required (spdep/igraph are not installed in this renv and
#' are not needed for a grid of rectangles where adjacency is directly
#' derivable from row/col indices).
#'
#' Rectangles whose same-tier neighbours (if any) are not reachable become
#' singleton zones (n=1 rectangle) — this is expected, not an error, and is
#' reported by the caller via zone-size distribution, not silently pooled
#' with other same-tier-but-disconnected rectangles.
build_contiguous_zones <- function(stat_recs, tier, adjacency_list) {
  stopifnot(length(stat_recs) == length(tier))
  tier_lookup <- setNames(tier, stat_recs)
  visited <- setNames(rep(FALSE, length(stat_recs)), stat_recs)
  zone_num <- setNames(rep(NA_integer_, length(stat_recs)), stat_recs)
  current_zone <- 0L

  for (s in stat_recs) {
    if (visited[[s]]) next
    current_zone <- current_zone + 1L
    queue <- s
    visited[[s]] <- TRUE
    zone_num[[s]] <- current_zone
    while (length(queue) > 0L) {
      cur <- queue[1]
      queue <- queue[-1]
      nbrs <- adjacency_list[[cur]]
      nbrs <- nbrs[nbrs %in% stat_recs]
      for (nb in nbrs) {
        if (!visited[[nb]] && identical(tier_lookup[[nb]], tier_lookup[[cur]])) {
          visited[[nb]] <- TRUE
          zone_num[[nb]] <- current_zone
          queue <- c(queue, nb)
        }
      }
    }
  }

  zone_num_ordered <- unname(zone_num[stat_recs])
  tibble(
    stat_rec = stat_recs,
    tier = tier,
    zone_num = zone_num_ordered,
    unit_id = sprintf("B_tier%d_zone%03d", tier, zone_num_ordered)
  )
}

# ---------------------------------------------------------------------------
# Temporal period assignment (Scheme-agnostic)
# ---------------------------------------------------------------------------

#' "pre" (year <= break_year - 1) / "post" (year >= break_year) label.
#' POLICY_BREAK_YEAR = 2003 means pre = 1985-2002, post = 2003-2015.
assign_period <- function(year, break_year) {
  ifelse(year < break_year, "pre", "post")
}

# ---------------------------------------------------------------------------
# Unit x period summaries
# ---------------------------------------------------------------------------

#' Haul count per spatial unit x period. `unit_map` is stat_rec -> unit_id.
summarise_unit_period_hauls <- function(haul_full, unit_map, break_year, year_min, year_max) {
  haul_full %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(year >= year_min, year <= year_max) %>%
    inner_join(unit_map, by = "stat_rec") %>%
    mutate(period = assign_period(year, break_year)) %>%
    group_by(unit_id, period) %>%
    summarise(n_hauls = dplyr::n(), n_rects_contributing = dplyr::n_distinct(stat_rec), .groups = "drop")
}

#' Mean residual (canonical sign convention: residual = log(B_obs) - log(B_pred),
#' i.e. haul_eeos$residual as-is) per spatial unit x period.
summarise_unit_period_residual <- function(haul_eeos, unit_map, break_year, year_min, year_max) {
  stopifnot("residual" %in% names(haul_eeos))
  haul_eeos %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(year >= year_min, year <= year_max, is.finite(residual)) %>%
    inner_join(unit_map, by = "stat_rec") %>%
    mutate(period = assign_period(year, break_year)) %>%
    group_by(unit_id, period) %>%
    summarise(n_hauls = dplyr::n(), mean_residual = mean(residual, na.rm = TRUE), .groups = "drop")
}

#' Mean Couce fishing hours per spatial unit x period.
summarise_unit_period_fishing <- function(couce_year, unit_map, break_year, year_min, year_max) {
  stopifnot("hours_total" %in% names(couce_year))
  couce_year %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(year >= year_min, year <= year_max) %>%
    inner_join(unit_map, by = "stat_rec") %>%
    mutate(period = assign_period(year, break_year)) %>%
    group_by(unit_id, period) %>%
    summarise(n_rect_years = dplyr::n(), mean_fishing_hours = mean(hours_total, na.rm = TRUE), .groups = "drop")
}

# ---------------------------------------------------------------------------
# Feasibility summary (the key deliverable)
# ---------------------------------------------------------------------------

#' Wide unit_id x {pre, post} haul-count table, zero-filled.
pivot_unit_period_hauls_wide <- function(unit_period_hauls) {
  wide <- unit_period_hauls %>%
    select(unit_id, period, n_hauls) %>%
    tidyr::pivot_wider(names_from = period, values_from = n_hauls, values_fill = 0L)
  if (!"pre" %in% names(wide)) wide$pre <- 0L
  if (!"post" %in% names(wide)) wide$post <- 0L
  wide
}

#' One feasibility-summary row for a scheme/parameter combination.
#' unit_period_hauls: output of summarise_unit_period_hauls(), long
#' (unit_id, period, n_hauls). min_hauls: MIN_HAULS_PER_CELL threshold.
#' NOTE: min_hauls is ONE point on a continuum, not a committed cutoff — see
#' feasibility_across_thresholds() / cell_variance_decomposition() for the
#' data-driven, reliability-based alternative thresholds explored alongside it.
feasibility_summary_row <- function(unit_period_hauls, scheme_label, param_label, min_hauls, n_units_total) {
  wide <- pivot_unit_period_hauls_wide(unit_period_hauls)
  n_units_usable_both <- sum(wide$pre >= min_hauls & wide$post >= min_hauls)

  tibble(
    scheme = scheme_label,
    parameter = param_label,
    n_units = n_units_total,
    n_units_with_any_data = nrow(wide),
    n_units_usable_both_periods = n_units_usable_both,
    min_hauls_per_cell = min_hauls,
    pre_mean_hauls = mean(wide$pre),
    pre_median_hauls = median(wide$pre),
    pre_min_hauls = min(wide$pre),
    post_mean_hauls = mean(wide$post),
    post_median_hauls = median(wide$post),
    post_min_hauls = min(wide$post),
    overall_mean_hauls = mean(c(wide$pre, wide$post)),
    overall_median_hauls = median(c(wide$pre, wide$post)),
    overall_min_hauls = min(c(wide$pre, wide$post))
  )
}

#' Feasibility (n_units_usable_both_periods) evaluated across a VECTOR of
#' candidate thresholds, tagged by method — the "explore beyond a single
#' fixed cutoff" sensitivity table. threshold_labels (same length as
#' thresholds) documents where each value came from (e.g. "fixed_5",
#' "reliability_0.8") so the table is self-explanatory without the run log.
feasibility_across_thresholds <- function(unit_period_hauls, scheme_label, param_label,
                                          thresholds, threshold_labels, n_units_total) {
  stopifnot(length(thresholds) == length(threshold_labels))
  wide <- pivot_unit_period_hauls_wide(unit_period_hauls)

  bind_rows(lapply(seq_along(thresholds), function(i) {
    thr <- thresholds[i]
    tibble(
      scheme = scheme_label,
      parameter = param_label,
      threshold_value = thr,
      threshold_method = threshold_labels[i],
      n_units = n_units_total,
      n_units_usable_both_periods = sum(wide$pre >= thr & wide$post >= thr),
      pct_units_usable_both_periods = round(100 * sum(wide$pre >= thr & wide$post >= thr) / n_units_total, 1)
    )
  }))
}

# ---------------------------------------------------------------------------
# Statistically defensible sample-size threshold — cell-mean reliability
# ---------------------------------------------------------------------------

#' Haul-level residual rows tagged with cell_id = unit_id x period, for
#' variance-components decomposition (how much of haul-to-haul residual
#' variance is between-cell "signal" vs within-cell noise).
label_haul_cells <- function(haul_eeos, unit_map, break_year, year_min, year_max) {
  stopifnot("residual" %in% names(haul_eeos))
  haul_eeos %>%
    mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
    filter(year >= year_min, year <= year_max, is.finite(residual)) %>%
    inner_join(unit_map, by = "stat_rec") %>%
    mutate(period = assign_period(year, break_year), cell_id = paste(unit_id, period, sep = "__"))
}

#' Required n (replicates per cell) for a target Spearman-Brown reliability
#' of the CELL MEAN, given the haul-level intraclass correlation (icc) —
#' i.e. how many hauls are needed so that sampling noise in the cell mean is
#' small relative to genuine between-cell signal. n = R(1-icc) / (icc(1-R)).
#' icc <= 0 (no detectable between-cell signal at all) -> Inf, not a typo:
#' no amount of additional replication fixes a zero-signal cell structure.
required_n_for_reliability <- function(icc, target_reliability) {
  vapply(target_reliability, function(target) {
    if (!is.finite(icc) || icc <= 0) {
      return(Inf)
    }
    n <- target * (1 - icc) / (icc * (1 - target))
    max(1, ceiling(n))
  }, numeric(1))
}

#' Full sample-size-adequacy row for one scheme/parameter combination:
#' variance decomposition of haul-level residual by cell (unit x period),
#' plus the reliability-derived required-n for each target in
#' reliability_targets. This is the "statistically defensible" alternative
#' to an arbitrary fixed MIN_HAULS_PER_CELL.
sample_size_adequacy_row <- function(haul_cells, scheme_label, param_label, reliability_targets) {
  vc <- variance_components_anova(haul_cells$residual, haul_cells$cell_id)
  required_n <- required_n_for_reliability(vc$icc, reliability_targets)

  out <- tibble(
    scheme = scheme_label,
    parameter = param_label,
    n_cells = vc$n_groups,
    n_hauls = vc$n_obs,
    min_hauls_per_cell_observed = vc$min_n_per_group,
    max_hauls_per_cell_observed = vc$max_n_per_group,
    var_between_cell = vc$var_between,
    var_within_cell = vc$var_within,
    icc = vc$icc
  )
  for (i in seq_along(reliability_targets)) {
    out[[sprintf("required_n_reliability_%.1f", reliability_targets[i])]] <- required_n[i]
  }
  out
}
