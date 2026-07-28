# Helpers for the "Missing Visuals" task: standalone haul map + composition
# figures. All per-haul quantities (S, D, size_CV, N, B_obs, B_pred, etc.) and
# the rectangle-level total_hauls / has_couce_coverage flags are read from
# existing outputs, never recomputed here.

#' Join the dominance table (D, size_CV, ...) with species richness (S) from
#' the main haul-level prediction table, by haul_id. Errors if the join is
#' not a lossless 1:1 match (no recomputation fallback).
build_composition_table <- function(dominance_tbl, haul_eeos_tbl) {
  s_lookup <- haul_eeos_tbl[, c("haul_id", "S")]
  out <- merge(dominance_tbl, s_lookup, by = "haul_id", all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(dominance_tbl)) {
    stop("Composition join changed row count: ", nrow(dominance_tbl), " -> ", nrow(out))
  }
  if (anyNA(out$S)) {
    stop("Composition join produced ", sum(is.na(out$S)), " unmatched haul_id (missing S).")
  }
  out
}

#' Read the pre-existing top-decile boundary for D from the H1 dominance
#' bin table (`h1_dominance_by_D_bins.csv`, produced by
#' `explore_h1_haul_dominance.R`), rather than recomputing a fresh quantile.
#' Returns list(threshold, n_at_or_above, source_n_bin).
read_D_top_decile_reference <- function(d_bins_path) {
  stopifnot(file.exists(d_bins_path))
  bins <- read.csv(d_bins_path, stringsAsFactors = FALSE)
  top <- bins[bins$metric == "D" & bins$n_bins == 10L & bins$bin == 10L, ]
  if (nrow(top) != 1L) {
    stop("Could not locate the D top-decile (n_bins=10, bin=10) row in ", d_bins_path)
  }
  list(threshold = top$metric_min[1], n_in_decile = top$n[1], source_n_bins = 10L)
}
