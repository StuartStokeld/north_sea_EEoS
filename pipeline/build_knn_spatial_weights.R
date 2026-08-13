# Build k-NN spatial weights (k=4) for the 158-rectangle H2 panel (Spec A).
#
# Centroids: st_centroid() on ICES shapefile subset (WGS84 geographic).
# Neighbours: great-circle distance; tie-break by (distance, stat_rec) ascending.
# Weights: row-standardised style W (1/k), left asymmetric, zero.policy=TRUE.
#
# Run from repo root:
#   Rscript pipeline/build_knn_spatial_weights.R
#
# Optional env: KNN_K=4  KNN_TIE_TOL_KM=1e-6

suppressPackageStartupMessages({
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/build_knn_spatial_weights.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_spatial_autocorr_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))

paths <- h2_output_paths(project_root)
path_panel <- paths$panel
path_out_listw <- file.path(project_root, "outputs", "knn_listw_k4.rds")
path_out_audit <- file.path(project_root, "outputs", "knn_listw_k4_audit.csv")
path_out_run_log <- file.path(project_root, "outputs", "knn_listw_k4_run_log.md")

K <- as.integer(Sys.getenv("KNN_K", unset = as.character(KNN_K_PRIMARY)))
tie_tol <- as.numeric(Sys.getenv("KNN_TIE_TOL_KM", unset = as.character(KNN_TIE_TOL_KM)))

stopifnot(file.exists(path_panel))

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# k-NN spatial weights build — run log")
logmsg("")
logmsg("Spec A: k-nearest-neighbour weights for lagged FP_between.")
logmsg(sprintf("k = %d; tie_tol_km = %.1e", K, tie_tol))

panel <- readRDS(path_panel)
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
if (nrow(panel) != 158L) {
  stop("Expected 158 rectangles in h2_rectangle_panel.rds; got ", nrow(panel))
}

logmsg("")
logmsg("## Inputs")
logmsg("Panel: ", path_panel)
logmsg("Centroid source: ICES DBF SOUTH/WEST + half cell (0.5° lat × 1° lon)")
logmsg("Row order: panel$stat_rec (", nrow(panel), " rectangles)")

pack <- build_knn_spatial_weights(panel, project_root, k = K, tie_tol_km = tie_tol)

logmsg("")
logmsg("## Weights summary")
logmsg(sprintf(
  "N = %d; n_isolated = %d; mean neighbours = %.3f; k = %d; style = W; zero.policy = TRUE",
  length(pack$listw$neighbours), pack$n_isolated,
  mean(lengths(pack$nb)), K
))
logmsg("CRS / distance: ", pack$crs)
logmsg("Symmetry: left asymmetric (k-NN directed)")

n_ties <- sum(pack$audit$tie_at_k)
logmsg(sprintf(
  "Tie-break audit: %d rectangle(s) with |dist_k - dist_{k+1}| <= tie_tol_km (%.1e)",
  n_ties, tie_tol
))
if (n_ties > 0L) {
  tie_rows <- pack$audit[pack$audit$tie_at_k, ]
  for (i in seq_len(nrow(tie_rows))) {
    r <- tie_rows[i, ]
    logmsg(sprintf(
      "  %s: dist_k=%.6f km, dist_k+1=%.6f km; neighbours=%s",
      r$stat_rec, r$dist_kth_km, r$dist_kplus1_km, r$neighbour_stat_rec
    ))
  }
}

neighbour_counts <- lengths(pack$nb)
if (!all(neighbour_counts == K)) {
  stop("Expected exactly k neighbours per rectangle; got range [",
       min(neighbour_counts), ", ", max(neighbour_counts), "]")
}

# Cross-check Couce rect_lon/lat if present
if (all(c("rect_lon", "rect_lat") %in% names(panel))) {
  d_couce <- sqrt(
    (pack$centroids$rect_lon - panel$rect_lon)^2 +
      (pack$centroids$rect_lat - panel$rect_lat)^2
  )
  logmsg(sprintf(
    "Centroid vs Couce mean coordinate offset (degrees): mean=%.4f, max=%.4f",
    mean(d_couce), max(d_couce)
  ))
}

saveRDS(
  list(
    listw = pack$listw,
    nb = pack$nb,
    n_isolated = pack$n_isolated,
    centroids = pack$centroids,
    k = pack$k,
    tie_tol_km = pack$tie_tol_km,
    crs = pack$crs,
    panel_n = nrow(panel),
    stat_rec_order = pack$centroids$stat_rec,
    built = Sys.time()
  ),
  path_out_listw
)
write_csv(pack$audit, path_out_audit)

logmsg("")
logmsg("## Outputs")
logmsg("Saved: ", path_out_listw)
logmsg("Saved: ", path_out_audit)

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== k-NN spatial weights build complete. ===\n")
