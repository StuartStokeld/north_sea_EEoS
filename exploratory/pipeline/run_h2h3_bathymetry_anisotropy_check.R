# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Bathymetry-informed spatial anisotropy diagnostic on primary H2/H3 residuals.
#
# Design: display_discussion/Design_bathymetry_spatial_anisotropy.md
# Run:    Rscript --vanilla pipeline/run_h2h3_bathymetry_anisotropy_check.R
#
# Dependencies: system Python3 + GDAL (osgeo) for zonal extraction;
#               circular + gstat (+ sp) in .R_libs or site library.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_bathymetry_anisotropy_check.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_bathymetry_anisotropy_helpers.R"))

ensure_bathymetry_libs(project_root)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
paths <- h2_output_paths(project_root)
gebco <- gebco_paths(project_root)
shp <- h2_ices_shapefile_path(project_root)

path_resids <- file.path(project_root, "outputs", "residuals_by_rectangle.csv")
path_blups <- file.path(project_root, "outputs", "blups_by_rectangle.csv")
path_panel <- paths$panel

path_bathy_csv <- file.path(project_root, "outputs", "bathymetry_by_rectangle.csv")
path_vgram_resid_png <- file.path(project_root, "outputs", "directional_variogram_resid.png")
path_vgram_blup_png <- file.path(project_root, "outputs", "directional_variogram_blup.png")
path_alignment <- file.path(project_root, "outputs", "bearing_alignment_test.csv")
path_verdict <- file.path(project_root, "outputs", "bathymetry_anisotropy_verdict.md")
path_run_log <- file.path(project_root, "outputs", "h2h3_bathymetry_anisotropy_run_log.md")
path_session <- file.path(project_root, "outputs", "h2h3_bathymetry_anisotropy_sessionInfo.txt")
path_tmp_ids <- file.path(project_root, "outputs", "_tmp_panel_stat_rec.txt")
py_script <- file.path(script_dir, "python", "extract_rectangle_bathymetry.py")  # exploratory/pipeline/python/

stopifnot(
  file.exists(path_resids),
  file.exists(path_blups),
  file.exists(path_panel),
  file.exists(shp),
  file.exists(gebco$bathy),
  file.exists(gebco$tid),
  file.exists(py_script)
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2/H3 bathymetry spatial anisotropy check — run log")
logmsg("")
logmsg(
  "Diagnostic test of whether primary-model residual spatial correlation is ",
  "anisotropic and aligned with local shelf (depth-gradient) geometry. ",
  "Design: `display_discussion/Design_bathymetry_spatial_anisotropy.md`."
)

# ---------------------------------------------------------------------------
# Session + locked conventions (BEFORE any confirmatory p-value)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session")
sink(path_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_session)

logmsg("")
logmsg("## Locked conventions (documented before unblinding confirmatory p)")
logmsg(
  "1. Bearings use gstat plane convention: degrees counterclockwise from +east ",
  "(0° = east, 90° = north)."
)
logmsg(
  "2. `bearing_depth` = along-shelf axis = local depth-gradient (cross-shelf) ",
  "bearing + 90°. Depth = −GEBCO elevation (positive deeper)."
)
logmsg(
  "3. `bearing_resid` / `bearing_blup` = direction bin (0/45/90/135 ±22.5°) with ",
  "lowest local semivariance among pairs in lag window [0.5°, 4°] involving that ",
  "rectangle (strongest residual similarity = along-shelf residual structure)."
)
logmsg(
  "4. Confirmatory test: fold both bearings to axial [0°, 180), convert to radians, ",
  "double the angle (axial→circular map), then ",
  "`circular::cor.circular(..., test = TRUE)` (Jammalamadaka–Sarma)."
)
logmsg(
  "5. TID: direct soundings = codes {10–17}; land TID=0 excluded from denominator. ",
  "Flag if pct_sounding < 50% (sensitivity threshold 25%). Flagged rectangles ",
  "excluded from the circular alignment test only."
)
logmsg("6. Primary object: signed rectangle-mean residual; BLUP = robustness only.")
logmsg("7. Alignment confirmed if primary p < 0.05.")

TID_FLAG_PCT <- 50
TID_SENS_PCT <- 25
ALPHA <- 0.05

# ---------------------------------------------------------------------------
# Load rectangle panel + residuals / BLUPs
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
panel <- readRDS(path_panel)
panel_ids <- sort(unique(as.character(panel$stat_rec)))
logmsg("Panel rectangles: ", length(panel_ids))

resids <- utils::read.csv(path_resids, stringsAsFactors = FALSE)
blups <- utils::read.csv(path_blups, stringsAsFactors = FALSE)
resids$stat_rec <- as.character(resids$stat_rec)
blups$stat_rec <- as.character(blups$stat_rec)
stopifnot(setequal(resids$stat_rec, panel_ids), setequal(blups$stat_rec, panel_ids))
logmsg("Loaded residuals: ", path_resids)
logmsg("Loaded BLUPs: ", path_blups)

# ---------------------------------------------------------------------------
# Step 1 — zonal bathymetry + TID (Python / GDAL)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1 — Zonal bathymetry and TID extraction")
writeLines(panel_ids, path_tmp_ids)
py <- Sys.which("python3")
if (!nzchar(py)) stop("python3 not found on PATH")
cmd <- paste(
  shQuote(py),
  shQuote(py_script),
  "--shapefile", shQuote(shp),
  "--panel-ids", shQuote(path_tmp_ids),
  "--bathy", shQuote(gebco$bathy),
  "--tid", shQuote(gebco$tid),
  "--out", shQuote(path_bathy_csv)
)
logmsg("Running: ", cmd)
status <- system(cmd)
if (!identical(as.integer(status), 0L)) {
  stop("Python bathymetry extraction failed with status ", status)
}
bathy <- utils::read.csv(path_bathy_csv, stringsAsFactors = FALSE)
bathy$stat_rec <- as.character(bathy$stat_rec)
stopifnot(nrow(bathy) == length(panel_ids), setequal(bathy$stat_rec, panel_ids))

bathy$tid_flag_lt50 <- is.finite(bathy$pct_sounding) & bathy$pct_sounding < TID_FLAG_PCT
bathy$tid_flag_lt25 <- is.finite(bathy$pct_sounding) & bathy$pct_sounding < TID_SENS_PCT
utils::write.csv(bathy, path_bathy_csv, row.names = FALSE)
logmsg("Wrote: ", path_bathy_csv)
logmsg(sprintf(
  "TID flags: n_lt50=%d (excluded from alignment); n_lt25=%d (sensitivity)",
  sum(bathy$tid_flag_lt50), sum(bathy$tid_flag_lt25)
))
logmsg(sprintf(
  "Depth summary: mean_depth median=%.1f m; depth_range median=%.1f m; grad magnitude median=%.5g",
  stats::median(bathy$mean_depth_m, na.rm = TRUE),
  stats::median(bathy$depth_range_m, na.rm = TRUE),
  stats::median(bathy$grad_magnitude_m_per_m, na.rm = TRUE)
))
logmsg(sprintf(
  "GEBCO edge distance: min=%.3f deg, median=%.3f deg",
  min(bathy$dist_to_gebco_edge_deg, na.rm = TRUE),
  stats::median(bathy$dist_to_gebco_edge_deg, na.rm = TRUE)
))

# Align join order to bathymetry table
ord <- match(bathy$stat_rec, resids$stat_rec)
z_resid <- resids$resid[ord]
z_blup <- blups$blup[match(bathy$stat_rec, blups$stat_rec)]
lon <- bathy$lon
lat <- bathy$lat

# ---------------------------------------------------------------------------
# Step 2 — directional variograms (descriptive)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 2 — Directional variograms (descriptive)")
vg_resid <- empirical_directional_variogram(z_resid, lon, lat)
vg_blup <- empirical_directional_variogram(z_blup, lon, lat)
fit_resid <- fit_directional_exponentials(vg_resid)
fit_blup <- fit_directional_exponentials(vg_blup)
logmsg("Residual directional fits:")
for (i in seq_len(nrow(fit_resid))) {
  r <- fit_resid[i, ]
  logmsg(sprintf(
    "  dir=%g°  nugget=%.4g  sill=%.4g  range=%.4g  fit_ok=%s",
    r$dir.hor, r$nugget, r$sill, r$range, r$fit_ok
  ))
}
logmsg("BLUP directional fits:")
for (i in seq_len(nrow(fit_blup))) {
  r <- fit_blup[i, ]
  logmsg(sprintf(
    "  dir=%g°  nugget=%.4g  sill=%.4g  range=%.4g  fit_ok=%s",
    r$dir.hor, r$nugget, r$sill, r$range, r$fit_ok
  ))
}
save_directional_variogram_plot(
  vg_resid, fit_resid, path_vgram_resid_png,
  "Directional variogram — rectangle-mean residual"
)
save_directional_variogram_plot(
  vg_blup, fit_blup, path_vgram_blup_png,
  "Directional variogram — BLUP"
)
logmsg("Wrote: ", path_vgram_resid_png)
logmsg("Wrote: ", path_vgram_blup_png)

# ---------------------------------------------------------------------------
# Step 3 — per-rectangle bearings + confirmatory circular correlation
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 3 — Per-rectangle bearings and Jammalamadaka–Sarma test")
loc_resid <- local_correlation_bearing(z_resid, lon, lat)
loc_blup <- local_correlation_bearing(z_blup, lon, lat)

bearing_tab <- data.frame(
  stat_rec = bathy$stat_rec,
  lon = lon,
  lat = lat,
  bearing_resid = loc_resid$bearing_corr_deg,
  n_local_pairs_resid = loc_resid$n_local_pairs,
  bearing_blup = loc_blup$bearing_corr_deg,
  n_local_pairs_blup = loc_blup$n_local_pairs,
  bearing_depth_cross = bathy$grad_bearing_cross_deg,
  bearing_depth = bathy$grad_bearing_along_deg,
  pct_sounding = bathy$pct_sounding,
  tid_flag_lt50 = bathy$tid_flag_lt50,
  tid_flag_lt25 = bathy$tid_flag_lt25,
  stringsAsFactors = FALSE
)

use_primary <- !bearing_tab$tid_flag_lt50
logmsg(sprintf(
  "Alignment sample after TID exclusion: %d / %d rectangles",
  sum(use_primary), nrow(bearing_tab)
))
logmsg(sprintf(
  "Finite bearing_resid: %d; bearing_blup: %d (within TID-ok set)",
  sum(use_primary & is.finite(bearing_tab$bearing_resid)),
  sum(use_primary & is.finite(bearing_tab$bearing_blup))
))

logmsg("")
logmsg("### Confirmatory test (primary residuals)")
test_resid <- jammalamadaka_sarma_test(
  bearing_tab$bearing_resid[use_primary],
  bearing_tab$bearing_depth[use_primary]
)
logmsg(sprintf(
  "n=%d  rho=%.4f  statistic=%.4g  p.value=%.6g  (%s)",
  test_resid$n, test_resid$rho, test_resid$statistic, test_resid$p.value, test_resid$note
))

logmsg("")
logmsg("### Robustness test (BLUPs)")
test_blup <- jammalamadaka_sarma_test(
  bearing_tab$bearing_blup[use_primary],
  bearing_tab$bearing_depth[use_primary]
)
logmsg(sprintf(
  "n=%d  rho=%.4f  statistic=%.4g  p.value=%.6g  (%s)",
  test_blup$n, test_blup$rho, test_blup$statistic, test_blup$p.value, test_blup$note
))

# Sensitivity: also at 25% TID threshold
use_sens <- !bearing_tab$tid_flag_lt25
test_resid_25 <- jammalamadaka_sarma_test(
  bearing_tab$bearing_resid[use_sens],
  bearing_tab$bearing_depth[use_sens]
)
logmsg(sprintf(
  "Sensitivity TID<25%% exclusion: n=%d  rho=%.4f  p=%.6g",
  test_resid_25$n, test_resid_25$rho, test_resid_25$p.value
))

verdict <- classify_anisotropy_verdict(test_resid$p.value, test_blup$p.value, ALPHA)
logmsg("")
logmsg("## Verdict")
logmsg("Primary p < 0.05? ", is.finite(test_resid$p.value) && test_resid$p.value < ALPHA)
logmsg("Robustness p < 0.05? ", is.finite(test_blup$p.value) && test_blup$p.value < ALPHA)
logmsg("Verdict: ", verdict$verdict)
logmsg("Next step: ", verdict$next_step)

# ---------------------------------------------------------------------------
# Write deliverables
# ---------------------------------------------------------------------------
align_summary <- data.frame(
  test = c("primary_resid", "robustness_blup", "sensitivity_resid_tid25"),
  n = c(test_resid$n, test_blup$n, test_resid_25$n),
  rho_js = c(test_resid$rho, test_blup$rho, test_resid_25$rho),
  statistic = c(test_resid$statistic, test_blup$statistic, test_resid_25$statistic),
  p_value = c(test_resid$p.value, test_blup$p.value, test_resid_25$p.value),
  alpha = ALPHA,
  significant = c(
    is.finite(test_resid$p.value) && test_resid$p.value < ALPHA,
    is.finite(test_blup$p.value) && test_blup$p.value < ALPHA,
    is.finite(test_resid_25$p.value) && test_resid_25$p.value < ALPHA
  ),
  verdict = c(verdict$verdict, NA_character_, NA_character_),
  stringsAsFactors = FALSE
)
# Per-rectangle bearings appended below summary via two-file approach: write
# alignment test summary, and also a detailed bearing table alongside.
utils::write.csv(align_summary, path_alignment, row.names = FALSE)
path_bearings <- file.path(project_root, "outputs", "bearings_by_rectangle.csv")
utils::write.csv(bearing_tab, path_bearings, row.names = FALSE)
logmsg("Wrote: ", path_alignment)
logmsg("Wrote: ", path_bearings)

# Fit tables for provenance
utils::write.csv(
  fit_resid,
  file.path(project_root, "outputs", "directional_variogram_resid_fits.csv"),
  row.names = FALSE
)
utils::write.csv(
  fit_blup,
  file.path(project_root, "outputs", "directional_variogram_blup_fits.csv"),
  row.names = FALSE
)

verdict_md <- c(
  "# Bathymetry spatial anisotropy — verdict",
  "",
  paste0("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste0("**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`"),
  paste0("**Run log:** `outputs/h2h3_bathymetry_anisotropy_run_log.md`"),
  "",
  "## Confirmatory test (Jammalamadaka–Sarma)",
  "",
  sprintf(
    "| Test | n | ρ_JS | p | Significant (α=%.2f) |",
    ALPHA
  ),
  "|---|---:|---:|---:|---|",
  sprintf(
    "| Primary (signed resid vs along-shelf depth) | %d | %.4f | %.4g | %s |",
    test_resid$n, test_resid$rho, test_resid$p.value,
    is.finite(test_resid$p.value) && test_resid$p.value < ALPHA
  ),
  sprintf(
    "| Robustness (BLUP vs along-shelf depth) | %d | %.4f | %.4g | %s |",
    test_blup$n, test_blup$rho, test_blup$p.value,
    is.finite(test_blup$p.value) && test_blup$p.value < ALPHA
  ),
  sprintf(
    "| Sensitivity (resid, TID flag at 25%%) | %d | %.4f | %.4g | %s |",
    test_resid_25$n, test_resid_25$rho, test_resid_25$p.value,
    is.finite(test_resid_25$p.value) && test_resid_25$p.value < ALPHA
  ),
  "",
  paste0("## Verdict: **", verdict$verdict, "**"),
  "",
  paste0("**Next step:** ", verdict$next_step),
  "",
  "## Notes",
  "",
  paste0(
    "- TID exclusion at 50%: ", sum(bathy$tid_flag_lt50),
    " rectangles flagged; retained in tables/maps but excluded from alignment test."
  ),
  paste0(
    "- Directional variograms (descriptive): `outputs/directional_variogram_resid.png`, ",
    "`outputs/directional_variogram_blup.png`."
  ),
  paste0(
    "- Per-rectangle bathymetry: `outputs/bathymetry_by_rectangle.csv`; ",
    "bearings: `outputs/bearings_by_rectangle.csv`."
  )
)
writeLines(verdict_md, path_verdict)
logmsg("Wrote: ", path_verdict)

writeLines(run_log, path_run_log)
logmsg("Wrote: ", path_run_log)

# Cleanup temp ids
if (file.exists(path_tmp_ids)) file.remove(path_tmp_ids)

invisible(NULL)
