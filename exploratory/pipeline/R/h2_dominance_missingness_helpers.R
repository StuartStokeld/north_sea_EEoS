# Step 0 robustness FOLLOW-UP helpers — characterising the 26 rectangles
# dropped between the post-haul-threshold (187) and post-Couce-coverage (161)
# panels. Does NOT modify Step 0 or the prior robustness-check outputs;
# reads only outputs/step0_rectangle_panel.csv, outputs/h2_couce_rectangle_effort.rds,
# and the ICES rectangle shapefile already in gis/ICES_rectangles/.
#
# Requires h2_common.R (STEP0_MIN_HAULS_PROVISIONAL is defined in
# h2_dominance_diagnostic_helpers.R; pass explicitly here to avoid a second
# source-order dependency).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# Follow-up 2: haul-count verification for the dropped rectangles
# ---------------------------------------------------------------------------

#' Min/median/max haul count for a set of rectangles, and whether all meet
#' the provisional haul-inclusion threshold (verification only — these
#' rectangles already passed the threshold to appear in the 187-set, so this
#' is not expected to change anything).
summarise_haul_counts <- function(df, min_hauls, n_col = "n_hauls") {
  tibble(
    n_rectangles = nrow(df),
    min_n_hauls = min(df[[n_col]], na.rm = TRUE),
    median_n_hauls = median(df[[n_col]], na.rm = TRUE),
    max_n_hauls = max(df[[n_col]], na.rm = TRUE),
    min_hauls_threshold = min_hauls,
    all_meet_threshold = all(df[[n_col]] >= min_hauls, na.rm = TRUE)
  )
}

# ---------------------------------------------------------------------------
# Follow-up 3: spatial / temporal characterisation
# ---------------------------------------------------------------------------

#' Load ICES rectangle geometry (bounding box + Ecoregion) as a plain
#' data.frame, one row per ICESNAME (rectangle code), keeping the
#' highest-PERCENTAGE Ecoregion for rectangles that straddle more than one
#' (per the shapefile's own overlap-percentage field).
load_ices_rectangle_geo <- function(project_root) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package required for load_ices_rectangle_geo()")
  }
  path <- file.path(project_root, "gis", "ICES_rectangles", "ICES_Statistical_Rectangles_Eco.shp")
  stopifnot(file.exists(path))
  shp <- sf::st_read(path, quiet = TRUE)
  geo <- as.data.frame(shp)[, c("ICESNAME", "SOUTH", "WEST", "NORTH", "EAST", "Ecoregion", "PERCENTAGE")]
  geo %>%
    group_by(ICESNAME) %>%
    slice_max(PERCENTAGE, n = 1L, with_ties = FALSE) %>%
    ungroup() %>%
    rename(stat_rec = ICESNAME)
}

#' Join rectangle-level geometry/Ecoregion onto a dominance/fishing-pressure panel.
join_rectangle_geo <- function(panel, geo) {
  panel %>% left_join(geo, by = "stat_rec")
}

#' Count / percent of rectangles per Ecoregion, labelled by group (e.g.
#' "dropped (n=26)" vs "retained (n=161)") for side-by-side reporting.
ecoregion_breakdown_table <- function(df, group_label) {
  df %>%
    count(Ecoregion, name = "n_rectangles") %>%
    mutate(
      pct = round(100 * n_rectangles / sum(n_rectangles), 1),
      group = group_label
    ) %>%
    relocate(group) %>%
    arrange(desc(n_rectangles))
}

#' Geographic bounding-box summary (lat/lon min/max) for a group of rectangles.
geo_bbox_summary <- function(df, group_label) {
  tibble(
    group = group_label,
    n_rectangles = nrow(df),
    lat_min = min(df$SOUTH, na.rm = TRUE),
    lat_max = max(df$NORTH, na.rm = TRUE),
    lon_min = min(df$WEST, na.rm = TRUE),
    lon_max = max(df$EAST, na.rm = TRUE)
  )
}
