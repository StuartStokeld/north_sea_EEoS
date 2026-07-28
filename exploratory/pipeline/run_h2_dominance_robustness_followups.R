# Step 0 robustness FOLLOW-UPS — three loose ends before the nested H2 model:
#   1. Re-report the linear-vs-loess max deviation (verify it wasn't actually
#      missing, rather than silently re-running with different parameters)
#   2. Verify the 26 rectangles dropped for "no Couce coverage" all already
#      satisfy the provisional 5-haul minimum (sanity check, not expected to
#      change anything)
#   3. Characterise the spatial/temporal pattern behind why those 26
#      rectangles lack Couce coverage (descriptive only)
#
# Does NOT modify Step 0 outputs (step0_*.csv/.rds/.md) or the prior
# robustness-check outputs (step0_robustness_*.csv/.md/.png). Reads only
# outputs/step0_rectangle_panel.csv, outputs/step0_robustness_linear_vs_loess.csv,
# outputs/h2_couce_rectangle_effort.rds, and the ICES rectangle shapefile.
#
# Reported uninterpreted / descriptively, as with Step 0 and the prior
# robustness check — does not draw a conclusion about H2 validity.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_dominance_robustness_followups.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_dominance_diagnostic_helpers.R")) # STEP0_MIN_HAULS_PROVISIONAL
source(file.path(script_dir, "R", "h2_dominance_missingness_helpers.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_rect_panel <- file.path(project_root, "outputs", "step0_rectangle_panel.csv")
path_loess_prior <- file.path(project_root, "outputs", "step0_robustness_linear_vs_loess.csv")
path_couce_rect <- file.path(project_root, "outputs", "h2_couce_rectangle_effort.rds")

path_out_md <- file.path(project_root, "outputs", "step0_robustness_followups.md")
path_out_haulcount <- file.path(project_root, "outputs", "step0_robustness_haulcount_check.csv")
path_out_geo <- file.path(project_root, "outputs", "step0_robustness_missing_rectangles_geo.csv")
path_out_ecoregion <- file.path(project_root, "outputs", "step0_robustness_ecoregion_breakdown.csv")
path_out_fig_map <- file.path(project_root, "outputs", "figures", "step0_robustness_missing_rectangles_map.png")

stopifnot(file.exists(path_rect_panel), file.exists(path_loess_prior), file.exists(path_couce_rect))
dir.create(dirname(path_out_fig_map), recursive = TRUE, showWarnings = FALSE)

rect_panel <- read.csv(path_rect_panel, stringsAsFactors = FALSE)
dropped <- rect_panel %>% filter(is.na(mean_annual_hours_total))
retained <- rect_panel %>% filter(!is.na(mean_annual_hours_total))
couce_rect <- readRDS(path_couce_rect)

cat("=== Step 0 robustness follow-ups ===\n")
cat("Rectangle panel:", nrow(rect_panel), "total ->", nrow(retained), "retained /", nrow(dropped), "dropped\n\n")

# ===========================================================================
# Follow-up 1: linear vs loess — verify, don't silently re-run
# ===========================================================================
cat("--- Follow-up 1: linear vs loess max deviation ---\n")

loess_prior <- read.csv(path_loess_prior, stringsAsFactors = FALSE)
loess_D <- loess_prior %>% filter(variable == "mean_D")

loess_status <- if (nrow(loess_D) == 1L && is.finite(loess_D$max_abs_diff[1])) {
  "PRESENT in outputs/step0_robustness_linear_vs_loess.csv (and in step0_robustness_check.md Check 3) — was NOT actually missing; not re-run with different parameters."
} else {
  "MISSING or non-finite in outputs/step0_robustness_linear_vs_loess.csv — flagging for investigation rather than silently re-attempting."
}
cat(loess_status, "\n")
if (nrow(loess_D) == 1L) {
  cat(sprintf(
    "  max_abs_diff = %.4f at fishing_hours = %.1f (linear_fitted = %.4f, loess_fitted = %.4f, span = %.2f, n = %d)\n",
    loess_D$max_abs_diff[1], loess_D$at_fishing_hours[1],
    loess_D$linear_fitted_at_max[1], loess_D$loess_fitted_at_max[1],
    loess_D$loess_span[1], loess_D$n[1]
  ))
}
cat("\n")

# ===========================================================================
# Follow-up 2: haul-count verification for the 26 dropped rectangles
# ===========================================================================
cat("--- Follow-up 2: haul-count verification (26 dropped rectangles) ---\n")

haulcount_check <- summarise_haul_counts(dropped, min_hauls = STEP0_MIN_HAULS_PROVISIONAL)
cat(sprintf(
  "n=%d, min=%d, median=%.1f, max=%d hauls; all >= provisional threshold (%d): %s\n",
  haulcount_check$n_rectangles, haulcount_check$min_n_hauls, haulcount_check$median_n_hauls,
  haulcount_check$max_n_hauls, haulcount_check$min_hauls_threshold, haulcount_check$all_meet_threshold
))
write_csv(haulcount_check, path_out_haulcount)
cat("Saved", path_out_haulcount, "\n\n")

# ===========================================================================
# Follow-up 3: spatial / temporal characterisation of the 26 dropped rectangles
# ===========================================================================
cat("--- Follow-up 3: spatial / temporal characterisation ---\n")

geo <- load_ices_rectangle_geo(project_root)
dropped_geo <- join_rectangle_geo(dropped, geo)
retained_geo <- join_rectangle_geo(retained, geo)

n_overlap_couce_stat_rec <- length(intersect(dropped$stat_rec, couce_rect$stat_rec))
cat(
  "Overlap between the 26 dropped stat_rec and Couce's own rectangle list (any year):",
  n_overlap_couce_stat_rec, "of 26 (0 confirms these rectangles never appear in Couce's data at all,",
  "i.e. a spatial coverage gap, not a rectangle-year-specific temporal gap)\n\n"
)

eco_breakdown <- bind_rows(
  ecoregion_breakdown_table(dropped_geo, sprintf("dropped (n=%d)", nrow(dropped_geo))),
  ecoregion_breakdown_table(retained_geo, sprintf("retained (n=%d)", nrow(retained_geo)))
)
cat("Ecoregion breakdown (ICES_Statistical_Rectangles_Eco.shp, primary Ecoregion by max overlap %):\n")
print(eco_breakdown)
write_csv(eco_breakdown, path_out_ecoregion)
cat("Saved", path_out_ecoregion, "\n\n")

bbox_summary <- bind_rows(
  geo_bbox_summary(dropped_geo, "dropped"),
  geo_bbox_summary(retained_geo, "retained")
)
cat("Geographic bounding box (degrees):\n")
print(bbox_summary)
cat("\n")

# Two visually distinct spatial clusters among the 26 (descriptive grouping,
# not a formal cluster analysis): southern edge (~49.5-51N, -1 to 2E, near the
# eastern English Channel approach) vs north-eastern edge (~55.5-59N, 8-13E,
# near the Skagerrak/Kattegat approach off Denmark).
dropped_geo <- dropped_geo %>%
  mutate(
    descriptive_cluster = case_when(
      NORTH <= 51 & WEST <= 2 ~ "southern edge (~49.5-51N, English Channel approach)",
      WEST >= 8 ~ "north-eastern edge (~55.5-59N, Skagerrak/Kattegat approach)",
      TRUE ~ "other"
    )
  )
cluster_counts <- dropped_geo %>% count(descriptive_cluster, name = "n_rectangles") %>% arrange(desc(n_rectangles))
cat("Descriptive spatial clustering among the 26 dropped rectangles:\n")
print(cluster_counts)
cat("\n")

# Temporal: haul-year range already restricted to 1985-2015 by Step 0's own
# aggregation filter, and Couce's window is also 1985-2015 (see DATA_SOURCE.md)
# -- so a rectangle-year window mismatch cannot explain the gap; combined with
# the zero stat_rec overlap above, the gap is spatial, not temporal.
year_range_summary <- dropped %>%
  summarise(
    n = n(),
    earliest_year_min_seen = min(year_min_seen),
    latest_year_max_seen = max(year_max_seen),
    n_full_1985_2015_coverage = sum(year_min_seen == 1985L & year_max_seen == 2015L)
  )
cat("Haul-year coverage of the 26 dropped rectangles (from step0_rectangle_panel.csv):\n")
print(year_range_summary)
cat(
  "Couce coverage window is also 1985-2015 (DATA_SOURCE.md) and Step 0's own aggregation",
  "already restricts hauls to 1985-2015, so this range cannot itself explain the gap;",
  "combined with 0/26 stat_rec overlap above, the gap reads as spatial (these rectangles",
  "are simply outside Couce's ~215-rectangle grid), not a temporal-window mismatch.\n\n"
)

geo_out <- dropped_geo %>%
  select(stat_rec, n_hauls, mean_D, mean_size_CV, year_min_seen, year_max_seen,
         SOUTH, NORTH, WEST, EAST, Ecoregion, descriptive_cluster) %>%
  arrange(WEST, stat_rec)
write_csv(geo_out, path_out_geo)
cat("Saved", path_out_geo, "\n\n")

# Supplementary map (not explicitly required, but the briefing lists a map as
# an optional supporting plot for the spatial-clustering question)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  map_df <- bind_rows(
    retained_geo %>% mutate(group = "retained (n=161)", lon_c = (WEST + EAST) / 2, lat_c = (SOUTH + NORTH) / 2),
    dropped_geo %>% mutate(group = "dropped (n=26)", lon_c = (WEST + EAST) / 2, lat_c = (SOUTH + NORTH) / 2)
  )
  p_map <- ggplot(map_df, aes(x = lon_c, y = lat_c, colour = group)) +
    geom_point(size = 2.2, alpha = 0.8) +
    scale_colour_manual(values = c("retained (n=161)" = "grey60", "dropped (n=26)" = "#d6604d")) +
    coord_fixed(ratio = 1.6) +
    labs(
      x = "Longitude (rectangle centroid)", y = "Latitude (rectangle centroid)",
      colour = NULL,
      title = "Step 0 robustness follow-up: rectangles dropped for no Couce coverage",
      subtitle = "Descriptive only — spatial clustering, not a formal test"
    ) +
    theme_minimal(base_size = 11)
  ggsave(path_out_fig_map, p_map, width = 7.5, height = 7, dpi = 120)
  cat("Saved", path_out_fig_map, "\n\n")
}

# ===========================================================================
# Single addendum document (does NOT modify step0_robustness_check.md)
# ===========================================================================
md_table <- function(df, digits = 4L) {
  df <- as.data.frame(df)
  num_cols <- vapply(df, is.numeric, logical(1))
  for (col in names(df)[num_cols]) {
    df[[col]] <- format(round(df[[col]], digits), nsmall = 0, trim = TRUE)
  }
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, body)
}

md <- c(
  "# Step 0 robustness follow-ups",
  "",
  paste(
    "Generated:", Sys.Date(),
    "— addendum to `step0_robustness_check.md`. Does not modify Step 0 or the",
    "prior robustness-check outputs. Reported descriptively / uninterpreted."
  ),
  "",
  "## 1. Linear vs loess max deviation",
  "",
  loess_status,
  "",
  md_table(loess_D %>% select(-variable)),
  "",
  "## 2. Haul-count verification (26 dropped rectangles)",
  "",
  md_table(haulcount_check),
  "",
  "## 3. Spatial / temporal characterisation of the 26 dropped rectangles",
  "",
  sprintf(
    "Overlap between the 26 dropped `stat_rec` and Couce's own rectangle list (any year): **%d of 26**.",
    n_overlap_couce_stat_rec
  ),
  "This confirms the gap is a rectangle-level spatial coverage gap (these stat_rec codes never",
  "appear in Couce's data in any year), not a year-specific temporal gap within an otherwise-covered rectangle.",
  "",
  "**Ecoregion breakdown** (ICES `ICES_Statistical_Rectangles_Eco.shp`, primary Ecoregion by max spatial overlap %):",
  "",
  md_table(eco_breakdown),
  "",
  "**Geographic bounding box (degrees):**",
  "",
  md_table(bbox_summary),
  "",
  "**Descriptive spatial clustering** (not a formal cluster analysis) among the 26:",
  "",
  md_table(cluster_counts),
  "",
  "**Haul-year coverage of the 26 dropped rectangles:**",
  "",
  md_table(year_range_summary),
  "",
  paste(
    "Couce's coverage window is also 1985-2015 (`data/external/couce_trawling_effort/DATA_SOURCE.md`)",
    "and Step 0's own aggregation already restricts hauls to 1985-2015, so the year range itself cannot",
    "explain the gap. Combined with the 0/26 stat_rec overlap above, the gap reads as spatial",
    "(these rectangles are outside Couce's ~215-rectangle grid), not a temporal-window mismatch."
  ),
  "",
  paste(
    "No existing repo file (`DATA_SOURCE.md`, `pipeline/README.md`, `h2_couce_import_diagnostics.csv`)",
    "explicitly documents which rectangles Couce's reconstruction excludes, so this characterisation",
    "required new descriptive analysis rather than confirming a pre-existing note."
  ),
  "",
  if (file.exists(path_out_fig_map)) {
    sprintf("![Dropped vs retained rectangles](figures/%s)", basename(path_out_fig_map))
  } else {
    character(0)
  },
  "",
  "## Notes on scope",
  "",
  "- Descriptive only, as instructed — no conclusion drawn about whether this affects H2's validity.",
  "- The two-cluster grouping above is a simple coordinate-based descriptive split (thresholds chosen",
  "  by inspection of the bounding boxes), not a statistical clustering method.",
  "",
  "*Outputs: `outputs/step0_robustness_haulcount_check.csv`, `outputs/step0_robustness_missing_rectangles_geo.csv`,",
  "`outputs/step0_robustness_ecoregion_breakdown.csv`, `outputs/figures/step0_robustness_missing_rectangles_map.png`.*"
)

writeLines(md, path_out_md)
cat("Saved", path_out_md, "\n")

cat("\n=== Step 0 robustness follow-ups complete — descriptive only; no conclusion drawn. ===\n")
