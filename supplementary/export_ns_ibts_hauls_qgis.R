# Export NS-IBTS haul points for QGIS (2000–2020)
# - One point per haul (haul_id + year; same haul_id in different years = separate points)
# - wgt_total = sum of catch weight (wgt) across species for that haul
# Output: gis/ns_ibts_hauls_2000_2020.gpkg (EPSG:4326)
#
# Requires: dplyr, sf
#   install.packages(c("dplyr", "sf"))
#
# Working directory: north_sea_eeos (folder containing FishGlob_data/)

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

find_ns_ibts_rdata <- function() {
  candidates <- c(
    file.path("FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData"),
    file.path("..", "FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("NS-IBTS_clean.RData not found. Set working directory to north_sea_eeos.")
}

y0 <- 2000L
y1 <- 2020L
out_dir <- "gis"
out_gpkg <- file.path(out_dir, sprintf("ns_ibts_hauls_%d_%d.gpkg", y0, y1))

path_rdata <- find_ns_ibts_rdata()
load(path_rdata)

d <- data |>
  filter(.data$year >= y0, .data$year <= y1) |>
  filter(!is.na(.data$longitude), !is.na(.data$latitude))

haul_pts <- d |>
  group_by(.data$haul_id, .data$year) |>
  summarise(
    wgt_total = sum(.data$wgt, na.rm = TRUE),
    latitude = dplyr::first(.data$latitude),
    longitude = dplyr::first(.data$longitude),
    survey = dplyr::first(.data$survey),
    .groups = "drop"
  )

haul_sf <- haul_pts |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write_sf(haul_sf, out_gpkg, layer = "hauls", delete_dsn = TRUE)

cat("Wrote:", normalizePath(out_gpkg), "\n")
cat("Features:", nrow(haul_sf), "  Year range:", min(haul_sf$year), "-", max(haul_sf$year), "\n")

# QGIS (brief):
# - Layer > Add Layer > Add Vector Layer → select this .gpkg, layer "hauls"
# - Symbology: graduated or proportional symbols on "wgt_total" (size and/or color)
# - Temporal: layer properties → Temporal → enable; assign "year" as datetime or range field
#   (or filter / categorize by year without temporal controller)
