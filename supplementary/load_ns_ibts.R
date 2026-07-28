# FishGlob: NS-IBTS (North Sea International Bottom Trawl Survey)
# Source: https://github.com/fishglob/FishGlob_data/tree/main/outputs/Cleaned_data
#
# In RStudio: Session > Set Working Directory > To Source File Location,
# or run: setwd("/path/to/north_sea_eeos")

cleaned_dir <- file.path("FishGlob_data", "outputs", "Cleaned_data")
rdata <- file.path(cleaned_dir, "NS-IBTS_clean.RData")

if (!file.exists(rdata)) {
  stop("Cannot find ", rdata, " — set working directory to north_sea_eeos (folder containing FishGlob_data/).")
}

load(rdata)
# Loaded objects: `data` (species-by-haul tibble), `readme` (metadata)

cat("Rows:", nrow(data), "  Columns:", ncol(data), "\n")
print(names(data))

if (interactive()) {
  utils::View(utils::head(data, 100))
}
