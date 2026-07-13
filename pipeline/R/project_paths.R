# Resolve repo workspace root for knitr reports outside pipeline/

#' Find north_sea_eeos workspace root (contains FishGlob_data/).
find_workspace_root <- function(start_dir = getwd()) {
  if (exists("get_project_root_from", mode = "function")) {
    return(get_project_root_from(start_dir))
  }
  dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(dir, "FishGlob_data"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      break
    }
    dir <- parent
  }
  stop("Could not find workspace root (FishGlob_data/).")
}

#' Load H1 helpers from pipeline/R/ given workspace root.
source_pipeline_helpers <- function(workspace_root) {
  r_dir <- file.path(workspace_root, "pipeline", "R")
  source(file.path(r_dir, "h1_common.R"))
  source(file.path(r_dir, "h1_lne_helpers.R"))
  source(file.path(r_dir, "h1_null_helpers.R"))
  invisible(r_dir)
}
