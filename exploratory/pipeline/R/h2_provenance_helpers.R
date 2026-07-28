# Provenance metadata for H2 pipeline outputs
# Requires h2_common.R

#' Collect input/version metadata for H2 reproducibility.
h2_collect_provenance <- function(project_root) {
  path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
  path_couce <- h2_couce_raw_path(project_root)

  haul_mtime <- if (file.exists(path_haul)) {
    format(file.info(path_haul)$mtime, "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
  } else {
    NA_character_
  }

  couce_md5 <- if (file.exists(path_couce)) {
    as.character(tools::md5sum(path_couce))
  } else {
    NA_character_
  }

  couce_mtime <- if (file.exists(path_couce)) {
    format(file.info(path_couce)$mtime, "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
  } else {
    NA_character_
  }

  git_commit <- tryCatch(
    {
      out <- system2(
        "git",
        c("-C", shQuote(project_root), "rev-parse", "HEAD"),
        stdout = TRUE,
        stderr = FALSE
      )
      if (length(out) == 0L) NA_character_ else trimws(out[[1L]])
    },
    error = function(e) NA_character_
  )

  spdep_ver <- if (requireNamespace("spdep", quietly = TRUE)) {
    as.character(utils::packageVersion("spdep"))
  } else {
    NA_character_
  }

  spatialreg_ver <- if (requireNamespace("spatialreg", quietly = TRUE)) {
    as.character(utils::packageVersion("spatialreg"))
  } else {
    NA_character_
  }

  list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
    haul_eeos_predictions_rds = path_haul,
    haul_eeos_predictions_mtime_utc = haul_mtime,
    couce_csv = path_couce,
    couce_csv_md5 = couce_md5,
    couce_csv_mtime_utc = couce_mtime,
    git_commit = git_commit,
    spdep_version = spdep_ver,
    spatialreg_version = spatialreg_ver,
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
}

#' Attach provenance list as attribute on a result object.
h2_stamp_result <- function(result, provenance) {
  attr(result, "h2_provenance") <- provenance
  result
}
