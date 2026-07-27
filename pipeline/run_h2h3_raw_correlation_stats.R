# Supplementary raw correlation & range statistics for H2/H3.
#
# Descriptive only — no model fitting. Reuses decomposed haul-level data from
# h2h3_wb_model_objects.rds (same panel / phase / FP_between / FP_within as the
# within-between models). Compares raw Pearson correlations to already-reported
# adjusted phase slopes for sign-disagreement flags only.
#
# Run: Rscript --vanilla pipeline/run_h2h3_raw_correlation_stats.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_raw_correlation_stats.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

PHASE_LEVELS <- c("1985-1988", "1989-2000", "2001-2007", "2008-2015")

path_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_slopes <- file.path(project_root, "outputs", "h2h3_wb_fp_slopes_by_phase.csv")

path_out_h2 <- file.path(project_root, "outputs", "h2h3_raw_correlation_H2.csv")
path_out_h3 <- file.path(project_root, "outputs", "h2h3_raw_correlation_H3.csv")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_raw_correlation_run_log.md")

stopifnot(file.exists(path_models), file.exists(path_slopes))

FRAMING_NOTE <- paste0(
  "These correlations are raw and unadjusted: they do not control for the other ",
  "decomposed fishing-pressure term, CAR spatial structure, or partial pooling the ",
  "way the fitted within-between models do. They are expected to be directionally ",
  "consistent with, but not numerically identical to, the adjusted phase slopes ",
  "already reported."
)

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

sign_disagreement <- function(raw_r, model_slope) {
  if (!is.finite(raw_r) || !is.finite(model_slope)) return(NA)
  (raw_r > 0 && model_slope < 0) || (raw_r < 0 && model_slope > 0)
}

extreme_row <- function(df, which = c("min", "max")) {
  which <- match.arg(which)
  idx <- if (which == "min") which.min(df$mean_residual) else which.max(df$mean_residual)
  df[idx, , drop = FALSE]
}

compute_h2_phase <- function(dat_phase, phase_label, model_slope) {
  rect <- dat_phase %>%
    group_by(stat_rec) %>%
    summarise(
      mean_residual = mean(residual, na.rm = TRUE),
      FP_between = FP_between[1L],
      .groups = "drop"
    )
  ct <- stats::cor.test(rect$FP_between, rect$mean_residual, method = "pearson")
  r <- unname(ct$estimate)
  r2 <- r * r
  lo <- extreme_row(rect, "min")
  hi <- extreme_row(rect, "max")
  tibble(
    phase = phase_label,
    n_rectangles = nrow(rect),
    pearson_r = r,
    r_squared = r2,
    cor_p_value = ct$p.value,
    model_fp_between_slope = model_slope,
    model_id_for_slope = "wb_car",
    sign_disagreement = sign_disagreement(r, model_slope),
    min_mean_residual = lo$mean_residual,
    min_rectangle = as.character(lo$stat_rec),
    min_FP_between = lo$FP_between,
    max_mean_residual = hi$mean_residual,
    max_rectangle = as.character(hi$stat_rec),
    max_FP_between = hi$FP_between
  )
}

compute_h3_phase <- function(dat_phase, phase_label, model_slope) {
  rect_year <- dat_phase %>%
    group_by(stat_rec, year) %>%
    summarise(
      mean_residual = mean(residual, na.rm = TRUE),
      FP_within = FP_within[1L],
      .groups = "drop"
    )
  ct <- stats::cor.test(rect_year$FP_within, rect_year$mean_residual, method = "pearson")
  r <- unname(ct$estimate)
  r2 <- r * r
  lo <- extreme_row(rect_year, "min")
  hi <- extreme_row(rect_year, "max")
  tibble(
    phase = phase_label,
    n_rectangle_years = nrow(rect_year),
    pearson_r = r,
    r_squared = r2,
    cor_p_value = ct$p.value,
    model_fp_within_slope = model_slope,
    model_id_for_slope = "wb_primary",
    sign_disagreement = sign_disagreement(r, model_slope),
    min_mean_residual = lo$mean_residual,
    min_rectangle = as.character(lo$stat_rec),
    min_year = lo$year,
    min_FP_within = lo$FP_within,
    max_mean_residual = hi$mean_residual,
    max_rectangle = as.character(hi$stat_rec),
    max_year = hi$year,
    max_FP_within = hi$FP_within
  )
}

# ---------------------------------------------------------------------------
# Load reused data (no new prep)
# ---------------------------------------------------------------------------
mod <- readRDS(path_models)
dat <- mod$data
stopifnot(all(c("residual", "FP_between", "FP_within", "phase", "stat_rec", "year") %in% names(dat)))

slopes <- read_csv(path_slopes, show_col_types = FALSE)
h2_model_slopes <- slopes %>%
  filter(model_id == "wb_car", component == "FP_between") %>%
  select(phase, model_slope = fp_slope)
h3_model_slopes <- slopes %>%
  filter(model_id == "wb_primary", component == "FP_within") %>%
  select(phase, model_slope = fp_slope)

if (nrow(h2_model_slopes) != 4L || nrow(h3_model_slopes) != 4L) {
  stop("Expected 4 CAR FP_between and 4 primary FP_within phase slopes.")
}

logmsg("# H2/H3 supplementary raw correlation & range statistics — run log")
logmsg("")
logmsg("## Framing (raw vs adjusted)")
logmsg(FRAMING_NOTE)

logmsg("")
logmsg("## Data sources (reused; no new data prep)")
logmsg("- Haul-level panel: ", path_models, " (", nrow(dat), " hauls, ", n_distinct(dat$stat_rec), " rectangles)")
logmsg("- Model slopes for sign comparison: ", path_slopes)
logmsg("  - H2 reference slopes: wb_car FP_between (same as primary H2 reporting)")
logmsg("  - H3 reference slopes: wb_primary FP_within (same as primary H3 reporting)")

# ---------------------------------------------------------------------------
# H2: rectangle-level, per phase
# ---------------------------------------------------------------------------
h2_out <- bind_rows(lapply(PHASE_LEVELS, function(ph) {
  ms <- h2_model_slopes$model_slope[h2_model_slopes$phase == ph]
  compute_h2_phase(dat %>% filter(phase == ph), ph, ms)
})) %>%
  mutate(phase = factor(phase, levels = PHASE_LEVELS))

write_csv(h2_out, path_out_h2)

logmsg("")
logmsg("## H2 — rectangle-level correlation (mean residual vs FP_between), per phase")
logmsg("Saved: ", path_out_h2)
for (i in seq_len(nrow(h2_out))) {
  r <- h2_out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: n=%d rectangles; r=%+.4f (r²=%.4f, p=%.4g); model slope=%+.4f; ",
      "min=%+.3f (rect %s, FP_between=%.3f); max=%+.3f (rect %s, FP_between=%.3f)%s"
    ),
    r$phase, r$n_rectangles, r$pearson_r, r$r_squared, r$cor_p_value, r$model_fp_between_slope,
    r$min_mean_residual, r$min_rectangle, r$min_FP_between,
    r$max_mean_residual, r$max_rectangle, r$max_FP_between,
    if (isTRUE(r$sign_disagreement)) "  ** SIGN DISAGREEMENT **" else ""
  ))
}

# ---------------------------------------------------------------------------
# H3: rectangle-year-level, per phase
# ---------------------------------------------------------------------------
h3_out <- bind_rows(lapply(PHASE_LEVELS, function(ph) {
  ms <- h3_model_slopes$model_slope[h3_model_slopes$phase == ph]
  compute_h3_phase(dat %>% filter(phase == ph), ph, ms)
})) %>%
  mutate(phase = factor(phase, levels = PHASE_LEVELS))

write_csv(h3_out, path_out_h3)

logmsg("")
logmsg("## H3 — rectangle-year-level correlation (mean residual vs FP_within), per phase")
logmsg("Saved: ", path_out_h3)
for (i in seq_len(nrow(h3_out))) {
  r <- h3_out[i, ]
  logmsg(sprintf(
    paste0(
      "  - %s: n=%d rectangle-years; r=%+.4f (r²=%.4f, p=%.4g); model slope=%+.4f; ",
      "min=%+.3f (rect %s, year %d, FP_within=%+.3f); ",
      "max=%+.3f (rect %s, year %d, FP_within=%+.3f)%s"
    ),
    r$phase, r$n_rectangle_years, r$pearson_r, r$r_squared, r$cor_p_value, r$model_fp_within_slope,
    r$min_mean_residual, r$min_rectangle, r$min_year, r$min_FP_within,
    r$max_mean_residual, r$max_rectangle, r$max_year, r$max_FP_within,
    if (isTRUE(r$sign_disagreement)) "  ** SIGN DISAGREEMENT **" else ""
  ))
}

# ---------------------------------------------------------------------------
# Sign-disagreement summary
# ---------------------------------------------------------------------------
h2_disagree <- h2_out %>% filter(sign_disagreement)
h3_disagree <- h3_out %>% filter(sign_disagreement)

logmsg("")
logmsg("## Sign-disagreement flags (raw r vs adjusted model slope)")
if (nrow(h2_disagree) == 0L && nrow(h3_disagree) == 0L) {
  logmsg("None: raw Pearson r agrees in sign with the corresponding adjusted slope in every phase.")
} else {
  if (nrow(h2_disagree) > 0L) {
    for (i in seq_len(nrow(h2_disagree))) {
      r <- h2_disagree[i, ]
      logmsg(sprintf(
        "  - H2 %s: raw r=%+.4f vs model slope=%+.4f",
        r$phase, r$pearson_r, r$model_fp_between_slope
      ))
    }
  }
  if (nrow(h3_disagree) > 0L) {
    for (i in seq_len(nrow(h3_disagree))) {
      r <- h3_disagree[i, ]
      logmsg(sprintf(
        "  - H3 %s: raw r=%+.4f vs model slope=%+.4f",
        r$phase, r$pearson_r, r$model_fp_within_slope
      ))
    }
  }
}

logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_h2)
logmsg("- ", path_out_h3)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Raw correlation statistics complete (no figures). ===\n")
