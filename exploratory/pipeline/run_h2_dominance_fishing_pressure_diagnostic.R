# Step 0 diagnostic — dominance (D) / size-homogeneity (size_CV) vs Couce et
# al. (2020) fishing pressure. See CURSOR_BRIEFING "Step 0 Diagnostic —
# Dominance/Size-Homogeneity vs Fishing Pressure" (chat-supplied, not a repo
# file) for the full spec this script implements.
#
# LOAD-BEARING DIAGNOSTIC, NOT A ROBUSTNESS CHECK: if D/size_CV correlate with
# fishing pressure, they become a required control in H2/H3 rather than an
# optional side-analysis. This script only computes and reports the specified
# tables — it does NOT interpret the result, fold D/size_CV into the H2
# model, or decide confound-control vs mediator framing. That is Stuart's
# call after reviewing this output.
#
# Two thresholds below are PROVISIONAL (flagged in code + console + output):
#   - STEP0_MIN_HAULS_PROVISIONAL (h2_dominance_diagnostic_helpers.R): minimum
#     hauls per rectangle for the cross-sectional (H2) case. The final H2/H3
#     haul-inclusion threshold is pending supervisor review; re-run this
#     script once it is decided.
#   - STEP0_MIN_YEARS_PROVISIONAL: minimum distinct years of haul coverage
#     for a rectangle to enter the temporal (H3) case.
# No H3 policy-change dates are on file yet in this repo, so the temporal
# case reports only the annual within-rectangle correlation (no before/after
# split against policy dates — flagged explicitly in output rather than
# guessing at dates).
#
# Does NOT modify any existing H1/H2 output (haul_eeos_predictions.rds,
# h1_dominance_haul_table.csv, h2_rectangle_panel.rds, h2_couce_*.rds, etc.).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_dominance_fishing_pressure_diagnostic.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_couce_helpers.R"))
source(file.path(script_dir, "R", "h2_provenance_helpers.R"))
source(file.path(script_dir, "R", "h2_dominance_diagnostic_helpers.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_haul_dom <- file.path(project_root, "outputs", "h1_dominance_haul_table.csv")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
path_couce_rect <- file.path(project_root, "outputs", "h2_couce_rectangle_effort.rds")

path_out_summary_csv <- file.path(project_root, "outputs", "step0_dominance_fishing_pressure_diagnostic.csv")
path_out_summary_rds <- file.path(project_root, "outputs", "step0_dominance_fishing_pressure_diagnostic.rds")
path_out_rect_panel <- file.path(project_root, "outputs", "step0_rectangle_panel.csv")
path_out_temporal_by_rect <- file.path(project_root, "outputs", "step0_temporal_by_rectangle.csv")
path_out_exclusions <- file.path(project_root, "outputs", "step0_exclusion_counts.csv")
path_out_fig_D <- file.path(project_root, "outputs", "figures", "step0_D_vs_fishing_hours.png")
path_out_fig_cv <- file.path(project_root, "outputs", "figures", "step0_sizeCV_vs_fishing_hours.png")

stopifnot(
  file.exists(path_haul_dom),
  file.exists(path_couce_year),
  file.exists(path_couce_rect)
)
dir.create(dirname(path_out_fig_D), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load inputs. Haul-level D and size_CV already exist from the H1 analysis
# (outputs/h1_dominance_haul_table.csv, built by explore_h1_haul_dominance.R
# using build_haul_dominance_table() / R/h1_dominance_helpers.R) — reused
# here as-is rather than recomputed, per the briefing.
# ---------------------------------------------------------------------------
haul_dom_raw <- read.csv(path_haul_dom, stringsAsFactors = FALSE)
couce_year <- readRDS(path_couce_year)
couce_rect <- readRDS(path_couce_rect)

cat("=== Step 0 diagnostic: D / size_CV vs Couce et al. (2020) fishing pressure ===\n")
cat("PROVISIONAL thresholds in force:\n")
cat("  - min hauls per rectangle (cross-sectional / H2 case):", STEP0_MIN_HAULS_PROVISIONAL, "\n")
cat("  - min distinct years of coverage (temporal / H3 case):", STEP0_MIN_YEARS_PROVISIONAL, "\n")
cat("  - analysis period:", H2_YEAR_MIN, "-", H2_YEAR_MAX, "(pending final confirmation)\n")
cat("  - H3 policy-change dates: NONE ON FILE in this repo; before/after split not computed.\n\n")

cat("Loaded", nrow(haul_dom_raw), "hauls from h1_dominance_haul_table.csv\n")

# ===========================================================================
# Step 1: aggregate to rectangle level (cross-sectional, H2 case)
# ===========================================================================
rect_panel_all <- build_rectangle_dominance_panel(haul_dom_raw, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX)
n_rect_raw <- nrow(rect_panel_all)

hauls_thresholded <- apply_min_hauls_threshold(rect_panel_all, min_hauls = STEP0_MIN_HAULS_PROVISIONAL)
rect_panel_thresholded <- hauls_thresholded$panel

cat("\n--- Step 1: cross-sectional rectangle aggregation ---\n")
cat("Rectangles (raw, any n_hauls, pooled", H2_YEAR_MIN, "-", H2_YEAR_MAX, "):", n_rect_raw, "\n")
cat(
  "Rectangles after n_hauls >=", STEP0_MIN_HAULS_PROVISIONAL, "threshold:",
  hauls_thresholded$n_after, "(dropped", hauls_thresholded$n_dropped, ")\n"
)

rect_panel <- rect_panel_thresholded %>%
  left_join(couce_rect %>% select(stat_rec, mean_annual_hours_total), by = "stat_rec")

n_rect_no_fishing <- sum(is.na(rect_panel$mean_annual_hours_total))
rect_panel_with_fishing <- rect_panel %>% filter(!is.na(mean_annual_hours_total))

cat(
  "Rectangles with Couce fishing-hours data joined:", nrow(rect_panel_with_fishing),
  "(", n_rect_no_fishing, "dropped for missing Couce coverage)\n"
)

write_csv(rect_panel, path_out_rect_panel)
cat("Saved", path_out_rect_panel, "\n")

# ===========================================================================
# Step 2: cross-sectional association test (H2 case)
# ===========================================================================
cat("\n--- Step 2: cross-sectional association (H2 case) ---\n")

cross_D <- cross_sectional_association(rect_panel_with_fishing, "mean_D")
cross_cv <- cross_sectional_association(rect_panel_with_fishing, "mean_size_CV")

cat(sprintf(
  "mean_D ~ fishing_hours:       r=%.3f, R2=%.3f, slope=%.3g, p=%.3g, n=%d\n",
  cross_D$correlation, cross_D$r_squared, cross_D$slope, cross_D$p_value, cross_D$n
))
cat(sprintf(
  "mean_size_CV ~ fishing_hours: r=%.3f, R2=%.3f, slope=%.3g, p=%.3g, n=%d\n",
  cross_cv$correlation, cross_cv$r_squared, cross_cv$slope, cross_cv$p_value, cross_cv$n
))

step0_scatterplot <- function(df, x_col, x_lab, title, colour) {
  ggplot(df, aes(x = mean_annual_hours_total, y = .data[[x_col]])) +
    geom_point(alpha = 0.55, size = 1.7, colour = "grey30") +
    geom_smooth(method = "lm", se = TRUE, colour = colour, linewidth = 0.9) +
    geom_smooth(
      method = "loess", se = FALSE, colour = "#1a1a1a",
      linewidth = 0.7, linetype = "dashed"
    ) +
    labs(
      x = "Mean annual fishing hours (Couce et al. 2020)",
      y = x_lab,
      title = title,
      subtitle = sprintf(
        "n = %d rectangles; min %d hauls/rectangle (provisional); solid = linear fit, dashed = loess",
        nrow(df), STEP0_MIN_HAULS_PROVISIONAL
      )
    ) +
    theme_minimal(base_size = 11)
}

p_D <- step0_scatterplot(
  rect_panel_with_fishing, "mean_D", "Rectangle-mean D (Berger-Parker dominance)",
  "Step 0: rectangle-mean D vs fishing hours", "#d6604d"
)
ggsave(path_out_fig_D, p_D, width = 8, height = 5.5, dpi = 120)

p_cv <- step0_scatterplot(
  rect_panel_with_fishing, "mean_size_CV", "Rectangle-mean size_CV (dominant-species mass CV)",
  "Step 0: rectangle-mean size_CV vs fishing hours", "#4393c3"
)
ggsave(path_out_fig_cv, p_cv, width = 8, height = 5.5, dpi = 120)
cat("Saved", path_out_fig_D, "and", path_out_fig_cv, "\n")

# ===========================================================================
# Step 3: temporal within-rectangle association test (H3 case)
# ===========================================================================
cat("\n--- Step 3: temporal within-rectangle association (H3 case) ---\n")

annual_panel <- build_annual_rectangle_dominance(haul_dom_raw, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX)
years_thresholded <- apply_min_years_threshold(annual_panel, min_years = STEP0_MIN_YEARS_PROVISIONAL)

cat("Rectangles with any haul-year data:", years_thresholded$n_before, "\n")
cat(
  "Rectangles with >=", STEP0_MIN_YEARS_PROVISIONAL, "distinct years (provisional coverage threshold):",
  years_thresholded$n_after, "(dropped", years_thresholded$n_dropped, "for insufficient temporal coverage)\n"
)
cat(
  "H3 policy-change dates: none finalised / on file in this repo — reporting only the",
  "annual within-rectangle correlation; before/after mean-difference NOT computed.\n"
)

temporal_D <- within_rectangle_temporal_correlation(
  annual_panel, couce_year, years_thresholded$kept_rectangles,
  metric_col = "mean_D", min_paired_years = 3L
)
temporal_cv <- within_rectangle_temporal_correlation(
  annual_panel, couce_year, years_thresholded$kept_rectangles,
  metric_col = "mean_size_CV", min_paired_years = 3L
)

temporal_D_summary <- summarise_within_rectangle_correlations(temporal_D)
temporal_cv_summary <- summarise_within_rectangle_correlations(temporal_cv)

cat(sprintf(
  "D vs fishing_hours (within-rectangle):       n_rectangles=%d, n_with_valid_corr=%d, median_r=%.3f, IQR=[%.3f, %.3f], pct_pos=%s%%, pct_neg=%s%%\n",
  temporal_D_summary$n_rectangles, temporal_D_summary$n_with_valid_corr,
  temporal_D_summary$median_corr, temporal_D_summary$iqr_low, temporal_D_summary$iqr_high,
  temporal_D_summary$pct_positive, temporal_D_summary$pct_negative
))
cat(sprintf(
  "size_CV vs fishing_hours (within-rectangle): n_rectangles=%d, n_with_valid_corr=%d, median_r=%.3f, IQR=[%.3f, %.3f], pct_pos=%s%%, pct_neg=%s%%\n",
  temporal_cv_summary$n_rectangles, temporal_cv_summary$n_with_valid_corr,
  temporal_cv_summary$median_corr, temporal_cv_summary$iqr_low, temporal_cv_summary$iqr_high,
  temporal_cv_summary$pct_positive, temporal_cv_summary$pct_negative
))

temporal_by_rect <- bind_rows(
  temporal_D %>% mutate(variable = "D vs fishing_hours"),
  temporal_cv %>% mutate(variable = "size_CV vs fishing_hours")
) %>%
  relocate(variable, stat_rec)
write_csv(temporal_by_rect, path_out_temporal_by_rect)
cat("Saved", path_out_temporal_by_rect, "(per-rectangle detail underlying the summarised temporal rows)\n")

# ===========================================================================
# Step 4: output summary table (deliverable 1)
# ===========================================================================
cat("\n--- Step 4: summary table ---\n")

summary_table <- build_step0_summary_table(cross_D, cross_cv, temporal_D_summary, temporal_cv_summary)
print(summary_table)

provenance <- h2_collect_provenance(project_root)
summary_stamped <- h2_stamp_result(
  list(
    summary_table = summary_table,
    rectangle_panel = rect_panel,
    temporal_by_rectangle = temporal_by_rect,
    provisional_min_hauls = STEP0_MIN_HAULS_PROVISIONAL,
    provisional_min_years = STEP0_MIN_YEARS_PROVISIONAL,
    policy_dates_on_file = FALSE
  ),
  provenance
)
saveRDS(summary_stamped, path_out_summary_rds)
write_csv(summary_table, path_out_summary_csv)
cat("Saved", path_out_summary_csv, "\n")
cat("Saved", path_out_summary_rds, "\n")

# ---------------------------------------------------------------------------
# Exclusion-count audit trail (every drop visible, per the briefing)
# ---------------------------------------------------------------------------
exclusion_counts <- tibble::tribble(
  ~stage, ~n_rectangles, ~n_dropped_this_stage, ~note,
  "raw (any n_hauls, pooled period)", n_rect_raw, 0L, "",
  "after min-hauls threshold (cross-sectional)", hauls_thresholded$n_after, hauls_thresholded$n_dropped,
    sprintf("provisional threshold = %d hauls/rectangle", STEP0_MIN_HAULS_PROVISIONAL),
  "after requiring Couce fishing-hours join (cross-sectional)", nrow(rect_panel_with_fishing), n_rect_no_fishing,
    "rectangles with no Couce coverage",
  "after min-years threshold (temporal / H3 case)", years_thresholded$n_after, years_thresholded$n_dropped,
    sprintf("provisional threshold = %d distinct years", STEP0_MIN_YEARS_PROVISIONAL)
)
write_csv(exclusion_counts, path_out_exclusions)
cat("Saved", path_out_exclusions, "\n")

cat("\n=== Step 0 diagnostic complete — results are for review; no interpretation applied. ===\n")
