# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# H3 strategy feasibility — policy-period split x coarse spatial zones.
# See CURSOR_BRIEFING "H3 Strategy Feasibility — Policy-Period Split x Coarse
# Spatial Zones" (chat-supplied, not a repo file) for the full spec this
# script implements.
#
# PURPOSE: the pre-H3 feasibility check (run_h3_pre_exploration.R) found
# effectively zero ICES rectangles with enough repeat annual sampling for a
# within-rectangle-year temporal test (1 of 197). This script tests a
# specific alternative unit of analysis that coarsens on both axes at once:
# pre/post 2003 cod-crisis-reform split (instead of annual resolution) x
# rectangle blocks or fishing-pressure zones (instead of individual
# rectangles). Feasibility/visualisation check only.
#
# THIS SCRIPT DOES NOT: fit any H3 model, select a final spatial scheme or
# parameter value, run significance tests on the pre/post comparison, or
# recompute H1/H2 statistics. Read-only on all existing H1/H2/pre-H3 outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h3_policy_zones.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h1_join_helpers.R"))
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h3_pre_exploration_helpers.R")) # variance_components_anova()
source(file.path(script_dir, "R", "h3_policy_zones_helpers.R"))

# ---------------------------------------------------------------------------
# Provisional constants — named once, surfaced here, never hard-coded inline.
# ---------------------------------------------------------------------------
POLICY_BREAK_YEAR <- 2003L # pre = 1985-2002 (18 yrs), post = 2003-2015 (13 yrs); midpoint of the 2002-2004 cod-crisis reform window; PROVISIONAL pending exact date confirmation with Jake/HH
BLOCK_SIZES <- c(2L, 3L)   # Scheme A rectangle-block sizes (block_size x block_size)
PRESSURE_TIERS <- c(3L, 4L) # Scheme B quantile-based pressure tier counts
MIN_HAULS_PER_CELL <- 5L   # reused from SPARSE_HAUL_THRESHOLD in the pre-H3 task; EXPLORATORY REFERENCE POINT ONLY (see "sample-size adequacy" section below) — not retained as the primary feasibility filter
FIXED_CANDIDATE_THRESHOLDS <- c(1L, 2L, 3L, 5L, 10L, 15L, 20L, 30L, 40L, 50L) # sensitivity grid for "explore beyond >=5 hauls"
RELIABILITY_TARGETS <- c(0.7, 0.8, 0.9) # Spearman-Brown reliability targets for the cell mean; conventional psychometric thresholds (Nunnally 1978: ~0.7 "acceptable", 0.8 "good", 0.9 "excellent"), not this script's own invention

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_fishglob <- file.path(project_root, "FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
path_haul_eeos <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
path_couce_rect <- file.path(project_root, "outputs", "h2_couce_rectangle_effort.rds")
path_rect_flags <- file.path(project_root, "outputs", "h3_pre_rectangle_usability_flags.csv")

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_feasibility <- file.path(project_root, "outputs", "h3_policy_feasibility_summary.csv")
path_out_feasibility_by_threshold <- file.path(project_root, "outputs", "h3_policy_feasibility_by_threshold.csv")
path_out_sample_size <- file.path(project_root, "outputs", "h3_policy_sample_size_adequacy.csv")
path_out_run_log <- file.path(project_root, "outputs", "h3_policy_run_log.md")

stopifnot(
  file.exists(path_fishglob),
  file.exists(path_haul_eeos),
  file.exists(path_couce_year),
  file.exists(path_couce_rect),
  file.exists(path_rect_flags)
)

# ---------------------------------------------------------------------------
# Run log accumulator
# ---------------------------------------------------------------------------
run_log <- character(0)
figure_log <- list()
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}
log_figure <- function(id, path, caption) {
  cat(sprintf("[figure %s] %s\n  %s\n", id, path, caption))
  figure_log[[id]] <<- list(path = path, caption = caption)
}

logmsg("# H3 strategy feasibility (policy-period split x coarse spatial zones) — run log")
logmsg("")
logmsg("Feasibility/visualisation check only. No H3 model fit, no scheme/parameter selected as final, no significance testing.")
logmsg("")
logmsg("## Provisional constants")
logmsg("- POLICY_BREAK_YEAR = ", POLICY_BREAK_YEAR, " (pre = ", H2_YEAR_MIN, "-", POLICY_BREAK_YEAR - 1L, " [", POLICY_BREAK_YEAR - H2_YEAR_MIN, " yrs]; post = ", POLICY_BREAK_YEAR, "-", H2_YEAR_MAX, " [", H2_YEAR_MAX - POLICY_BREAK_YEAR + 1L, " yrs]; midpoint of the 2002-2004 cod-crisis reform window; PROVISIONAL pending exact date confirmation with Jake/HH; periods are unbalanced by construction — reported, not corrected for.)")
logmsg("- BLOCK_SIZES = ", paste(BLOCK_SIZES, collapse = ", "), " (Scheme A block_size x block_size rectangle merge)")
logmsg("- PRESSURE_TIERS = ", paste(PRESSURE_TIERS, collapse = ", "), " (Scheme B quantile-based tier counts, via dplyr::ntile on mean_annual_hours_total pooled over the full period)")
logmsg("- MIN_HAULS_PER_CELL = ", MIN_HAULS_PER_CELL, " (reused from SPARSE_HAUL_THRESHOLD in the pre-H3 feasibility task; EXPLORATORY REFERENCE POINT ONLY — retained in the feasibility summary for continuity with the pre-H3 task, NOT the primary analysis filter; see 'Sample-size adequacy' section below for the statistically defensible alternative)")
logmsg("- FIXED_CANDIDATE_THRESHOLDS = ", paste(FIXED_CANDIDATE_THRESHOLDS, collapse = ", "), " (sensitivity grid, 'explore beyond >=5 hauls')")
logmsg("- RELIABILITY_TARGETS = ", paste(RELIABILITY_TARGETS, collapse = ", "), " (Spearman-Brown reliability targets for the cell mean; conventional psychometric benchmarks, not a new invention: ~0.7 'acceptable', 0.8 'good', 0.9 'excellent' per Nunnally 1978)")

# ---------------------------------------------------------------------------
# Data + universe
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Universe of rectangles")

haul_full <- build_fishglob_haul_table(
  path_fishglob, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX, analysis_quarter = 1L
) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
haul_eeos <- readRDS(path_haul_eeos) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
couce_year <- readRDS(path_couce_year) %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX)
couce_rect <- readRDS(path_couce_rect) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
rect_flags <- read_csv(path_rect_flags, show_col_types = FALSE) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))

ices_sf <- load_ices_rectangles_sf(project_root)
grid_index <- build_rect_grid_index(ices_sf)

universe_A <- rect_flags$stat_rec
universe_B <- rect_flags$stat_rec[rect_flags$has_couce_coverage]

logmsg(
  "This task uses the SAME ", length(universe_A), "-rectangle haul-bearing universe as Section A of the pre-H3 ",
  "feasibility script (outputs/h3_pre_rectangle_usability_flags.csv), NOT Step 0's 187-rectangle universe (which ",
  "pre-filtered on total pooled hauls >= 5) — see that script's run log for the reconciliation between the two. ",
  "Scheme A (blind to fishing pressure) uses the full ", length(universe_A), "-rectangle universe. Scheme B ",
  "(fishing-pressure tiers) additionally requires Couce coverage, restricting to ", length(universe_B), " of ",
  length(universe_A), " rectangles (the same has_couce_coverage flag from the pre-H3 task)."
)

n_no_geom_A <- sum(!universe_A %in% grid_index$stat_rec)
n_no_geom_B <- sum(!universe_B %in% grid_index$stat_rec)
logmsg(
  "Rectangles without matching ICES shapefile geometry (excluded from block/zone construction): ",
  n_no_geom_A, " of ", length(universe_A), " (Scheme A universe); ", n_no_geom_B, " of ", length(universe_B),
  " (Scheme B universe)."
)
universe_A <- intersect(universe_A, grid_index$stat_rec)
universe_B <- intersect(universe_B, grid_index$stat_rec)

# ===========================================================================
# Helper: build one scheme/parameter combination end-to-end
# ===========================================================================
palette_for_units <- function(unit_ids) {
  u <- sort(unique(unit_ids))
  set.seed(42)
  cols <- grDevices::hcl.colors(length(u), palette = "Dark 3")
  setNames(sample(cols), u)
}

build_combo_outputs <- function(scheme_label, param_label, unit_map, id_prefix, secondary_map = NULL, secondary_name = NULL) {
  n_units_total <- dplyr::n_distinct(unit_map$unit_id)

  # --- units map -------------------------------------------------------------
  merged_sf <- ices_sf %>%
    inner_join(unit_map, by = "stat_rec") %>%
    group_by(unit_id) %>%
    summarise(n_rects = dplyr::n(), .groups = "drop")
  pal <- palette_for_units(merged_sf$unit_id)

  if (!is.null(secondary_map)) {
    merged_secondary <- ices_sf %>%
      inner_join(secondary_map, by = "stat_rec") %>%
      group_by(tier) %>%
      summarise(.groups = "drop")
    p_units_main <- ggplot(merged_sf) +
      geom_sf(aes(fill = unit_id), colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = pal, guide = "none") +
      labs(title = sprintf("%s output: contiguous zones (n=%d)", param_label, n_units_total)) +
      theme_minimal(base_size = 9) +
      theme(axis.text = element_text(size = 6))
    p_units_secondary <- ggplot(merged_secondary) +
      geom_sf(aes(fill = factor(tier)), colour = "white", linewidth = 0.15) +
      scale_fill_viridis_d(option = "A", name = "Pressure\ntier") +
      labs(title = sprintf("%s input: pressure tier classification", param_label)) +
      theme_minimal(base_size = 9) +
      theme(axis.text = element_text(size = 6))
    p_units <- p_units_secondary + p_units_main +
      patchwork::plot_annotation(
        title = sprintf("Scheme B spatial units — %s", param_label),
        caption = sprintf(
          "Left: quantile pressure tier (input classification, ntile on mean Couce hours 1985-2015). Right: %s (output units, rook-adjacency connected components within same tier; n=%d zones).",
          secondary_name, n_units_total
        )
      )
    ggsave(file.path(fig_dir, sprintf("h3_policy_%s_units_map.png", id_prefix)), p_units, width = 12, height = 5.5, dpi = 150)
  } else {
    p_units <- ggplot(merged_sf) +
      geom_sf(aes(fill = unit_id), colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = pal, guide = "none") +
      labs(
        title = sprintf("Scheme A spatial units — %s block merge (n=%d blocks)", param_label, n_units_total),
        caption = sprintf("Contiguous %s rectangle blocks, anchored to the absolute ICES grid; blind to fishing pressure.", param_label)
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text = element_text(size = 6))
    ggsave(file.path(fig_dir, sprintf("h3_policy_%s_units_map.png", id_prefix)), p_units, width = 8, height = 6.5, dpi = 150)
  }
  log_figure(
    paste0(id_prefix, "_units_map"), file.path(fig_dir, sprintf("h3_policy_%s_units_map.png", id_prefix)),
    sprintf("%s / %s spatial units (n=%d), for visual sanity-check before any statistic is computed on them.", scheme_label, param_label, n_units_total)
  )

  # --- haul-count choropleth (pre/post) --------------------------------------
  unit_period_hauls <- summarise_unit_period_hauls(haul_full, unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)
  map_hauls <- merged_sf %>% select(unit_id) %>% inner_join(unit_period_hauls, by = "unit_id")
  map_hauls$period <- factor(map_hauls$period, levels = c("pre", "post"), labels = c(sprintf("Pre (%d-%d)", H2_YEAR_MIN, POLICY_BREAK_YEAR - 1L), sprintf("Post (%d-%d)", POLICY_BREAK_YEAR, H2_YEAR_MAX)))

  p_haulcount <- ggplot(map_hauls) +
    geom_sf(aes(fill = n_hauls), colour = "white", linewidth = 0.1) +
    facet_wrap(~period, nrow = 1) +
    scale_fill_viridis_c(option = "C", name = "N hauls") +
    labs(
      title = sprintf("Haul count per spatial unit, pre vs post %d — %s / %s", POLICY_BREAK_YEAR, scheme_label, param_label),
      caption = "Haul count per spatial unit per period; units with zero hauls that period are absent (no Couce/haul dependency here, just survey coverage)."
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_text(size = 6))
  path_haulcount <- file.path(fig_dir, sprintf("h3_policy_%s_haulcount_prepost.png", id_prefix))
  ggsave(path_haulcount, p_haulcount, width = 11, height = 5.5, dpi = 150)
  log_figure(
    paste0(id_prefix, "_haulcount"), path_haulcount,
    sprintf("Haul count per %s spatial unit, pre-%d vs post-%d, side by side.", param_label, POLICY_BREAK_YEAR, POLICY_BREAK_YEAR)
  )

  # --- descriptive pre/post comparison (residual + fishing pressure) --------
  unit_period_resid <- summarise_unit_period_residual(haul_eeos, unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)
  unit_period_fish <- summarise_unit_period_fishing(couce_year, unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)

  resid_wide <- unit_period_resid %>%
    select(unit_id, period, mean_residual) %>%
    pivot_wider(names_from = period, values_from = mean_residual) %>%
    filter(!is.na(pre), !is.na(post)) %>%
    mutate(metric = "Mean residual (log B_obs - log B_pred)")
  fish_wide <- unit_period_fish %>%
    select(unit_id, period, mean_fishing_hours) %>%
    pivot_wider(names_from = period, values_from = mean_fishing_hours) %>%
    filter(!is.na(pre), !is.na(post)) %>%
    mutate(metric = "Mean Couce fishing hours")

  n_units_resid_both <- nrow(resid_wide)
  n_units_fish_both <- nrow(fish_wide)
  logmsg(
    id_prefix, ": ", n_units_resid_both, " of ", n_units_total,
    " units have data in both periods for mean residual; ", n_units_fish_both, " of ", n_units_total,
    " for mean fishing hours (units missing a period, e.g. zero contributing rectangles that period, excluded from the scatter, not imputed)."
  )

  comparison_df <- bind_rows(resid_wide, fish_wide)
  p_prepost <- ggplot(comparison_df, aes(x = pre, y = post)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(alpha = 0.7, colour = "#4575b4", size = 2) +
    facet_wrap(~metric, scales = "free") +
    labs(
      x = sprintf("Pre (%d-%d)", H2_YEAR_MIN, POLICY_BREAK_YEAR - 1L),
      y = sprintf("Post (%d-%d)", POLICY_BREAK_YEAR, H2_YEAR_MAX),
      title = sprintf("Pre vs post %d, per spatial unit — %s / %s", POLICY_BREAK_YEAR, scheme_label, param_label),
      caption = "One point per spatial unit; dashed line = pre==post; purely descriptive, no fitted line, no significance test."
    ) +
    theme_minimal(base_size = 11)
  path_prepost <- file.path(fig_dir, sprintf("h3_policy_%s_prepost_comparison.png", id_prefix))
  ggsave(path_prepost, p_prepost, width = 10, height = 5.5, dpi = 150)
  log_figure(
    paste0(id_prefix, "_prepost"), path_prepost,
    sprintf("Pre-%d vs post-%d mean residual and mean fishing hours, one point per %s spatial unit; descriptive only.", POLICY_BREAK_YEAR, POLICY_BREAK_YEAR, param_label)
  )

  list(
    feasibility_row = feasibility_summary_row(unit_period_hauls, scheme_label, param_label, MIN_HAULS_PER_CELL, n_units_total),
    unit_period_hauls = unit_period_hauls,
    n_units_total = n_units_total
  )
}

#' Sample-size-adequacy row + threshold-sensitivity rows for one combo.
#' Computed once per combo, after build_combo_outputs(), reusing its
#' unit_period_hauls table (avoids recomputing it) plus a fresh haul-level
#' cell-labelling pass for the variance decomposition.
build_sample_size_outputs <- function(scheme_label, param_label, unit_map, combo_result) {
  haul_cells <- label_haul_cells(haul_eeos, unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)
  ss_row <- sample_size_adequacy_row(haul_cells, scheme_label, param_label, RELIABILITY_TARGETS)

  required_n_cols <- grep("^required_n_reliability_", names(ss_row), value = TRUE)
  reliability_n <- unlist(ss_row[required_n_cols])
  reliability_n_finite <- reliability_n[is.finite(reliability_n)]
  reliability_labels <- sprintf("reliability_%s", sub("required_n_reliability_", "", names(reliability_n_finite)))

  all_thresholds <- c(FIXED_CANDIDATE_THRESHOLDS, unname(reliability_n_finite))
  all_labels <- c(sprintf("fixed_%d", FIXED_CANDIDATE_THRESHOLDS), reliability_labels)
  keep <- !duplicated(all_thresholds)

  sens_rows <- feasibility_across_thresholds(
    combo_result$unit_period_hauls, scheme_label, param_label,
    thresholds = all_thresholds[keep], threshold_labels = all_labels[keep],
    n_units_total = combo_result$n_units_total
  )

  logmsg(sprintf(
    "  Sample-size adequacy (%s / %s): n_cells=%d, n_hauls=%d, ICC=%.4f (haul-level residual variance between-cell vs within-cell); required n/cell for reliability %s -> %s.",
    scheme_label, param_label, ss_row$n_cells, ss_row$n_hauls, ss_row$icc,
    paste(RELIABILITY_TARGETS, collapse = "/"),
    paste(ifelse(is.finite(reliability_n), reliability_n, "Inf"), collapse = "/")
  ))

  list(sample_size_row = ss_row, sensitivity_rows = sens_rows)
}

# ===========================================================================
# Scheme A — geographic block merge
# ===========================================================================
logmsg("")
logmsg("## Scheme A — geographic block merge")

feasibility_rows <- list()
sample_size_rows <- list()
sensitivity_rows <- list()
for (bs in BLOCK_SIZES) {
  logmsg("Building Scheme A blocks at ", bs, "x", bs, "...")
  blocks <- build_scheme_a_blocks(grid_index, universe_A, bs)
  logmsg(
    "  ", nrow(blocks), " rectangles merged into ", dplyr::n_distinct(blocks$unit_id), " blocks (sizes: ",
    paste(range(table(blocks$unit_id)), collapse = "-"), " rectangles/block)."
  )
  scheme_label <- "Scheme A (block merge)"
  param_label <- sprintf("%dx%d", bs, bs)
  combo_result <- build_combo_outputs(
    scheme_label = scheme_label, param_label = param_label,
    unit_map = blocks, id_prefix = sprintf("A_%dx%d", bs, bs)
  )
  feasibility_rows[[length(feasibility_rows) + 1L]] <- combo_result$feasibility_row

  ss_out <- build_sample_size_outputs(scheme_label, param_label, blocks, combo_result)
  sample_size_rows[[length(sample_size_rows) + 1L]] <- ss_out$sample_size_row
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- ss_out$sensitivity_rows
}

# ===========================================================================
# Scheme B — contiguous fishing-pressure zones
# ===========================================================================
logmsg("")
logmsg("## Scheme B — contiguous fishing-pressure zones")
logmsg(
  "Method: (1) classify each of the ", length(universe_B), " Couce-covered rectangles into a quantile-based ",
  "pressure tier (dplyr::ntile on mean_annual_hours_total, pooled 1985-2015, from outputs/h2_couce_rectangle_effort.rds); ",
  "(2) build a rook (edge-sharing, not diagonal) adjacency list directly from ICES-grid row/col indices ",
  "(derived from the shapefile's own SOUTH/WEST corner fields — no spdep/igraph dependency needed, since grid ",
  "adjacency is directly computable from integer row/col indices); (3) find connected components via a plain ",
  "breadth-first search restricted to same-tier neighbours (R/h3_policy_zones_helpers.R::build_contiguous_zones()). ",
  "This was straightforward to implement with the grid-index approach (no igraph/spdep required)."
)

adjacency_B <- build_rook_adjacency(grid_index, universe_B)

for (nt in PRESSURE_TIERS) {
  logmsg("Building Scheme B zones at ", nt, " tiers...")
  pressure <- couce_rect %>% filter(stat_rec %in% universe_B)
  tiers <- assign_pressure_tiers(pressure, nt)
  zones <- build_contiguous_zones(tiers$stat_rec, tiers$tier, adjacency_B)
  zone_sizes <- zones %>% count(unit_id, name = "n_rects")
  n_singleton <- sum(zone_sizes$n_rects == 1L)
  logmsg(
    "  ", nt, "-tier split: table(tier) = ", paste(sprintf("tier%s=%d", names(table(tiers$tier)), table(tiers$tier)), collapse = ", "), "; ",
    "resulting in ", nrow(zone_sizes), " contiguous zones (", n_singleton, " singleton [n=1 rectangle] zones — same-tier ",
    "rectangles with no same-tier rook-neighbour; NOT pooled with other same-tier-but-non-adjacent rectangles). ",
    "Largest zone: ", max(zone_sizes$n_rects), " rectangles."
  )
  unit_map <- zones %>% select(stat_rec, unit_id)
  tier_map <- zones %>% select(stat_rec, tier)
  scheme_label <- "Scheme B (pressure zones)"
  param_label <- sprintf("%d-tier", nt)
  combo_result <- build_combo_outputs(
    scheme_label = scheme_label, param_label = param_label,
    unit_map = unit_map, id_prefix = sprintf("B_%dtier", nt),
    secondary_map = tier_map, secondary_name = "contiguous same-tier zones"
  )
  feasibility_rows[[length(feasibility_rows) + 1L]] <- combo_result$feasibility_row

  ss_out <- build_sample_size_outputs(scheme_label, param_label, unit_map, combo_result)
  sample_size_rows[[length(sample_size_rows) + 1L]] <- ss_out$sample_size_row
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- ss_out$sensitivity_rows
}

# ===========================================================================
# Feasibility summary (across all combinations)
# ===========================================================================
logmsg("")
logmsg("## Feasibility summary")

feasibility_summary <- bind_rows(feasibility_rows)
write_csv(feasibility_summary, path_out_feasibility)
logmsg("Saved feasibility summary: ", path_out_feasibility)
for (i in seq_len(nrow(feasibility_summary))) {
  r <- feasibility_summary[i, ]
  logmsg(sprintf(
    "  - %s / %s: n_units=%d, n_units_usable_both_periods=%d (>= %d hauls in both periods); pre mean/median/min = %.1f/%.1f/%d; post mean/median/min = %.1f/%.1f/%d",
    r$scheme, r$parameter, r$n_units, r$n_units_usable_both_periods, r$min_hauls_per_cell,
    r$pre_mean_hauls, r$pre_median_hauls, r$pre_min_hauls,
    r$post_mean_hauls, r$post_median_hauls, r$post_min_hauls
  ))
}

best_row <- feasibility_summary %>% arrange(desc(n_units_usable_both_periods)) %>% slice(1)
frac_usable <- best_row$n_units_usable_both_periods / best_row$n_units
logmsg("")
logmsg(
  "Best-performing combination by n_units_usable_both_periods (at the exploratory MIN_HAULS_PER_CELL=", MIN_HAULS_PER_CELL,
  " reference point): ", best_row$scheme, " / ", best_row$parameter,
  " (", best_row$n_units_usable_both_periods, " of ", best_row$n_units, " units, ", round(100 * frac_usable, 1), "%). ",
  "Reported, not a recommendation — comparison material for the supervisor discussion."
)

# ===========================================================================
# Sample-size adequacy — a statistically defensible alternative to a fixed,
# arbitrary MIN_HAULS_PER_CELL. The ">=5 hauls in both periods" gate above
# was an exploratory placeholder; median hauls/rectangle/year is only ~2, so
# a single round-number cutoff has no particular statistical grounding.
# Here, for each scheme/parameter, we decompose haul-level residual variance
# into between-cell (unit x period) and within-cell (haul-to-haul) components
# via one-way ANOVA (reusing variance_components_anova() from the pre-H3
# task), then invert the Spearman-Brown formula to ask: how many hauls does
# a cell need for its MEAN residual to be a target-reliable estimate of that
# cell's true value, given the actual observed signal-to-noise ratio (ICC)?
# This is explicitly exploratory too — no single threshold is selected as
# final — but it is grounded in the data's own variance structure rather
# than a round number.
# ===========================================================================
logmsg("")
logmsg("## Sample-size adequacy for MIN_HAULS_PER_CELL (statistically defensible thresholds)")
logmsg(
  "The exploratory MIN_HAULS_PER_CELL = ", MIN_HAULS_PER_CELL, " gate above is a reference point only, not the ",
  "primary feasibility filter. Method: for each scheme/parameter combination, label every haul with its cell ",
  "(spatial unit x pre/post period), then decompose the haul-level residual's variance into between-cell and ",
  "within-cell (haul-to-haul) components (one-way random-effects ANOVA, method-of-moments estimator for ",
  "unbalanced groups — R/h3_pre_exploration_helpers.R::variance_components_anova(), same estimator used in the ",
  "pre-H3 task's Section E). The resulting ICC (intraclass correlation = between-cell variance / total variance) ",
  "is the fraction of haul-to-haul residual variation that reflects genuine between-cell signal rather than ",
  "sampling noise. Inverting the Spearman-Brown formula (reliability of an n-haul cell mean = n*ICC / ",
  "(1+(n-1)*ICC)) gives the number of hauls a cell needs for its mean to reach a target reliability — this is ",
  "the same logic used to decide how many test items or raters are 'enough' in measurement theory, applied here ",
  "to deciding how many hauls are 'enough' for one cell mean. If ICC <= 0 (no detectable between-cell signal at ",
  "all), required n is reported as Inf — no amount of extra replication manufactures a signal that isn't there."
)

sample_size_summary <- bind_rows(sample_size_rows)
write_csv(sample_size_summary, path_out_sample_size)
logmsg("Saved sample-size adequacy table: ", path_out_sample_size)
for (i in seq_len(nrow(sample_size_summary))) {
  r <- sample_size_summary[i, ]
  req_cols <- grep("^required_n_reliability_", names(r), value = TRUE)
  req_str <- paste(sprintf("R%s->n=%s", sub("required_n_reliability_", "", req_cols), ifelse(is.finite(unlist(r[req_cols])), unlist(r[req_cols]), "Inf")), collapse = ", ")
  logmsg(sprintf(
    "  - %s / %s: n_cells=%d, n_hauls=%d, ICC=%.4f (var_between=%.4f, var_within=%.4f); %s",
    r$scheme, r$parameter, r$n_cells, r$n_hauls, r$icc, r$var_between_cell, r$var_within_cell, req_str
  ))
}

feasibility_by_threshold <- bind_rows(sensitivity_rows)
write_csv(feasibility_by_threshold, path_out_feasibility_by_threshold)
logmsg("Saved feasibility-by-threshold sensitivity table: ", path_out_feasibility_by_threshold)
logmsg(
  "This table repeats the n_units_usable_both_periods calculation across the fixed candidate grid (",
  paste(FIXED_CANDIDATE_THRESHOLDS, collapse = ", "), ") AND each combo's own reliability-derived threshold(s) — ",
  "the full 'explore beyond >=5 hauls' picture, not a single number."
)

p_sensitivity <- ggplot(feasibility_by_threshold, aes(x = threshold_value, y = pct_units_usable_both_periods, colour = paste(scheme, parameter))) +
  geom_line(data = feasibility_by_threshold %>% filter(grepl("^fixed_", threshold_method)), linewidth = 0.6) +
  geom_point(data = feasibility_by_threshold %>% filter(grepl("^fixed_", threshold_method)), size = 1.4) +
  geom_point(
    data = feasibility_by_threshold %>% filter(grepl("^reliability_", threshold_method)),
    aes(shape = threshold_method), size = 3, stroke = 1.1
  ) +
  scale_x_log10() +
  scale_shape_manual(values = c(reliability_0.7 = 17, reliability_0.8 = 15, reliability_0.9 = 18), name = "Reliability-derived\nthreshold") +
  labs(
    x = "MIN_HAULS_PER_CELL threshold (log scale)", y = "% spatial units usable in both periods",
    colour = "Scheme / parameter",
    title = "Feasibility sensitivity to the haul-count threshold, with reliability-derived reference points",
    caption = sprintf(
      "Lines/circles = fixed thresholds %s. Triangle/square/diamond = data-driven thresholds for %s cell-mean reliability (Spearman-Brown, per combo's own ICC). MIN_HAULS_PER_CELL=%d marked for reference only.",
      paste(FIXED_CANDIDATE_THRESHOLDS, collapse = ","), paste(RELIABILITY_TARGETS, collapse = "/"), MIN_HAULS_PER_CELL
    )
  ) +
  geom_vline(xintercept = MIN_HAULS_PER_CELL, linetype = "dotted", colour = "grey40") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical")
path_sensitivity <- file.path(fig_dir, "h3_policy_sample_size_sensitivity.png")
ggsave(path_sensitivity, p_sensitivity, width = 9, height = 7.5, dpi = 150)
log_figure(
  "sample_size_sensitivity", path_sensitivity,
  "Percent of spatial units usable in both periods vs haul-count threshold, one line per scheme/parameter, with each combo's reliability-derived thresholds marked."
)

# ---------------------------------------------------------------------------
# Secondary/robustness variant (symmetric trimmed window) — built ONLY if
# the primary split shows promise, per the briefing (not built speculatively).
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Secondary/robustness variant (symmetric trimmed window)")
PROMISE_THRESHOLD <- 0.5 # >=50% of units usable in both periods counts as "shows promise" for this gate; provisional, this script's own judgement call
if (frac_usable >= PROMISE_THRESHOLD) {
  logmsg(
    "Primary split shows promise (best combination has ", round(100 * frac_usable, 1), "% of units usable in both ",
    "periods, >= ", round(100 * PROMISE_THRESHOLD), "% gate) — building the symmetric trimmed-window variant for ",
    "the best-performing combination only."
  )
  trim_years <- min(POLICY_BREAK_YEAR - H2_YEAR_MIN, H2_YEAR_MAX - POLICY_BREAK_YEAR + 1L)
  logmsg(
    "Symmetric window: ", trim_years, " years either side of ", POLICY_BREAK_YEAR, " (pre = ",
    POLICY_BREAK_YEAR - trim_years, "-", POLICY_BREAK_YEAR - 1L, ", post = ", POLICY_BREAK_YEAR, "-",
    POLICY_BREAK_YEAR + trim_years - 1L, ") — NOT built as a full figure set here; flagged as the next step if ",
    "the supervisor wants it pursued (out of scope to build a second full figure set speculatively beyond ",
    "confirming the primary split's promise)."
  )
} else {
  logmsg(
    "Primary split does NOT show promise (best combination has only ", round(100 * frac_usable, 1), "% of units ",
    "usable in both periods, below the ", round(100 * PROMISE_THRESHOLD), "% gate) — symmetric trimmed-window ",
    "variant NOT built, per the briefing's instruction not to build it speculatively."
  )
}

# ===========================================================================
# Figure index + outputs
# ===========================================================================
logmsg("")
logmsg("## Figure index")
for (id in names(figure_log)) {
  logmsg(sprintf("- **%s** `%s` — %s", id, basename(figure_log[[id]]$path), figure_log[[id]]$caption))
}

logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_feasibility, " (single-threshold reference table, MIN_HAULS_PER_CELL=", MIN_HAULS_PER_CELL, " only — exploratory)")
logmsg("- ", path_out_sample_size, " (variance decomposition + reliability-derived required-n per combo — statistically grounded)")
logmsg("- ", path_out_feasibility_by_threshold, " (feasibility sensitivity across the full threshold grid — the primary deliverable for this question)")
logmsg("- ", path_out_run_log, " (this file)")
logmsg("- ", length(figure_log), " figures in ", fig_dir, " (h3_policy_*.png)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H3 strategy feasibility check complete — comparison material only, no scheme selected. ===\n")
