# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Quantitative support for the H2/H3 design decision (temporal phasing +
# spatial-unit definition). See CURSOR_BRIEFING "Quantitative Support for
# H2/H3 Design Decision" (chat-supplied, not a repo file) for the full spec
# this script implements.
#
# PURPOSE: extract specific numbers behind patterns already visible in
# existing figures (h3_pre_D1/D2/D3), and extend the existing
# variance-decomposition/ICC method to compare five candidate spatial-unit
# definitions side by side. REPORTS NUMBERS ONLY — no temporal phase
# boundary, spatial scheme, or resolution is recommended anywhere in this
# script; that decision is made separately, informed by this output.
#
# THIS SCRIPT DOES NOT: fit any H3 model, run new significance tests beyond
# the specified correlations/linear trends, or run structural-break
# detection (flagged as a possible follow-up only). Read-only on all
# existing H1/H2/pre-H3/policy-zone outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(sf)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_design_support.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h1_join_helpers.R"))
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h3_pre_exploration_helpers.R")) # variance_components_anova, assign_decade_bin, build_rect_decade_summary, add_resid_signed
source(file.path(script_dir, "R", "h3_policy_zones_helpers.R"))    # grid index, Scheme A/B builders, required_n_for_reliability, unit x period summaries
source(file.path(script_dir, "R", "h2h3_design_support_helpers.R"))

# ---------------------------------------------------------------------------
# Provisional constants — named once, surfaced here and in the run log,
# never hard-coded inline. NONE of these constitute a design recommendation.
# ---------------------------------------------------------------------------
PHASE_BOUNDARIES <- list(
  list(label = "1985-1990", year_from = 1985L, year_to = 1990L),
  list(label = "1990-2000", year_from = 1990L, year_to = 2000L),
  list(label = "2000-2010", year_from = 2000L, year_to = 2010L),
  list(label = "2010-2015", year_from = 2010L, year_to = 2015L)
) # visual/approximate boundaries from h3_pre_D1/D3, per the briefing; boundary years (1990/2000/2010) are shared between adjacent phases, not split — provisional
DECADE_BINS <- list(c(1985L, 1994L), c(1995L, 2004L), c(2005L, 2015L)) # reused as-is from run_h3_pre_exploration.R — same bins, not re-chosen here
RELIABILITY_TARGETS <- c(0.7, 0.8, 0.9) # reused as-is from run_h3_policy_zones.R (Spearman-Brown targets, Nunnally 1978 benchmarks)
BLOCK_SIZES <- c(2L, 3L)   # Scheme A, reused as-is from run_h3_policy_zones.R
PRESSURE_TIERS <- c(3L, 4L) # Scheme B, reused as-is from run_h3_policy_zones.R
POLICY_BREAK_YEAR <- 2003L  # reused as-is from run_h3_policy_zones.R (pre/post split for Section D only)
MIN_YEARS_RELAXED <- 2L     # Section C individual-rectangle relaxed inclusion rule: has_couce_coverage AND n_years_present >= this (replaces the stricter usable_for_fishing_analysis flag, which returned 0 rectangles)

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

path_out_A1 <- file.path(project_root, "outputs", "h2h3_designA1_year_fishing_summary.csv")
path_out_A2 <- file.path(project_root, "outputs", "h2h3_designA2_phase_boundary_changes.csv")
path_out_A3 <- file.path(project_root, "outputs", "h2h3_designA3_phase_trends.csv")
path_out_B1 <- file.path(project_root, "outputs", "h2h3_designB1_rect_decade_wide.csv")
path_out_B2 <- file.path(project_root, "outputs", "h2h3_designB2_rect_decade_change.csv")
path_out_B3 <- file.path(project_root, "outputs", "h2h3_designB3_change_classification.csv")
path_out_B4 <- file.path(project_root, "outputs", "h2h3_designB4_decade_persistence_correlation.csv")
path_out_C <- file.path(project_root, "outputs", "h2h3_designC_icc_by_spatial_unit.csv")
path_out_D <- file.path(project_root, "outputs", "h2h3_designD_haul_unevenness.csv")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_design_support_run_log.md")

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
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# Quantitative support for the H2/H3 design decision — run log")
logmsg("")
logmsg("Reports numbers only. No temporal phase boundary, spatial scheme, or resolution is recommended anywhere in this script or its outputs — that decision is made separately.")
logmsg("")
logmsg("## Deviation from the brief's output format")
logmsg(
  "The brief requests 'one CSV per section (A, B, C, D)'. Sections A, B and D each specify multiple genuinely ",
  "different table shapes (e.g. A: a year x metric series, a phase-boundary-change table, a within-phase-trend ",
  "table). Rather than force heterogeneous tables into one wide/long CSV (which would be harder to read and use ",
  "than the source tables), each sub-table is saved as its own CSV, prefixed by section (A1/A2/A3, B1-B4), plus a ",
  "single CSV for C (which genuinely is one table) and one for D. This mirrors the existing repo convention (see ",
  "outputs/h3_pre_*.csv, outputs/h3_policy_*.csv — always one file per logical table, never force-merged). Flagged ",
  "here as a deviation from the literal instruction, not silently done."
)

# ---------------------------------------------------------------------------
# Data + rectangle universes
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Rectangle universes (reconciliation)")

couce_year <- readRDS(path_couce_year) %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX)
couce_rect <- readRDS(path_couce_rect) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
rect_flags <- read_csv(path_rect_flags, show_col_types = FALSE) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
haul_full <- build_fishglob_haul_table(
  path_fishglob, year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX, analysis_quarter = 1L
) %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
haul_eeos <- readRDS(path_haul_eeos) %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
  add_resid_signed()

ices_sf <- load_ices_rectangles_sf(project_root)
grid_index <- build_rect_grid_index(ices_sf)

n_couce_full <- dplyr::n_distinct(couce_year$stat_rec)
universe_197 <- rect_flags$stat_rec
universe_165 <- rect_flags$stat_rec[rect_flags$has_couce_coverage]

logmsg(
  "THREE distinct rectangle universes appear across this task's tables — none is more 'correct' than another, ",
  "they answer different questions:"
)
logmsg(
  "  (1) **", n_couce_full, "-rectangle full Couce universe** (`outputs/h2_couce_year_effort.rds`'s own distinct ",
  "stat_rec set) — every rectangle Couce et al. (2020) provide a fishing-hours reconstruction for, REGARDLESS of ",
  "whether NS-IBTS Q1 ever hauled there. This is the SAME universe used by the existing h3_pre_D1/D2 figures ",
  "('across all rectangles with coverage' — see outputs/exploratory_review/index.md Task 1 D.1/D.2 rows) and is ",
  "the universe used in **Section A (A1-A3) and Section B (B1-B4)** below, since those sections describe the ",
  "Couce fishing-pressure dataset on its own terms, extending those two existing figures directly."
)
logmsg(
  "  (2) **197-rectangle NS-IBTS haul-bearing universe** (`outputs/h3_pre_rectangle_usability_flags.csv`, all ",
  "rows) — every rectangle with >=1 Q1 haul, 1985-2015, regardless of Couce coverage. Used to build Scheme A ",
  "(2x2/3x3 block-merge) spatial units in **Section C and Section D** (Scheme A is blind to fishing pressure by ",
  "construction, per run_h3_policy_zones.R, so it is built on the full haul-bearing universe, not just the ",
  "Couce-covered subset)."
)
logmsg(
  "  (3) **", length(universe_165), "-rectangle Couce-covered subset of the 197** (`has_couce_coverage` flag in ",
  "the same CSV) — rectangles with BOTH >=1 Q1 haul AND >=1 Couce fishing-hours record. Used to build Scheme B ",
  "(3-/4-tier pressure zones) in **Section C and Section D** (Scheme B needs a fishing-pressure value per ",
  "rectangle to assign tiers, per run_h3_policy_zones.R), AND as the relaxed-rule individual-rectangle universe ",
  "for **Section C item 1** (see below — the relaxed rule turns out to select exactly this same 165-rectangle set)."
)
logmsg(
  n_couce_full - length(universe_165), " rectangles are in the full Couce universe (1) but NOT in the 197-rectangle ",
  "haul-bearing universe (2) at all, i.e. Couce provides a fishing-hours reconstruction for water NS-IBTS Q1 never ",
  "sampled in this period — consistent with the existing h3_pre run log's note that D uses outputs/h2_couce_year_effort.rds directly, unrestricted by the haul universe."
)

# ===========================================================================
# Section A — whole-study-area temporal trend in fishing pressure
# ===========================================================================
logmsg("")
logmsg("## Section A — whole-study-area temporal trend in fishing pressure")
logmsg(
  "Universe: full Couce dataset (", n_couce_full, " rectangles), matching h3_pre_D1 exactly. Phase boundaries ",
  "(", paste(vapply(PHASE_BOUNDARIES, function(p) p$label, character(1)), collapse = ", "), ") are the visual/",
  "approximate boundaries named in the briefing (from h3_pre_D1/D3), NOT statistically derived — treated as ",
  "provisional throughout this section. Boundary years are shared between adjacent phases (e.g. 1990 is the last ",
  "year of phase 1 AND the first year of phase 2), a provisional modelling choice, not an error."
)

year_summary <- build_year_fishing_summary(couce_year, H2_YEAR_MIN, H2_YEAR_MAX)
write_csv(year_summary, path_out_A1)
logmsg(
  "A1 saved: ", path_out_A1, " (", nrow(year_summary), " year rows, 1985-2015; n_rect_contributing ranges ",
  min(year_summary$n_rect_contributing), "-", max(year_summary$n_rect_contributing), " rectangles/year)."
)

a2_transitions <- bind_rows(
  lapply(seq_len(length(PHASE_BOUNDARIES)), function(i) {
    p <- PHASE_BOUNDARIES[[i]]
    phase_boundary_change(year_summary, p$year_from, p$year_to, p$label)
  })
)
a2_full <- phase_boundary_change(year_summary, H2_YEAR_MIN, H2_YEAR_MAX, sprintf("%d-%d (full period)", H2_YEAR_MIN, H2_YEAR_MAX))
a2_table <- bind_rows(a2_transitions, a2_full)
write_csv(a2_table, path_out_A2)
logmsg("A2 saved: ", path_out_A2, " (", nrow(a2_table), " rows: 4 phase-boundary transitions + full 1985 vs 2015).")
for (i in seq_len(nrow(a2_table))) {
  r <- a2_table[i, ]
  logmsg(sprintf(
    "  - %s: mean hours %.2f -> %.2f (abs change %+.2f, %+.1f%%)",
    r$transition, r$mean_hours_from, r$mean_hours_to, r$abs_change, r$pct_change
  ))
}

a3_table <- bind_rows(
  lapply(PHASE_BOUNDARIES, function(p) fit_phase_trend(year_summary, p$year_from, p$year_to, p$label))
)
write_csv(a3_table, path_out_A3)
logmsg("A3 saved: ", path_out_A3, " (within-phase linear trends, one row per phase; NOT one 30-year trend).")
for (i in seq_len(nrow(a3_table))) {
  r <- a3_table[i, ]
  logmsg(sprintf(
    "  - %s (n=%d yrs): slope=%.3f hours/yr (%s, %.2f%%/yr of phase mean), R^2=%.3f",
    r$phase, r$n_years, r$slope_hours_per_year, r$direction, r$slope_pct_of_mean_per_year, r$r_squared
  ))
}
logmsg(
  "Follow-up flag (NOT run here, per the briefing): a formal structural-break detection (e.g. `strucchange::breakpoints()` ",
  "on the annual mean-hours series) could sharpen these four visual phase boundaries into statistically located ",
  "break years, rather than the eyeballed 1990/2000/2010 boundaries used above. Flagged as a possible follow-up only."
)

# ===========================================================================
# Section B — spatial heterogeneity of fishing pressure by decade
# ===========================================================================
logmsg("")
logmsg("## Section B — spatial heterogeneity of fishing pressure by decade")
logmsg(
  "Universe: full Couce dataset (", n_couce_full, " rectangles), matching h3_pre_D2 exactly. DECADE_BINS reused ",
  "as-is from run_h3_pre_exploration.R: ", paste(vapply(DECADE_BINS, function(b) sprintf("%d-%d", b[1], b[2]), character(1)), collapse = ", "), "."
)

rect_decade_long <- build_rect_decade_summary(couce_year, "hours_total", DECADE_BINS)
decade_cols <- levels(rect_decade_long$decade)
rect_decade_wide <- pivot_rect_decade_wide(rect_decade_long, "mean_val")
write_csv(rect_decade_wide, path_out_B1)
logmsg(
  "B1 saved: ", path_out_B1, " (", nrow(rect_decade_wide), " rectangles x ", length(decade_cols),
  " decade columns [", paste(decade_cols, collapse = ", "), "]; NA = no Couce record for that rectangle in that decade)."
)

b2_comparisons <- list(
  list(from = decade_cols[2], to = decade_cols[1], label = "decade2_minus_decade1"),
  list(from = decade_cols[3], to = decade_cols[2], label = "decade3_minus_decade2"),
  list(from = decade_cols[3], to = decade_cols[1], label = "decade3_minus_decade1")
)
b2_list <- list(
  rect_decade_change(rect_decade_wide, decade_cols[1], decade_cols[2], "decade2_minus_decade1"),
  rect_decade_change(rect_decade_wide, decade_cols[2], decade_cols[3], "decade3_minus_decade2"),
  rect_decade_change(rect_decade_wide, decade_cols[1], decade_cols[3], "decade3_minus_decade1")
)
b2_table <- bind_rows(b2_list)
write_csv(b2_table, path_out_B2)
for (i in seq_along(b2_list)) {
  logmsg(sprintf("B2 comparison '%s': n_rectangles_used = %d (both decades non-NA).", b2_comparisons[[i]]$label, attr(b2_list[[i]], "n_used")))
}
logmsg("B2 saved: ", path_out_B2, " (", nrow(b2_table), " rows across 3 decade-pair comparisons, absolute + percent change).")

# --- B3: tercile classification on the full-period change (decade3 - decade1) ---
b2_full_period <- b2_list[[3]] # decade3_minus_decade1
b3_table <- classify_change_tercile(b2_full_period, "abs_change")
tercile_low <- attr(b3_table, "tercile_boundary_low")
tercile_high <- attr(b3_table, "tercile_boundary_high")
write_csv(b3_table, path_out_B3)
logmsg(
  "B3 saved: ", path_out_B3, " — tercile classification of |decade3-decade1 absolute change| via dplyr::ntile(., 3) ",
  "on ", nrow(b3_table), " rectangles with both decades present. PROVISIONAL exact boundary values: 'stable' <= ",
  sprintf("%.2f", tercile_low), " hours; 'large_change' >= ", sprintf("%.2f", tercile_high),
  " hours (bottom/top tercile of |absolute change|, NOT a formal cutoff)."
)
logmsg(sprintf(
  "  table(change_class): stable=%d, moderate=%d, large_change=%d",
  sum(b3_table$change_class == "stable"), sum(b3_table$change_class == "moderate"), sum(b3_table$change_class == "large_change")
))

# --- B4: spatial persistence correlations ---
b4_table <- bind_rows(
  decade_persistence_correlation(rect_decade_wide, decade_cols[1], decade_cols[2], "decade1_vs_decade2"),
  decade_persistence_correlation(rect_decade_wide, decade_cols[2], decade_cols[3], "decade2_vs_decade3"),
  decade_persistence_correlation(rect_decade_wide, decade_cols[1], decade_cols[3], "decade1_vs_decade3")
)
write_csv(b4_table, path_out_B4)
logmsg("B4 saved: ", path_out_B4, " — spatial persistence (Pearson + Spearman correlation of rectangle-level fishing pressure across decade pairs):")
for (i in seq_len(nrow(b4_table))) {
  r <- b4_table[i, ]
  logmsg(sprintf("  - %s: n=%d, Pearson r=%.3f, Spearman rho=%.3f", r$comparison, r$n_rectangles, r$pearson_r, r$spearman_rho))
}

# ===========================================================================
# Section C — ICC by spatial-unit definition (5 resolutions)
# ===========================================================================
logmsg("")
logmsg("## Section C — ICC by spatial-unit definition (5 resolutions)")

relaxed_rects <- rect_flags %>%
  filter(has_couce_coverage, n_years_present >= MIN_YEARS_RELAXED) %>%
  pull(stat_rec)
n_has_couce <- sum(rect_flags$has_couce_coverage)
logmsg(
  "Relaxed inclusion rule for individual-rectangle resolution (#1): has_couce_coverage AND n_years_present >= ",
  MIN_YEARS_RELAXED, " (n_years_present = years with >=1 haul, ANY count — NOT the stricter usable_for_fishing_analysis ",
  "flag, which required >=10 years with >=5 hauls/year and returned 0 rectangles; see h3_pre_exploration_run_log.md). ",
  length(relaxed_rects), " of ", nrow(rect_flags), " rectangles qualify under this relaxed rule."
)
if (length(relaxed_rects) == n_has_couce) {
  logmsg(
    "This relaxed rule turns out to select EXACTLY the same set as has_couce_coverage alone (", n_has_couce,
    " of ", nrow(rect_flags), ") — i.e. essentially every Couce-covered rectangle already has >=2 years of Q1 ",
    "haul data (median n_years_present across all 197 rectangles = ", median(rect_flags$n_years_present),
    " of a possible 31), so the >=2-year bar adds no further restriction beyond Couce coverage itself. Reported ",
    "plainly, not treated as a coincidence requiring adjustment."
  )
}

fishing_df_C <- couce_year %>% filter(stat_rec %in% relaxed_rects)
resid_df_C <- haul_eeos %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX, stat_rec %in% relaxed_rects) %>%
  group_by(stat_rec, year) %>%
  summarise(n_hauls = dplyr::n(), mean_signed = mean(resid_signed, na.rm = TRUE), mean_mag = mean(resid_magnitude, na.rm = TRUE), .groups = "drop")
logmsg(
  "Input data for ALL FIVE resolutions below is the SAME: ", nrow(fishing_df_C), " rectangle-year fishing-pressure ",
  "rows and ", nrow(resid_df_C), " rectangle-year residual-mean rows (", dplyr::n_distinct(resid_df_C$stat_rec),
  " rectangles), both restricted to the ", length(relaxed_rects), "-rectangle relaxed universe. NOTE: unlike the ",
  "pre-H3 script's Section E, residual rectangle-year means here use ALL years with >=1 haul (no SPARSE_HAUL_THRESHOLD ",
  ">=5-hauls/year gate) — a deliberate further relaxation, consistent with the brief's relaxed-rule instruction. Only ",
  "the GROUPING (unit_id) changes across the five rows below — this isolates the effect of spatial coarsening from ",
  "any change in the underlying data."
)

universe_A_geo <- intersect(universe_197, grid_index$stat_rec)
universe_B_geo <- intersect(universe_165, grid_index$stat_rec)
logmsg(
  "Scheme A spatial units are constructed over the full ", length(universe_A_geo), "-rectangle haul-bearing universe ",
  "(blind to fishing pressure, matching run_h3_policy_zones.R); Scheme B over the ", length(universe_B_geo),
  "-rectangle Couce-covered subset (needs a pressure value to assign tiers). Both then joined against the SAME ",
  length(relaxed_rects), "-rectangle relaxed input data above — a block/zone that happens to span rectangles ",
  "outside the relaxed universe simply contributes fewer observations, not zero (reported via n_units_with_*_data)."
)

unit_map_individual <- tibble(stat_rec = relaxed_rects, unit_id = relaxed_rects)

blocks_2x2 <- build_scheme_a_blocks(grid_index, universe_A_geo, 2L)
blocks_3x3 <- build_scheme_a_blocks(grid_index, universe_A_geo, 3L)

adjacency_B <- build_rook_adjacency(grid_index, universe_B_geo)
pressure_B <- couce_rect %>% filter(stat_rec %in% universe_B_geo)
tiers_3 <- assign_pressure_tiers(pressure_B, 3L)
zones_3 <- build_contiguous_zones(tiers_3$stat_rec, tiers_3$tier, adjacency_B)
unit_map_b3 <- zones_3 %>% select(stat_rec, unit_id)
tiers_4 <- assign_pressure_tiers(pressure_B, 4L)
zones_4 <- build_contiguous_zones(tiers_4$stat_rec, tiers_4$tier, adjacency_B)
unit_map_b4 <- zones_4 %>% select(stat_rec, unit_id)

c_table <- bind_rows(
  icc_comparison_row("1. Individual ICES rectangle (relaxed rule)", 1L, length(relaxed_rects), unit_map_individual, fishing_df_C, resid_df_C, RELIABILITY_TARGETS),
  icc_comparison_row("2. Scheme A, 2x2 block", 2L, dplyr::n_distinct(blocks_2x2$unit_id), blocks_2x2, fishing_df_C, resid_df_C, RELIABILITY_TARGETS),
  icc_comparison_row("3. Scheme A, 3x3 block", 3L, dplyr::n_distinct(blocks_3x3$unit_id), blocks_3x3, fishing_df_C, resid_df_C, RELIABILITY_TARGETS),
  icc_comparison_row("4. Scheme B, 3-tier zone", 4L, dplyr::n_distinct(unit_map_b3$unit_id), unit_map_b3, fishing_df_C, resid_df_C, RELIABILITY_TARGETS),
  icc_comparison_row("5. Scheme B, 4-tier zone", 5L, dplyr::n_distinct(unit_map_b4$unit_id), unit_map_b4, fishing_df_C, resid_df_C, RELIABILITY_TARGETS)
) %>% arrange(resolution_order)
write_csv(c_table, path_out_C)
logmsg("")
logmsg("C saved: ", path_out_C, " — ordered finest to coarsest:")
for (i in seq_len(nrow(c_table))) {
  r <- c_table[i, ]
  logmsg(sprintf(
    "  - %s: n_units_constructed=%d, n_units_with_data(fish/resid)=%d/%d | ICC fishing=%s, resid_mag=%s, resid_signed=%s",
    r$resolution, r$n_units_constructed, r$n_units_with_fishing_data, r$n_units_with_residual_data,
    ifelse(is.na(r$icc_fishing_pressure), "NA", sprintf("%.4f", r$icc_fishing_pressure)),
    ifelse(is.na(r$icc_residual_magnitude), "NA", sprintf("%.4f", r$icc_residual_magnitude)),
    ifelse(is.na(r$icc_residual_signed), "NA", sprintf("%.4f", r$icc_residual_signed))
  ))
}

# ===========================================================================
# Section D — haul-effort unevenness and its relationship to fishing pressure
# ===========================================================================
logmsg("")
logmsg("## Section D — haul-effort unevenness vs fishing pressure")
logmsg(
  "Reuses the SAME Scheme A (2x2/3x3, built over the ", length(universe_A_geo), "-rectangle haul-bearing universe) ",
  "and Scheme B (3-/4-tier, built over the ", length(universe_B_geo), "-rectangle Couce-covered universe) spatial ",
  "units as Section C above (identical construction to run_h3_policy_zones.R), with the same POLICY_BREAK_YEAR = ",
  POLICY_BREAK_YEAR, " pre/post split reused as-is from that script."
)

d_combos <- list(
  list(scheme = "Scheme A (block merge)", param = "2x2", unit_map = blocks_2x2, id_prefix = "A_2x2"),
  list(scheme = "Scheme A (block merge)", param = "3x3", unit_map = blocks_3x3, id_prefix = "A_3x3"),
  list(scheme = "Scheme B (pressure zones)", param = "3-tier", unit_map = unit_map_b3, id_prefix = "B_3tier"),
  list(scheme = "Scheme B (pressure zones)", param = "4-tier", unit_map = unit_map_b4, id_prefix = "B_4tier")
)

d_rows <- list()
scatter_dfs <- list()
for (combo in d_combos) {
  unit_period_hauls <- summarise_unit_period_hauls(haul_full, combo$unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)
  wide_hauls <- pivot_unit_period_hauls_wide(unit_period_hauls)
  unit_period_fish <- summarise_unit_period_fishing(couce_year, combo$unit_map, POLICY_BREAK_YEAR, H2_YEAR_MIN, H2_YEAR_MAX)

  for (per in c("pre", "post")) {
    fish_per <- unit_period_fish %>%
      filter(period == per) %>%
      select(unit_id, mean_fishing_hours) %>%
      inner_join(wide_hauls %>% select(unit_id, n_hauls = all_of(per)), by = "unit_id")

    row <- haul_unevenness_row(combo$scheme, combo$param, per, wide_hauls[[per]], fish_per)
    d_rows[[length(d_rows) + 1L]] <- row

    scatter_dfs[[length(scatter_dfs) + 1L]] <- fish_per %>%
      mutate(scheme = combo$scheme, parameter = combo$param, period = per)
  }
}
d_table <- bind_rows(d_rows)
write_csv(d_table, path_out_D)
logmsg("D saved: ", path_out_D, ":")
for (i in seq_len(nrow(d_table))) {
  r <- d_table[i, ]
  logmsg(sprintf(
    "  - %s / %s / %s: n_units=%d, CV=%.1f%%, Gini=%.3f, max/min ratio (strict/excl-zero)=%s/%s, cor(fishing,haulcount) Pearson=%.3f Spearman=%.3f (n=%d)",
    r$scheme, r$parameter, r$period, r$n_units, r$cv_pct_haul_count, r$gini_haul_count,
    ifelse(is.infinite(r$max_min_ratio_strict), "Inf", sprintf("%.1f", r$max_min_ratio_strict)),
    ifelse(is.na(r$max_min_ratio_excl_zero_units), "NA", sprintf("%.1f", r$max_min_ratio_excl_zero_units)),
    r$pearson_r_fishing_vs_haulcount, r$spearman_rho_fishing_vs_haulcount, r$n_units_used_for_correlation
  ))
}

# --- D.4 (optional byproduct): scatter figures, zone haul count vs mean fishing pressure, faceted pre/post ---
scatter_all <- bind_rows(scatter_dfs)
fig_paths_D <- character(0)
for (combo in d_combos) {
  df <- scatter_all %>% filter(scheme == combo$scheme, parameter == combo$param)
  df$period <- factor(df$period, levels = c("pre", "post"), labels = c(sprintf("Pre (%d-%d)", H2_YEAR_MIN, POLICY_BREAK_YEAR - 1L), sprintf("Post (%d-%d)", POLICY_BREAK_YEAR, H2_YEAR_MAX)))
  p <- ggplot(df, aes(x = mean_fishing_hours, y = n_hauls)) +
    geom_point(alpha = 0.7, colour = "#4575b4", size = 2) +
    facet_wrap(~period, scales = "free") +
    labs(
      x = "Zone mean Couce fishing hours", y = "Zone haul count",
      title = sprintf("Zone haul count vs zone mean fishing pressure — %s / %s", combo$scheme, combo$param),
      caption = "One point per spatial unit; descriptive only, no fitted line, no significance test."
    ) +
    theme_minimal(base_size = 11)
  path_fig <- file.path(fig_dir, sprintf("h2h3_designD4_%s_haulcount_vs_pressure.png", combo$id_prefix))
  ggsave(path_fig, p, width = 9, height = 5, dpi = 150)
  fig_paths_D <- c(fig_paths_D, path_fig)
}
logmsg("D4 (optional byproduct) figures saved: ", paste(basename(fig_paths_D), collapse = ", "))

# ===========================================================================
# Outputs index
# ===========================================================================
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_A1)
logmsg("- ", path_out_A2)
logmsg("- ", path_out_A3)
logmsg("- ", path_out_B1)
logmsg("- ", path_out_B2)
logmsg("- ", path_out_B3)
logmsg("- ", path_out_B4)
logmsg("- ", path_out_C)
logmsg("- ", path_out_D)
logmsg("- ", paste(fig_paths_D, collapse = "\n  - "), " (optional D4 byproduct figures)")
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== H2/H3 design-support quantitative extraction complete — numbers only, no design recommendation made. ===\n")
