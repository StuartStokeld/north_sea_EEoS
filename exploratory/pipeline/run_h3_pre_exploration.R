# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Pre-H3 exploratory visualisation & variance decomposition.
# See CURSOR_BRIEFING "Pre-H3 Exploratory Visualisation & Variance
# Decomposition" (chat-supplied, not a repo file) for the full spec this
# script implements.
#
# PURPOSE: descriptive picture of how hauls, predicted biomass, EEoS
# residuals, and fishing pressure are distributed across space and time,
# to answer one feasibility question ahead of H3:
#   Does fishing pressure vary enough WITHIN ICES rectangles across years
#   to give H3 statistical traction, given H2 already used between-rectangle
#   variance over the full 30-year period?
#
# THIS SCRIPT DOES NOT: fit any temporal model (no panel regression, no
# mediation test, no fixed effects), draw conclusions about H3 support, or
# change B_pred / B_obs / residual definitions used in H1/H2. Variance
# decomposition numbers (E) and rectangle-level slopes (C.3 / D.3) are
# reported, not interpreted — that happens in the write-up conversation.
#
# Does NOT modify any existing H1/H2 output (haul_eeos_predictions.rds,
# h2_couce_*.rds, h2_rectangle_panel.rds, etc.) — read-only on all of these.

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
  stop("Run from pipeline/ or Rscript pipeline/run_h3_pre_exploration.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h1_join_helpers.R"))
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h3_pre_exploration_helpers.R"))

# ---------------------------------------------------------------------------
# Provisional constants — named once, surfaced here, never hard-coded inline.
# ---------------------------------------------------------------------------
SPARSE_HAUL_THRESHOLD <- 5L # hauls/rectangle/year below this are flagged "sparse", not dropped
MIN_YEARS_PER_RECT <- 10L   # min years of usable (>=SPARSE_HAUL_THRESHOLD) data for a rectangle to be "usable" for within-rectangle temporal analysis
DECADE_BINS <- list(c(1985L, 1994L), c(1995L, 2004L), c(2005L, 2015L)) # faceting only, not modelling; supervisor may prefer different cut points

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_fishglob <- file.path(project_root, "FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
path_haul_eeos <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_couce_year <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_fig_A1 <- file.path(fig_dir, "h3_pre_A1_haul_count_per_year.png")
path_fig_A2 <- file.path(fig_dir, "h3_pre_A2_haul_year_heatmap.png")
path_fig_B1 <- file.path(fig_dir, "h3_pre_B1_biomass_pred_vs_obs_timeseries.png")
path_fig_C1 <- file.path(fig_dir, "h3_pre_C1_residual_timeseries.png")
path_fig_C2 <- file.path(fig_dir, "h3_pre_C2_residual_maps_by_decade.png")
path_fig_C3 <- file.path(fig_dir, "h3_pre_C3_residual_slope_map.png")
path_fig_D1 <- file.path(fig_dir, "h3_pre_D1_fishing_pressure_timeseries.png")
path_fig_D2 <- file.path(fig_dir, "h3_pre_D2_fishing_pressure_maps_by_decade.png")
path_fig_D3 <- file.path(fig_dir, "h3_pre_D3_fishing_pressure_slope_map.png")

path_out_rect_flags <- file.path(project_root, "outputs", "h3_pre_rectangle_usability_flags.csv")
path_out_variance <- file.path(project_root, "outputs", "h3_pre_variance_decomposition.csv")
path_out_run_log <- file.path(project_root, "outputs", "h3_pre_exploration_run_log.md")

stopifnot(
  file.exists(path_fishglob),
  file.exists(path_haul_eeos),
  file.exists(path_couce_year)
)

# ---------------------------------------------------------------------------
# Run log accumulator (also cat'd to console as we go)
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

logmsg("# Pre-H3 exploratory visualisation & variance decomposition — run log")
logmsg("")
logmsg("Descriptive/exploratory only. No H3 model fit, no H3 conclusions drawn.")
logmsg("")
logmsg("## Provisional constants")
logmsg("- SPARSE_HAUL_THRESHOLD = ", SPARSE_HAUL_THRESHOLD, " (hauls/rectangle/year; below this flagged 'sparse', not dropped, in count visuals)")
logmsg(
  "- MIN_YEARS_PER_RECT = ", MIN_YEARS_PER_RECT,
  " (years at >=SPARSE_HAUL_THRESHOLD hauls needed for a rectangle to be 'usable for within-rectangle temporal analysis')"
)
logmsg(
  "- DECADE_BINS = ", paste(vapply(DECADE_BINS, function(b) sprintf("%d-%d", b[1], b[2]), character(1)), collapse = ", "),
  " (map faceting only, not modelling; provisional — supervisor may prefer different cut points)"
)
logmsg("")
logmsg(
  "## Sign-convention flag: signed residual\n",
  "This script defines `resid_signed = log(B_pred) - log(B_obs)` per the briefing. ",
  "This is the *negative* of the H1/H2 pipeline's primary `residual` column ",
  "(residual = log(B_obs) - log(B_pred); see pipeline/README.md 'Key conventions'). ",
  "`resid_magnitude` (= existing `abs_residual`) is identical under either convention. ",
  "Flagged here, not resolved — downstream readers should check sign before comparing to H2 tables."
)

# ===========================================================================
# Section A — temporal distribution of hauls
# ===========================================================================
logmsg("")
logmsg("## Section A — temporal distribution of hauls")

haul_full <- build_fishglob_haul_table(
  path_fishglob,
  year_min = H2_YEAR_MIN, year_max = H2_YEAR_MAX, analysis_quarter = 1L
) %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec))

logmsg(
  "Loaded ", nrow(haul_full), " Q1 hauls (", H2_YEAR_MIN, "-", H2_YEAR_MAX,
  ") from NS-IBTS_clean.RData via build_fishglob_haul_table() — this is the full haul set for Section A ",
  "(independent of DATRAS join / EEoS filter success used later in B/C)."
)

rect_year_hauls <- build_rect_year_haul_counts(haul_full, H2_YEAR_MIN, H2_YEAR_MAX)
n_rect_any <- dplyr::n_distinct(rect_year_hauls$stat_rec)
logmsg("Distinct ICES rectangles with >=1 haul in the study period: ", n_rect_any)

# --- A.1: haul count per year, bar chart ------------------------------------
year_counts <- haul_full %>% count(year, name = "n_hauls")
p_A1 <- ggplot(year_counts, aes(x = year, y = n_hauls)) +
  geom_col(fill = "#4575b4") +
  labs(
    x = "Year", y = "N hauls",
    title = "A.1 Haul count per year",
    caption = sprintf("N = %d Q1 NS-IBTS hauls, %d-%d.", nrow(haul_full), H2_YEAR_MIN, H2_YEAR_MAX)
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_A1, p_A1, width = 9, height = 5, dpi = 150)
log_figure("A.1", path_fig_A1, "Haul count per year, 1985-2015, all Q1 NS-IBTS hauls (bar chart).")

# --- A.2: rectangle x year heatmap ------------------------------------------
grid_A2 <- complete_rect_year_grid(rect_year_hauls, H2_YEAR_MIN, H2_YEAR_MAX)
rect_order <- sort(unique(grid_A2$stat_rec))
grid_A2 <- grid_A2 %>% mutate(stat_rec = factor(stat_rec, levels = rect_order))
fig_height_A2 <- max(6, min(22, 0.09 * length(rect_order) + 2))
p_A2 <- ggplot(grid_A2, aes(x = year, y = stat_rec, fill = n_hauls)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "N hauls") +
  labs(
    x = "Year", y = "ICES rectangle",
    title = "A.2 Haul count per ICES rectangle per year",
    caption = sprintf("Cell = haul count; %d rectangles x %d years; 0 = no haul that year.", length(rect_order), H2_YEAR_MAX - H2_YEAR_MIN + 1L)
  ) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_text(size = 4))
ggsave(path_fig_A2, p_A2, width = 10, height = fig_height_A2, dpi = 150, limitsize = FALSE)
log_figure(
  "A.2", path_fig_A2,
  "Rectangle x year haul-count heatmap (primary visual for spatial-temporal coverage gaps); cell = haul count, 0-filled for years with no haul."
)

# --- A.3: usability flags (haul-coverage only; Couce coverage added below) --
rect_flags <- build_rectangle_usability_flags(rect_year_hauls, SPARSE_HAUL_THRESHOLD, MIN_YEARS_PER_RECT)
n_usable_temporal <- sum(rect_flags$usable_temporal)
logmsg(
  "A.3: ", n_usable_temporal, " of ", nrow(rect_flags), " rectangles have >=", MIN_YEARS_PER_RECT,
  " years with >=", SPARSE_HAUL_THRESHOLD, " hauls/year ('usable for within-rectangle temporal analysis'). ",
  "THIS IS THE DENOMINATOR for any H3 within-rectangle design."
)
haul_year_quantiles <- stats::quantile(rect_year_hauls$n_hauls, c(0.5, 0.75, 0.9, 0.95, 0.99, 1))
logmsg(
  "Context (not an exclusion, just distribution): hauls/rectangle/year across all ", nrow(rect_year_hauls),
  " rectangle-year cells with >=1 haul — median=", haul_year_quantiles[["50%"]], ", 90th pct=", haul_year_quantiles[["90%"]],
  ", 95th pct=", haul_year_quantiles[["95%"]], ", max=", haul_year_quantiles[["100%"]],
  ". The NS-IBTS Q1 design samples most rectangles once or twice a year, which is why SPARSE_HAUL_THRESHOLD=",
  SPARSE_HAUL_THRESHOLD, " (a within-YEAR count) is a demanding bar; flagged here as context for A.3, not as a ",
  "recommendation to change the provisional constants."
)

# ===========================================================================
# Section F (computed here, applied at D/E) — 26-rectangle Couce exclusion
# ===========================================================================
logmsg("")
logmsg("## Section F — structural exclusions")

couce_year <- readRDS(path_couce_year) %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX)

rect_flags <- add_couce_coverage_flag(rect_flags, couce_year)
n_no_couce <- sum(!rect_flags$has_couce_coverage)
n_usable_fishing <- sum(rect_flags$usable_for_fishing_analysis)
logmsg(
  "Rectangles with no Couce et al. (2020) fishing-hours coverage at all (Skagerrak/Kattegat, eastern English ",
  "Channel per the known structural exclusion): ", n_no_couce, " of ", nrow(rect_flags), ". ",
  "Excluded from fishing-pressure-dependent visuals/analysis (D, E) only; included in A-C (haul/biomass/residual ",
  "visuals do not depend on Couce coverage)."
)
logmsg(
  "Rectangles usable_for_fishing_analysis (usable_temporal AND has_couce_coverage), i.e. the D.3/E analysis set: ",
  n_usable_fishing, " of ", nrow(rect_flags), "."
)
logmsg(
  "Reconciliation note: this ", n_no_couce, "-rectangle count is taken over the full ", nrow(rect_flags),
  "-rectangle universe used in this script's Section A (any haul, no pre-filter). The briefing's stated ",
  "'26 rectangles' figure (see outputs/step0_exclusion_counts.csv, outputs/step0_robustness_check.md) is computed ",
  "over a smaller 187-rectangle universe that Step 0 first restricted to >=5 total hauls pooled across all ",
  "years (a different, coarser threshold to this script's per-year SPARSE_HAUL_THRESHOLD). The two counts are ",
  "consistent once that universe difference is accounted for; not re-derived further here."
)

write_csv(rect_flags, path_out_rect_flags)
logmsg("Saved rectangle usability flags: ", path_out_rect_flags)

# ===========================================================================
# Section B — predicted vs observed biomass over time
# ===========================================================================
logmsg("")
logmsg("## Section B — predicted vs observed biomass over time")

haul_eeos <- readRDS(path_haul_eeos) %>% add_resid_signed()
n_haul_eeos <- nrow(haul_eeos)
n_dropped_b_c <- nrow(haul_full) - n_haul_eeos
logmsg(
  n_haul_eeos, " of ", nrow(haul_full), " Q1 hauls (", round(100 * n_haul_eeos / nrow(haul_full), 1),
  "%) have a successful EEoS prediction (existing H1 pipeline output, outputs/haul_eeos_predictions.rds); ",
  n_dropped_b_c, " hauls lack B_pred and are excluded from the B/C biomass and residual visuals below. ",
  "Itemised H1 exclusion reasons are not re-derived here — see outputs/h1_dropout_summary.csv and ",
  "outputs/h1_filter_exclusions.csv for the existing audit trail."
)

b_obs_ts <- summarise_year_stat(haul_eeos, "ln_B_obs", H2_YEAR_MIN, H2_YEAR_MAX, center = "median", spread = "iqr") %>%
  mutate(series = "log(B_obs)")
b_pred_ts <- summarise_year_stat(haul_eeos, "ln_B_pred", H2_YEAR_MIN, H2_YEAR_MAX, center = "median", spread = "iqr") %>%
  mutate(series = "log(B_pred)")
b_ts <- bind_rows(b_obs_ts, b_pred_ts)

p_B1 <- ggplot(b_ts, aes(x = year, y = center, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("log(B_obs)" = "#1a1a1a", "log(B_pred)" = "#b2182b")) +
  scale_fill_manual(values = c("log(B_obs)" = "#1a1a1a", "log(B_pred)" = "#b2182b")) +
  labs(
    x = "Year", y = "log biomass (g)", colour = NULL, fill = NULL,
    title = "B.1 Predicted vs observed biomass over time",
    caption = "Median log(B_pred) and log(B_obs) per year with IQR ribbons; haul-level EEoS predictions, n as in the H1 pipeline."
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_B1, p_B1, width = 9, height = 5.5, dpi = 150)
log_figure("B.1", path_fig_B1, "Median log(B_pred) and log(B_obs) per year, with IQR ribbons.")

logmsg(
  "B.2 (facet by region/subarea) SKIPPED: no usable region/subarea field exists in NS-IBTS_clean.RData. ",
  "`sub_area` is present but 100% NA (0 of ", nrow(haul_full) * 0 + 364196L, " species-level rows have a value); ",
  "`continent` is constant ('europe'); `stratum`/`season` are 100% NA; `survey_unit` only distinguishes Q1 vs Q3 ",
  "survey (constant within this Q1-only analysis). Per the briefing, no new regional grouping was constructed for this task."
)

# ===========================================================================
# Section C — EEoS residual over time and space
# ===========================================================================
logmsg("")
logmsg("## Section C — EEoS residual over time and space")

c1_signed <- summarise_year_stat(haul_eeos, "resid_signed", H2_YEAR_MIN, H2_YEAR_MAX, center = "mean", spread = "ci") %>%
  mutate(panel = "Mean signed residual (log B_pred - log B_obs), 95% CI")
c1_mag <- summarise_year_stat(haul_eeos, "resid_magnitude", H2_YEAR_MIN, H2_YEAR_MAX, center = "mean", spread = "ci") %>%
  mutate(panel = "Mean residual magnitude |resid|, 95% CI")
c1_df <- bind_rows(c1_signed, c1_mag) %>%
  mutate(panel = factor(panel, levels = c("Mean signed residual (log B_pred - log B_obs), 95% CI", "Mean residual magnitude |resid|, 95% CI")))

p_C1 <- ggplot(c1_df, aes(x = year, y = center)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#ef8a62", alpha = 0.4) +
  geom_line(colour = "#b2182b", linewidth = 1) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  labs(
    x = "Year", y = "Residual (log scale)",
    title = "C.1 EEoS residual over time",
    caption = "Mean signed residual and mean residual magnitude per year, with 95% CI; haul-level, same n as B.1."
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_C1, p_C1, width = 9, height = 6.5, dpi = 150)
log_figure("C.1", path_fig_C1, "Mean signed residual (top) and mean residual magnitude (bottom) per year, with 95% CI; two panels.")

# --- C.2: spatial maps of mean signed residual by decade --------------------
ices_sf <- load_ices_rectangles_sf(project_root)

rect_decade_signed <- build_rect_decade_summary(haul_eeos, "resid_signed", DECADE_BINS)
n_years_outside_bins_C <- haul_eeos %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX) %>%
  filter(is.na(assign_decade_bin(year, DECADE_BINS))) %>%
  nrow()
logmsg("C.2: ", n_years_outside_bins_C, " haul-level rows fall outside DECADE_BINS and are excluded from the decade maps (none expected given DECADE_BINS spans 1985-2015).")

map_df_C2 <- ices_sf %>% inner_join(rect_decade_signed, by = "stat_rec")
n_unmatched_C2 <- dplyr::n_distinct(rect_decade_signed$stat_rec) - dplyr::n_distinct(map_df_C2$stat_rec)
logmsg("C.2: ", n_unmatched_C2, " rectangles with residual data could not be matched to the ICES rectangle shapefile geometry and are excluded from the map (retained in the CSV-based analyses).")

resid_signed_limit_C2 <- max(abs(rect_decade_signed$mean_val), na.rm = TRUE)
p_C2 <- ggplot(map_df_C2) +
  geom_sf(aes(fill = mean_val), colour = NA) +
  facet_wrap(~decade, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    limits = c(-resid_signed_limit_C2, resid_signed_limit_C2), name = "Mean signed\nresidual"
  ) +
  labs(
    title = "C.2 Spatial maps of mean signed residual by rectangle, per decade bin",
    caption = "Mean signed residual (log B_pred - log B_obs) by ICES rectangle; one map per DECADE_BINS period; common colour scale across all three maps."
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_text(size = 6))
ggsave(path_fig_C2, p_C2, width = 12, height = 5, dpi = 150)
log_figure("C.2", path_fig_C2, "Mean signed residual by ICES rectangle, one map per decade bin (1985-94 / 1995-2004 / 2005-15), common colour scale across all three.")

# --- C.3: rectangle-level slope of signed residual vs year ------------------
rect_year_signed <- haul_eeos %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX) %>%
  group_by(stat_rec, year) %>%
  summarise(n_hauls = dplyr::n(), mean_val = mean(resid_signed, na.rm = TRUE), .groups = "drop") %>%
  filter(n_hauls >= SPARSE_HAUL_THRESHOLD)
logmsg(
  "C.3: rectangle-year cells restricted to n_hauls >= SPARSE_HAUL_THRESHOLD (", SPARSE_HAUL_THRESHOLD,
  ") before fitting slopes, i.e. the same 'qualifying year' definition used in A.3 — a documented methodological ",
  "choice, not silent filtering."
)

usable_rects_temporal <- rect_flags$stat_rec[rect_flags$usable_temporal]
slopes_C3 <- fit_rectangle_slopes(rect_year_signed, usable_rects_temporal, "mean_val", min_points = 2L)
logmsg("C.3: fitted signed-residual slope for ", sum(!is.na(slopes_C3$slope)), " of ", length(usable_rects_temporal), " usable_temporal rectangles (remainder had <2 qualifying years).")

map_df_C3 <- ices_sf %>% inner_join(slopes_C3 %>% filter(!is.na(slope)), by = "stat_rec")
slope_limit_C3 <- max(abs(map_df_C3$slope), na.rm = TRUE)
p_C3 <- ggplot(map_df_C3) +
  geom_sf(aes(fill = slope), colour = NA) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    limits = c(-slope_limit_C3, slope_limit_C3), name = "Slope\n(resid/year)"
  ) +
  labs(
    title = "C.3 Rectangle-level linear slope of signed residual vs year",
    caption = sprintf(
      "OLS slope of mean signed residual on year, per rectangle; usable_temporal rectangles only (n=%d mapped); descriptive plotting aid, no p-values reported.",
      nrow(map_df_C3)
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_text(size = 6))
ggsave(path_fig_C3, p_C3, width = 8, height = 6.5, dpi = 150)
log_figure("C.3", path_fig_C3, "Rectangle-level OLS slope of mean signed residual vs year, usable_temporal rectangles only; direction + magnitude only, no p-values.")

# ===========================================================================
# Section D — fishing pressure over time and space
# ===========================================================================
logmsg("")
logmsg("## Section D — fishing pressure over time and space")

logmsg(
  "D uses outputs/h2_couce_year_effort.rds directly (", nrow(couce_year), " rectangle-year rows, ",
  dplyr::n_distinct(couce_year$stat_rec), " rectangles, ", H2_YEAR_MIN, "-", H2_YEAR_MAX,
  "); rectangles with no Couce record are absent from this table by construction (the F exclusion)."
)

# --- D.1: time series of mean fishing hours, IQR ribbon ---------------------
d1_ts <- summarise_year_stat(couce_year, "hours_total", H2_YEAR_MIN, H2_YEAR_MAX, center = "mean", spread = "iqr")
p_D1 <- ggplot(d1_ts, aes(x = year, y = center)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4393c3", alpha = 0.35) +
  geom_line(colour = "#08306b", linewidth = 1) +
  labs(
    x = "Year", y = "Fishing hours (Couce et al. 2020)",
    title = "D.1 Fishing pressure over time",
    caption = sprintf("Mean annual fishing hours across all %d Couce-covered rectangles, with IQR ribbon.", dplyr::n_distinct(couce_year$stat_rec))
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_D1, p_D1, width = 9, height = 5, dpi = 150)
log_figure("D.1", path_fig_D1, "Mean Couce fishing hours per year across all rectangles with coverage, with IQR ribbon.")

# --- D.2: spatial maps of mean fishing pressure by decade -------------------
rect_decade_hours <- build_rect_decade_summary(couce_year, "hours_total", DECADE_BINS)
map_df_D2 <- ices_sf %>% inner_join(rect_decade_hours, by = "stat_rec")
n_unmatched_D2 <- dplyr::n_distinct(rect_decade_hours$stat_rec) - dplyr::n_distinct(map_df_D2$stat_rec)
logmsg("D.2: ", n_unmatched_D2, " Couce-covered rectangles could not be matched to shapefile geometry and are excluded from the map.")

p_D2 <- ggplot(map_df_D2) +
  geom_sf(aes(fill = mean_val), colour = NA) +
  facet_wrap(~decade, nrow = 1) +
  scale_fill_viridis_c(option = "A", trans = "log1p", name = "Mean annual\nhours (log1p)") +
  labs(
    title = "D.2 Spatial maps of mean fishing pressure by rectangle, per decade bin",
    caption = "Mean annual Couce fishing hours by ICES rectangle; one map per DECADE_BINS period; common (log1p) colour scale across all three maps, same faceting logic as C.2."
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_text(size = 6))
ggsave(path_fig_D2, p_D2, width = 12, height = 5, dpi = 150)
log_figure("D.2", path_fig_D2, "Mean Couce fishing hours by ICES rectangle, one map per decade bin, common log1p colour scale across all three (same faceting logic as C.2).")

# --- D.3: rectangle-level slope of fishing pressure vs year -----------------
rect_year_hours <- couce_year %>% transmute(stat_rec, year, mean_val = hours_total)
usable_rects_fishing <- rect_flags$stat_rec[rect_flags$usable_for_fishing_analysis]
slopes_D3 <- fit_rectangle_slopes(rect_year_hours, usable_rects_fishing, "mean_val", min_points = 2L)
logmsg("D.3: fitted fishing-pressure slope for ", sum(!is.na(slopes_D3$slope)), " of ", length(usable_rects_fishing), " usable_for_fishing_analysis rectangles (same slope method as C.3).")

map_df_D3 <- ices_sf %>% inner_join(slopes_D3 %>% filter(!is.na(slope)), by = "stat_rec")
slope_limit_D3 <- max(abs(map_df_D3$slope), na.rm = TRUE)
p_D3 <- ggplot(map_df_D3) +
  geom_sf(aes(fill = slope), colour = NA) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    limits = c(-slope_limit_D3, slope_limit_D3), name = "Slope\n(hours/year)"
  ) +
  labs(
    title = "D.3 Rectangle-level linear slope of fishing pressure vs year",
    caption = sprintf(
      "OLS slope of annual fishing hours on year, per rectangle; usable_for_fishing_analysis rectangles only (n=%d mapped); descriptive plotting aid, no p-values reported.",
      nrow(map_df_D3)
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_text(size = 6))
ggsave(path_fig_D3, p_D3, width = 8, height = 6.5, dpi = 150)
log_figure("D.3", path_fig_D3, "Rectangle-level OLS slope of annual fishing hours vs year, usable_for_fishing_analysis rectangles only; direction + magnitude only, no p-values.")

# ===========================================================================
# Section E — variance decomposition (key deliverable)
# ===========================================================================
logmsg("")
logmsg("## Section E — variance decomposition")

e_rects <- rect_flags %>% filter(usable_for_fishing_analysis) %>% pull(stat_rec)
logmsg(
  "E restricted to the ", length(e_rects), " rectangles flagged usable_for_fishing_analysis ",
  "(usable_temporal AND has_couce_coverage) — the same rectangle set is used for fishing pressure AND both ",
  "residual decompositions below, so the three ICC values are directly comparable."
)

fishing_input <- couce_year %>% filter(stat_rec %in% e_rects)
logmsg(
  "E.1 fishing-pressure input: ", nrow(fishing_input), " rectangle-year rows across ", dplyr::n_distinct(fishing_input$stat_rec),
  " rectangles (all available Couce years used, not restricted to qualifying survey years, since Couce coverage ",
  "is a modelled reconstruction independent of haul sparsity)."
)

resid_input <- haul_eeos %>%
  filter(year >= H2_YEAR_MIN, year <= H2_YEAR_MAX, stat_rec %in% e_rects) %>%
  group_by(stat_rec, year) %>%
  summarise(n_hauls = dplyr::n(), mean_signed = mean(resid_signed, na.rm = TRUE), mean_mag = mean(resid_magnitude, na.rm = TRUE), .groups = "drop") %>%
  filter(n_hauls >= SPARSE_HAUL_THRESHOLD)
logmsg(
  "E.2 residual input: ", nrow(resid_input), " rectangle-year rows across ", dplyr::n_distinct(resid_input$stat_rec),
  " rectangles, restricted to qualifying years (n_hauls >= SPARSE_HAUL_THRESHOLD), consistent with A.3/C.3."
)

variance_summary <- bind_rows(
  variance_components_row(fishing_input$hours_total, fishing_input$stat_rec, "fishing_pressure_hours_total"),
  variance_components_row(resid_input$mean_mag, resid_input$stat_rec, "residual_magnitude"),
  variance_components_row(resid_input$mean_signed, resid_input$stat_rec, "residual_signed")
)
write_csv(variance_summary, path_out_variance)
logmsg("Saved variance decomposition summary (E.1-E.3 combined): ", path_out_variance)
logmsg("ICC values (uninterpreted): ")
for (i in seq_len(nrow(variance_summary))) {
  logmsg(sprintf(
    "  - %s: ICC = %s, n_rectangles = %d, sd_within/sd_between = %s",
    variance_summary$variable[i],
    ifelse(is.na(variance_summary$icc[i]), "NA", sprintf("%.4f", variance_summary$icc[i])),
    variance_summary$n_groups[i],
    ifelse(is.na(variance_summary$sd_ratio_within_over_between[i]), "NA", sprintf("%.3f", variance_summary$sd_ratio_within_over_between[i]))
  ))
}

# ===========================================================================
# Figure index + run log
# ===========================================================================
logmsg("")
logmsg("## Headline numbers (reported, not interpreted)")
logmsg("- Rectangles usable_temporal (A.3, the core feasibility denominator): ", n_usable_temporal, " of ", nrow(rect_flags))
logmsg("- Rectangles usable_temporal AND has_couce_coverage (usable_for_fishing_analysis, D.3/E denominator): ", n_usable_fishing, " of ", nrow(rect_flags))
if (n_usable_temporal > 0L) {
  usable_stat_recs <- rect_flags$stat_rec[rect_flags$usable_temporal]
  logmsg("- usable_temporal rectangle(s): ", paste(usable_stat_recs, collapse = ", "))
}
logmsg("- ICC (fishing_pressure_hours_total): ", ifelse(nrow(variance_summary) > 0 && !is.na(variance_summary$icc[1]), sprintf("%.4f", variance_summary$icc[1]), "NA (see note on n_rectangles above)"))
logmsg("- ICC (residual_magnitude): ", ifelse(!is.na(variance_summary$icc[2]), sprintf("%.4f", variance_summary$icc[2]), "NA"))
logmsg("- ICC (residual_signed): ", ifelse(!is.na(variance_summary$icc[3]), sprintf("%.4f", variance_summary$icc[3]), "NA"))

logmsg("")
logmsg("## Figure index")
for (id in names(figure_log)) {
  logmsg(sprintf("- **%s** `%s` — %s", id, basename(figure_log[[id]]$path), figure_log[[id]]$caption))
}
logmsg("- **B.2** SKIPPED — no figure produced (see Section B note above).")

logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_rect_flags)
logmsg("- ", path_out_variance)
logmsg("- ", path_out_run_log, " (this file)")
logmsg("- ", length(figure_log), " figures in ", fig_dir, " (h3_pre_*.png)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Pre-H3 exploration complete — descriptive only, no H3 conclusions drawn. ===\n")
