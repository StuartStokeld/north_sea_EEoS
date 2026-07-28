# Missing Visuals — Haul Map & Composition Figures.
# See CURSOR_BRIEFING "Missing Visuals — Haul Map & Composition Figures"
# (chat-supplied, not a repo file) for the full spec this script implements.
#
# PURPOSE: fill two gaps flagged by the exploratory-outputs collation
# (`outputs/exploratory_review/index.md`, Task 3): (1) a standalone
# geographic haul map (not a rectangle x year matrix), and (2) the plain
# composition/big-shoal distribution figures — histograms of S/D/size_CV,
# a plain D vs size_CV scatter, and a "big shoal" reference-share summary.
#
# THIS SCRIPT DOES NOT: recompute D, size_CV, S, N, B_obs, B_pred, or any
# other per-haul quantity (all read from existing outputs); define a new
# formal "big shoal" classification (falls back to a reference-only top
# decile of D, per the briefing, because no combined D x size_CV cutoff was
# found anywhere in the H1 dominance analysis); or produce the H2 residual
# map (explicitly deferred).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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
  stop("Run from pipeline/ or Rscript pipeline/run_missing_visuals.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "missing_visuals_helpers.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_rect_flags <- file.path(project_root, "outputs", "h3_pre_rectangle_usability_flags.csv")
path_dominance <- file.path(project_root, "outputs", "h1_dominance_haul_table.csv")
path_haul_eeos <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_D_bins <- file.path(project_root, "outputs", "h1_dominance_by_D_bins.csv")

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_fig_haul_map_1 <- file.path(fig_dir, "haul_map_1.png")
path_fig_composition_1 <- file.path(fig_dir, "composition_1.png")
path_fig_composition_2 <- file.path(fig_dir, "composition_2.png")
path_fig_composition_3 <- file.path(fig_dir, "composition_3.png")

path_out_run_log <- file.path(project_root, "outputs", "missing_visuals_run_log.md")

stopifnot(
  file.exists(path_rect_flags),
  file.exists(path_dominance),
  file.exists(path_haul_eeos),
  file.exists(path_D_bins)
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

logmsg("# Missing visuals — haul map & composition figures — run log")
logmsg("")
logmsg("Fills two gaps flagged in `outputs/exploratory_review/index.md` (Task 3, ",
       "\"core visual suite\"): a standalone geographic haul map, and the plain ",
       "S/D/size_CV composition distributions (as opposed to the residual-ratio ",
       "overlays already produced for Task 2). H2 residual map deliberately out of scope.")
logmsg("")

# ---------------------------------------------------------------------------
# Load data (read-only; no per-haul quantity is recomputed)
# ---------------------------------------------------------------------------
rect_flags <- read.csv(path_rect_flags, stringsAsFactors = FALSE)
dominance <- read.csv(path_dominance, stringsAsFactors = FALSE)
haul_eeos <- readRDS(path_haul_eeos)

logmsg("Loaded `h3_pre_rectangle_usability_flags.csv` (", nrow(rect_flags), " rectangles), ",
       "`h1_dominance_haul_table.csv` (", nrow(dominance), " hauls), ",
       "`haul_eeos_predictions.rds` (", nrow(haul_eeos), " hauls).")

comp <- build_composition_table(dominance, haul_eeos)
logmsg("Composition table: joined dominance table to `haul_eeos_predictions.rds` on `haul_id` ",
       "for `S` only; ", nrow(comp), " hauls, zero-loss join (verified 1:1 match, no recomputation).")

# ---------------------------------------------------------------------------
# haul_map_1: total haul count per ICES rectangle, full period, with the
# 32 no-Couce-coverage rectangles marked distinctly.
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Haul map")
logmsg("")
logmsg("Shapefile source: `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp`, loaded via ",
       "`load_ices_rectangles_sf()` (`pipeline/R/h2_spatial_helpers.R`), the same ICES StatRec grid ",
       "already used for the H2/H3 spatial maps (e.g. `h3_pre_C2`/`D2`, `h3_policy_*` zone maps). ",
       "No new shapefile was needed or added.")

ices_sf <- load_ices_rectangles_sf(project_root)
map_df_haul <- ices_sf %>% inner_join(rect_flags, by = "stat_rec")
n_unmatched_map <- dplyr::n_distinct(rect_flags$stat_rec) - dplyr::n_distinct(map_df_haul$stat_rec)
logmsg("haul_map_1: ", n_unmatched_map, " of ", nrow(rect_flags),
       " rectangles in the usability-flags table could not be matched to shapefile geometry ",
       "and are excluded from the map (total_hauls / coverage counts below use the full ", nrow(rect_flags), ").")

n_no_couce <- sum(!rect_flags$has_couce_coverage)
logmsg("Rectangles without Couce fishing-pressure coverage: ", n_no_couce, " of ", nrow(rect_flags),
       " (matches the confirmed universe from Task 1's run log).")

map_df_haul <- map_df_haul %>%
  mutate(couce_coverage = ifelse(has_couce_coverage, "Covered", paste0("Not covered (", n_no_couce, " rect.)")))

p_haul_map_1 <- ggplot(map_df_haul) +
  geom_sf(aes(fill = total_hauls, colour = couce_coverage, linewidth = couce_coverage)) +
  scale_fill_viridis_c(option = "C", name = "Total hauls\n(1985-2015)") +
  scale_colour_manual(values = setNames(c("grey50", "#d6604d"), c("Covered", paste0("Not covered (", n_no_couce, " rect.)"))),
                       name = "Couce fishing-pressure\ncoverage") +
  scale_linewidth_manual(values = setNames(c(0.15, 0.7), c("Covered", paste0("Not covered (", n_no_couce, " rect.)"))),
                          guide = "none") +
  labs(
    title = "haul_map_1: Total haul count per ICES rectangle, 1985-2015",
    caption = sprintf(
      "Choropleth of total Q1 NS-IBTS haul count per ICES rectangle, full period (%d rectangles mapped of %d in the usability-flags universe); rectangles outlined in red (n=%d) lack Couce et al. (2020) fishing-pressure coverage.",
      dplyr::n_distinct(map_df_haul$stat_rec), nrow(rect_flags), n_no_couce
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_text(size = 7), legend.key.size = unit(0.8, "lines"))
ggsave(path_fig_haul_map_1, p_haul_map_1, width = 9, height = 7.5, dpi = 150)
log_figure("haul_map_1", path_fig_haul_map_1,
           sprintf("Total haul count per ICES rectangle, 1985-2015 (choropleth, N=%d rectangles); the %d rectangles without Couce fishing-pressure coverage are outlined in red.",
                   nrow(rect_flags), n_no_couce))

# ---------------------------------------------------------------------------
# composition_1: histograms of S, D, size_CV
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Composition figures")
logmsg("")
logmsg("All values (S, D, size_CV) read as-is from `haul_eeos_predictions.rds` / ",
       "`h1_dominance_haul_table.csv`; none recomputed.")

comp_long <- comp %>%
  select(S, D, size_CV) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("S", "D", "size_CV"),
                          labels = c("Species richness (S)", "Berger-Parker dominance (D)", "Dominant-species size_CV")))

p_composition_1 <- ggplot(comp_long, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#4393c3", colour = "white", linewidth = 0.15) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  labs(
    x = NULL, y = "N hauls",
    title = "composition_1: Distribution of S, D, and size_CV across hauls",
    caption = sprintf("Histograms of species richness (S), Berger-Parker dominance (D), and dominant-species size_CV; N=%d hauls, `h1_dominance_haul_table.csv` joined to `haul_eeos_predictions.rds` for S.", nrow(comp))
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_composition_1, p_composition_1, width = 11, height = 4.2, dpi = 150)
log_figure("composition_1", path_fig_composition_1,
           sprintf("Distribution histograms of S, D, and size_CV per haul (N=%d).", nrow(comp)))

# ---------------------------------------------------------------------------
# composition_2: plain scatter of D vs size_CV (no ratio/residual overlay)
# ---------------------------------------------------------------------------
p_composition_2 <- ggplot(comp, aes(x = D, y = size_CV)) +
  geom_point(alpha = 0.12, size = 0.6, colour = "#2166ac") +
  geom_density2d(colour = "grey30", linewidth = 0.25) +
  labs(
    x = "D (Berger-Parker numerical dominance)", y = "size_CV (dominant-species mass CV)",
    title = "composition_2: D vs size_CV, plain distribution",
    caption = sprintf("Plain scatter of Berger-Parker dominance (D) vs dominant-species size_CV, N=%d hauls; density contours added for readability. No ratio/residual overlay (cf. Task 2's h1_dominance_ratio_vs_D.png / h1_dominance_ratio_vs_sizeCV.png, which plot log(B_pred/B_obs) against these axes instead).", nrow(comp))
  ) +
  theme_minimal(base_size = 11)
ggsave(path_fig_composition_2, p_composition_2, width = 7.5, height = 6, dpi = 150)
log_figure("composition_2", path_fig_composition_2,
           sprintf("Plain scatter of D vs size_CV, N=%d hauls, density contours only (no ratio/residual colouring).", nrow(comp)))

# ---------------------------------------------------------------------------
# composition_3: "big shoal" reference-share summary
# ---------------------------------------------------------------------------
logmsg("")
logmsg("### \"Big shoal\" cutoff check")
logmsg("Searched `pipeline/*.R`, `pipeline/R/*.R`, and `display_discussion/*.md` for any existing ",
       "formal cutoff jointly combining D and size_CV into a single \"big shoal\" classification. ",
       "None found: the H1 dominance analysis (`explore_h1_haul_dominance.R`, ",
       "`display_discussion/H1_dominance_results_summary.md`) treats the top decile of D and the ",
       "bottom decile of size_CV as separate marginal extremes for descriptive/taxonomic breakdown ",
       "only (e.g. \"D, top decile\" and \"size_CV, bottom decile\" reported as distinct rows), never ",
       "combined into one joint threshold. `display_discussion/Methods_note_for_discussion.md` uses ",
       "only qualitative language (\"high D and low size_CV\") with no numeric cutoff.")
logmsg("Per the briefing's fallback, this script marks the top decile of D as a visual reference ",
       "line only — explicitly not a formal \"big shoal\" classification. Reused the existing top-decile ",
       "boundary already computed in `outputs/h1_dominance_by_D_bins.csv` (produced by ",
       "`explore_h1_haul_dominance.R`) rather than computing a fresh quantile.")

d_ref <- read_D_top_decile_reference(path_D_bins)
comp <- comp %>%
  mutate(D_top_decile_ref = D >= d_ref$threshold)
n_ref <- sum(comp$D_top_decile_ref)
pct_ref <- round(100 * n_ref / nrow(comp), 1)
logmsg(sprintf("D top-decile reference threshold = %.4f (from h1_dominance_by_D_bins.csv, n_bins=10 bin=10, n=%d there vs n=%d here after the S-join; %.1f%% of hauls).",
               d_ref$threshold, d_ref$n_in_decile, n_ref, pct_ref))

ref_label <- sprintf("Top decile of D (D >= %.3f); n=%d, %.1f%% of hauls\n[reference line only, not a formal classification]", d_ref$threshold, n_ref, pct_ref)
comp <- comp %>%
  mutate(ref_group = ifelse(D_top_decile_ref, ref_label, "Remaining 9 deciles"))

p_composition_3 <- ggplot(comp, aes(x = D, y = size_CV, colour = ref_group)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_vline(xintercept = d_ref$threshold, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = setNames(c("#b2182b", "grey70"), c(ref_label, "Remaining 9 deciles")), name = NULL) +
  labs(
    x = "D (Berger-Parker numerical dominance)", y = "size_CV (dominant-species mass CV)",
    title = "composition_3: Share of hauls above the D top-decile reference line",
    caption = sprintf(
      "%d of %d hauls (%.1f%%) have D at or above the top-decile boundary (D >= %.3f, reused from `h1_dominance_by_D_bins.csv`); dashed line marks the boundary. No combined D x size_CV \"big shoal\" cutoff exists elsewhere in the H1 dominance analysis, so this top-decile-of-D line is a visual reference only, not a formal classification.",
      n_ref, nrow(comp), pct_ref, d_ref$threshold
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(path_fig_composition_3, p_composition_3, width = 7.5, height = 6.5, dpi = 150)
log_figure("composition_3", path_fig_composition_3,
           sprintf("%d of %d hauls (%.1f%%) fall at/above the D top-decile reference line (D >= %.3f); reference only, not a formal \"big shoal\" classification (none exists elsewhere in the H1 dominance analysis).",
                   n_ref, nrow(comp), pct_ref, d_ref$threshold))

# ---------------------------------------------------------------------------
# Wrap-up
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
for (id in names(figure_log)) {
  logmsg("- `", figure_log[[id]]$path, "` — ", figure_log[[id]]$caption)
}
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Missing visuals complete — haul map + composition figures only, H2 residual map deferred. ===\n")
