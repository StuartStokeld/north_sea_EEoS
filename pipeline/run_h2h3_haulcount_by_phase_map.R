# Haul-count choropleth by analysis phase — individual ICES rectangles.
#
# Replaces the exploratory Scheme A 2x2 pre/post 2003 map for presentation in
# display_discussion/H2_H3_methods_justification. Uses the live H2/H3 analysis
# universe and the four structural-break phases (1985-1988, 1989-2000,
# 2001-2007, 2008-2015).
#
# Run: Rscript --vanilla pipeline/run_h2h3_haulcount_by_phase_map.R

suppressPackageStartupMessages({
  library(dplyr)
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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_haulcount_by_phase_map.R")
}

source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2_spatial_helpers.R"))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))

haul <- readRDS(file.path(project_root, "outputs", "haul_eeos_predictions.rds"))
panel <- readRDS(file.path(project_root, "outputs", "h2_rectangle_panel.rds"))
couce_year <- readRDS(file.path(project_root, "outputs", "h2_couce_year_effort.rds"))

# Same universe filter as the live within-between model, without needing glmmTMB
# (build_feasibility_data adds a numFactor spatial column we do not need here).
panel <- panel %>% mutate(stat_rec = normalize_stat_rec(stat_rec))
dat <- haul %>%
  mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
  filter(
    stat_rec %in% panel$stat_rec,
    is.finite(residual),
    year >= H2_YEAR_MIN,
    year <= H2_YEAR_MAX
  ) %>%
  inner_join(
    couce_year %>%
      mutate(stat_rec = normalize_stat_rec(stat_rec)) %>%
      select(stat_rec, year, hours_total),
    by = c("stat_rec", "year")
  ) %>%
  mutate(phase = build_phase_factor(year))

phase_counts <- dat %>%
  count(stat_rec, phase, name = "n_hauls")

ices_sf <- load_ices_rectangles_sf(project_root) %>%
  filter(stat_rec %in% panel$stat_rec) %>%
  select(stat_rec)

map_df <- ices_sf %>%
  inner_join(phase_counts, by = "stat_rec")

map_df$phase <- factor(
  as.character(map_df$phase),
  levels = c("1985-1988", "1989-2000", "2001-2007", "2008-2015")
)

p <- ggplot(map_df) +
  geom_sf(aes(fill = n_hauls), colour = "white", linewidth = 0.08) +
  facet_wrap(~phase, nrow = 2) +
  scale_fill_viridis_c(option = "C", name = "N hauls") +
  labs(
    title = "Haul count per ICES rectangle by analysis phase",
    subtitle = sprintf(
      "Live H2/H3 analysis universe: %d rectangles, %d hauls (1985–2015)",
      dplyr::n_distinct(dat$stat_rec), nrow(dat)
    ),
    caption = paste0(
      "Phases from structural breaks 1989 / 2001 / 2008. ",
      "Individual rectangles (no block merge). ",
      "Rectangles with zero hauls in a phase are absent from that panel."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(size = 8, hjust = 0)
  )

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(fig_dir, "h2h3_haulcount_by_phase_rectangles.png")
ggsave(out_path, p, width = 10, height = 9, dpi = 150)
message("Wrote: ", out_path)

dd_fig <- file.path(project_root, "display_discussion", "figures")
dir.create(dd_fig, recursive = TRUE, showWarnings = FALSE)
file.copy(out_path, file.path(dd_fig, basename(out_path)), overwrite = TRUE)
message("Copied to display_discussion/figures/")
