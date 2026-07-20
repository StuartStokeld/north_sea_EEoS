# H1 haul-level numerical dominance (D) and dominant-species size-homogeneity
# (size_CV) diagnostic — parallel track, exploratory only.
#
# Does NOT modify the primary H1 pipeline (log_r2 / cor2, null model, dropout
# diagnostics, catchability scaling). Reuses the already-built EEoS
# predictions (haul_eeos_predictions.rds) and DATRAS HL bins
# (datras_hl_raw.rds) exactly as the production pipeline built them.
#
# Question: does EEoS prediction failure vary systematically with how
# "typical" a haul's community structure is (numerically dominated by one
# species; that species caught at a near-uniform size), independent of / on
# top of the already-reported biomass-magnitude bias?

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/explore_h1_haul_dominance.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_hl_helpers.R"))
source(file.path(r_dir, "h1_dominance_helpers.R"))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
path_hl_raw <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_lw_lookup <- file.path(project_root, "outputs", "fishbase_lw_lookup_v2.csv")
path_haul_state <- file.path(project_root, "outputs", "haul_state_variables.rds")
path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_dropout_by_rec <- file.path(project_root, "outputs", "h1_dropout_by_stat_rec.csv")

path_out_haul_tbl <- file.path(project_root, "outputs", "h1_dominance_haul_table.csv")
path_out_corr <- file.path(project_root, "outputs", "h1_dominance_correlation_matrix.csv")
path_out_cond_corr <- file.path(project_root, "outputs", "h1_dominance_conditional_correlation.csv")
path_out_dq_by_year <- file.path(project_root, "outputs", "h1_dominance_dataquality_by_year.csv")
path_out_dq_by_rec <- file.path(project_root, "outputs", "h1_dominance_dataquality_by_stat_rec.csv")
path_out_dq_summary <- file.path(project_root, "outputs", "h1_dominance_dataquality_summary.csv")
path_out_tax <- file.path(project_root, "outputs", "h1_dominance_taxonomic_breakdown.csv")
path_out_D_bins <- file.path(project_root, "outputs", "h1_dominance_by_D_bins.csv")
path_out_cv_bins <- file.path(project_root, "outputs", "h1_dominance_by_sizeCV_bins.csv")
path_out_D_cond <- file.path(project_root, "outputs", "h1_dominance_by_D_x_bobs_quartile.csv")
path_out_cv_cond <- file.path(project_root, "outputs", "h1_dominance_by_sizeCV_x_bobs_quartile.csv")
path_out_fig_D <- file.path(project_root, "outputs", "figures", "h1_dominance_ratio_vs_D.png")
path_out_fig_cv <- file.path(project_root, "outputs", "figures", "h1_dominance_ratio_vs_sizeCV.png")
path_out_md <- file.path(
  project_root, "display_discussion", "H1_dominance_size_homogeneity_exploration.md"
)

stopifnot(
  file.exists(path_hl_raw),
  file.exists(path_lw_lookup),
  file.exists(path_haul_state),
  file.exists(path_preds)
)
dir.create(dirname(path_out_fig_D), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load inputs (all already-built pipeline artefacts; nothing recomputed
#    that would change N / E / B_pred for the primary H1 result)
# ---------------------------------------------------------------------------
hl_raw <- readRDS(path_hl_raw)
lw_lookup <- read.csv(path_lw_lookup, stringsAsFactors = FALSE)
haul_state <- readRDS(path_haul_state)
preds <- readRDS(path_preds)

# ---------------------------------------------------------------------------
# 2. Rebuild the same LW-matched bin set the production pipeline used for
#    N / E_raw / m_min (clean_hl_raw -> add_lw_mass_to_bins), so per-species
#    abundance and dominant-species bin mass are on an identical footing.
# ---------------------------------------------------------------------------
hl_bins <- clean_hl_raw(hl_raw)
hl_mass <- add_lw_mass_to_bins(hl_bins, lw_lookup)

# ---------------------------------------------------------------------------
# 3. Outcome table: haul_key <-> haul_id <-> N / B_obs / B_pred / residual
# ---------------------------------------------------------------------------
haul_key_map <- haul_state %>% select(haul_id, haul_key)
haul_outcomes <- preds %>%
  inner_join(haul_key_map, by = "haul_id") %>%
  select(haul_key, haul_id, year, stat_rec, N, B_obs, B_pred, residual)

n_preds <- nrow(preds)
n_outcomes <- nrow(haul_outcomes)
cat("EEoS predictions:", n_preds, "-> matched to haul_key via haul_state_variables:", n_outcomes, "\n")
if (n_outcomes < n_preds) {
  warning(
    n_preds - n_outcomes,
    " prediction hauls could not be mapped back to haul_key ",
    "(excluded from this diagnostic only; primary H1 result unaffected)."
  )
}

# ---------------------------------------------------------------------------
# 4. Per-haul dominance / size-CV table (deliverable 1)
# ---------------------------------------------------------------------------
haul_tbl <- build_haul_dominance_table(hl_mass, lw_lookup, haul_outcomes)

n_mismatch <- sum(abs(haul_tbl$n_recomputed_check - haul_tbl$N) > 1e-6)
cat(
  "Validation — haul-level abundance recomputed from LW-matched bins vs ",
  "production N: ", n_mismatch, " mismatches of ", nrow(haul_tbl), " hauls\n",
  sep = ""
)

haul_tbl_out <- haul_tbl %>%
  select(
    haul_id, haul_key, year, stat_rec,
    D, dominant_aphia_id, dominant_species,
    size_CV, n_bins_dominant_species,
    N, B_obs, B_pred, pred_obs_ratio, residual, ln_ratio, B_obs_quartile
  )
write.csv(haul_tbl_out, path_out_haul_tbl, row.names = FALSE)
cat("Saved", path_out_haul_tbl, "(", nrow(haul_tbl_out), "hauls )\n")

# ---------------------------------------------------------------------------
# 5. Confound-check step 1: correlation matrix (deliverable 2)
# ---------------------------------------------------------------------------
corr_vars <- c("D", "size_CV", "N", "B_obs", "n_bins_dominant_species")
corr_tbl <- dominance_correlation_matrix(haul_tbl, corr_vars)
write.csv(corr_tbl, path_out_corr, row.names = FALSE)
cat("Saved", path_out_corr, "\n")

# ---------------------------------------------------------------------------
# 6. Confound-check step 2: conditional (within B_obs quartile) (deliverable 3)
# ---------------------------------------------------------------------------
cond_corr_tbl <- conditional_correlation_by_bobs_quartile(haul_tbl, c("D", "size_CV"))
write.csv(cond_corr_tbl, path_out_cond_corr, row.names = FALSE)
cat("Saved", path_out_cond_corr, "\n")

# ---------------------------------------------------------------------------
# 7. Confound-check step 3: data-quality cross-check (deliverable 4)
#    Year flags: 1998 / 2013-2014 dropout spikes + 1985-only discovery-file
#    era (all historical artefacts already resolved; see H1_dropout_diagnosis.md).
#    Rectangle flags: any ICES rectangle with nonzero dropout in
#    h1_dropout_by_stat_rec.csv (bonus check beyond the year-based ask).
# ---------------------------------------------------------------------------
dq_year <- dataquality_crosstab(
  haul_tbl, group_col = "year", flagged_values = c(1985L, 1998L, 2013L, 2014L)
)
write.csv(dq_year$by_group, path_out_dq_by_year, row.names = FALSE)

dq_summary <- bind_rows(
  dq_year$flagged_summary %>% mutate(grouping = "year"),
  {
    if (file.exists(path_dropout_by_rec)) {
      dropout_rec <- read.csv(path_dropout_by_rec, stringsAsFactors = FALSE)
      flagged_rec <- dropout_rec$stat_rec[dropout_rec$pct_dropped > 0]
      dq_rec <- dataquality_crosstab(
        haul_tbl, group_col = "stat_rec", flagged_values = flagged_rec
      )
      write.csv(dq_rec$by_group, path_out_dq_by_rec, row.names = FALSE)
      cat("Saved", path_out_dq_by_rec, "\n")
      dq_rec$flagged_summary %>% mutate(grouping = "stat_rec")
    } else {
      tibble()
    }
  }
)
write.csv(dq_summary, path_out_dq_summary, row.names = FALSE)
cat("Saved", path_out_dq_by_year, "and", path_out_dq_summary, "\n")

# ---------------------------------------------------------------------------
# 8. Confound-check step 4: taxonomic breakdown, descriptive only (deliverable 5)
# ---------------------------------------------------------------------------
tax_tbl <- bind_rows(
  taxonomic_breakdown(haul_tbl, "D", extreme = "top"),
  taxonomic_breakdown(haul_tbl, "size_CV", extreme = "top"),
  taxonomic_breakdown(haul_tbl, "size_CV", extreme = "bottom")
)
write.csv(tax_tbl, path_out_tax, row.names = FALSE)
cat("Saved", path_out_tax, "\n")

# ---------------------------------------------------------------------------
# 9. Binned reporting tables — D and size_CV, analogous to the existing
#    B_obs-quartile magnitude-bias table (deliverable 6)
# ---------------------------------------------------------------------------
D_bins <- bind_rows(
  bin_ratio_table(haul_tbl, "D", n_bins = 4L),
  bin_ratio_table(haul_tbl, "D", n_bins = 10L)
)
cv_bins <- bind_rows(
  bin_ratio_table(haul_tbl, "size_CV", n_bins = 4L),
  bin_ratio_table(haul_tbl, "size_CV", n_bins = 10L)
)
write.csv(D_bins, path_out_D_bins, row.names = FALSE)
write.csv(cv_bins, path_out_cv_bins, row.names = FALSE)
cat("Saved", path_out_D_bins, "and", path_out_cv_bins, "\n")

D_cond <- bin_ratio_table_conditional(haul_tbl, "D", n_bins = 4L)
cv_cond <- bin_ratio_table_conditional(haul_tbl, "size_CV", n_bins = 4L)
write.csv(D_cond, path_out_D_cond, row.names = FALSE)
write.csv(cv_cond, path_out_cv_cond, row.names = FALSE)
cat("Saved", path_out_D_cond, "and", path_out_cv_cond, "\n")

# ---------------------------------------------------------------------------
# 10. Figures
# ---------------------------------------------------------------------------
p_D <- ggplot(haul_tbl, aes(x = D, y = ln_ratio)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "grey35") +
  geom_smooth(method = "loess", se = FALSE, colour = "#d6604d", linewidth = 0.9) +
  labs(
    x = "D (Berger-Parker numerical dominance)",
    y = "log(B_pred / B_obs)",
    title = "EEoS residual ratio vs haul numerical dominance",
    subtitle = "Diagnostic only — parallel to the B_obs-quartile magnitude-bias result"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_D, p_D, width = 8, height = 5.5, dpi = 120)

p_cv <- ggplot(haul_tbl, aes(x = size_CV, y = ln_ratio)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "grey35") +
  geom_smooth(method = "loess", se = FALSE, colour = "#4393c3", linewidth = 0.9) +
  labs(
    x = "size_CV (dominant-species mass CV)",
    y = "log(B_pred / B_obs)",
    title = "EEoS residual ratio vs dominant-species size homogeneity",
    subtitle = "Diagnostic only — parallel to the B_obs-quartile magnitude-bias result"
  ) +
  theme_minimal(base_size = 11)
ggsave(path_out_fig_cv, p_cv, width = 8, height = 5.5, dpi = 120)
cat("Saved", path_out_fig_D, "and", path_out_fig_cv, "\n")

# ---------------------------------------------------------------------------
# 11. Markdown write-up
# ---------------------------------------------------------------------------
md_table <- function(df, digits = 3L) {
  df <- as.data.frame(df)
  num_cols <- vapply(df, is.numeric, logical(1))
  for (col in names(df)[num_cols]) {
    df[[col]] <- format(round(df[[col]], digits), nsmall = 0, trim = TRUE)
  }
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(header, sep, body)
}

n_flag_corr <- sum(corr_tbl$flag_gt_0.4, na.rm = TRUE)

md <- c(
  "# H1 haul-level dominance and size-homogeneity exploration",
  "",
  paste(
    "Generated:", Sys.Date(),
    "— **parallel diagnostic track; does not modify the primary H1 pipeline",
    "(log_r2/cor2, null model, dropout diagnostics, catchability scaling).**"
  ),
  "",
  "## Question",
  "",
  "Does EEoS prediction failure vary systematically with how \"typical\" a haul's",
  "community structure is (numerically dominated by one species; that species",
  "caught at a near-uniform size), independent of / conditional on the already-",
  "reported biomass-magnitude bias (median `B_pred/B_obs` rising from ~3x in the",
  "lowest `B_obs` quartile to >5x in the highest)?",
  "",
  sprintf("N hauls in this diagnostic: **%d** (of %d EEoS predictions).", nrow(haul_tbl), n_preds),
  "",
  "## Metric definitions",
  "",
  "- **D** (Berger-Parker dominance) = `n_max_species / N_haul`, using the same",
  "  SubFactor-raised DATRAS HL abundance already used for `N` in the EEoS call",
  "  (LW-matched bins only — a species with no FishBase LW parameters is not",
  "  N-eligible under the production pipeline's own convention, so it cannot be",
  "  numerically dominant here either).",
  "- **size_CV** = raised-count-weighted coefficient of variation of bin-level",
  "  LW mass (`W = a*L^b` per 1cm bin, not on a single mean length) for the",
  "  dominant species' own length bins in that haul. `n_bins_dominant_species`",
  "  is reported alongside so mechanically-low CV from sparse bins can be told",
  "  apart from genuine single-cohort catches (see data-quality check below).",
  "",
  "## Confound check 1 — correlation matrix",
  "",
  md_table(corr_tbl),
  "",
  sprintf(
    "**%d of %d pairs** exceed |r| > 0.4 (either Pearson or Spearman) and are",
    n_flag_corr, nrow(corr_tbl)
  ),
  " therefore read via the conditional (within-B_obs-quartile) check below rather than marginally.",
  "",
  "## Confound check 2 — conditional on B_obs quartile",
  "",
  "Correlation of `D` and `size_CV` with `ln_ratio` (`log(B_pred/B_obs)`),",
  "recomputed separately within each existing `B_obs` quartile:",
  "",
  md_table(cond_corr_tbl),
  "",
  "If these within-quartile correlations are materially weaker than the marginal",
  "correlations above, dominance/size-homogeneity are largely re-describing the",
  "biomass-magnitude bias rather than an independent axis of failure.",
  "",
  "## Confound check 3 — data-quality cross-check",
  "",
  sprintf(
    "Bottom decile of `size_CV` (n = %d of %d hauls). Fraction falling in",
    dq_year$n_low_cv_total, dq_year$n_hauls_total
  ),
  " flagged years (1985 discovery-file era, 1998, 2013-2014 dropout spikes):",
  "",
  md_table(dq_year$flagged_summary),
  "",
  if (file.exists(path_dropout_by_rec)) {
    c(
      "Same check against ICES rectangles with any recorded EEoS-filter dropout",
      "(`h1_dropout_by_stat_rec.csv`, bonus check beyond the year-based ask):",
      "",
      md_table(dq_summary %>% filter(grouping == "stat_rec") %>% select(-grouping))
    )
  } else {
    character(0)
  },
  "",
  "Per `H1_dropout_diagnosis.md`, the 1998/2013-2014 spikes and 1985-only HL",
  "coverage are historical artefacts already fixed upstream (NA-LW propagation);",
  "the current `datras_hl_raw.rds` build (`hl_bins_full` mode) has 0 mean-length-",
  "fallback hauls. A concentration of low-CV hauls in these years would still be",
  "worth flagging as a residual data-quality signal rather than pure ecology.",
  "",
  "## Confound check 4 — taxonomic breakdown (descriptive only)",
  "",
  "Dominant-species identity in the extreme decile of each metric — descriptive",
  "only, **not** used to construct a species-based flag or filter:",
  "",
  md_table(tax_tbl),
  "",
  "## Binned reporting — D (Berger-Parker dominance)",
  "",
  "Quartile table, formatted analogously to the existing `B_obs`-quartile",
  "magnitude-bias result (`h1_catchability_by_quartile.csv`), with IQR added:",
  "",
  md_table(D_bins %>% filter(n_bins == 4L) %>% select(-metric, -n_bins)),
  "",
  "Decile table (finer resolution; sample size per bin ≈ n/10):",
  "",
  md_table(D_bins %>% filter(n_bins == 10L) %>% select(-metric, -n_bins)),
  "",
  "Broken out within each existing B_obs quartile (D quartile x B_obs quartile):",
  "",
  md_table(D_cond %>% select(-metric, -n_bins)),
  "",
  "## Binned reporting — size_CV (dominant-species size homogeneity)",
  "",
  "Quartile table:",
  "",
  md_table(cv_bins %>% filter(n_bins == 4L) %>% select(-metric, -n_bins)),
  "",
  "Decile table:",
  "",
  md_table(cv_bins %>% filter(n_bins == 10L) %>% select(-metric, -n_bins)),
  "",
  "Broken out within each existing B_obs quartile (size_CV quartile x B_obs quartile):",
  "",
  md_table(cv_cond %>% select(-metric, -n_bins)),
  "",
  "## Notes",
  "",
  "- `D` and `size_CV` are reported and binned **separately** throughout — no",
  "  composite dominance score is constructed at this stage.",
  "- This diagnostic reuses `outputs/datras_hl_raw.rds`, `outputs/fishbase_lw_lookup_v2.csv`,",
  "  `outputs/haul_state_variables.rds`, and `outputs/haul_eeos_predictions.rds` exactly",
  "  as built by the production pipeline; it does not recompute or alter N, E, B_pred,",
  "  or the primary residual.",
  "",
  sprintf("![Ratio vs D](../outputs/figures/%s)", basename(path_out_fig_D)),
  "",
  sprintf("![Ratio vs size_CV](../outputs/figures/%s)", basename(path_out_fig_cv)),
  "",
  "*Outputs: `outputs/h1_dominance_haul_table.csv`, `outputs/h1_dominance_correlation_matrix.csv`,",
  "`outputs/h1_dominance_conditional_correlation.csv`, `outputs/h1_dominance_dataquality_*.csv`,",
  "`outputs/h1_dominance_taxonomic_breakdown.csv`, `outputs/h1_dominance_by_D_bins.csv`,",
  "`outputs/h1_dominance_by_sizeCV_bins.csv`, `outputs/h1_dominance_by_D_x_bobs_quartile.csv`,",
  "`outputs/h1_dominance_by_sizeCV_x_bobs_quartile.csv`.*"
)

writeLines(md, path_out_md)
cat("Saved", path_out_md, "\n")

cat("\n=== H1 dominance / size-homogeneity diagnostic complete ===\n")
