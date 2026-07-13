# H1 dropout diagnosis — year funnel and 1998 / 2013–14 spike check
# Run after build_datras_state_variables.R and build_eeos_predictions.R.

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h1_dropout_diagnosis.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
source(file.path(r_dir, "h1_dropout_helpers.R"))
project_root <- get_project_root_from(script_dir)

path_hl <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_datras <- file.path(project_root, "outputs", "datras_haul_state.rds")
path_hs <- file.path(project_root, "outputs", "haul_state_variables.rds")
path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_out_funnel <- file.path(project_root, "outputs", "h1_dropout_funnel_by_year.csv")
path_out_spike <- file.path(project_root, "outputs", "h1_dropout_spike_years.csv")
path_out_reasons <- file.path(project_root, "outputs", "h1_dropout_by_reason.csv")
path_out_filter_excl <- file.path(project_root, "outputs", "h1_filter_exclusions.csv")
path_out_md <- file.path(project_root, "display_discussion", "H1_dropout_diagnosis.md")

stopifnot(
  file.exists(path_hl),
  file.exists(path_datras),
  file.exists(path_hs),
  file.exists(path_preds)
)

hl <- readRDS(path_hl)
datras <- readRDS(path_datras)
hs <- readRDS(path_hs)
preds <- readRDS(path_preds)

funnel <- dropout_funnel_by_year(hl, datras, hs, preds)
spike <- diagnose_spike_years(funnel)
verdict <- dropout_diagnosis_verdict(funnel)

if (file.exists(path_out_filter_excl)) {
  reasons <- read.csv(path_out_filter_excl, stringsAsFactors = FALSE) %>%
    count(exclusion_reason, name = "n_hauls") %>%
    arrange(desc(n_hauls))
} else {
  reasons <- dropout_reason_table(hs, preds)$summary
}
write.csv(funnel, path_out_funnel, row.names = FALSE)
write.csv(spike, path_out_spike, row.names = FALSE)
write.csv(reasons, path_out_reasons, row.names = FALSE)

# Largest year-level filter drop (expected: 1991 zero-biomass hauls)
top_filter_year <- funnel %>%
  arrange(desc(join_to_pred_drop)) %>%
  slice(1)

md_lines <- c(
  "# H1 dropout diagnosis",
  "",
  paste("Generated:", Sys.Date(), "from `pipeline/run_h1_dropout_diagnosis.R`."),
  "",
  "## Root cause of historical 1998 / 2013–14 spikes",
  "",
  "Earlier drafts showed large year-level dropout in 1998 and 2013–14. That was **not** missing DATRAS HL years.",
  "Bins with missing FishBase LW parameters previously propagated `NA` through `sum()` / `min()`, producing hauls with `NA E_raw` / `NA m_min` that survived the join and appeared as mass dropout.",
  "",
  "**Fix (implemented):** `add_lw_mass_to_bins()` drops bad LW bins; `filter_valid_haul_state()` excludes incomplete hauls before `datras_haul_state.rds` is saved; `prepare_datras_for_join()` enforces valid normalised E before the FishGlob join (`stopifnot(null E == 0)`).",
  "",
  "## Spike-year check (1998, 2013, 2014)",
  "",
  "| Year | HL hauls | DATRAS state | FishGlob joined | EEoS predictions | HL→DATRAS drop | DATRAS→join drop | Join→pred drop | Status |",
  "|------|----------|--------------|-----------------|------------------|----------------|------------------|----------------|--------|"
)

for (i in seq_len(nrow(spike))) {
  row <- spike[i, ]
  md_lines <- c(
    md_lines,
    sprintf(
      "| %d | %d | %d | %d | %d | %d | %d | %d | %s |",
      row$Year, row$n_hl_hauls, row$n_datras_state, row$n_fishglob_joined,
      row$n_eeos_predictions, row$hl_to_datras_drop, row$datras_to_join_drop,
      row$join_to_pred_drop, row$spike_status
    )
  )
}

md_lines <- c(
  md_lines,
  "",
  "## Current largest filter-year drop",
  "",
  sprintf(
    "Year **%d**: %d joined hauls → %d predictions (%d excluded at EEoS filters).",
    top_filter_year$Year,
    top_filter_year$n_fishglob_joined,
    top_filter_year$n_eeos_predictions,
    top_filter_year$join_to_pred_drop
  ),
  "",
  "The dominant exclusion reason is **`bad_B_obs`** (FishGlob `sum(wgt) == 0`), concentrated in **1991** (40 hauls). This is intentional: `B_obs > 0` is required for log-scale comparison.",
  "",
  "## Filter exclusions (all years)",
  "",
  "| Reason | Hauls |",
  "|--------|------:|"
)

for (i in seq_len(nrow(reasons))) {
  md_lines <- c(
    md_lines,
    sprintf("| %s | %d |", reasons$exclusion_reason[i], reasons$n_hauls[i])
  )
}

md_lines <- c(
  md_lines,
  "",
  "## Audit verdict",
  "",
  "| Check | Status | Detail |",
  "|-------|--------|--------|"
)

for (i in seq_len(nrow(verdict))) {
  md_lines <- c(
    md_lines,
    sprintf("| %s | %s | %s |", verdict$check[i], verdict$status[i], verdict$detail[i])
  )
}

md_lines <- c(
  md_lines,
  "",
  sprintf(
    "**Panel:** %s hauls with EEoS predictions (%d DATRAS state hauls; %d FishGlob-matched).",
    format(nrow(preds), big.mark = ","),
    nrow(datras),
    nrow(hs)
  ),
  "",
  "*Outputs: `outputs/h1_dropout_funnel_by_year.csv`, `outputs/h1_filter_exclusions.csv`, `outputs/h1_join_gaps.csv`.*"
)

writeLines(md_lines, path_out_md)

cat("=== H1 dropout diagnosis ===\n\n")
print(verdict)
cat("\nSpike years:\n")
print(spike)
cat("\nSaved:\n")
cat(" ", path_out_funnel, "\n")
cat(" ", path_out_spike, "\n")
cat(" ", path_out_reasons, "\n")
cat(" ", path_out_md, "\n")

if (any(verdict$status == "WARN")) {
  message("Dropout diagnosis completed with warnings — review spike-year table.")
} else {
  cat("\nNo dropout spikes in 1998 / 2013 / 2014.\n")
}
