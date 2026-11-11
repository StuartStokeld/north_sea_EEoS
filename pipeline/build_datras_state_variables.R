# Build DATRAS haul-level state variables (S, N, E_raw, m_min) from HL data.
#
# Audit requirement: S, N, and E must come from the same filtered DATRAS HL rows.
# outputs/datras_hl_raw.rds must cover the full Q1 1985-2015 period; all variables
# are computed from length bins (including m_min from the smallest individual bin).
#
# Effort alignment (before build_state_from_hl): HL is joined to the HH haul records
# (outputs/datras_hh_raw.rds, from pipeline/import_datras_hh.R) and restricted to
# valid hauls (HaulVal == "V"), DataType C/R/P/S, and HaulDur within the window
# (primary: exactly 30 min; the 26-34 min sensitivity sets EEOS_HAUL_DUR_MIN=26
# EEOS_HAUL_DUR_MAX=34).
# DataType C counts are per hour and are converted to raw per-haul counts:
#   n_raw = HLNoAtLngt * SubFactor * f,  f = HaulDur / 60 (C), 1 (R/P/S)
# No duration rescaling. The former mean-length hybrid build is disabled: it cannot
# be effort-aligned.

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
  stop("Run from pipeline/ or Rscript pipeline/build_datras_state_variables.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_hl_helpers.R"))
source(file.path(r_dir, "datras_state_helpers.R"))
source(file.path(r_dir, "datras_csv_import.R"))
source(file.path(r_dir, "datras_hh_helpers.R"))

path_hl_raw <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_hh_raw <- file.path(project_root, "outputs", "datras_hh_raw.rds")
path_lw_lookup <- file.path(project_root, "outputs", "fishbase_lw_lookup_v2.csv")
path_out_state <- file.path(project_root, "outputs", "datras_haul_state.rds")
path_out_E <- file.path(project_root, "outputs", "datras_haul_E.rds")
path_out_mean_length <- file.path(project_root, "outputs", "datras_haul_mean_length.rds")
path_out_lngt_audit <- file.path(project_root, "outputs", "datras_lngt_code_audit.csv")
path_out_diag <- file.path(project_root, "outputs", "datras_build_diagnostics.csv")
path_out_excluded <- file.path(project_root, "outputs", "datras_excluded_hauls.csv")
path_out_effort_summary <- file.path(project_root, "outputs", "datras_effort_alignment_summary.csv")
path_out_raw_count_check <- file.path(project_root, "outputs", "datras_effort_raw_count_check.csv")

stopifnot(file.exists(path_lw_lookup))

dir.create(dirname(path_out_state), recursive = TRUE, showWarnings = FALSE)

lw_lookup <- read.csv(path_lw_lookup, stringsAsFactors = FALSE)

# Import full HL from ICES CSV when RDS missing or incomplete
ensure_datras_hl_raw(project_root, path_hl_raw)
# Import HH (haul duration, validity, DataType) when RDS missing
ensure_datras_hh_raw(project_root, path_hh_raw)

if (!file.exists(path_hl_raw)) {
  stop("Missing ", path_hl_raw, " — place the ICES unaggregated CSV under the project root.")
}

hl <- readRDS(path_hl_raw)
coverage <- assess_hl_coverage(hl)
lngt_audit <- lngt_code_audit_table(hl)
write.csv(lngt_audit, path_out_lngt_audit, row.names = FALSE)

cat("=== DATRAS HL coverage ===\n")
cat(coverage$message, "\n")
cat("LngtCode audit saved:", path_out_lngt_audit, "\n\n")

if (!coverage$complete) {
  stop(
    "datras_hl_raw.rds is incomplete (", coverage$message, ").\n",
    "The mean-length hybrid build is disabled: it cannot be aligned to haul effort ",
    "(DataType C counts are per hour; HL bins are needed to convert them per haul).\n",
    "Re-import the full HL: Rscript pipeline/import_datras_hl_from_csv.R --force"
  )
}

# --- Effort alignment: HL x HH ---------------------------------------------------
window <- haul_dur_window()
hh <- readRDS(path_hh_raw)
aligned <- align_hl_to_hh(hl, hh, dur_min = window[["min"]], dur_max = window[["max"]])

write.csv(aligned$summary, path_out_effort_summary, row.names = FALSE)
write.csv(aligned$raw_count_check, path_out_raw_count_check, row.names = FALSE)

cat("\n=== Effort alignment (hauls by reason x DataType) ===\n")
print(as.data.frame(aligned$summary))
cat("\n=== Raw-count check (after conversion) ===\n")
print(as.data.frame(aligned$raw_count_check))
cat("\n")

built <- build_state_from_hl(aligned$hl, lw_lookup)
haul_state <- built$haul_state %>%
  left_join(
    aligned$hauls %>% select(haul_key, DataType, HaulDur, effort_f),
    by = "haul_key"
  )
haul_mean_length_out <- built$haul_mean_length
build_mode <- "hl_bins_effort_aligned"
cat("Using effort-aligned HL bin pipeline (S, N, E_raw, m_min from raw per-haul counts).\n\n")

n_hauls_before_filter <- nrow(haul_state)
excluded_hauls <- haul_state %>%
  anti_join(filter_valid_haul_state(haul_state), by = "haul_key")
haul_state <- filter_valid_haul_state(haul_state)
n_excluded <- nrow(excluded_hauls)

if (n_excluded > 0L) {
  write.csv(
    excluded_hauls %>%
      select(haul_key, Survey, Year, Quarter, Country, Platform, HaulNumber, S, N, E_raw, m_min),
    path_out_excluded,
    row.names = FALSE
  )
  cat("Excluded hauls (no valid LW mass / E_raw):", n_excluded, "\n")
  cat("Saved:", path_out_excluded, "\n\n")
}

haul_state <- normalize_haul_E(haul_state)

validation <- validate_haul_state(haul_state)
cat("=== Internal validation ===\n")
print(as.data.frame(validation))
if (!all(validation$pass[validation$check != "e_gt_n"])) {
  warning("Structural haul-state validation checks failed — review before EEoS run.")
}
if (!validation$pass[validation$check == "e_gt_n"]) {
  n_bad <- sum(haul_state$E <= haul_state$N, na.rm = TRUE)
  message(
    n_bad, " hauls have E <= N (excluded by build_eeos_predictions.R filters)."
  )
}

# Haul_E table (E_raw + metadata; E normalised applied downstream in build_eeos)
haul_E <- haul_state %>%
  select(
    haul_key, Survey, Year, Quarter, Country, Platform, HaulNumber,
    E = E_raw, n_species_with_lw, S, N, m_min, m_min_method
  )

saveRDS(haul_state, path_out_state)
saveRDS(haul_E, path_out_E)
saveRDS(haul_mean_length_out, path_out_mean_length)

effort_counts <- table(factor(aligned$hauls$effort_reason, levels = c(
  "kept", "no_hh", "ambiguous_hh_key", "invalid", "datatype", "duration"
)))

diag <- tibble(
  metric = c(
    "build_mode",
    "haul_dur_min", "haul_dur_max",
    "n_hl_hauls", "hl_hh_match_rate",
    paste0("n_effort_", names(effort_counts)),
    "n_hauls_before_filter",
    "n_hauls_excluded",
    "n_hauls",
    "S_min", "S_max",
    "N_min", "N_max",
    "E_raw_min", "E_raw_max",
    "m_min_min", "m_min_max",
    "n_min_length_bin", "n_mean_length_fallback",
    "n_species_gt_s",
    "n_na_E"
  ),
  value = c(
    build_mode,
    window[["min"]], window[["max"]],
    nrow(aligned$hauls), round(aligned$match_rate, 5),
    as.integer(effort_counts),
    n_hauls_before_filter,
    n_excluded,
    nrow(haul_state),
    min(haul_state$S, na.rm = TRUE), max(haul_state$S, na.rm = TRUE),
    min(haul_state$N, na.rm = TRUE), max(haul_state$N, na.rm = TRUE),
    min(haul_state$E_raw, na.rm = TRUE), max(haul_state$E_raw, na.rm = TRUE),
    min(haul_state$m_min, na.rm = TRUE), max(haul_state$m_min, na.rm = TRUE),
    sum(haul_state$m_min_method == "min_length_bin", na.rm = TRUE),
    sum(haul_state$m_min_method == "mean_length_fallback", na.rm = TRUE),
    sum(haul_state$n_species_with_lw > haul_state$S, na.rm = TRUE),
    sum(is.na(haul_state$E))
  )
)
write.csv(diag, path_out_diag, row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Build mode:", build_mode, "\n")
cat("Haul-duration window:", window[["min"]], "-", window[["max"]], "min\n")
cat("HL hauls:", nrow(aligned$hauls), " effort-kept:", effort_counts[["kept"]], "\n")
cat("Hauls:", nrow(haul_state), "\n")
cat("m_min method — min_length_bin:", sum(haul_state$m_min_method == "min_length_bin"),
    " mean_length_fallback:", sum(haul_state$m_min_method == "mean_length_fallback"), "\n")
cat("n_species_with_lw > S:", sum(haul_state$n_species_with_lw > haul_state$S), "\n\n")

cat("Saved:\n")
cat(" ", path_out_state, "\n")
cat(" ", path_out_E, "\n")
cat(" ", path_out_mean_length, "\n")
cat(" ", path_out_diag, "\n")
cat(" ", path_out_effort_summary, "\n")
cat(" ", path_out_raw_count_check, "\n")
