# Build DATRAS haul-level state variables (S, N, E_raw, m_min) from HL data.
#
# Audit requirement: S, N, and E must come from the same filtered DATRAS HL rows.
# When outputs/datras_hl_raw.rds covers the full Q1 1985-2015 period, all variables
# are computed from length bins (including m_min from the smallest individual bin).
# When HL coverage is incomplete, S/N/E_raw are derived from datras_haul_mean_length.rds
# and m_min uses HL bins where available (mean-length fallback otherwise).

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

path_hl_raw <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_mean_length <- file.path(project_root, "outputs", "datras_haul_mean_length.rds")
path_lw_lookup <- file.path(project_root, "outputs", "fishbase_lw_lookup_v2.csv")
path_out_state <- file.path(project_root, "outputs", "datras_haul_state.rds")
path_out_E <- file.path(project_root, "outputs", "datras_haul_E.rds")
path_out_mean_length <- file.path(project_root, "outputs", "datras_haul_mean_length.rds")
path_out_lngt_audit <- file.path(project_root, "outputs", "datras_lngt_code_audit.csv")
path_out_diag <- file.path(project_root, "outputs", "datras_build_diagnostics.csv")
path_out_excluded <- file.path(project_root, "outputs", "datras_excluded_hauls.csv")

stopifnot(file.exists(path_lw_lookup))

dir.create(dirname(path_out_state), recursive = TRUE, showWarnings = FALSE)

lw_lookup <- read.csv(path_lw_lookup, stringsAsFactors = FALSE)

# Import full HL from ICES CSV when RDS missing or incomplete
ensure_datras_hl_raw(project_root, path_hl_raw)

build_mode <- "mean_length_hybrid"
haul_state <- NULL
haul_mean_length_out <- NULL

if (file.exists(path_hl_raw)) {
  hl <- readRDS(path_hl_raw)
  coverage <- assess_hl_coverage(hl)
  lngt_audit <- lngt_code_audit_table(hl)
  write.csv(lngt_audit, path_out_lngt_audit, row.names = FALSE)

  cat("=== DATRAS HL coverage ===\n")
  cat(coverage$message, "\n")
  cat("LngtCode audit saved:", path_out_lngt_audit, "\n\n")

  if (coverage$complete) {
    built <- build_state_from_hl(hl, lw_lookup)
    haul_state <- built$haul_state
    haul_mean_length_out <- built$haul_mean_length
    build_mode <- "hl_bins_full"
    cat("Using full HL bin pipeline (S, N, E_raw, m_min from length bins).\n\n")
  } else {
    cat("HL file incomplete — hybrid mode:\n")
    cat("  S, N, E_raw from datras_haul_mean_length.rds + LW lookup\n")
    cat("  m_min from HL bins where available; mean-length fallback otherwise\n\n")

    if (!file.exists(path_mean_length)) {
      stop("Hybrid mode requires ", path_mean_length, " or a complete datras_hl_raw.rds")
    }
    haul_mean_length <- readRDS(path_mean_length)
    haul_state <- derive_state_from_mean_length(haul_mean_length, lw_lookup)
    m_min_hl <- derive_m_min_from_hl(hl, lw_lookup)
    haul_state <- apply_m_min_with_fallback(haul_state, m_min_hl)
    haul_mean_length_out <- haul_mean_length
    build_mode <- "mean_length_hybrid"
  }
} else {
  stop(
    "Missing ", path_hl_raw,
    " — upload NS-IBTS HL data or run with existing mean-length only via derive_state_from_mean_length."
  )
}

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

diag <- tibble(
  metric = c(
    "build_mode",
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
cat("Hauls:", nrow(haul_state), "\n")
cat("m_min method — min_length_bin:", sum(haul_state$m_min_method == "min_length_bin"),
    " mean_length_fallback:", sum(haul_state$m_min_method == "mean_length_fallback"), "\n")
cat("n_species_with_lw > S:", sum(haul_state$n_species_with_lw > haul_state$S), "\n\n")

cat("Saved:\n")
cat(" ", path_out_state, "\n")
cat(" ", path_out_E, "\n")
cat(" ", path_out_mean_length, "\n")
cat(" ", path_out_diag, "\n")
