# Pipeline audit diagnostics (cursor_pipeline_audit.md checks)
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
  stop("Run from pipeline/ or Rscript pipeline/run_pipeline_diagnostics.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "h1_dropout_helpers.R"))
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_hl_helpers.R"))
source(file.path(r_dir, "datras_state_helpers.R"))

path_state <- file.path(project_root, "outputs", "datras_haul_state.rds")
path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_failed <- file.path(project_root, "outputs", "eeos_failed_hauls.rds")
path_haul_state <- file.path(project_root, "outputs", "haul_state_variables.rds")
path_hl_raw <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_out <- file.path(project_root, "outputs", "pipeline_audit_results.csv")
path_dropout_year <- file.path(project_root, "outputs", "h1_dropout_by_year.csv")
path_dropout_rect <- file.path(project_root, "outputs", "h1_dropout_by_stat_rec.csv")

stopifnot(file.exists(path_state))

haul_state <- readRDS(path_state)
validation <- validate_haul_state(haul_state)

# Structural checks only (E > N is enforced at EEoS filter step, not on full DATRAS table)
structural <- validation %>%
  filter(check != "e_gt_n")

results <- structural %>%
  transmute(
    check_id = check,
    description = detail,
    status = if_else(pass, "PASS", "FAIL"),
    notes = NA_character_
  )

n_e_le_n_state <- sum(haul_state$E <= haul_state$N, na.rm = TRUE)
results <- bind_rows(
  results,
  tibble(
    check_id = "e_gt_n_pool",
    description = "Hauls with E <= N in full DATRAS pool (excluded before EEoS)",
    status = "INFO",
    notes = as.character(n_e_le_n_state)
  )
)

add_check <- function(id, desc, status, notes = NA_character_) {
  tibble(check_id = id, description = desc, status = status, notes = notes)
}

# Check 1 & 2: S/N from DATRAS (not FishGlob)
results <- bind_rows(
  results,
  add_check(
    "1_N_source",
    "N derived from DATRAS HL pipeline (not FishGlob num/num_cpue)",
    if ("N" %in% names(haul_state) && all(haul_state$N > 0, na.rm = TRUE)) "PASS" else "FAIL",
    "See build_datras_state_variables.R"
  ),
  add_check(
    "2_S_source",
    "S = n_distinct(AphiaID) in same LW-filtered DATRAS rows as E",
    if (all(haul_state$n_species_with_lw == haul_state$S, na.rm = TRUE)) "PASS" else "FAIL",
    sprintf("%d hauls with n_species_with_lw != S",
            sum(haul_state$n_species_with_lw != haul_state$S, na.rm = TRUE))
  )
)

# Check 7: n_species_with_lw > S
results <- bind_rows(
  results,
  add_check(
    "7_n_lw_vs_S",
    "n_species_with_lw <= S for all hauls",
    if (all(haul_state$n_species_with_lw <= haul_state$S, na.rm = TRUE)) "PASS" else "FAIL",
    sprintf("%d violations", sum(haul_state$n_species_with_lw > haul_state$S, na.rm = TRUE))
  )
)

# Check 4: m_min method
n_fallback <- sum(haul_state$m_min_method == "mean_length_fallback", na.rm = TRUE)
results <- bind_rows(
  results,
  add_check(
    "4_m_min_method",
    "m_min from smallest length bin (not species mean length)",
    if (n_fallback == 0L) "PASS" else "WARN",
    sprintf("%d hauls still on mean_length_fallback — upload full datras_hl_raw.rds", n_fallback)
  )
)

# Check 3a: LngtCode audit file
lngt_path <- file.path(project_root, "outputs", "datras_lngt_code_audit.csv")
if (file.exists(lngt_path)) {
  lngt <- read.csv(lngt_path, stringsAsFactors = FALSE)
  unknown <- lngt$lngt_code[!lngt$known_code & !is.na(lngt$lngt_code)]
  results <- bind_rows(
    results,
    add_check(
      "3a_LngtCode",
      "All LngtCode values in HL have known conversion",
      if (length(unknown) == 0L) "PASS" else "WARN",
      if (length(unknown)) paste("Unknown codes:", paste(unique(unknown), collapse = ", ")) else NA_character_
    )
  )
}

# Check 6 & 5: EEoS join funnel
if (file.exists(path_haul_state) && file.exists(path_preds)) {
  hs <- readRDS(path_haul_state)
  preds <- readRDS(path_preds)
  failed <- if (file.exists(path_failed)) readRDS(path_failed) else tibble()

  n_null_e <- sum(is.na(hs$E))
  n_s_gt_n <- sum(hs$S > hs$N, na.rm = TRUE)
  n_eeos <- nrow(preds)

  results <- bind_rows(
    results,
    add_check(
      "R1_null_E_joined",
      "No joined hauls with NA E in haul_state_variables",
      if (n_null_e == 0L) "PASS" else "FAIL",
      as.character(n_null_e)
    ),
    add_check(
      "6_E_gt_N_filtered",
      "All EEoS prediction hauls satisfy E > N",
      if (n_eeos > 0L && all(preds$E > preds$N, na.rm = TRUE)) "PASS" else "FAIL",
      sprintf("joined=%d; filtered predictions=%d", nrow(hs), n_eeos)
    ),
    add_check("6_S_le_N", "S <= N on joined hauls", if (n_s_gt_n == 0L) "PASS" else "FAIL", as.character(n_s_gt_n)),
    add_check(
      "5_join_funnel",
      "Documented haul counts at each pipeline step",
      "INFO",
      sprintf(
        "DATRAS valid=%d; FishGlob joined=%d; null E=%d; predictions=%d; EEoS failed=%d",
        nrow(haul_state), nrow(hs), n_null_e, nrow(preds), nrow(failed)
      )
    )
  )

  write.csv(dropout_by_year(hs, preds), path_dropout_year, row.names = FALSE)
  write.csv(dropout_by_stat_rec(hs, preds), path_dropout_rect, row.names = FALSE)

  if (file.exists(path_hl_raw)) {
    hl <- readRDS(path_hl_raw)
    funnel <- dropout_funnel_by_year(hl, haul_state, hs, preds)
    spike_verdict <- dropout_diagnosis_verdict(funnel)
    for (i in seq_len(nrow(spike_verdict))) {
      results <- bind_rows(
        results,
        add_check(
          spike_verdict$check[i],
          paste("1998/2013-14 dropout diagnosis:", spike_verdict$check[i]),
          spike_verdict$status[i],
          spike_verdict$detail[i]
        )
      )
    }
  }
}

# Check 10: residual convention
if (file.exists(path_preds)) {
  preds <- readRDS(path_preds)
  resid_ok <- all(abs(preds$residual - (log(preds$B_obs) - log(preds$B_pred))) < 1e-10, na.rm = TRUE)
  results <- bind_rows(
    results,
    add_check(
      "10_residual",
      "residual = log(B_obs) - log(B_pred); abs_residual primary metric",
      if (resid_ok) "PASS" else "FAIL",
      "residual_norm is diagnostic only"
    )
  )

  ratio_med <- median(preds$B_pred / preds$B_obs, na.rm = TRUE)
  r2_log <- log_r2(log(preds$B_obs), log(preds$B_pred))
  cor2_log <- log_cor2(log(preds$B_obs), log(preds$B_pred))
  results <- bind_rows(
    results,
    add_check(
      "4_B_scale",
      "B_obs in grams; median B_pred/B_obs documented (scale offset, not units)",
      "INFO",
      sprintf("median B_pred/B_obs = %.3f", ratio_med)
    ),
    add_check(
      "H1_r2_metric",
      "Primary H1 R² is log_r2 (SS-based), not cor(log B)²",
      "INFO",
      sprintf("log_r2 = %.3f; cor² = %.3f", r2_log, cor2_log)
    )
  )
}

# Check 8: quarter filter on HL
if (file.exists(path_hl_raw)) {
  hl <- readRDS(path_hl_raw)
  cov <- assess_hl_coverage(hl)
  q_ok <- all(hl$Quarter[hl$Year >= ANALYSIS_YEAR_MIN & hl$Year <= ANALYSIS_YEAR_MAX] == ANALYSIS_QUARTER, na.rm = TRUE)
  results <- bind_rows(
    results,
    add_check(
      "8_quarter_filter",
      "Analysis restricted to Q1 1985-2015 in HL processing",
      if (q_ok) "PASS" else "WARN",
      cov$message
    )
  )
}

write.csv(results, path_out, row.names = FALSE)

cat("=== Pipeline audit diagnostics ===\n\n")
print(results, n = Inf)
cat("\nSaved:", path_out, "\n")

n_fail <- sum(results$status == "FAIL")
n_warn <- sum(results$status == "WARN")
if (n_fail > 0L) {
  stop("Audit failed ", n_fail, " check(s). See ", path_out)
}
if (n_warn > 0L) {
  message("Audit completed with ", n_warn, " warning(s).")
} else {
  cat("All checks passed.\n")
}
