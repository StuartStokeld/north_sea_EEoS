# Build EEoS haul-level B predictions (Hypothesis 1)
# B_obs from FishGlob; S, N, E (normalised) from DATRAS haul state (see build_datras_state_variables.R).

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(reticulate)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/build_eeos_predictions.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
source(file.path(script_dir, "R", "h1_join_helpers.R"))
project_root <- get_project_root_from(script_dir)

path_fishglob <- file.path(
  project_root, "FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData"
)
path_datras_state <- file.path(project_root, "outputs", "datras_haul_state.rds")
path_biomass_py <- file.path(project_root, "equation_of_state", "biomass.py")
path_venv_python <- file.path(project_root, ".venv", "bin", "python")
path_out_state <- file.path(project_root, "outputs", "haul_state_variables.rds")
path_out_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_out_failed <- file.path(project_root, "outputs", "eeos_failed_hauls.rds")
path_out_dropout <- file.path(project_root, "outputs", "h1_dropout_summary.csv")
path_out_join_gaps <- file.path(project_root, "outputs", "h1_join_gaps.csv")
path_out_filter_excl <- file.path(project_root, "outputs", "h1_filter_exclusions.csv")

stopifnot(
  file.exists(path_fishglob),
  file.exists(path_datras_state),
  file.exists(path_biomass_py),
  file.exists(path_venv_python)
)

# ---------------------------------------------------------------------------
# 1. Load FishGlob — B_obs only (grams) + join key
# ---------------------------------------------------------------------------
fishglob_haul <- build_fishglob_haul_table(path_fishglob)
n_fishglob_hauls <- nrow(fishglob_haul)

# ---------------------------------------------------------------------------
# 2. Prepare DATRAS state (valid E only) and join — null-E fix
# ---------------------------------------------------------------------------
datras_state <- readRDS(path_datras_state)
n_datras_total <- nrow(datras_state)
datras_ready <- prepare_datras_for_join(datras_state)
n_datras_ready <- nrow(datras_ready)
n_datras_excluded_prejoin <- n_datras_total - n_datras_ready

haul_joined <- join_fishglob_datras(fishglob_haul, datras_ready)
n_joined <- nrow(haul_joined)
n_unmatched_fishglob <- n_fishglob_hauls - n_joined
n_null_e <- sum(is.na(haul_joined$E))

join_gaps <- summarise_join_gaps(fishglob_haul, datras_ready)
write.csv(join_gaps, path_out_join_gaps, row.names = FALSE)

cat("FishGlob hauls (Q1, 1985-2015):", n_fishglob_hauls, "\n")
cat("DATRAS hauls with valid E (pre-join):", n_datras_ready, "\n")
cat("DATRAS excluded before join (invalid E/m_min):", n_datras_excluded_prejoin, "\n")
cat("Matched to DATRAS state:", n_joined, "\n")
cat("Unmatched FishGlob hauls:", n_unmatched_fishglob, "\n")
cat("Null E after join:", n_null_e, "\n")
stopifnot(n_null_e == 0L)

haul_state_variables <- select_haul_state_columns(haul_joined)

saveRDS(haul_state_variables, path_out_state)
cat("Saved", path_out_state, "\n")

# ---------------------------------------------------------------------------
# 4. Apply filters (EEoS input validation)
# ---------------------------------------------------------------------------
haul_joined$exclusion_reason <- classify_eeos_exclusion(haul_joined)
excluded_hauls <- haul_joined %>%
  filter(exclusion_reason != "passed") %>%
  select(haul_id, join_key, year, stat_rec, S, N, E, B_obs, exclusion_reason)
write.csv(excluded_hauls, path_out_filter_excl, row.names = FALSE)

haul_filtered <- filter_eeos_inputs(haul_joined)

n_filtered <- nrow(haul_filtered)
n_dropped <- n_joined - n_filtered

dropout_summary <- tibble(
  stage = c(
    "fishglob_q1_hauls",
    "datras_valid_E",
    "matched_datras_state",
    "null_E_after_join",
    "failed_eeos_filters",
    "eeos_predictions"
  ),
  n_hauls = c(
    n_fishglob_hauls,
    n_datras_ready,
    n_joined,
    n_null_e,
    n_dropped,
    NA_integer_
  )
)
write.csv(dropout_summary, path_out_dropout, row.names = FALSE)

cat("Saved", path_out_join_gaps, "(", nrow(join_gaps), "unmatched hauls )\n")
cat("Saved", path_out_filter_excl, "(", nrow(excluded_hauls), "excluded )\n")

cat("Hauls passing all filters:", n_filtered, "\n")
if (n_dropped > 0L) {
  cat("Hauls dropped at filter step:", n_dropped, "\n")
  if (n_null_e > 0L) {
    cat("  (includes ", n_null_e, " with NA E from join)\n", sep = "")
  }
}

# ---------------------------------------------------------------------------
# 5. Configure Python and load biomass.py
# ---------------------------------------------------------------------------
use_python(path_venv_python, required = TRUE)
source_python(path_biomass_py)

py_run_string("
def predict_biomass(S, N, E):
    import pandas as pd
    s = pd.Series({'S': float(S), 'N': float(N), 'E': float(E)})
    return biomass(s)
")

predict_biomass <- function(S, N, E) {
  py$predict_biomass(S, N, E)
}

# ---------------------------------------------------------------------------
# 6. Call EEoS per haul
# ---------------------------------------------------------------------------
predict_one_haul <- function(row) {
  tryCatch(
    {
      B_pred_norm <- predict_biomass(row$S, row$N, row$E)
      if (!is.finite(B_pred_norm) || B_pred_norm <= 0) {
        stop("Non-finite or non-positive B_pred_norm")
      }
      B_pred <- B_pred_norm * row$m_min
      tibble(
        haul_id = row$haul_id,
        haul_key = row$haul_key,
        year = row$year,
        stat_rec = row$stat_rec,
        S = row$S,
        N = row$N,
        E = row$E,
        E_raw = row$E_raw,
        B_obs = row$B_obs,
        B_pred_norm = B_pred_norm,
        B_pred = B_pred,
        m_min = row$m_min,
        m_min_method = row$m_min_method,
        min_epsilon = row$min_epsilon,
        n_species_with_lw = row$n_species_with_lw,
        ln_B_obs = log(B_obs),
        ln_B_pred = log(B_pred),
        ln_B_pred_norm = log(B_pred_norm),
        # Primary H1 residual (log grams); negative => EEoS over-predicts biomass.
        residual = log(B_obs) - log(B_pred),
        abs_residual = abs(log(B_obs) - log(B_pred)),
        residual_norm = log(B_obs) - log(B_pred_norm),
        abs_residual_norm = abs(log(B_obs) - log(B_pred_norm))
      )
    },
    error = function(e) {
      tibble(
        haul_id = row$haul_id,
        haul_key = row$haul_key,
        year = row$year,
        S = row$S,
        N = row$N,
        E = row$E,
        error_msg = conditionMessage(e)
      )
    }
  )
}

cat("Running EEoS predictions on", n_filtered, "hauls...\n")

results <- map(
  seq_len(n_filtered),
  function(i) {
    if (i %% 500 == 0) cat("  processed", i, "of", n_filtered, "\n")
    predict_one_haul(haul_filtered[i, ])
  }
)

results_df <- bind_rows(results)

if ("error_msg" %in% names(results_df)) {
  failed_hauls <- results_df %>% filter(!is.na(error_msg))
  haul_predictions <- results_df %>% filter(is.na(error_msg))
} else {
  failed_hauls <- tibble()
  haul_predictions <- results_df
}

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------
saveRDS(haul_predictions, path_out_preds)
saveRDS(failed_hauls, path_out_failed)
dropout_summary$n_hauls[dropout_summary$stage == "eeos_predictions"] <- nrow(haul_predictions)
write.csv(dropout_summary, path_out_dropout, row.names = FALSE)
cat("Saved", path_out_preds, "\n")
cat("Saved", path_out_failed, "\n")
cat("Saved", path_out_dropout, "\n")

# ---------------------------------------------------------------------------
# 8. Summary statistics
# ---------------------------------------------------------------------------
cat("\n========== Summary ==========\n")
cat("Total hauls assembled (joined):", n_joined, "\n")
cat("Hauls passing all filters:", n_filtered, "\n")
cat("Hauls with successful EEoS prediction:", nrow(haul_predictions), "\n")
cat("Hauls failed to converge:", nrow(failed_hauls), "\n")
if (n_filtered > 0) {
  cat(
    "Convergence failure rate:",
    round(100 * nrow(failed_hauls) / n_filtered, 2), "%\n"
  )
}

if (nrow(haul_predictions) > 0) {
  norm_ok <- all(abs(haul_predictions$min_epsilon - 1) < 1e-10, na.rm = TRUE)
  cat("E normalisation check — min epsilon == 1:", norm_ok, "\n")
  cat(
    "Median B_pred/B_obs (systematic scale offset; not a unit bug):",
    round(median(haul_predictions$B_pred / haul_predictions$B_obs, na.rm = TRUE), 3), "\n"
  )
  cat(
    "Mean log residual (negative => EEoS over-predicts):",
    round(mean(haul_predictions$residual, na.rm = TRUE), 3), "\n"
  )
  cat(
    "OLS slope log(B_obs) ~ log(B_pred):",
    round(coef(lm(ln_B_obs ~ ln_B_pred, data = haul_predictions))[2], 3), "\n"
  )

  r2_scaled <- log_r2(haul_predictions$ln_B_obs, haul_predictions$ln_B_pred)
  r2_norm <- log_r2(haul_predictions$ln_B_obs, haul_predictions$ln_B_pred_norm)
  cor2_scaled <- log_cor2(haul_predictions$ln_B_obs, haul_predictions$ln_B_pred)

  cat("R² log scale (primary H1 metric, scaled B_pred):", round(r2_scaled, 4), "\n")
  cat("cor(log B)² (association only, not H1 R²):", round(cor2_scaled, 4), "\n")
  cat("R² (log scale, unscaled B_pred_norm):", round(r2_norm, 4), "\n")
  cat(
    "Median absolute residual (scaled):",
    round(median(haul_predictions$abs_residual, na.rm = TRUE), 4), "\n"
  )
  cat(
    "Median absolute residual (unscaled):",
    round(median(haul_predictions$abs_residual_norm, na.rm = TRUE), 4), "\n"
  )
}

cat("=============================\n")
