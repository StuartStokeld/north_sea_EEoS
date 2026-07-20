# Haul-level numerical dominance (D) and dominant-species size-homogeneity (size_CV)
# helpers — parallel diagnostic track, does not modify N/E/B_pred used by the
# primary H1 pipeline (log_r2 / cor2, null model, dropout diagnostics, catchability
# scaling all proceed unchanged).
#
# Source R/datras_hl_helpers.R (clean_hl_raw, add_lw_mass_to_bins) before this file.

suppressPackageStartupMessages({
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Metric 1: Berger-Parker numerical dominance D = n_max_species / N_haul
# ---------------------------------------------------------------------------

#' Per-haul, per-species raised abundance from LW-matched length bins.
#'
#' Uses `hl_mass` (the output of `add_lw_mass_to_bins()`), not the raw
#' `clean_hl_raw()` bins, so the abundance pool exactly matches the N already
#' used in the EEoS call: `aggregate_haul_state_from_bins()` sums
#' `NumberAtLength` over this same LW-matched bin set. A species with no
#' FishBase LW match can never be "N-eligible" under the production pipeline's
#' own convention, so it is excluded here too for internal consistency.
species_haul_abundance <- function(hl_mass) {
  hl_mass %>%
    group_by(haul_key, AphiaID) %>%
    summarise(n_species = sum(NumberAtLength, na.rm = TRUE), .groups = "drop")
}

#' Identify the numerically dominant species per haul.
#'
#' Ties (equal raised abundance) broken by lowest AphiaID for determinism —
#' this is a tie-break convention only, not an ecological claim.
dominant_species_per_haul <- function(species_abundance) {
  species_abundance %>%
    group_by(haul_key) %>%
    arrange(desc(n_species), AphiaID, .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup() %>%
    rename(dominant_aphia_id = AphiaID, n_max_species = n_species)
}

#' Berger-Parker dominance D = n_max_species / N_haul.
#'
#' `haul_N` must supply the *production* N per haul (haul_key, N) — the same
#' N column already used in the EEoS call (haul_state_variables.rds /
#' haul_eeos_predictions.rds) — so D is guaranteed consistent with N rather
#' than recomputed from a differently-filtered bin set. `n_recomputed_check`
#' is the independently summed raised abundance across the same LW-matched
#' bins, returned purely for validation (should equal `N` after the join).
compute_haul_dominance <- function(hl_mass, haul_N) {
  stopifnot(all(c("haul_key", "N") %in% names(haul_N)))

  species_abundance <- species_haul_abundance(hl_mass)
  dominant <- dominant_species_per_haul(species_abundance)

  n_check <- species_abundance %>%
    group_by(haul_key) %>%
    summarise(n_recomputed_check = sum(n_species, na.rm = TRUE), .groups = "drop")

  dominant %>%
    left_join(n_check, by = "haul_key") %>%
    inner_join(haul_N, by = "haul_key") %>%
    mutate(D = n_max_species / N)
}

# ---------------------------------------------------------------------------
# Metric 2: dominant-species size homogeneity (raised-count-weighted CV of mass)
# ---------------------------------------------------------------------------

#' Raised-count-weighted CV of mass for the dominant species' length bins.
#'
#' count_i = NumberAtLength (SubFactor-raised); mass_i = bin-level LW mass (g),
#' i.e. W = a * L^b applied per 1cm bin (not to a single mean length), matching
#' the existing haul-mean-length pipeline convention. Uses the population
#' (not sample) weighted variance per the spec: denominator = sum(count_i).
#'
#' `n_bins_dominant_species` is reported alongside `size_CV` so low-bin-count
#' hauls (mechanically low CV) can be flagged rather than silently trusted —
#' see the data-quality cross-check.
dominant_species_size_cv <- function(hl_mass, dominant_tbl) {
  dom_key <- dominant_tbl %>%
    select(haul_key, AphiaID = dominant_aphia_id)

  bins <- hl_mass %>%
    inner_join(dom_key, by = c("haul_key", "AphiaID"))

  bins %>%
    group_by(haul_key) %>%
    summarise(
      n_bins_dominant_species = n(),
      mean_mass = sum(NumberAtLength * mass_g) / sum(NumberAtLength),
      var_mass = sum(NumberAtLength * (mass_g - mean_mass)^2) / sum(NumberAtLength),
      size_CV = sqrt(var_mass) / mean_mass,
      .groups = "drop"
    ) %>%
    select(haul_key, size_CV, n_bins_dominant_species, mean_mass_dominant_species = mean_mass)
}

# ---------------------------------------------------------------------------
# Assembly: per-haul diagnostic table
# ---------------------------------------------------------------------------

#' Assemble the full per-haul dominance / size-homogeneity diagnostic table.
#'
#' `haul_outcomes` must supply haul_key, haul_id, year, stat_rec, N, B_obs,
#' B_pred, residual (log-scale, log(B_obs) - log(B_pred)) for the hauls that
#' passed the primary EEoS pipeline (haul_eeos_predictions.rds).
build_haul_dominance_table <- function(hl_mass, lw_lookup, haul_outcomes) {
  required <- c("haul_key", "N", "B_obs", "B_pred", "residual")
  missing <- setdiff(required, names(haul_outcomes))
  if (length(missing) > 0L) {
    stop("haul_outcomes missing columns: ", paste(missing, collapse = ", "))
  }

  dominant <- compute_haul_dominance(hl_mass, haul_outcomes %>% select(haul_key, N))
  size_cv <- dominant_species_size_cv(hl_mass, dominant)

  species_names <- lw_lookup %>%
    select(dominant_aphia_id = aphia_id, dominant_species = accepted_name) %>%
    distinct(dominant_aphia_id, .keep_all = TRUE)

  haul_outcomes_no_n <- haul_outcomes %>% select(-N)

  dominant %>%
    left_join(size_cv, by = "haul_key") %>%
    left_join(species_names, by = "dominant_aphia_id") %>%
    inner_join(haul_outcomes_no_n, by = "haul_key") %>%
    mutate(
      pred_obs_ratio = B_pred / B_obs,
      ln_ratio = log(B_pred) - log(B_obs),
      B_obs_quartile = ntile(log(B_obs), 4L)
    )
}

# ---------------------------------------------------------------------------
# Confound-check step 1: correlation matrix
# ---------------------------------------------------------------------------

#' Pearson + Spearman correlation matrix for a set of numeric columns.
#'
#' Flags |r| > 0.4 (either method) as needing the conditional (step 2) check
#' rather than a marginal reading, per the briefing.
dominance_correlation_matrix <- function(df, vars) {
  combos <- utils::combn(vars, 2L, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    x <- df[[pair[1L]]]
    y <- df[[pair[2L]]]
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)
    tibble(
      var1 = pair[1L],
      var2 = pair[2L],
      n = n_ok,
      pearson_r = if (n_ok >= 3L) suppressWarnings(cor(x[ok], y[ok], method = "pearson")) else NA_real_,
      spearman_r = if (n_ok >= 3L) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
    )
  })
  bind_rows(rows) %>%
    mutate(flag_gt_0.4 = pmax(abs(pearson_r), abs(spearman_r), na.rm = TRUE) > 0.4)
}

# ---------------------------------------------------------------------------
# Confound-check step 2: conditional (within B_obs-quartile) relationships
# ---------------------------------------------------------------------------

#' Correlation of each predictor with ln_ratio, computed separately within
#' each existing B_obs quartile (tests whether D / size_CV explain residual
#' variation beyond biomass magnitude, not merely tracking it).
conditional_correlation_by_bobs_quartile <- function(
    df, vars, target = "ln_ratio", quartile_col = "B_obs_quartile") {
  df %>%
    group_by(.data[[quartile_col]]) %>%
    group_modify(function(chunk, key) {
      rows <- lapply(vars, function(v) {
        x <- chunk[[v]]
        y <- chunk[[target]]
        ok <- is.finite(x) & is.finite(y)
        n_ok <- sum(ok)
        tibble(
          predictor = v,
          n = n_ok,
          pearson_r = if (n_ok >= 3L) suppressWarnings(cor(x[ok], y[ok])) else NA_real_,
          spearman_r = if (n_ok >= 3L) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
        )
      })
      bind_rows(rows)
    }) %>%
    ungroup() %>%
    rename(b_obs_quartile = !!quartile_col)
}

# ---------------------------------------------------------------------------
# Binning / reporting: median + IQR of pred_obs_ratio by metric bin
# ---------------------------------------------------------------------------

#' Bin hauls by a continuous metric and report median/IQR of pred_obs_ratio.
#'
#' Mirrors the existing biomass-quartile reporting style (median ratio, n)
#' with IQR added as an extension (not present in the original quartile
#' table format).
bin_ratio_table <- function(df, metric_col, n_bins = 4L, ratio_col = "pred_obs_ratio") {
  df <- df %>% filter(is.finite(.data[[metric_col]]), is.finite(.data[[ratio_col]]))
  df$.bin <- ntile(df[[metric_col]], n_bins)
  df %>%
    filter(!is.na(.bin)) %>%
    group_by(.bin) %>%
    summarise(
      n = n(),
      metric_min = min(.data[[metric_col]], na.rm = TRUE),
      metric_max = max(.data[[metric_col]], na.rm = TRUE),
      median_ratio = median(.data[[ratio_col]], na.rm = TRUE),
      iqr_ratio_low = quantile(.data[[ratio_col]], 0.25, na.rm = TRUE),
      iqr_ratio_high = quantile(.data[[ratio_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(bin = .bin) %>%
    mutate(metric = metric_col, n_bins = n_bins) %>%
    relocate(metric, n_bins, bin)
}

#' Nested binning: metric quartile/decile within each existing B_obs quartile.
bin_ratio_table_conditional <- function(
    df, metric_col, n_bins = 4L, ratio_col = "pred_obs_ratio",
    quartile_col = "B_obs_quartile") {
  df %>%
    group_by(.data[[quartile_col]]) %>%
    group_modify(function(chunk, key) {
      bin_ratio_table(chunk, metric_col, n_bins = n_bins, ratio_col = ratio_col)
    }) %>%
    ungroup() %>%
    rename(b_obs_quartile = !!quartile_col)
}

# ---------------------------------------------------------------------------
# Confound-check step 3: data-quality cross-check
# ---------------------------------------------------------------------------

#' Cross-tabulate bottom-decile size_CV hauls against a known problem grouping
#' (year or ICES rectangle) already flagged in the pipeline audits.
#'
#' `group_col` defaults to "year"; `flagged_values` defaults to the 1998 /
#' 2013-2014 dropout spikes and the 1985-only discovery-file era already
#' documented in run_h1_dropout_diagnosis.R / H1_dropout_diagnosis.md. Pass
#' group_col = "stat_rec" with flagged_values = known high-dropout rectangles
#' for the rectangle-level variant of the same check.
dataquality_crosstab <- function(
    df, size_cv_col = "size_CV", group_col = "year",
    flagged_values = c(1985L, 1998L, 2013L, 2014L), decile = 1L) {
  df <- df %>% filter(is.finite(.data[[size_cv_col]]))
  df$.cv_decile <- ntile(df[[size_cv_col]], 10L)
  low <- df %>% filter(.cv_decile == decile)

  by_group <- df %>%
    count(.data[[group_col]], name = "n_total_hauls") %>%
    left_join(
      low %>% count(.data[[group_col]], name = "n_low_cv_hauls"),
      by = group_col
    ) %>%
    mutate(
      n_low_cv_hauls = coalesce(n_low_cv_hauls, 0L),
      pct_of_group_in_low_cv_decile = round(100 * n_low_cv_hauls / n_total_hauls, 1),
      flagged = .data[[group_col]] %in% flagged_values
    ) %>%
    arrange(.data[[group_col]])

  flagged_summary <- tibble(
    flagged = c(TRUE, FALSE),
    n_low_cv_hauls = c(
      sum(low[[group_col]] %in% flagged_values),
      sum(!low[[group_col]] %in% flagged_values)
    )
  ) %>%
    mutate(pct_of_low_cv_total = round(100 * n_low_cv_hauls / nrow(low), 1))

  list(
    by_group = by_group,
    flagged_summary = flagged_summary,
    n_low_cv_total = nrow(low),
    n_hauls_total = nrow(df),
    group_col = group_col
  )
}

# ---------------------------------------------------------------------------
# Confound-check step 4: taxonomic breakdown (descriptive only)
# ---------------------------------------------------------------------------

#' Distribution of dominant-species identity in the top/bottom decile of a metric.
#'
#' Descriptive only — per the briefing this must not be used to construct a
#' species-based flag or filter.
taxonomic_breakdown <- function(
    df, extreme_col, extreme = c("top", "bottom"), species_col = "dominant_species") {
  extreme <- match.arg(extreme)
  df <- df %>% filter(is.finite(.data[[extreme_col]]))
  df$.decile <- ntile(df[[extreme_col]], 10L)
  target_decile <- if (extreme == "top") 10L else 1L

  df %>%
    filter(.decile == target_decile) %>%
    count(.data[[species_col]], name = "n_hauls") %>%
    arrange(desc(n_hauls)) %>%
    mutate(
      pct_of_decile = round(100 * n_hauls / sum(n_hauls), 1),
      metric = extreme_col,
      extreme = extreme
    ) %>%
    relocate(metric, extreme)
}
