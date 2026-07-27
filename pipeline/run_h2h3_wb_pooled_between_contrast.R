# H2 contrast: pooled FP_between (no phase interaction), CAR
#
# PURPOSE: additional discussion contrast — what H2 looks like if the
# between-rectangle fishing-pressure effect is assumed stable over time.
# NOT a replacement for the primary phase-specific H2 specification.
#
# Fit (CAR): residual ~ FP_between + FP_between * phase  [NO]
#            residual ~ FP_between + FP_within * phase + adjacency(1 | stat_rec)
# Only change vs primary within-between CAR: drop FP_between × phase.
# FP_within × phase unchanged. Same decomposed data / 158-rectangle panel.
#
# NO FIGURES — numbers and tables only.
#
# Run: Rscript --vanilla pipeline/run_h2h3_wb_pooled_between_contrast.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_wb_pooled_between_contrast.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop("spaMM required. Run with: Rscript --vanilla pipeline/run_h2h3_wb_pooled_between_contrast.R")
}
suppressPackageStartupMessages(library(spaMM))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))

path_wb <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_slopes <- file.path(project_root, "outputs", "h2h3_wb_fp_slopes_by_phase.csv")
path_round2 <- file.path(project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds")
stopifnot(file.exists(path_wb), file.exists(path_slopes), file.exists(path_round2))

path_out_pooled <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_coef.csv")
path_out_compare <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_comparison.csv")
path_out_conv <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_convergence.csv")
path_out_models <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_model_objects.rds")
path_out_run_log <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_run_log.md")
path_out_session <- file.path(project_root, "outputs", "h2h3_wb_pooled_between_sessionInfo.txt")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# H2 contrast: pooled FP_between (no phase interaction) — run log")
logmsg("")
logmsg(
  "ADDITIONAL CONTRAST ONLY — not a replacement for the primary phase-specific H2 ",
  "model (`FP_between * phase`). Quantifies a single time-stable between-rectangle ",
  "fishing-pressure slope for discussion. No figures produced."
)

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg("spaMM ", as.character(utils::packageVersion("spaMM")))

# ---------------------------------------------------------------------------
# Data reuse
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data and specification")

wb <- readRDS(path_wb)
dat <- wb$data
round2 <- readRDS(path_round2)
adjMatrix <- round2$adjMatrix
stopifnot(!is.null(adjMatrix))
stopifnot(all(c("FP_between", "FP_within", "phase", "residual", "stat_rec") %in% names(dat)))

dat_car <- dat
dat_car$stat_rec <- factor(as.character(dat_car$stat_rec), levels = rownames(adjMatrix))

formula_pooled <- residual ~ FP_between + FP_within * phase + adjacency(1 | stat_rec)
formula_primary_car <- residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)

logmsg(sprintf(
  "Data reused from h2h3_wb_model_objects.rds: %d hauls, %d rectangles.",
  nrow(dat_car), dplyr::n_distinct(dat_car$stat_rec)
))
logmsg("adjMatrix reused from Round 2 RDS (same queen adjacency as prior CAR fits).")
logmsg("")
logmsg("Contrast formula (this run):")
logmsg("  residual ~ FP_between + FP_within * phase + adjacency(1 | stat_rec)  [CAR REML]")
logmsg("Primary within-between CAR (unchanged; for reference only):")
logmsg("  residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)")
logmsg("")
logmsg("Confirmed changes: ONLY FP_between × phase dropped (FP_between enters as a single pooled main effect).")
logmsg("Confirmed unchanged: FP_within * phase; phase factor; CAR adjacency; no biomass; same panel.")

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Model fit")
time_fit <- system.time({
  fit_pooled <- spaMM::fitme(
    formula_pooled, data = dat_car, adjMatrix = adjMatrix, method = "REML"
  )
})
logmsg(sprintf("Fit time: %.2f sec.", time_fit["elapsed"]))

saveRDS(
  list(
    fit_pooled_between = fit_pooled,
    data = dat_car,
    adjMatrix = adjMatrix,
    formula = formula_pooled,
    note = "Contrast only; primary remains FP_between * phase + FP_within * phase + CAR"
  ),
  path_out_models
)
logmsg("Saved model object: ", path_out_models)

# ---------------------------------------------------------------------------
# Convergence
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Convergence")

rhorange <- car_rho_admissible_range(adjMatrix)
conv <- spamm_car_convergence_report(
  fit_pooled, "wb_pooled_between_car", formula_pooled, dat_car, adjMatrix, rhorange
)
write_csv(conv, path_out_conv)
logmsg(sprintf(
  paste0(
    "warnings=%d (%s); robust_to_starting_value=%s; fitted_rho=%.6f; ",
    "any non-finite FE SE=%s"
  ),
  conv$n_warnings_during_fit, conv$warnings_text, conv$robust_to_starting_value,
  conv$fitted_rho, conv$any_na_or_infinite_fixed_effect_se
))
logmsg("Saved: ", path_out_conv)

# ---------------------------------------------------------------------------
# 1. Pooled FP_between coefficient
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 1. Pooled FP_between coefficient")

fe <- tidy_fixed_effects_spamm(fit_pooled, "wb_pooled_between_car")
pooled_row <- fe %>% filter(term == "FP_between")
if (nrow(pooled_row) != 1L) {
  stop("Expected exactly one FP_between coefficient; found: ", paste(fe$term, collapse = ", "))
}

pooled_out <- pooled_row %>%
  mutate(
    ci_lo = estimate - 1.96 * std_error,
    ci_hi = estimate + 1.96 * std_error,
    term_plain = "Pooled FP_between (no phase interaction)",
    hypothesis = "H2 contrast — time-stable between-rectangle FP slope",
    model_note = "residual ~ FP_between + FP_within * phase + adjacency(1|stat_rec)"
  ) %>%
  select(
    model_id, term, term_plain, hypothesis,
    estimate, std_error, ci_lo, ci_hi, statistic, p_value, model_note
  )
write_csv(pooled_out, path_out_pooled)

logmsg(sprintf(
  "FP_between (pooled): estimate=%+.6f, SE=%.6f, 95%% CI=[%+.6f, %+.6f], z=%+.3f, p=%.6g",
  pooled_out$estimate, pooled_out$std_error, pooled_out$ci_lo, pooled_out$ci_hi,
  pooled_out$statistic, pooled_out$p_value
))
logmsg("Saved: ", path_out_pooled)

# Confirm FP_within × phase still present
within_terms <- fe$term[grepl("FP_within", fe$term)]
logmsg("FP_within-related terms still in this contrast fit: ", paste(within_terms, collapse = "; "))
between_int <- fe$term[grepl("FP_between", fe$term) & grepl("phase", fe$term)]
if (length(between_int) > 0L) {
  stop("Unexpected FP_between × phase terms found in pooled contrast: ", paste(between_int, collapse = ", "))
}
logmsg("Confirmed: no FP_between × phase interaction terms in the contrast fit.")

# ---------------------------------------------------------------------------
# 2. Side-by-side vs phase-specific CAR FP_between slopes
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 2. Pooled vs phase-specific FP_between (CAR)")

phase_specific <- read_csv(path_slopes, show_col_types = FALSE) %>%
  filter(model_id == "wb_car", component == "FP_between") %>%
  transmute(
    row_type = "phase_specific",
    phase = as.character(phase),
    estimate = fp_slope,
    std_error = fp_slope_se,
    ci_lo = fp_slope_lo,
    ci_hi = fp_slope_hi,
    p_value = p_value,
    source = "primary within-between CAR (FP_between * phase)"
  )

pooled_cmp <- data.frame(
  row_type = "pooled",
  phase = "all_years_pooled",
  estimate = pooled_out$estimate,
  std_error = pooled_out$std_error,
  ci_lo = pooled_out$ci_lo,
  ci_hi = pooled_out$ci_hi,
  p_value = pooled_out$p_value,
  source = "contrast CAR (FP_between main effect only)",
  stringsAsFactors = FALSE
)

compare <- bind_rows(pooled_cmp, phase_specific) %>%
  mutate(
    phase = factor(
      phase,
      levels = c("all_years_pooled", "1985-1988", "1989-2000", "2001-2007", "2008-2015")
    )
  ) %>%
  arrange(phase)

write_csv(compare, path_out_compare)
logmsg("Side-by-side:")
for (i in seq_len(nrow(compare))) {
  r <- compare[i, ]
  logmsg(sprintf(
    "  - %-18s  est=%+.6f  SE=%.6f  CI=[%+.6f, %+.6f]  p=%.6g",
    as.character(r$phase), r$estimate, r$std_error, r$ci_lo, r$ci_hi, r$p_value
  ))
}
logmsg("Saved: ", path_out_compare)

# ---------------------------------------------------------------------------
# 3. One sentence — plain numeric observation only
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## 3. Plain numeric observation (not an interpretation)")

ps <- phase_specific %>% arrange(match(phase, c("1985-1988", "1989-2000", "2001-2007", "2008-2015")))
est_pooled <- pooled_out$estimate
est_phases <- ps$estimate
min_ph <- min(est_phases)
max_ph <- max(est_phases)
# "sits between opposing phase-specific effects" — phases include both signs
signs <- sign(est_phases)
has_opposing <- any(signs < 0) && any(signs > 0)
between_range <- isTRUE(est_pooled > min_ph && est_pooled < max_ph)
# partial cancellation: pooled closer to zero than the largest |phase| effect
max_abs_phase <- max(abs(est_phases))
closer_to_zero <- abs(est_pooled) < max_abs_phase

obs <- sprintf(
  paste0(
    "The pooled FP_between estimate (%+.4f) lies between the most negative phase-specific ",
    "CAR slope (%+.4f in %s) and the most positive (%+.4f in %s)%s; |pooled| (%.4f) is ",
    "smaller than the largest |phase-specific| slope (%.4f), so pooling partially cancels ",
    "the phase-level reversal in magnitude."
  ),
  est_pooled,
  min_ph, ps$phase[which.min(est_phases)],
  max_ph, ps$phase[which.max(est_phases)],
  if (has_opposing && between_range) "" else " [FLAG: check range/sign pattern]",
  abs(est_pooled), max_abs_phase
)
logmsg(obs)

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs (no figures)")
logmsg("- ", path_out_pooled)
logmsg("- ", path_out_compare)
logmsg("- ", path_out_conv)
logmsg("- ", path_out_models)
logmsg("- ", path_out_session)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Pooled FP_between contrast complete (no figures). ===\n")
