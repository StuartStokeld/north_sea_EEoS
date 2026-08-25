# H2 / H3 multiplicity correction (Part A)
#
# PURPOSE:
#   Bonferroni family-wise correction on four phase-specific slopes, run as
#   TWO SEPARATE families (each m = 4; α_adj = 0.05/4 = 0.0125):
#     - H2: FP_between × phase_v2 from wb_car_v2 (CAR)
#     - H3: FP_within  × phase_v2 from wb_car_v2 (CAR)
#   H2 and H3 are FWER-corrected as separate families because they test
#   distinct, independently motivated disturbance pathways (spatial fishing
#   pressure vs temporal fishing pressure), not because of any data-driven
#   reason to keep them apart. They are NOT pooled into a joint m = 8 family.
#
# Part B (CAR rectangle sub-sampling) lives in:
#   pipeline/run_h2h3_rectangle_subsampling_car_refit.R
# Companion RE rectangle sub-sampling:
#   pipeline/run_h2h3_rectangle_subsampling_refit.R
# OLS residual-proxy Part B in this file is superseded.
#
# Prerequisite: outputs/phase_v2_fp_slopes_by_phase.csv
#   (+ phase_v2_proportional_effects_H{2,3}.csv).
#
# Run: Rscript --vanilla pipeline/run_h2h3_h2_multiplicity_subsampling.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_h2_multiplicity_subsampling.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root

# ---------------------------------------------------------------------------
# Paths / settings
# ---------------------------------------------------------------------------
path_models_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_slopes <- file.path(project_root, "outputs", "phase_v2_fp_slopes_by_phase.csv")
path_prop_h2 <- file.path(project_root, "outputs", "phase_v2_proportional_effects_H2.csv")
path_prop_h3 <- file.path(project_root, "outputs", "phase_v2_proportional_effects_H3.csv")

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_mult_h2_csv <- file.path(
  project_root, "outputs", "h2_multiplicity_correction_by_phase.csv"
)
path_out_mult_h2_md <- file.path(
  project_root, "outputs", "h2_multiplicity_correction_summary.md"
)
path_out_mult_h3_csv <- file.path(
  project_root, "outputs", "h3_multiplicity_correction_by_phase.csv"
)
path_out_mult_h3_md <- file.path(
  project_root, "outputs", "h3_multiplicity_correction_summary.md"
)
path_out_iter <- file.path(
  project_root, "outputs", "h2h3_residual_subsampling_iterations.csv"
)
path_out_sum_h2 <- file.path(
  project_root, "outputs", "h2h3_residual_subsampling_summary_H2.csv"
)
path_out_sum_h3 <- file.path(
  project_root, "outputs", "h2h3_residual_subsampling_summary_H3.csv"
)
path_out_sum_md <- file.path(
  project_root, "outputs", "h2h3_residual_subsampling_summary.md"
)
path_out_fig_h2 <- file.path(
  fig_dir, "h2_residual_subsampling_contrast_histograms.png"
)
path_out_fig_h3 <- file.path(
  fig_dir, "h3_residual_subsampling_contrast_histograms.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "h2_multiplicity_subsampling_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "h2_multiplicity_subsampling_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "h2h3_residual_subsampling_objects.rds"
)

PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
ALPHA <- 0.05
M_FAMILY <- 4L
ALPHA_ADJ <- ALPHA / M_FAMILY

N_ITER <- as.integer(Sys.getenv("N_ITER", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
KEEP_FRAC <- as.numeric(Sys.getenv("KEEP_FRAC", unset = "0.80"))
if (is.na(KEEP_FRAC) || KEEP_FRAC <= 0 || KEEP_FRAC >= 1) {
  stop("KEEP_FRAC must be in (0, 1); got ", Sys.getenv("KEEP_FRAC"))
}
DROP_FRAC <- 1 - KEEP_FRAC

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, sprintf("%.2e", p), sprintf("%.4f", p))
  )
}

yn <- function(x) ifelse(isTRUE(x), "Y", "N")

pct_gap_change_from_slope <- function(slope, dx, r0) {
  if (!is.finite(slope) || !is.finite(dx) || !is.finite(r0) || r0 <= 0) {
    return(NA_real_)
  }
  gap0 <- 1 - r0
  if (!is.finite(gap0) || abs(gap0) < .Machine$double.eps) return(NA_real_)
  r1 <- r0 * exp(slope * dx)
  gap1 <- 1 - r1
  ((gap0 - gap1) / gap0) * 100
}

# OLS IQR contrast on haul-level EEoS residuals (no mixed-model refit).
# H2: residual ~ FP_between; dx = IQR(FP_between) across unique rectangles.
# H3: residual ~ FP_within;  dx = IQR(FP_within) among hauls in the subset.
ols_iqr_contrast <- function(dat_sub, fp_col) {
  empty <- list(
    ok = FALSE, slope = NA_real_, dx = NA_real_, r0 = NA_real_,
    contrast = NA_real_, p_value = NA_real_, n_hauls = 0L, n_rect = 0L
  )
  if (nrow(dat_sub) < 10L) return(empty)
  if (!fp_col %in% names(dat_sub)) return(empty)

  n_rect <- dplyr::n_distinct(dat_sub$stat_rec)
  if (identical(fp_col, "FP_between")) {
    rect <- dat_sub %>% distinct(stat_rec, FP_between)
    if (nrow(rect) < 5L) return(empty)
    dx <- as.numeric(IQR(rect$FP_between, na.rm = TRUE))
  } else if (identical(fp_col, "FP_within")) {
    dx <- as.numeric(IQR(dat_sub$FP_within, na.rm = TRUE))
  } else {
    stop("fp_col must be FP_between or FP_within")
  }
  if (!is.finite(dx) || dx <= 0) return(empty)

  r0 <- exp(stats::median(dat_sub$residual, na.rm = TRUE))
  fml <- stats::as.formula(paste("residual ~", fp_col))
  fit <- tryCatch(
    stats::lm(fml, data = dat_sub),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(empty)
  sm <- summary(fit)$coefficients
  if (!fp_col %in% rownames(sm)) return(empty)

  slope <- as.numeric(sm[fp_col, "Estimate"])
  p_value <- as.numeric(sm[fp_col, "Pr(>|t|)"])
  contrast <- pct_gap_change_from_slope(slope, dx, r0)
  list(
    ok = is.finite(contrast) && is.finite(p_value),
    slope = slope,
    dx = dx,
    r0 = r0,
    contrast = contrast,
    p_value = p_value,
    n_hauls = nrow(dat_sub),
    n_rect = as.integer(n_rect)
  )
}

resolve_h2_model_id <- function(slopes) {
  if (any(slopes$model_id == "wb_car_v2" & slopes$component == "FP_between")) {
    "wb_car_v2"
  } else {
    "wb_primary_v2"
  }
}

resolve_h3_model_id <- function(slopes) {
  if (any(slopes$model_id == "wb_car_v2" & slopes$component == "FP_within")) {
    "wb_car_v2"
  } else {
    "wb_primary_v2"
  }
}

build_multiplicity_table <- function(phase_slopes, prop_effects) {
  raw_p <- phase_slopes$p_value
  stopifnot(length(raw_p) == M_FAMILY)
  bonf <- pmin(1, raw_p * M_FAMILY)
  prop <- prop_effects %>%
    transmute(
      phase = as.character(phase),
      pct_gap_change_iqr,
      gap_direction_iqr
    )
  tibble::tibble(
    phase = as.character(phase_slopes$phase),
    fp_slope = phase_slopes$fp_slope,
    fp_slope_se = phase_slopes$fp_slope_se,
    pct_gap_change_iqr = prop$pct_gap_change_iqr[
      match(as.character(phase_slopes$phase), prop$phase)
    ],
    gap_direction_iqr = prop$gap_direction_iqr[
      match(as.character(phase_slopes$phase), prop$phase)
    ],
    raw_p = raw_p,
    bonferroni_p = bonf,
    sig_alpha_0.05 = raw_p < ALPHA,
    sig_alpha_adj_0.0125 = bonf < ALPHA
  )
}

write_multiplicity_summary_md <- function(mult, path_md, hypothesis_label,
                                          component_label, model_id,
                                          family_note, paragraph) {
  md <- c(
    paste0(
      "# ", hypothesis_label,
      " multiplicity correction — four phase-specific contrasts"
    ),
    "",
    paste0(
      "Source model: `", model_id, "` (`", component_label,
      "` × `phase_v2` slopes). ",
      "Family size m = 4; α = 0.05; Bonferroni α_adj = 0.05/4 = 0.0125. ",
      family_note
    ),
    "",
    "## Results table",
    "",
    paste0(
      "| phase | slope | SE | gap change (IQR) | raw p | Bonferroni p | ",
      "sig α=0.05 | sig α_adj=0.0125 |"
    ),
    "|---|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(mult)), function(i) {
      r <- mult[i, ]
      sprintf(
        "| %s | %+.4f | %.4f | %s %.2f%% | %s | %s | %s | %s |",
        r$phase, r$fp_slope, r$fp_slope_se,
        r$gap_direction_iqr, abs(r$pct_gap_change_iqr),
        fmt_p(r$raw_p), fmt_p(r$bonferroni_p),
        yn(r$sig_alpha_0.05), yn(r$sig_alpha_adj_0.0125)
      )
    }, character(1)),
    "",
    "## Interpretation",
    "",
    paragraph,
    "",
    "### Notes",
    "",
    "- Bonferroni-adjusted p = min(1, p_raw × 4).",
    "- Significance under α_adj = 0.0125 is equivalent to Bonferroni p < 0.05.",
    paste0(
      "- H2 and H3 are corrected as **separate** m = 4 families (spatial vs ",
      "temporal disturbance pathways); not pooled to m = 8."
    ),
    "- Additive check only; does not replace spatial-null / permutation / KNN work.",
    ""
  )
  writeLines(md, path_md)
}

summarise_subsamples <- function(ok_df, ref_df, phases) {
  lapply(phases, function(ph) {
    sub <- ok_df %>% filter(as.character(phase) == ph)
    ref <- ref_df %>% filter(as.character(phase) == ph)
    stopifnot(nrow(ref) == 1L)
    ref_sign <- sign(ref$full_sample_contrast[[1]])
    if (ref_sign == 0) ref_sign <- 1
    ref_sig <- isTRUE(ref$full_sample_sig[[1]])
    sub_sign <- sign(sub$contrast)
    sub_sign[sub_sign == 0] <- 1
    same_sign <- sub_sign == ref_sign
    same_sig_call <- (sub$p_value < ALPHA) == ref_sig
    match_sign_and_sig <- same_sign & same_sig_call
    tibble::tibble(
      phase = ph,
      n_ok = nrow(sub),
      full_sample_contrast = ref$full_sample_contrast[[1]],
      full_sample_gap_direction = ref$full_sample_gap_direction[[1]],
      full_sample_p = ref$full_sample_p[[1]],
      full_sample_sig = ref_sig,
      subsample_median = stats::median(sub$contrast, na.rm = TRUE),
      subsample_p025 = as.numeric(stats::quantile(sub$contrast, 0.025, na.rm = TRUE)),
      subsample_p975 = as.numeric(stats::quantile(sub$contrast, 0.975, na.rm = TRUE)),
      pct_sign_match = 100 * mean(same_sign),
      pct_sign_and_sig_match = 100 * mean(match_sign_and_sig)
    )
  }) %>%
    bind_rows() %>%
    mutate(phase = factor(phase, levels = phases))
}

make_hist <- function(ok_df, ref_df, title, subtitle, outfile) {
  obs <- ref_df %>%
    transmute(
      phase = factor(as.character(phase), levels = PHASE_V2),
      contrast = full_sample_contrast
    )
  p <- ggplot(ok_df, aes(x = contrast)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 40, fill = "grey75", colour = "white"
    ) +
    geom_density(colour = "#2166ac", linewidth = 0.8) +
    geom_vline(
      data = obs, aes(xintercept = contrast),
      colour = "#b2182b", linetype = "dashed", linewidth = 0.7
    ) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~phase, ncol = 2, scales = "free_y") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "IQR contrast (% gap change; + = closed, − = widened)",
      y = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(outfile, p, width = 9, height = 7, dpi = 150)
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# H2 / H3 multiplicity + H2/H3 residual sub-sampling — run log")
logmsg("")
logmsg(
  "Part A: Bonferroni correction on two separate m = 4 families — H2 ",
  "(`FP_between` from CAR) and H3 (`FP_within` from CAR). Families are kept ",
  "apart because they test distinct disturbance pathways (spatial vs temporal), ",
  "not pooled to m = 8. Part B: keep ",
  sprintf("%.0f%%", 100 * KEEP_FRAC),
  " of rectangles per phase × ", N_ITER,
  " iterations; recompute IQR contrast via OLS on existing EEoS residuals ",
  "(no mixed-model refit)."
)
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf(
  "N_ITER = %d; SEED = %d; KEEP_FRAC = %.2f (drop %.0f%%)",
  N_ITER, SEED, KEEP_FRAC, 100 * DROP_FRAC
))
logmsg(sprintf(
  "Family size m = %d; α = %.2f; Bonferroni α_adj = %.4f",
  M_FAMILY, ALPHA, ALPHA_ADJ
))

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
for (p in c(path_models_v2, path_slopes, path_prop_h2, path_prop_h3)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

obj <- readRDS(path_models_v2)
dat <- obj$data
stopifnot(all(
  c("residual", "FP_between", "FP_within", "phase_v2", "stat_rec") %in% names(dat)
))
dat$stat_rec <- as.character(dat$stat_rec)
dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
if (anyNA(dat$phase_v2)) stop("NA in phase_v2 after relevel.")

slopes <- read_csv(path_slopes, show_col_types = FALSE)
prop_h2 <- read_csv(path_prop_h2, show_col_types = FALSE)
prop_h3 <- read_csv(path_prop_h3, show_col_types = FALSE)

h2_model_id <- resolve_h2_model_id(slopes)
h3_model_id <- resolve_h3_model_id(slopes)
h2_slopes <- slopes %>%
  filter(model_id == h2_model_id, component == "FP_between") %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
h3_slopes <- slopes %>%
  filter(model_id == h3_model_id, component == "FP_within") %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
if (nrow(h2_slopes) != 4L || nrow(h3_slopes) != 4L) {
  stop("Expected 4 H2 and 4 H3 phase slopes.")
}

prop_h2 <- prop_h2 %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
prop_h3 <- prop_h3 %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)

n_hauls <- nrow(dat)
n_rect_all <- dplyr::n_distinct(dat$stat_rec)
logmsg("Loaded: ", path_models_v2)
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  n_hauls, n_rect_all, min(dat$year), max(dat$year)
))
logmsg(sprintf(
  "Part A H2 source: `%s` FP_between slopes from %s",
  h2_model_id, path_slopes
))
logmsg(
  "Part A H3 source: `", h3_model_id, "` FP_within slopes — own independent ",
  "m = 4 Bonferroni family (spatial vs temporal pathways kept separate; ",
  "not pooled with H2)."
)
logmsg(
  "Part B statistic: OLS `residual ~ FP_*` × phase IQR → pct gap change ",
  "(same transform as proportional-effects tables)."
)

# ---------------------------------------------------------------------------
# Part A — Bonferroni (H2 and H3 as separate m = 4 families)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Part A — Bonferroni correction (separate H2 and H3 families)")
logmsg(
  "Each family = four phase-specific slopes (m = 4). H2 and H3 are corrected ",
  "separately because they test distinct, independently motivated disturbance ",
  "pathways (spatial fishing pressure vs temporal fishing pressure). ",
  "Existing spatial-null / permutation / KNN checks are unaffected."
)

# --- H2 family ---
mult_h2 <- build_multiplicity_table(h2_slopes, prop_h2)
write_csv(mult_h2, path_out_mult_h2_csv)
logmsg("Saved: ", path_out_mult_h2_csv)

logmsg("")
logmsg("### H2 (`FP_between`)")
logmsg(
  "| phase | gap change | raw p | Bonferroni p | sig α=0.05 | sig α_adj=0.0125 |"
)
logmsg("|---|---|---|---|---|---|")
for (i in seq_len(nrow(mult_h2))) {
  r <- mult_h2[i, ]
  gap_lab <- sprintf(
    "%s %.1f%%",
    r$gap_direction_iqr,
    abs(r$pct_gap_change_iqr)
  )
  logmsg(sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    r$phase, gap_lab, fmt_p(r$raw_p), fmt_p(r$bonferroni_p),
    yn(r$sig_alpha_0.05), yn(r$sig_alpha_adj_0.0125)
  ))
}

sig_raw_h2 <- mult_h2$phase[mult_h2$sig_alpha_0.05]
sig_bonf_h2 <- mult_h2$phase[mult_h2$sig_alpha_adj_0.0125]
paragraph_h2 <- paste0(
  "Family of four phase-specific H2 (`FP_between`) contrasts from `",
  h2_model_id, "` (α = 0.05; Bonferroni α_adj = 0.0125). ",
  "Uncorrected significant phases: ",
  if (length(sig_raw_h2)) paste(sig_raw_h2, collapse = ", ") else "none",
  ". After Bonferroni: ",
  if (length(sig_bonf_h2)) paste(sig_bonf_h2, collapse = ", ") else "none",
  ". H3 is corrected separately (see `h3_multiplicity_correction_summary.md`); ",
  "H2 and H3 are not pooled into a joint m = 8 family."
)

write_multiplicity_summary_md(
  mult_h2, path_out_mult_h2_md,
  hypothesis_label = "H2",
  component_label = "FP_between",
  model_id = h2_model_id,
  family_note = paste0(
    "H3 is corrected separately as its own m = 4 family ",
    "(`h3_multiplicity_correction_summary.md`); not pooled with H2."
  ),
  paragraph = paragraph_h2
)
logmsg("Saved: ", path_out_mult_h2_md)
logmsg("")
logmsg(paragraph_h2)

# --- H3 family ---
mult_h3 <- build_multiplicity_table(h3_slopes, prop_h3)
write_csv(mult_h3, path_out_mult_h3_csv)
logmsg("Saved: ", path_out_mult_h3_csv)

logmsg("")
logmsg("### H3 (`FP_within`)")
logmsg(
  "| phase | gap change | raw p | Bonferroni p | sig α=0.05 | sig α_adj=0.0125 |"
)
logmsg("|---|---|---|---|---|---|")
for (i in seq_len(nrow(mult_h3))) {
  r <- mult_h3[i, ]
  gap_lab <- sprintf(
    "%s %.1f%%",
    r$gap_direction_iqr,
    abs(r$pct_gap_change_iqr)
  )
  logmsg(sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    r$phase, gap_lab, fmt_p(r$raw_p), fmt_p(r$bonferroni_p),
    yn(r$sig_alpha_0.05), yn(r$sig_alpha_adj_0.0125)
  ))
}

sig_raw_h3 <- mult_h3$phase[mult_h3$sig_alpha_0.05]
sig_bonf_h3 <- mult_h3$phase[mult_h3$sig_alpha_adj_0.0125]
paragraph_h3 <- paste0(
  "Family of four phase-specific H3 (`FP_within`) contrasts from ",
  "`", h3_model_id, "` (α = 0.05; Bonferroni α_adj = 0.0125). ",
  "Uncorrected significant phases: ",
  if (length(sig_raw_h3)) paste(sig_raw_h3, collapse = ", ") else "none",
  ". After Bonferroni: ",
  if (length(sig_bonf_h3)) paste(sig_bonf_h3, collapse = ", ") else "none",
  ". H2 is corrected separately (see `h2_multiplicity_correction_summary.md`); ",
  "H2 and H3 are FWER-corrected as separate families because they test ",
  "distinct, independently motivated disturbance pathways (spatial vs temporal ",
  "fishing pressure), not because of any data-driven reason to keep them apart."
)

write_multiplicity_summary_md(
  mult_h3, path_out_mult_h3_md,
  hypothesis_label = "H3",
  component_label = "FP_within",
  model_id = h3_model_id,
  family_note = paste0(
    "H2 is corrected separately as its own m = 4 family ",
    "(`h2_multiplicity_correction_summary.md`); not pooled with H3. ",
    "Separate families reflect distinct spatial vs temporal disturbance ",
    "pathways, not a data-driven split."
  ),
  paragraph = paragraph_h3
)
logmsg("Saved: ", path_out_mult_h3_md)
logmsg("")
logmsg(paragraph_h3)

# ---------------------------------------------------------------------------
# Part B — residual sub-sampling (H2 + H3, no refit)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Part B — Rectangle residual sub-sampling (H2 and H3)")
logmsg(
  "Unit of resampling: rectangles present in each phase. Each iteration keeps ",
  sprintf("%.0f%%", 100 * KEEP_FRAC),
  " without replacement (uniform; no stratification), then recomputes the ",
  "IQR contrast via OLS on haul-level EEoS residuals from those rectangles. ",
  "Full-sample point estimates / significance calls come from the presented ",
  "model slopes (H2: `", h2_model_id, "`; H3: `", h3_model_id, "`) and ",
  "proportional-effects IQR gap-change. Subsample p-values are OLS ",
  "(derived-statistic check, not a mixed-model refit)."
)

ref_h2 <- tibble::tibble(
  phase = as.character(prop_h2$phase),
  full_sample_contrast = prop_h2$pct_gap_change_iqr,
  full_sample_gap_direction = prop_h2$gap_direction_iqr,
  full_sample_p = h2_slopes$p_value[match(as.character(prop_h2$phase), as.character(h2_slopes$phase))],
  full_sample_sig = h2_slopes$p_value[match(as.character(prop_h2$phase), as.character(h2_slopes$phase))] < ALPHA
)
ref_h3 <- tibble::tibble(
  phase = as.character(prop_h3$phase),
  full_sample_contrast = prop_h3$pct_gap_change_iqr,
  full_sample_gap_direction = prop_h3$gap_direction_iqr,
  full_sample_p = h3_slopes$p_value[match(as.character(prop_h3$phase), as.character(h3_slopes$phase))],
  full_sample_sig = h3_slopes$p_value[match(as.character(prop_h3$phase), as.character(h3_slopes$phase))] < ALPHA
)

logmsg("")
logmsg("### Full-sample reference (headline IQR gap-change)")
for (i in seq_len(nrow(ref_h2))) {
  r <- ref_h2[i, ]
  logmsg(sprintf(
    "  H2 %s: contrast=%+.2f%% (%s) p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}
for (i in seq_len(nrow(ref_h3))) {
  r <- ref_h3[i, ]
  logmsg(sprintf(
    "  H3 %s: contrast=%+.2f%% (%s) p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}

rects_by_phase <- lapply(PHASE_V2, function(ph) {
  sort(unique(dat$stat_rec[as.character(dat$phase_v2) == ph]))
})
names(rects_by_phase) <- PHASE_V2
for (ph in PHASE_V2) {
  logmsg(sprintf(
    "  phase %s: %d rectangles available for sampling",
    ph, length(rects_by_phase[[ph]])
  ))
}

t_start <- proc.time()[[3]]
set.seed(SEED)
seeds <- SEED + seq_len(N_ITER)

iter_rows <- vector("list", N_ITER * length(PHASE_V2) * 2L)
k <- 0L

for (i in seq_len(N_ITER)) {
  if (i == 1L || i %% 100L == 0L || i == N_ITER) {
    logmsg(sprintf("  iteration %d / %d ...", i, N_ITER))
  }
  set.seed(seeds[[i]])
  for (ph in PHASE_V2) {
    all_r <- rects_by_phase[[ph]]
    n_rect <- length(all_r)
    n_keep <- max(2L, as.integer(round(n_rect * KEEP_FRAC)))
    n_drop <- n_rect - n_keep
    keep <- sample(all_r, size = n_keep, replace = FALSE)
    dat_sub <- dat %>%
      filter(as.character(phase_v2) == ph, stat_rec %in% keep)

    for (hyp in c("H2", "H3")) {
      fp_col <- if (hyp == "H2") "FP_between" else "FP_within"
      res <- ols_iqr_contrast(dat_sub, fp_col)
      k <- k + 1L
      iter_rows[[k]] <- data.frame(
        iter = i,
        seed = seeds[[i]],
        phase = ph,
        hypothesis = hyp,
        n_rect_available = n_rect,
        n_keep = n_keep,
        n_drop = n_drop,
        n_hauls = res$n_hauls,
        n_rect_in_subset = res$n_rect,
        ok = res$ok,
        slope = res$slope,
        dx_iqr = res$dx,
        baseline_ratio = res$r0,
        contrast = res$contrast,
        p_value = res$p_value,
        stringsAsFactors = FALSE
      )
    }
  }
}

iter_df <- dplyr::bind_rows(iter_rows)
runtime_sec <- proc.time()[[3]] - t_start

n_fail <- iter_df %>%
  filter(!ok) %>%
  summarise(n = dplyr::n()) %>%
  pull(n)
logmsg(sprintf(
  "Completed %d iteration×phase×hypothesis rows (%d not-ok); runtime %.1f sec.",
  nrow(iter_df), n_fail, runtime_sec
))

write_csv(iter_df, path_out_iter)
logmsg("Saved: ", path_out_iter)

ok_h2 <- iter_df %>%
  filter(hypothesis == "H2", ok, is.finite(contrast), is.finite(p_value)) %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2))
ok_h3 <- iter_df %>%
  filter(hypothesis == "H3", ok, is.finite(contrast), is.finite(p_value)) %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2))

sum_h2 <- summarise_subsamples(ok_h2, ref_h2, PHASE_V2)
sum_h3 <- summarise_subsamples(ok_h3, ref_h3, PHASE_V2)
write_csv(sum_h2, path_out_sum_h2)
write_csv(sum_h3, path_out_sum_h3)
logmsg("Saved: ", path_out_sum_h2)
logmsg("Saved: ", path_out_sum_h3)

logmsg("")
logmsg("### H2 summary")
logmsg(
  "| phase | full contrast | subsample median [2.5, 97.5] | % sign match | % sign+sig match | n_ok |"
)
logmsg("|---|---|---|---|---|---|")
for (i in seq_len(nrow(sum_h2))) {
  r <- sum_h2[i, ]
  logmsg(sprintf(
    "| %s | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %d |",
    as.character(r$phase), r$full_sample_contrast,
    r$subsample_median, r$subsample_p025, r$subsample_p975,
    r$pct_sign_match, r$pct_sign_and_sig_match, r$n_ok
  ))
}

logmsg("")
logmsg("### H3 summary")
logmsg(
  "| phase | full contrast | subsample median [2.5, 97.5] | % sign match | % sign+sig match | n_ok |"
)
logmsg("|---|---|---|---|---|---|")
for (i in seq_len(nrow(sum_h3))) {
  r <- sum_h3[i, ]
  logmsg(sprintf(
    "| %s | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %d |",
    as.character(r$phase), r$full_sample_contrast,
    r$subsample_median, r$subsample_p025, r$subsample_p975,
    r$pct_sign_match, r$pct_sign_and_sig_match, r$n_ok
  ))
}

subtitle_common <- sprintf(
  "%d iterations; keep %.0f%% of phase rectangles; dashed red = full-sample headline IQR contrast; seed = %d",
  N_ITER, 100 * KEEP_FRAC, SEED
)
make_hist(
  ok_h2, ref_h2,
  "H2 residual sub-sampling: IQR gap-change contrast by phase",
  subtitle_common,
  path_out_fig_h2
)
make_hist(
  ok_h3, ref_h3,
  "H3 residual sub-sampling: IQR gap-change contrast by phase",
  subtitle_common,
  path_out_fig_h3
)
logmsg("Saved: ", path_out_fig_h2)
logmsg("Saved: ", path_out_fig_h3)

md_table <- function(sum_df) {
  c(
    "| phase | full-sample contrast | direction | full p | full sig | subsample median | p2.5 | p97.5 | % sign match | % sign+sig match | n_ok |",
    "|---|---|---|---|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(sum_df)), function(i) {
      r <- sum_df[i, ]
      sprintf(
        "| %s | %+.2f%% | %s | %s | %s | %+.2f%% | %+.2f%% | %+.2f%% | %.1f%% | %.1f%% | %d |",
        as.character(r$phase), r$full_sample_contrast, r$full_sample_gap_direction,
        fmt_p(r$full_sample_p), yn(r$full_sample_sig),
        r$subsample_median, r$subsample_p025, r$subsample_p975,
        r$pct_sign_match, r$pct_sign_and_sig_match, r$n_ok
      )
    }, character(1))
  )
}

sub_md <- c(
  "# H2/H3 residual sub-sampling sensitivity",
  "",
  paste0(
    "Each iteration keeps ", sprintf("%.0f%%", 100 * KEEP_FRAC),
    " of rectangles present in a phase (uniform, without replacement), then ",
    "**recomputes** the headline IQR contrast (25th-vs-75th FP shift → % gap change) ",
    "via OLS on haul-level EEoS residuals from those rectangles. ",
    "**No mixed-model refit.**"
  ),
  "",
  sprintf("- Iterations: %d", N_ITER),
  sprintf("- Seed: %d", SEED),
  sprintf("- Keep fraction: %.2f (drop %.0f%%)", KEEP_FRAC, 100 * DROP_FRAC),
  sprintf("- Runtime: %.1f sec", runtime_sec),
  "",
  "## H2 (`FP_between`)",
  "",
  paste0(
    "Full-sample reference: `", h2_model_id,
    "` slopes × IQR(`FP_between`) gap-change from ",
    "`phase_v2_proportional_effects_H2.csv`."
  ),
  "",
  md_table(sum_h2),
  "",
  "## H3 (`FP_within`)",
  "",
  paste0(
    "Full-sample reference: `", h3_model_id, "` slopes × phase IQR(`FP_within`) ",
    "from `phase_v2_proportional_effects_H3.csv`. Part B match rule uses ",
    "nominal α = 0.05 on raw model slopes (stability check). Part A FWER ",
    "for H3 is in `h3_multiplicity_correction_summary.md` (separate m = 4 family)."
  ),
  "",
  md_table(sum_h3),
  "",
  "## Figures",
  "",
  paste0("- H2: `", path_out_fig_h2, "`"),
  paste0("- H3: `", path_out_fig_h3, "`"),
  "",
  "### Notes",
  "",
  "- Match rule: same sign as full-sample contrast **and** same significance call (p < 0.05 vs not) as the presented model slope.",
  "- Subsample p-values are from OLS on EEoS residuals — a derived-statistic check, not a re-estimation of the mixed model.",
  "- Prior Part B (glmmTMB refit after rectangle drop) is superseded; archive files may remain as `outputs/h2_rectangle_subsampling_*`.",
  ""
)
writeLines(sub_md, path_out_sum_md)
logmsg("Saved: ", path_out_sum_md)

# ---------------------------------------------------------------------------
# Persist + finish
# ---------------------------------------------------------------------------
saveRDS(
  list(
    multiplicity_h2 = mult_h2,
    multiplicity_h3 = mult_h3,
    multiplicity_paragraph_h2 = paragraph_h2,
    multiplicity_paragraph_h3 = paragraph_h3,
    h2_model_id = h2_model_id,
    ref_h2 = ref_h2,
    ref_h3 = ref_h3,
    iterations = iter_df,
    summary_h2 = sum_h2,
    summary_h3 = sum_h3,
    settings = list(
      n_iter = N_ITER,
      seed = SEED,
      keep_frac = KEEP_FRAC,
      drop_frac = DROP_FRAC,
      alpha = ALPHA,
      alpha_adj = ALPHA_ADJ,
      m_family = M_FAMILY,
      m_family_note = paste0(
        "H2 and H3 each corrected as separate m = 4 families ",
        "(spatial vs temporal pathways); not pooled to m = 8"
      ),
      phases = PHASE_V2,
      method = "ols_iqr_contrast_on_eeos_residuals_no_mixed_refit"
    ),
    runtime_sec = runtime_sec
  ),
  path_out_rds
)
logmsg("Saved: ", path_out_rds)

writeLines(run_log, path_out_run_log)
logmsg("Saved: ", path_out_run_log)
cat("=== H2 multiplicity + H2/H3 residual sub-sampling complete. ===\n")
