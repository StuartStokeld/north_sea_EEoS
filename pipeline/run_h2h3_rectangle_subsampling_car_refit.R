# H2/H3 CAR rectangle sub-sampling sensitivity (canonical shared CAR model)
#
# PURPOSE: refit the presented CAR specification on the same global rectangle
# subsamples used by the companion RE script, and check sign / significance
# stability of the four phase-specific FP_between (H2) and FP_within (H3)
# IQR contrasts.
#
# Model (matches run_h2h3_phase_v2_reporting.R):
#   residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
#   spaMM::fitme(..., method = "REML")
#
# Paired with RE companion: same SEED blocks and rectangle order as
#   pipeline/run_h2h3_rectangle_subsampling_refit.R
#
# IQR convention: fixed full-sample IQR(FP_between) + fixed phase-specific
# IQR(FP_within) + fixed phase baselines (same as phase_v2 reporting).
#
# Prerequisite: outputs/primary_model_v2.rds,
#   outputs/h2h3_feasibility_round2_model_objects.rds$adjMatrix,
#   outputs/phase_v2_proportional_effects_H{2,3}.csv (wb_car_v2)
#
# Run: Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_car_refit.R
# Optional: N_ITER=1000 SEED=42 SKIP_REFIT=1
#   SKIP_REFIT=1 summarises existing iterations CSV (no new CAR fits).

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_rectangle_subsampling_car_refit.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))
source(file.path(script_dir, "R", "h2h3_feasibility_helpers.R"))

# Prefer project-local spaMM
local_libs <- file.path(project_root, ".R_libs")
if (dir.exists(local_libs)) {
  .libPaths(c(normalizePath(local_libs), .libPaths()))
}
if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop(
    "Package 'spaMM' required. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_car_refit.R"
  )
}
suppressPackageStartupMessages(library(spaMM))

# ---------------------------------------------------------------------------
# Paths / settings
# ---------------------------------------------------------------------------
path_models_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_round2 <- file.path(
  project_root, "outputs", "h2h3_feasibility_round2_model_objects.rds"
)
path_prop_h2 <- file.path(
  project_root, "outputs", "phase_v2_proportional_effects_H2.csv"
)
path_prop_h3 <- file.path(
  project_root, "outputs", "phase_v2_proportional_effects_H3.csv"
)
path_slopes <- file.path(
  project_root, "outputs", "phase_v2_fp_slopes_by_phase.csv"
)

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_iter <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_iterations.csv"
)
path_out_sum_h2 <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_summary_H2.csv"
)
path_out_sum_h3 <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_summary_H3.csv"
)
path_out_conv <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_convergence_by_retention.csv"
)
path_out_sum_md <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_summary.md"
)
path_out_fig_h2 <- file.path(
  fig_dir, "h2_car_subsampling_contrast_histograms.png"
)
path_out_fig_h3 <- file.path(
  fig_dir, "h3_car_subsampling_contrast_histograms.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "h2h3_car_subsampling_objects.rds"
)

PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
RETENTION <- c(0.80, 0.60, 0.40)
ALPHA <- 0.05

N_ITER <- as.integer(Sys.getenv("N_ITER", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
SKIP_REFIT <- identical(Sys.getenv("SKIP_REFIT", unset = "0"), "1")

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
  ((gap0 - (1 - r1)) / gap0) * 100
}

gap_direction <- function(g) {
  if (!is.finite(g)) return(NA_character_)
  if (g > 0) "closed" else if (g < 0) "widened" else "unchanged"
}

# ---------------------------------------------------------------------------
# One iteration: same seed → same rectangle draw as RE Part B
# ---------------------------------------------------------------------------
one_iteration_car <- function(
  dat, form, all_rects, adj_full, phases, keep_frac, seed_i,
  iqr_between, phase_baseline
) {
  warns <- character(0)
  set.seed(seed_i)
  n_rect <- length(all_rects)
  n_keep <- max(2L, as.integer(round(n_rect * keep_frac)))
  n_drop <- n_rect - n_keep
  keep <- sample(all_rects, size = n_keep, replace = FALSE)

  adj_sub <- adj_full[keep, keep, drop = FALSE]
  deg <- as.numeric(rowSums(adj_sub))
  n_isolates <- as.integer(sum(deg == 0))
  min_degree <- as.integer(min(deg))

  dat_sub <- dat[as.character(dat$stat_rec) %in% keep, , drop = FALSE]
  dat_sub$stat_rec <- factor(as.character(dat_sub$stat_rec), levels = keep)
  dat_sub$phase_v2 <- factor(as.character(dat_sub$phase_v2), levels = phases)

  thin_counts <- vapply(phases, function(ph) {
    length(unique(as.character(
      dat_sub$stat_rec[as.character(dat_sub$phase_v2) == ph]
    )))
  }, integer(1))
  names(thin_counts) <- phases
  thinnest_phase <- names(which.min(thin_counts))[1]
  n_rect_thinnest_phase <- as.integer(min(thin_counts))

  fit <- withCallingHandlers(
    tryCatch(
      spaMM::fitme(form, data = dat_sub, adjMatrix = adj_sub, method = "REML"),
      error = function(e) e
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  warn_txt <- if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_
  singular_or_boundary_warn <- length(warns) > 0L && any(grepl(
    "singular|boundary|false convergence|non-finite|diverg",
    warns, ignore.case = TRUE
  ))
  # Isolate-related warnings are logged via n_isolates; they do NOT alone
  # set converged = FALSE (see revision: separate logging from exclusion).

  fail_rows <- function(reason) {
    dplyr::bind_rows(lapply(c("H2", "H3"), function(hyp) {
      dplyr::bind_rows(lapply(phases, function(ph) {
        data.frame(
          seed = seed_i,
          keep_frac = keep_frac,
          n_keep = n_keep,
          n_drop = n_drop,
          n_hauls = nrow(dat_sub),
          n_rect_thinnest_phase = n_rect_thinnest_phase,
          thinnest_phase = thinnest_phase,
          n_isolates = n_isolates,
          min_degree = min_degree,
          converged = FALSE,
          singular_or_boundary_warn = singular_or_boundary_warn,
          warning_text = warn_txt,
          fail_reason = reason,
          hypothesis = hyp,
          h3_car_diagnostic_only = identical(hyp, "H3"),
          phase = ph,
          fp_slope = NA_real_,
          p_value = NA_real_,
          contrast = NA_real_,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }

  if (inherits(fit, "error") || inherits(fit, "try-error")) {
    reason <- if (inherits(fit, "error")) conditionMessage(fit) else "try-error"
    return(fail_rows(reason))
  }

  # spaMM does not expose glmmTMB-style convergence codes; treat successful
  # return + extractable finite slopes as converged (warnings logged separately).
  h2 <- tryCatch(
    extract_wb_phase_slopes_spamm(
      fit, "FP_between", "wb_car_v2_sub", "H2_spatial_between",
      phases = phases
    ),
    error = function(e) e
  )
  h3 <- tryCatch(
    extract_wb_phase_slopes_spamm(
      fit, "FP_within", "wb_car_v2_sub", "H3_temporal_within",
      phases = phases
    ),
    error = function(e) e
  )

  if (inherits(h2, "error") || inherits(h3, "error") ||
      !is.data.frame(h2) || !is.data.frame(h3) ||
      nrow(h2) != length(phases) || nrow(h3) != length(phases)) {
    reason <- "contrast_extract_failed"
    if (inherits(h2, "error")) reason <- paste0("H2: ", conditionMessage(h2))
    if (inherits(h3, "error")) {
      reason <- paste0(reason, " H3: ", conditionMessage(h3))
    }
    return(fail_rows(reason))
  }

  if (any(!is.finite(h2$fp_slope)) || any(!is.finite(h2$p_value)) ||
      any(!is.finite(h3$fp_slope)) || any(!is.finite(h3$p_value))) {
    return(fail_rows("non_finite_slopes_or_p"))
  }

  make_ok <- function(slopes, hyp, dx_fun) {
    dplyr::bind_rows(lapply(seq_len(nrow(slopes)), function(i) {
      ph <- as.character(slopes$phase[i])
      r0 <- phase_baseline$r0[phase_baseline$phase == ph]
      dx <- dx_fun(ph)
      slope <- slopes$fp_slope[i]
      p <- slopes$p_value[i]
      contrast <- pct_gap_change_from_slope(slope, dx, r0)
      data.frame(
        seed = seed_i,
        keep_frac = keep_frac,
        n_keep = n_keep,
        n_drop = n_drop,
        n_hauls = nrow(dat_sub),
        n_rect_thinnest_phase = n_rect_thinnest_phase,
        thinnest_phase = thinnest_phase,
        n_isolates = n_isolates,
        min_degree = min_degree,
        converged = TRUE,
        singular_or_boundary_warn = singular_or_boundary_warn,
        warning_text = warn_txt,
        fail_reason = NA_character_,
        hypothesis = hyp,
        h3_car_diagnostic_only = identical(hyp, "H3"),
        phase = ph,
        fp_slope = slope,
        p_value = p,
        contrast = contrast,
        stringsAsFactors = FALSE
      )
    }))
  }

  dplyr::bind_rows(
    make_ok(h2, "H2", function(ph) iqr_between),
    make_ok(h3, "H3", function(ph) {
      phase_baseline$iqr_within[phase_baseline$phase == ph]
    })
  )
}

summarise_h2 <- function(ok_df, ref_df, phases, keep_frac, n_attempted, n_converged) {
  lapply(phases, function(ph) {
    sub <- ok_df %>% filter(as.character(phase) == ph)
    ref <- ref_df %>% filter(as.character(phase) == ph)
    stopifnot(nrow(ref) == 1L)
    ref_sign <- sign(ref$full_sample_contrast[[1]])
    if (ref_sign == 0) ref_sign <- 1
    ref_sig <- isTRUE(ref$full_sample_sig[[1]])
    if (nrow(sub) == 0L) {
      return(tibble::tibble(
        phase = ph,
        retention = keep_frac,
        full_sample_contrast = ref$full_sample_contrast[[1]],
        full_sample_gap_direction = ref$full_sample_gap_direction[[1]],
        full_sample_p = ref$full_sample_p[[1]],
        full_sample_sig = ref_sig,
        subsample_median = NA_real_,
        subsample_p025 = NA_real_,
        subsample_p975 = NA_real_,
        pct_sign_match = NA_real_,
        pct_sign_and_sig_match = NA_real_,
        n_ok = 0L,
        n_attempted = n_attempted,
        n_converged = n_converged,
        pct_converged = 100 * n_converged / n_attempted
      ))
    }
    sub_sign <- sign(sub$contrast)
    sub_sign[sub_sign == 0] <- 1
    same_sign <- sub_sign == ref_sign
    same_sig <- (sub$p_value < ALPHA) == ref_sig
    tibble::tibble(
      phase = ph,
      retention = keep_frac,
      full_sample_contrast = ref$full_sample_contrast[[1]],
      full_sample_gap_direction = ref$full_sample_gap_direction[[1]],
      full_sample_p = ref$full_sample_p[[1]],
      full_sample_sig = ref_sig,
      subsample_median = stats::median(sub$contrast, na.rm = TRUE),
      subsample_p025 = as.numeric(stats::quantile(sub$contrast, 0.025, na.rm = TRUE)),
      subsample_p975 = as.numeric(stats::quantile(sub$contrast, 0.975, na.rm = TRUE)),
      pct_sign_match = 100 * mean(same_sign),
      pct_sign_and_sig_match = 100 * mean(same_sign & same_sig),
      n_ok = nrow(sub),
      n_attempted = n_attempted,
      n_converged = n_converged,
      pct_converged = 100 * n_converged / n_attempted
    )
  }) %>%
    bind_rows()
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# H2 CAR rectangle sub-sampling — run log")
logmsg("")
logmsg(
  "Refit `adjacency(1|stat_rec)` CAR on global rectangle subsamples ",
  "(retention 80/60/40%; ", N_ITER, " iterations each). ",
  "Paired seeds with RE companion. Deliverable = H2 and H3 from the shared CAR fit."
)
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf("spaMM %s", as.character(utils::packageVersion("spaMM"))))
logmsg(sprintf(
  "N_ITER = %d; SEED = %d; retention = %s",
  N_ITER, SEED, paste(sprintf("%.0f%%", 100 * RETENTION), collapse = ", ")
))
logmsg(
  "IQR convention: fixed full-sample IQR(FP_between) + fixed phase baselines ",
  "(matches phase_v2 reporting and RE Part B)."
)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
for (p in c(path_models_v2, path_round2, path_prop_h2, path_slopes)) {
  if (!file.exists(p)) stop("Missing required input: ", p)
}

obj <- readRDS(path_models_v2)
dat <- obj$data
stopifnot(all(
  c("residual", "FP_between", "FP_within", "phase_v2", "stat_rec") %in% names(dat)
))
dat$stat_rec <- factor(as.character(dat$stat_rec))
dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
if (anyNA(dat$phase_v2)) stop("NA in phase_v2 after relevel.")

round2 <- readRDS(path_round2)
if (is.null(round2$adjMatrix)) stop("Round 2 RDS missing adjMatrix.")
adjMatrix <- round2$adjMatrix

all_rects <- levels(dat$stat_rec)
if (!identical(all_rects, rownames(adjMatrix))) {
  stop(
    "Rectangle order mismatch vs adjMatrix rownames — paired RE draws would break. ",
    "RE levels and adj rownames must be identical."
  )
}
# Align factor levels to adjMatrix for full-sample reference fits if needed
dat$stat_rec <- factor(as.character(dat$stat_rec), levels = rownames(adjMatrix))

form <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)

# Fixed full-sample IQR + baselines (canonical convention)
rect_between <- dat %>% distinct(stat_rec, FP_between)
iqr_between <- as.numeric(IQR(rect_between$FP_between, na.rm = TRUE))
phase_baseline <- dat %>%
  group_by(phase_v2) %>%
  summarise(
    r0 = exp(stats::median(residual, na.rm = TRUE)),
    iqr_within = as.numeric(IQR(FP_within, na.rm = TRUE)),
    n_rect = dplyr::n_distinct(stat_rec),
    .groups = "drop"
  ) %>%
  transmute(phase = as.character(phase_v2), r0, iqr_within, n_rect)

# Canonical H2/H3 reference from reporting artifacts (wb_car_v2)
prop_h2 <- read_csv(path_prop_h2, show_col_types = FALSE) %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
prop_h3 <- read_csv(path_prop_h3, show_col_types = FALSE) %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
slopes <- read_csv(path_slopes, show_col_types = FALSE)
h2_car_slopes <- slopes %>%
  filter(model_id == "wb_car_v2", component == "FP_between") %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
h3_car_slopes <- slopes %>%
  filter(model_id == "wb_car_v2", component == "FP_within") %>%
  mutate(phase = factor(as.character(phase), levels = PHASE_V2)) %>%
  arrange(phase)
if (nrow(h2_car_slopes) != 4L || nrow(prop_h2) != 4L) {
  stop("Expected 4 CAR H2 slopes and 4 proportional H2 rows.")
}
if (nrow(h3_car_slopes) != 4L || nrow(prop_h3) != 4L) {
  stop("Expected 4 CAR H3 slopes and 4 proportional H3 rows.")
}
if (!identical(unique(as.character(prop_h2$model_id)), "wb_car_v2")) {
  stop("phase_v2_proportional_effects_H2.csv is not wb_car_v2 — refusing.")
}
if (!identical(unique(as.character(prop_h3$model_id)), "wb_car_v2")) {
  stop("phase_v2_proportional_effects_H3.csv is not wb_car_v2 — refusing.")
}

ref_h2 <- tibble::tibble(
  phase = as.character(prop_h2$phase),
  fp_slope = h2_car_slopes$fp_slope,
  full_sample_p = h2_car_slopes$p_value,
  full_sample_sig = h2_car_slopes$p_value < ALPHA,
  full_sample_contrast = prop_h2$pct_gap_change_iqr,
  full_sample_gap_direction = prop_h2$gap_direction_iqr
)

ref_h3 <- tibble::tibble(
  phase = as.character(prop_h3$phase),
  fp_slope = h3_car_slopes$fp_slope,
  full_sample_p = h3_car_slopes$p_value,
  full_sample_sig = h3_car_slopes$p_value < ALPHA,
  full_sample_contrast = prop_h3$pct_gap_change_iqr,
  full_sample_gap_direction = prop_h3$gap_direction_iqr
)

# Sanity: recompute H3 gap from slope × phase IQR should match prop table
recheck_h3 <- vapply(seq_len(4L), function(i) {
  ph <- as.character(h3_car_slopes$phase[i])
  pct_gap_change_from_slope(
    h3_car_slopes$fp_slope[i],
    phase_baseline$iqr_within[phase_baseline$phase == ph],
    phase_baseline$r0[phase_baseline$phase == ph]
  )
}, numeric(1))
if (max(abs(recheck_h3 - ref_h3$full_sample_contrast)) > 1e-6) {
  stop(
    "Recomputed CAR IQR gap-change does not match proportional_effects_H3.csv; ",
    "check IQR / baseline convention."
  )
}
recheck <- vapply(seq_len(4L), function(i) {
  pct_gap_change_from_slope(
    h2_car_slopes$fp_slope[i],
    iqr_between,
    phase_baseline$r0[phase_baseline$phase == as.character(h2_car_slopes$phase[i])]
  )
}, numeric(1))
if (max(abs(recheck - ref_h2$full_sample_contrast)) > 1e-6) {
  stop(
    "Recomputed CAR IQR gap-change does not match proportional_effects_H2.csv; ",
    "check IQR / baseline convention."
  )
}

logmsg("Loaded: ", path_models_v2)
logmsg("adjMatrix: ", path_round2, sprintf(" (%d x %d)", nrow(adjMatrix), ncol(adjMatrix)))
logmsg(sprintf(
  "Full sample: %d hauls, %d rectangles; IQR(FP_between)=%.4f (fixed)",
  nrow(dat), length(all_rects), iqr_between
))
logmsg("Formula: ", paste(deparse(form), collapse = " "))
logmsg("H2 reference: wb_car_v2 from phase_v2_proportional_effects_H2.csv")
for (i in seq_len(nrow(ref_h2))) {
  r <- ref_h2[i, ]
  logmsg(sprintf(
    "  H2 %s: contrast=%+.2f%% (%s) slope=%+.4f p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    r$fp_slope, fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}
logmsg("H3 reference: wb_car_v2 from phase_v2_proportional_effects_H3.csv")
for (i in seq_len(nrow(ref_h3))) {
  r <- ref_h3[i, ]
  logmsg(sprintf(
    "  H3 %s: contrast=%+.2f%% (%s) slope=%+.4f p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    r$fp_slope, fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}

# ---------------------------------------------------------------------------
# Benchmark + batch, or reuse existing iterations
# ---------------------------------------------------------------------------
if (SKIP_REFIT) {
  if (!file.exists(path_out_iter)) {
    stop("SKIP_REFIT=1 but missing ", path_out_iter)
  }
  logmsg("")
  logmsg("## SKIP_REFIT=1: loading existing iterations (no new CAR fits)")
  iter_df <- read_csv(path_out_iter, show_col_types = FALSE)
  prev <- if (file.exists(path_out_rds)) readRDS(path_out_rds) else NULL
  t_bench <- if (!is.null(prev$settings$single_fit_benchmark_sec_80)) {
    prev$settings$single_fit_benchmark_sec_80
  } else {
    NA_real_
  }
  t40 <- if (!is.null(prev$settings$single_fit_benchmark_sec_40)) {
    prev$settings$single_fit_benchmark_sec_40
  } else {
    NA_real_
  }
  runtime_sec <- if (!is.null(prev$runtime_sec)) prev$runtime_sec else NA_real_
  n_ok_h3 <- iter_df %>%
    filter(hypothesis == "H3") %>%
    distinct(keep_frac, seed) %>%
    nrow()
  N_ITER <- iter_df %>%
    filter(abs(keep_frac - 0.80) < 1e-9) %>%
    distinct(seed) %>%
    nrow()
  logmsg(sprintf(
    "Loaded %d iteration rows (%d H3 seed×retention cells; N_ITER=%d) from %s",
    nrow(iter_df), n_ok_h3, N_ITER, path_out_iter
  ))
} else {
  logmsg("")
  logmsg("## Benchmark")
  set.seed(SEED)
  keep_bench <- sample(all_rects, size = as.integer(round(length(all_rects) * 0.80)), replace = FALSE)
  adj_b <- adjMatrix[keep_bench, keep_bench, drop = FALSE]
  dat_b <- dat[as.character(dat$stat_rec) %in% keep_bench, , drop = FALSE]
  dat_b$stat_rec <- factor(as.character(dat_b$stat_rec), levels = keep_bench)
  t_bench0 <- proc.time()[[3]]
  fit_b <- spaMM::fitme(form, data = dat_b, adjMatrix = adj_b, method = "REML")
  t_bench <- proc.time()[[3]] - t_bench0
  set.seed(SEED + 200000L + 1L)
  keep40 <- sample(all_rects, size = as.integer(round(length(all_rects) * 0.40)), replace = FALSE)
  adj40 <- adjMatrix[keep40, keep40, drop = FALSE]
  dat40 <- dat[as.character(dat$stat_rec) %in% keep40, , drop = FALSE]
  dat40$stat_rec <- factor(as.character(dat40$stat_rec), levels = keep40)
  t40_0 <- proc.time()[[3]]
  fit40 <- spaMM::fitme(form, data = dat40, adjMatrix = adj40, method = "REML")
  t40 <- proc.time()[[3]] - t40_0
  logmsg(sprintf(
    "Single-fit 80%%: %.3f sec (isolates=%d); 40%%: %.3f sec (isolates=%d)",
    t_bench, sum(rowSums(adj_b) == 0), t40, sum(rowSums(adj40) == 0)
  ))
  logmsg(sprintf(
    "Projected %d fits ≈ %.1f min sequential (using 80%% time)",
    length(RETENTION) * N_ITER, length(RETENTION) * N_ITER * t_bench / 60
  ))
  rm(fit_b, dat_b, adj_b, keep_bench, fit40, dat40, adj40, keep40)

  logmsg("")
  logmsg("## Batch")
  seed_blocks <- lapply(seq_along(RETENTION), function(j) {
    SEED + (j - 1L) * 100000L + seq_len(N_ITER)
  })
  names(seed_blocks) <- sprintf("%.2f", RETENTION)

  t_batch0 <- proc.time()[[3]]
  iter_list <- list()

  for (j in seq_along(RETENTION)) {
    kf <- RETENTION[[j]]
    seeds <- seed_blocks[[j]]
    logmsg(sprintf("### Retention %.0f%% (%d iterations)", 100 * kf, N_ITER))
    level_list <- vector("list", N_ITER)
    for (i in seq_len(N_ITER)) {
      if (i == 1L || i %% 100L == 0L || i == N_ITER) {
        logmsg(sprintf("  iteration %d / %d ...", i, N_ITER))
      }
      level_list[[i]] <- one_iteration_car(
        dat, form, all_rects, adjMatrix, PHASE_V2, kf, seeds[[i]],
        iqr_between, phase_baseline
      )
    }
    iter_list[[as.character(kf)]] <- dplyr::bind_rows(level_list)
  }

  iter_df <- dplyr::bind_rows(iter_list)
  runtime_sec <- proc.time()[[3]] - t_batch0
  logmsg(sprintf(
    "Batch runtime: %.1f sec (%.2f min; %.2f sec/fit average)",
    runtime_sec, runtime_sec / 60, runtime_sec / (length(RETENTION) * N_ITER)
  ))

  write_csv(iter_df, path_out_iter)
  logmsg("Saved: ", path_out_iter)
}

# Convergence summary (one row per iteration seed) — isolates logged separately
conv_by_seed <- iter_df %>%
  distinct(
    keep_frac, seed, converged, singular_or_boundary_warn,
    n_rect_thinnest_phase, thinnest_phase, n_isolates, min_degree
  )

conv_summary <- conv_by_seed %>%
  group_by(keep_frac) %>%
  summarise(
    n_attempted = dplyr::n(),
    n_converged = sum(converged),
    n_excluded = sum(!converged),
    pct_converged = 100 * mean(converged),
    pct_excluded = 100 * mean(!converged),
    n_singular_or_boundary_warn = sum(singular_or_boundary_warn, na.rm = TRUE),
    n_iters_with_isolates = sum(n_isolates > 0L),
    max_n_isolates = max(n_isolates, na.rm = TRUE),
    min_n_rect_thinnest_phase = min(n_rect_thinnest_phase, na.rm = TRUE),
    median_n_rect_thinnest_phase = stats::median(n_rect_thinnest_phase, na.rm = TRUE),
    min_min_degree = min(min_degree, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(keep_frac))

write_csv(conv_summary, path_out_conv)
logmsg("Saved: ", path_out_conv)

logmsg("")
logmsg("### Convergence by retention (isolates logged, not auto-excluded)")
for (i in seq_len(nrow(conv_summary))) {
  r <- conv_summary[i, ]
  logmsg(sprintf(
    "  %.0f%%: attempted=%d; converged=%d (%.1f%%); excluded=%d; singular/boundary warns=%d; iters with isolates=%d (max isolates=%d); min thinnest-phase n_rect=%d; min degree=%d",
    100 * r$keep_frac, r$n_attempted, r$n_converged, r$pct_converged,
    r$n_excluded, r$n_singular_or_boundary_warn,
    r$n_iters_with_isolates, r$max_n_isolates,
    r$min_n_rect_thinnest_phase, r$min_min_degree
  ))
}

min_thin_40 <- conv_by_seed %>%
  filter(abs(keep_frac - 0.40) < 1e-9) %>%
  pull(n_rect_thinnest_phase) %>%
  min(na.rm = TRUE)
logmsg(sprintf(
  "Minimum realized thinnest-phase rectangle count at 40%% retention: %d",
  min_thin_40
))

# ---------------------------------------------------------------------------
# H2 summary only
# ---------------------------------------------------------------------------
sum_h2_list <- list()
for (kf in RETENTION) {
  level <- iter_df %>% filter(abs(keep_frac - kf) < 1e-9)
  n_att <- N_ITER
  n_conv <- level %>%
    distinct(seed, converged) %>%
    summarise(n = sum(converged)) %>%
    pull(n)
  ok_h2 <- level %>%
    filter(
      hypothesis == "H2", converged,
      is.finite(contrast), is.finite(p_value)
    )
  sum_h2_list[[as.character(kf)]] <- summarise_h2(
    ok_h2, ref_h2, PHASE_V2, kf, n_att, n_conv
  )
}
sum_h2_out <- bind_rows(sum_h2_list) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase, desc(retention))

write_csv(sum_h2_out, path_out_sum_h2)
logmsg("Saved: ", path_out_sum_h2)

logmsg("")
logmsg("### H2 summary")
logmsg(
  "| Phase | Retention | Full-sample GLMM contrast | Subsample median [2.5, 97.5] | % sign match | % sign+sig match | % converged |"
)
logmsg("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(sum_h2_out))) {
  r <- sum_h2_out[i, ]
  logmsg(sprintf(
    "| %s | %.0f%% | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %.1f%% |",
    as.character(r$phase), 100 * r$retention, r$full_sample_contrast,
    r$subsample_median, r$subsample_p025, r$subsample_p975,
    r$pct_sign_match, r$pct_sign_and_sig_match, r$pct_converged
  ))
}

sum_h3_list <- list()
for (kf in RETENTION) {
  level <- iter_df %>% filter(abs(keep_frac - kf) < 1e-9)
  n_att <- N_ITER
  n_conv <- level %>%
    distinct(seed, converged) %>%
    summarise(n = sum(converged)) %>%
    pull(n)
  ok_h3 <- level %>%
    filter(
      hypothesis == "H3", converged,
      is.finite(contrast), is.finite(p_value)
    )
  sum_h3_list[[as.character(kf)]] <- summarise_h2(
    ok_h3, ref_h3, PHASE_V2, kf, n_att, n_conv
  )
}
sum_h3_out <- bind_rows(sum_h3_list) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase, desc(retention))

write_csv(sum_h3_out, path_out_sum_h3)
logmsg("Saved: ", path_out_sum_h3)

logmsg("")
logmsg("### H3 summary")
logmsg(
  "| Phase | Retention | Full-sample GLMM contrast | Subsample median [2.5, 97.5] | % sign match | % sign+sig match | % converged |"
)
logmsg("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(sum_h3_out))) {
  r <- sum_h3_out[i, ]
  logmsg(sprintf(
    "| %s | %.0f%% | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %.1f%% |",
    as.character(r$phase), 100 * r$retention, r$full_sample_contrast,
    r$subsample_median, r$subsample_p025, r$subsample_p975,
    r$pct_sign_match, r$pct_sign_and_sig_match, r$pct_converged
  ))
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
make_car_hist <- function(hyp, ref_df, title, outfile) {
  ok_plot <- iter_df %>%
    filter(hypothesis == hyp, converged, is.finite(contrast)) %>%
    mutate(
      phase = factor(as.character(phase), levels = PHASE_V2),
      retention = factor(
        sprintf("%.0f%%", 100 * keep_frac),
        levels = sprintf("%.0f%%", 100 * RETENTION)
      )
    )
  obs <- ref_df %>%
    transmute(
      phase = factor(phase, levels = PHASE_V2),
      contrast = full_sample_contrast
    )
  p <- ggplot(ok_plot, aes(x = contrast, colour = retention, fill = retention)) +
    geom_density(alpha = 0.18, linewidth = 0.7) +
    geom_vline(
      data = obs, aes(xintercept = contrast),
      colour = "#b2182b", linetype = "dashed", linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~phase, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = c("80%" = "#2166ac", "60%" = "#67a9cf", "40%" = "#d1e5f0")) +
    scale_fill_manual(values = c("80%" = "#2166ac", "60%" = "#67a9cf", "40%" = "#d1e5f0")) +
    labs(
      title = title,
      subtitle = sprintf(
        "spaMM adjacency CAR refits; dashed red = canonical wb_car_v2 full-sample contrast; seed = %d",
        SEED
      ),
      x = "IQR contrast (% gap change; + = closed, − = widened)",
      y = "Density",
      colour = "Retention",
      fill = "Retention"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(outfile, p, width = 9, height = 7, dpi = 150)
}

make_car_hist("H2", ref_h2, "H2 CAR rectangle sub-sampling: IQR gap-change by phase", path_out_fig_h2)
logmsg("Saved: ", path_out_fig_h2)
make_car_hist("H3", ref_h3, "H3 CAR rectangle sub-sampling: IQR gap-change by phase", path_out_fig_h3)
logmsg("Saved: ", path_out_fig_h3)

fmt_md_table <- function(sum_out) {
  c(
    "| Phase | Retention | Full-sample GLMM contrast | Subsample median [2.5, 97.5] | % sign match | % sign+sig match | % converged |",
    "|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(sum_out)), function(i) {
      r <- sum_out[i, ]
      sprintf(
        "| %s | %.0f%% | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %.1f%% |",
        as.character(r$phase), 100 * r$retention, r$full_sample_contrast,
        r$subsample_median, r$subsample_p025, r$subsample_p975,
        r$pct_sign_match, r$pct_sign_and_sig_match, r$pct_converged
      )
    }, character(1))
  )
}

fmt_sec <- function(x) {
  if (!is.finite(x)) "NA (reused existing batch)" else sprintf("%.3f", x)
}
fmt_runtime <- function(sec) {
  if (!is.finite(sec)) {
    "NA (reused existing batch)"
  } else {
    sprintf("%.1f sec (%.2f min)", sec, sec / 60)
  }
}

sub_md <- c(
  "# H2/H3 CAR rectangle sub-sampling sensitivity",
  "",
  "Refit of the canonical shared CAR model",
  "`residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)`",
  "on global rectangle subsamples (adjacency matrix subset to retained rectangles).",
  "Full-sample reference = `wb_car_v2` IQR gap-change from",
  "`phase_v2_proportional_effects_H{2,3}.csv`.",
  "",
  if (SKIP_REFIT) {
    paste0(
      "This run used `SKIP_REFIT=1`: H2 and H3 summaries were built from ",
      "existing `h2h3_car_subsampling_iterations.csv` (no new CAR fits)."
    )
  } else {
    "H2 and H3 slopes are extracted from the same CAR refit on each subsample."
  },
  "",
  "## Batch header",
  "",
  sprintf("- Seed: %d (paired with RE companion seed blocks)", SEED),
  sprintf("- Iterations per retention level: %d", N_ITER),
  sprintf("- Retention levels: %s", paste(sprintf("%.0f%%", 100 * RETENTION), collapse = ", ")),
  paste0("- Single-fit benchmark (80%): ", if (is.finite(t_bench)) sprintf("%.3f sec", t_bench) else "NA (reused existing batch)"),
  paste0("- Total batch runtime: ", fmt_runtime(runtime_sec)),
  sprintf("- Minimum realized thinnest-phase rectangle count at 40%% retention: %d", min_thin_40),
  "- IQR convention: fixed full-sample `IQR(FP_between)` + phase-specific `IQR(FP_within)` + fixed phase baselines",
  "",
  "### Convergence / exclusions by retention",
  "",
  "| Retention | Attempted | Converged | Excluded | % excluded | Singular/boundary warns | Iters with isolates | Max isolates | Min thinnest-phase n_rect |",
  "|---|---|---|---|---|---|---|---|---|",
  vapply(seq_len(nrow(conv_summary)), function(i) {
    r <- conv_summary[i, ]
    sprintf(
      "| %.0f%% | %d | %d | %d | %.1f%% | %d | %d | %d | %d |",
      100 * r$keep_frac, r$n_attempted, r$n_converged, r$n_excluded,
      r$pct_excluded, r$n_singular_or_boundary_warn,
      r$n_iters_with_isolates, r$max_n_isolates, r$min_n_rect_thinnest_phase
    )
  }, character(1)),
  "",
  "Isolates are logged (`n_isolates`, `min_degree`) but do **not** alone exclude",
  "an iteration; exclusion is optimizer/extract failure only.",
  "",
  "## H2 (`FP_between`, `wb_car_v2`)",
  "",
  fmt_md_table(sum_h2_out),
  "",
  "## H3 (`FP_within`, `wb_car_v2`)",
  "",
  fmt_md_table(sum_h3_out),
  "",
  "## Figures",
  "",
  paste0("- H2: ", path_out_fig_h2),
  paste0("- H3: ", path_out_fig_h3),
  "",
  "### Notes",
  "",
  "- Match rule: same sign as canonical full-sample IQR contrast **and** same significance call (p < 0.05) as the full-sample CAR phase slope.",
  "- Paired rectangle draws with `run_h2h3_rectangle_subsampling_refit.R` (same SEED blocks).",
  ""
)
writeLines(sub_md, path_out_sum_md)
logmsg("Saved: ", path_out_sum_md)

saveRDS(
  list(
    ref_h2 = ref_h2,
    ref_h3 = ref_h3,
    iterations = iter_df,
    summary_h2 = sum_h2_out,
    summary_h3 = sum_h3_out,
    convergence = conv_summary,
    settings = list(
      n_iter = N_ITER,
      seed = SEED,
      retention = RETENTION,
      alpha = ALPHA,
      phases = PHASE_V2,
      formula = form,
      method = "spaMM_CAR_refit_global_rectangle_subsample",
      paired_with_re_companion = TRUE,
      skip_refit = SKIP_REFIT,
      iqr_convention = "fixed_full_sample",
      iqr_between = iqr_between,
      single_fit_benchmark_sec_80 = t_bench,
      single_fit_benchmark_sec_40 = t40,
      min_thinnest_phase_n_rect_at_40 = min_thin_40
    ),
    runtime_sec = runtime_sec
  ),
  path_out_rds
)
logmsg("Saved: ", path_out_rds)

writeLines(run_log, path_out_run_log)
logmsg("Saved: ", path_out_run_log)
cat("=== H2/H3 CAR rectangle sub-sampling complete. ===\n")
