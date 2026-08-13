# H2/H3 rectangle sub-sampling sensitivity — GLMM refit (primary Part B)
#
# PURPOSE: refit primary_model_v2 on random rectangle subsamples and check
# whether phase-specific H2 (FP_between) and H3 (FP_within) IQR contrasts are
# stable. Same model form on subsample and full sample — mismatches are from
# rectangle removal only.
#
# Method:
#   - Global rectangle subsample (one draw per iteration; all phases)
#   - Retention levels: 80%, 60%, 40%
#   - 1,000 iterations per level (3,000 refits)
#   - Extract eight phase contrasts per refit (4 H2 + 4 H3)
#   - Reference: primary_model_v2 full-sample RE contrasts
#   - Seed = 42
#
# Supersedes OLS residual-proxy Part B. Archived earlier H2-only 25%-drop run:
#   outputs/h2_rectangle_subsampling_*
#
# Prerequisite: outputs/primary_model_v2.rds
#
# Run: Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_refit.R
# Optional: N_ITER=1000 SEED=42 N_CORES=1

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_rectangle_subsampling_refit.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_refit.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

# ---------------------------------------------------------------------------
# Paths / settings
# ---------------------------------------------------------------------------
path_models_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_iter <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_iterations.csv"
)
path_out_sum_h2 <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_summary_H2.csv"
)
path_out_sum_h3 <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_summary_H3.csv"
)
path_out_sum_md <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_summary.md"
)
path_out_conv <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_convergence_by_retention.csv"
)
path_out_fig_h2 <- file.path(
  fig_dir, "h2_glmm_subsampling_contrast_histograms.png"
)
path_out_fig_h3 <- file.path(
  fig_dir, "h3_glmm_subsampling_contrast_histograms.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "h2h3_glmm_subsampling_objects.rds"
)

PHASE_V2 <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
RETENTION <- c(0.80, 0.60, 0.40)
ALPHA <- 0.05

N_ITER <- as.integer(Sys.getenv("N_ITER", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
N_CORES <- as.integer(Sys.getenv("N_CORES", unset = "1"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

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
# One iteration: global rectangle subsample → one GLMM refit → 8 contrasts
# ---------------------------------------------------------------------------
one_iteration <- function(
  dat, form, all_rects, phases, keep_frac, seed_i,
  iqr_between, phase_baseline
) {
  warns <- character(0)
  set.seed(seed_i)
  n_rect <- length(all_rects)
  n_keep <- max(2L, as.integer(round(n_rect * keep_frac)))
  n_drop <- n_rect - n_keep
  keep <- sample(all_rects, size = n_keep, replace = FALSE)

  dat_sub <- dat[as.character(dat$stat_rec) %in% keep, , drop = FALSE]
  dat_sub$stat_rec <- factor(as.character(dat_sub$stat_rec))
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
      glmmTMB::glmmTMB(form, data = dat_sub, REML = TRUE),
      error = function(e) e
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  warn_txt <- if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_
  singular_warn <- length(warns) > 0L && any(grepl(
    "singular|boundary|false convergence|non-finite",
    warns, ignore.case = TRUE
  ))

  fail_row <- function(hyp, ph, reason) {
    data.frame(
      seed = seed_i,
      keep_frac = keep_frac,
      n_keep = n_keep,
      n_drop = n_drop,
      n_hauls = nrow(dat_sub),
      n_rect_thinnest_phase = n_rect_thinnest_phase,
      thinnest_phase = thinnest_phase,
      converged = FALSE,
      singular_or_boundary_warn = singular_warn,
      warning_text = warn_txt,
      fail_reason = reason,
      hypothesis = hyp,
      phase = ph,
      fp_slope = NA_real_,
      p_value = NA_real_,
      contrast = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (inherits(fit, "error") || inherits(fit, "try-error")) {
    reason <- if (inherits(fit, "error")) conditionMessage(fit) else "try-error"
    return(dplyr::bind_rows(lapply(c("H2", "H3"), function(hyp) {
      dplyr::bind_rows(lapply(phases, function(ph) fail_row(hyp, ph, reason)))
    })))
  }

  conv_code <- tryCatch(fit$fit$convergence, error = function(e) NA_integer_)
  converged <- isTRUE(conv_code == 0)

  h2 <- tryCatch(
    extract_wb_phase_slopes(
      fit, "FP_between", "wb_primary_v2_sub", "H2_spatial_between",
      phases = phases
    ),
    error = function(e) e
  )
  h3 <- tryCatch(
    extract_wb_phase_slopes(
      fit, "FP_within", "wb_primary_v2_sub", "H3_temporal_within",
      phases = phases
    ),
    error = function(e) e
  )

  if (inherits(h2, "error") || inherits(h3, "error") ||
      !is.data.frame(h2) || !is.data.frame(h3) ||
      nrow(h2) != length(phases) || nrow(h3) != length(phases)) {
    reason <- "contrast_extract_failed"
    if (inherits(h2, "error")) reason <- paste0("H2: ", conditionMessage(h2))
    if (inherits(h3, "error")) reason <- paste0(reason, " H3: ", conditionMessage(h3))
    return(dplyr::bind_rows(lapply(c("H2", "H3"), function(hyp) {
      dplyr::bind_rows(lapply(phases, function(ph) fail_row(hyp, ph, reason)))
    })))
  }

  if (!converged) {
    return(dplyr::bind_rows(lapply(c("H2", "H3"), function(hyp) {
      dplyr::bind_rows(lapply(phases, function(ph) {
        fail_row(hyp, ph, sprintf("optimizer_convergence=%s", conv_code))
      }))
    })))
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
        converged = TRUE,
        singular_or_boundary_warn = singular_warn,
        warning_text = warn_txt,
        fail_reason = NA_character_,
        hypothesis = hyp,
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

summarise_level_phase <- function(ok_df, ref_df, phases, keep_frac, n_attempted, n_converged) {
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
    bind_rows() %>%
    mutate(
      phase = factor(phase, levels = phases),
      retention = keep_frac
    )
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# H2/H3 GLMM rectangle sub-sampling — run log")
logmsg("")
logmsg(
  "Primary Part B: refit `primary_model_v2` on global rectangle subsamples ",
  "(retention 80/60/40%; ", N_ITER, " iterations each). ",
  "Contrast = IQR gap-change % from phase slopes (fixed full-sample IQR / baselines). ",
  "OLS residual proxy superseded."
)
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf("glmmTMB %s", as.character(utils::packageVersion("glmmTMB"))))
logmsg(sprintf(
  "N_ITER = %d; SEED = %d; N_CORES = %d; retention = %s",
  N_ITER, SEED, N_CORES, paste(sprintf("%.0f%%", 100 * RETENTION), collapse = ", ")
))

# ---------------------------------------------------------------------------
# Load + full-sample reference
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
if (!file.exists(path_models_v2)) {
  stop("Missing ", path_models_v2, " — run pipeline/run_h2h3_phase_v2_refit.R first.")
}
obj <- readRDS(path_models_v2)
if (!inherits(obj$primary_model_v2, "glmmTMB")) {
  stop("primary_model_v2.rds missing glmmTMB fit under $primary_model_v2")
}
dat <- obj$data
stopifnot(all(
  c("residual", "FP_between", "FP_within", "phase_v2", "stat_rec") %in% names(dat)
))
dat$stat_rec <- factor(as.character(dat$stat_rec))
dat$phase_v2 <- factor(as.character(dat$phase_v2), levels = PHASE_V2)
if (anyNA(dat$phase_v2)) stop("NA in phase_v2 after relevel.")

form <- if (!is.null(obj$formula_v2)) {
  obj$formula_v2
} else {
  residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)
}

fit_full <- obj$primary_model_v2
h2_full <- extract_wb_phase_slopes(
  fit_full, "FP_between", "wb_primary_v2", "H2_spatial_between",
  phases = PHASE_V2
)
h3_full <- extract_wb_phase_slopes(
  fit_full, "FP_within", "wb_primary_v2", "H3_temporal_within",
  phases = PHASE_V2
)
h2_full$phase <- factor(as.character(h2_full$phase), levels = PHASE_V2)
h3_full$phase <- factor(as.character(h3_full$phase), levels = PHASE_V2)

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
  transmute(
    phase = as.character(phase_v2),
    r0, iqr_within, n_rect
  )

ref_from_slopes <- function(slopes, dx_map) {
  tibble::tibble(
    phase = as.character(slopes$phase),
    fp_slope = slopes$fp_slope,
    full_sample_p = slopes$p_value,
    full_sample_sig = slopes$p_value < ALPHA,
    full_sample_contrast = vapply(seq_len(nrow(slopes)), function(i) {
      ph <- as.character(slopes$phase[i])
      pct_gap_change_from_slope(
        slopes$fp_slope[i],
        dx_map[[ph]],
        phase_baseline$r0[phase_baseline$phase == ph]
      )
    }, numeric(1))
  ) %>%
    mutate(full_sample_gap_direction = vapply(full_sample_contrast, gap_direction, character(1)))
}

dx_h2 <- setNames(rep(iqr_between, 4L), PHASE_V2)
dx_h3 <- setNames(phase_baseline$iqr_within, phase_baseline$phase)
ref_h2 <- ref_from_slopes(h2_full, dx_h2)
ref_h3 <- ref_from_slopes(h3_full, dx_h3)

all_rects <- levels(dat$stat_rec)
if (length(all_rects) == 0L) {
  all_rects <- sort(unique(as.character(dat$stat_rec)))
}
n_rect_all <- length(all_rects)

logmsg("Loaded: ", path_models_v2)
logmsg(sprintf(
  "Full sample: %d hauls, %d rectangles, years %d–%d",
  nrow(dat), n_rect_all, min(dat$year), max(dat$year)
))
logmsg("Formula: ", paste(deparse(form), collapse = " "))
logmsg(sprintf("IQR(FP_between) = %.4f (fixed for H2 contrast transform)", iqr_between))
for (i in seq_len(nrow(phase_baseline))) {
  r <- phase_baseline[i, ]
  logmsg(sprintf(
    "  phase %s: n_rect=%d; baseline_ratio=%.4f; IQR(FP_within)=%.3f",
    r$phase, r$n_rect, r$r0, r$iqr_within
  ))
}

logmsg("")
logmsg("### Full-sample reference (primary_model_v2 RE → IQR gap-change %)")
for (i in seq_len(nrow(ref_h2))) {
  r <- ref_h2[i, ]
  logmsg(sprintf(
    "  H2 %s: contrast=%+.2f%% (%s) slope=%+.4f p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    r$fp_slope, fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}
for (i in seq_len(nrow(ref_h3))) {
  r <- ref_h3[i, ]
  logmsg(sprintf(
    "  H3 %s: contrast=%+.2f%% (%s) slope=%+.4f p=%s sig=%s",
    r$phase, r$full_sample_contrast, r$full_sample_gap_direction,
    h3_full$fp_slope[as.character(h3_full$phase) == r$phase],
    fmt_p(r$full_sample_p), yn(r$full_sample_sig)
  ))
}

# ---------------------------------------------------------------------------
# Benchmark one refit
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Benchmark (single 80% refit)")
set.seed(SEED)
keep_bench <- sample(all_rects, size = as.integer(round(n_rect_all * 0.80)), replace = FALSE)
dat_bench <- dat[as.character(dat$stat_rec) %in% keep_bench, , drop = FALSE]
dat_bench$stat_rec <- factor(as.character(dat_bench$stat_rec))
t_bench0 <- proc.time()[[3]]
fit_bench <- glmmTMB::glmmTMB(form, data = dat_bench, REML = TRUE)
t_bench <- proc.time()[[3]] - t_bench0
logmsg(sprintf(
  "Single-fit runtime: %.3f sec (converged=%s); projected %d×%d = %d fits ≈ %.1f min sequential",
  t_bench, yn(isTRUE(fit_bench$fit$convergence == 0)),
  length(RETENTION), N_ITER, length(RETENTION) * N_ITER,
  length(RETENTION) * N_ITER * t_bench / 60
))
rm(fit_bench, dat_bench, keep_bench)

# ---------------------------------------------------------------------------
# Batch
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Batch")
# Distinct seed blocks per retention level; shared SEED base
seed_blocks <- lapply(seq_along(RETENTION), function(j) {
  SEED + (j - 1L) * 100000L + seq_len(N_ITER)
})
names(seed_blocks) <- sprintf("%.2f", RETENTION)

t_batch0 <- proc.time()[[3]]
iter_list <- list()

for (j in seq_along(RETENTION)) {
  kf <- RETENTION[[j]]
  seeds <- seed_blocks[[j]]
  logmsg(sprintf(
    "### Retention %.0f%% (%d iterations)", 100 * kf, N_ITER
  ))

  if (N_CORES > 1L && .Platform$OS.type == "unix") {
    logmsg(sprintf(
      "Parallel mclapply mc.cores=%d (TMB/fork can be fragile; prefer N_CORES=1 if issues)",
      N_CORES
    ))
    level_list <- parallel::mclapply(seq_len(N_ITER), function(i) {
      one_iteration(
        dat, form, all_rects, PHASE_V2, kf, seeds[[i]],
        iqr_between, phase_baseline
      )
    }, mc.cores = N_CORES)
  } else {
    level_list <- vector("list", N_ITER)
    for (i in seq_len(N_ITER)) {
      if (i == 1L || i %% 100L == 0L || i == N_ITER) {
        logmsg(sprintf("  iteration %d / %d ...", i, N_ITER))
      }
      level_list[[i]] <- one_iteration(
        dat, form, all_rects, PHASE_V2, kf, seeds[[i]],
        iqr_between, phase_baseline
      )
    }
  }
  iter_list[[as.character(kf)]] <- dplyr::bind_rows(level_list)
}

iter_df <- dplyr::bind_rows(iter_list)
runtime_sec <- proc.time()[[3]] - t_batch0
logmsg(sprintf("Batch runtime: %.1f sec (%.2f sec/fit average)", runtime_sec, runtime_sec / (length(RETENTION) * N_ITER)))

write_csv(iter_df, path_out_iter)
logmsg("Saved: ", path_out_iter)

# Convergence summary (per retention; one row per iteration seed)
conv_by_seed <- iter_df %>%
  distinct(keep_frac, seed, converged, singular_or_boundary_warn, n_rect_thinnest_phase, thinnest_phase)

conv_summary <- conv_by_seed %>%
  group_by(keep_frac) %>%
  summarise(
    n_attempted = dplyr::n(),
    n_converged = sum(converged),
    n_excluded = sum(!converged),
    pct_converged = 100 * mean(converged),
    pct_excluded = 100 * mean(!converged),
    n_singular_or_boundary_warn = sum(singular_or_boundary_warn, na.rm = TRUE),
    min_n_rect_thinnest_phase = min(n_rect_thinnest_phase, na.rm = TRUE),
    median_n_rect_thinnest_phase = stats::median(n_rect_thinnest_phase, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(keep_frac))

write_csv(conv_summary, path_out_conv)
logmsg("Saved: ", path_out_conv)

logmsg("")
logmsg("### Convergence by retention")
for (i in seq_len(nrow(conv_summary))) {
  r <- conv_summary[i, ]
  logmsg(sprintf(
    "  %.0f%%: attempted=%d; converged=%d (%.1f%%); excluded=%d (%.1f%%); singular/boundary warns=%d; min thinnest-phase n_rect=%d",
    100 * r$keep_frac, r$n_attempted, r$n_converged, r$pct_converged,
    r$n_excluded, r$pct_excluded, r$n_singular_or_boundary_warn,
    r$min_n_rect_thinnest_phase
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
# Summaries H2 / H3
# ---------------------------------------------------------------------------
sum_h2_list <- list()
sum_h3_list <- list()
for (kf in RETENTION) {
  level <- iter_df %>% filter(abs(keep_frac - kf) < 1e-9)
  conv_seeds <- level %>%
    distinct(seed, converged) %>%
    filter(converged) %>%
    pull(seed)
  n_att <- N_ITER
  n_conv <- length(unique(conv_seeds))
  ok_h2 <- level %>%
    filter(hypothesis == "H2", converged, is.finite(contrast), is.finite(p_value))
  ok_h3 <- level %>%
    filter(hypothesis == "H3", converged, is.finite(contrast), is.finite(p_value))
  sum_h2_list[[as.character(kf)]] <- summarise_level_phase(
    ok_h2, ref_h2, PHASE_V2, kf, n_att, n_conv
  )
  sum_h3_list[[as.character(kf)]] <- summarise_level_phase(
    ok_h3, ref_h3, PHASE_V2, kf, n_att, n_conv
  )
}
sum_h2 <- bind_rows(sum_h2_list) %>%
  mutate(
    phase = factor(phase, levels = PHASE_V2),
    retention = factor(
      sprintf("%.0f%%", 100 * retention),
      levels = sprintf("%.0f%%", 100 * RETENTION)
    )
  ) %>%
  arrange(phase, retention)
sum_h3 <- bind_rows(sum_h3_list) %>%
  mutate(
    phase = factor(phase, levels = PHASE_V2),
    retention = factor(
      sprintf("%.0f%%", 100 * retention),
      levels = sprintf("%.0f%%", 100 * RETENTION)
    )
  ) %>%
  arrange(phase, retention)

# Restore numeric retention for CSV clarity
sum_h2_out <- bind_rows(sum_h2_list) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase, desc(retention))
sum_h3_out <- bind_rows(sum_h3_list) %>%
  mutate(phase = factor(phase, levels = PHASE_V2)) %>%
  arrange(phase, desc(retention))

write_csv(sum_h2_out, path_out_sum_h2)
write_csv(sum_h3_out, path_out_sum_h3)
logmsg("Saved: ", path_out_sum_h2)
logmsg("Saved: ", path_out_sum_h3)

print_sum_table <- function(sum_df, label) {
  logmsg("")
  logmsg(sprintf("### %s summary", label))
  logmsg(
    "| Phase | Retention | Full-sample GLMM contrast | Subsample median [2.5, 97.5] | % sign match | % sign+sig match | % converged |"
  )
  logmsg("|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(sum_df))) {
    r <- sum_df[i, ]
    logmsg(sprintf(
      "| %s | %.0f%% | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %.1f%% |",
      as.character(r$phase), 100 * r$retention, r$full_sample_contrast,
      r$subsample_median, r$subsample_p025, r$subsample_p975,
      r$pct_sign_match, r$pct_sign_and_sig_match, r$pct_converged
    ))
  }
}
print_sum_table(sum_h2_out, "H2")
print_sum_table(sum_h3_out, "H3")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
ok_plot <- iter_df %>%
  filter(converged, is.finite(contrast)) %>%
  mutate(
    phase = factor(as.character(phase), levels = PHASE_V2),
    retention = factor(
      sprintf("%.0f%%", 100 * keep_frac),
      levels = sprintf("%.0f%%", 100 * RETENTION)
    )
  )

make_fig <- function(hyp, ref_df, outfile, title) {
  d <- ok_plot %>% filter(hypothesis == hyp)
  obs <- ref_df %>%
    transmute(
      phase = factor(phase, levels = PHASE_V2),
      contrast = full_sample_contrast
    )
  p <- ggplot(d, aes(x = contrast, colour = retention, fill = retention)) +
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
        "GLMM refits; dashed red = full-sample primary_model_v2 IQR contrast; seed = %d",
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

make_fig("H2", ref_h2, path_out_fig_h2, "H2 GLMM rectangle sub-sampling: IQR gap-change by phase")
make_fig("H3", ref_h3, path_out_fig_h3, "H3 GLMM rectangle sub-sampling: IQR gap-change by phase")
logmsg("Saved: ", path_out_fig_h2)
logmsg("Saved: ", path_out_fig_h3)

# ---------------------------------------------------------------------------
# Markdown summary
# ---------------------------------------------------------------------------
md_table <- function(sum_df) {
  c(
    "| Phase | Retention | Full-sample GLMM contrast | Subsample median [2.5, 97.5] | % sign match | % sign+sig match | % converged |",
    "|---|---|---|---|---|---|---|",
    vapply(seq_len(nrow(sum_df)), function(i) {
      r <- sum_df[i, ]
      sprintf(
        "| %s | %.0f%% | %+.2f%% | %+.2f%% [%+.2f, %+.2f] | %.1f%% | %.1f%% | %.1f%% |",
        as.character(r$phase), 100 * r$retention, r$full_sample_contrast,
        r$subsample_median, r$subsample_p025, r$subsample_p975,
        r$pct_sign_match, r$pct_sign_and_sig_match, r$pct_converged
      )
    }, character(1))
  )
}

sub_md <- c(
  "# H2/H3 GLMM rectangle sub-sampling sensitivity",
  "",
  "Refit of `primary_model_v2` (`residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)`) on global rectangle subsamples. Contrast = phase slope × fixed full-sample IQR → % gap change. Same model form as the full-sample reference — no OLS proxy.",
  "",
  "## Batch header",
  "",
  sprintf("- Seed: %d", SEED),
  sprintf("- Iterations per retention level: %d", N_ITER),
  sprintf("- Retention levels: %s", paste(sprintf("%.0f%%", 100 * RETENTION), collapse = ", ")),
  sprintf("- Single-fit benchmark: %.3f sec", t_bench),
  sprintf("- Total batch runtime: %.1f sec (%.2f min)", runtime_sec, runtime_sec / 60),
  sprintf("- Minimum realized thinnest-phase rectangle count at 40%% retention: %d", min_thin_40),
  "",
  "### Convergence / exclusions by retention",
  "",
  "| Retention | Attempted | Converged | Excluded | % excluded | Singular/boundary warns | Min thinnest-phase n_rect |",
  "|---|---|---|---|---|---|---|",
  vapply(seq_len(nrow(conv_summary)), function(i) {
    r <- conv_summary[i, ]
    sprintf(
      "| %.0f%% | %d | %d | %d | %.1f%% | %d | %d |",
      100 * r$keep_frac, r$n_attempted, r$n_converged, r$n_excluded,
      r$pct_excluded, r$n_singular_or_boundary_warn, r$min_n_rect_thinnest_phase
    )
  }, character(1)),
  "",
  "## H2 (`FP_between`)",
  "",
  md_table(sum_h2_out),
  "",
  "## H3 (`FP_within`)",
  "",
  md_table(sum_h3_out),
  "",
  "## Figures",
  "",
  paste0("- ", path_out_fig_h2),
  paste0("- ", path_out_fig_h3),
  "",
  "### Notes",
  "",
  "- Match rule: same sign as full-sample IQR contrast **and** same significance call (p < 0.05) as the full-sample phase slope.",
  "- Non-converged iterations excluded from medians / percentiles / match rates; counted in % converged.",
  "- IQR transforms use full-sample `IQR(FP_between)` and phase-specific `IQR(FP_within)` so subsample variation tracks the slope, not a shifting Δx.",
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
      n_cores = N_CORES,
      retention = RETENTION,
      alpha = ALPHA,
      phases = PHASE_V2,
      formula = form,
      method = "glmmTMB_refit_global_rectangle_subsample",
      single_fit_benchmark_sec = t_bench,
      min_thinnest_phase_n_rect_at_40 = min_thin_40
    ),
    runtime_sec = runtime_sec
  ),
  path_out_rds
)
logmsg("Saved: ", path_out_rds)

writeLines(run_log, path_out_run_log)
logmsg("Saved: ", path_out_run_log)
cat("=== H2/H3 GLMM rectangle sub-sampling complete. ===\n")
