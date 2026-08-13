# Spatial confounding bootstrap on FP_between (H2) — CAR specification
#
# PURPOSE: test whether the reported H2 effect depends on the spatial arrangement
# of rectangle-level FP_between under the primary CAR model used for H2 spatial
# contrasts. Shuffles FP_between once per rectangle (preserving the set of
# values), leaves residual / FP_within / phase_v2 / adjacency / rectangle IDs
# untouched, and refits via spaMM::fitme.
#
# Primary model (policy-anchored phase_v2):
#   residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
#   phases: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015
#
# Adjacency matrix: loaded from the SAME source used by the primary CAR H2 fit
# in run_h2h3_phase_v2_reporting.R —
#   outputs/h2h3_feasibility_round2_model_objects.rds$adjMatrix
# (queen-contiguity binary matrix). Confirmed identical to the copy stored in
# outputs/phase_v2_reporting_model_objects.rds when that artifact exists.
#
# Supersedes the RE (1 | stat_rec) version archived at:
#   exploratory/pipeline/permutation_bootstrap_FP_between_RE.R
#
# Run: Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR.R
#
# Optional env overrides:
#   N_BOOT=1000  SEED=42  N_CORES=1

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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
  stop("Run from pipeline/ or Rscript pipeline/permutation_bootstrap_FP_between_CAR.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

# Prefer project-local library (e.g. .R_libs/spaMM) when present
local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop(
    "Package 'spaMM' required. Run with: ",
    "Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR.R ",
    "(ensure .R_libs/spaMM or ambient spaMM is available)"
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
path_reporting <- file.path(
  project_root, "outputs", "phase_v2_reporting_model_objects.rds"
)
path_re_summary <- file.path(
  project_root, "exploratory", "outputs",
  "permutation_bootstrap_FP_between_RE_summary.md"
)
path_re_objects <- file.path(
  project_root, "exploratory", "outputs",
  "permutation_bootstrap_FP_between_RE_objects.rds"
)

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_results <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_results.csv"
)
path_out_summary <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_summary.md"
)
path_out_fig <- file.path(
  fig_dir, "permutation_bootstrap_FP_between_CAR_null.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_objects.rds"
)
path_out_compare <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_RE_vs_CAR_null_width.csv"
)

N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
n_cores_env <- Sys.getenv("N_CORES", unset = "1")
N_CORES <- as.integer(n_cores_env)
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

PHASE_V2_LEVELS <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
slope_col <- function(phase_label) {
  paste0("slope_", gsub("-", "_", phase_label, fixed = TRUE))
}

na_coef_template <- function(phases) {
  stats::setNames(
    rep(NA_real_, length(phases) + 1L),
    c("FP_between", vapply(phases, slope_col, character(1)))
  )
}

extract_fp_between_coefs_spamm <- function(fit, phases) {
  out <- na_coef_template(phases)
  b <- tryCatch(spaMM::fixef(fit), error = function(e) NULL)
  if (is.null(b) || !"FP_between" %in% names(b)) return(out)

  main <- unname(b[["FP_between"]])
  out[["FP_between"]] <- main
  out[[slope_col(phases[[1]])]] <- main
  if (length(phases) > 1L) {
    for (ph in phases[-1]) {
      int_nm <- tryCatch(
        find_fp_phase_interaction_name(names(b), "FP_between", ph),
        error = function(e) NA_character_
      )
      if (is.na(int_nm) || !int_nm %in% names(b)) {
        out[[slope_col(ph)]] <- NA_real_
      } else {
        out[[slope_col(ph)]] <- main + unname(b[[int_nm]])
      }
    }
  }
  out
}

permute_fp_between_data <- function(data, rectangle_col = "stat_rec",
                                    fp_col = "FP_between") {
  rect_fp_map <- unique(data[, c(rectangle_col, fp_col), drop = FALSE])
  if (nrow(rect_fp_map) != length(unique(data[[rectangle_col]]))) {
    stop("FP_between is not unique per rectangle; cannot permute at rectangle level.")
  }
  shuffled_fp <- rect_fp_map
  shuffled_fp[[fp_col]] <- sample(shuffled_fp[[fp_col]])
  data_perm <- data
  idx <- match(data_perm[[rectangle_col]], shuffled_fp[[rectangle_col]])
  if (anyNA(idx)) stop("Rectangle match failed during FP_between permutation.")
  data_perm[[fp_col]] <- shuffled_fp[[fp_col]][idx]
  data_perm
}

fit_car_on_data <- function(data, model_formula, adjMatrix) {
  tryCatch(
    spaMM::fitme(
      model_formula,
      data = data,
      adjMatrix = adjMatrix,
      method = "REML"
    ),
    error = function(e) structure(list(error = conditionMessage(e)), class = "car_fit_fail")
  )
}

is_car_fit_ok <- function(fit) {
  if (inherits(fit, "car_fit_fail") || is.null(fit)) return(FALSE)
  b <- tryCatch(spaMM::fixef(fit), error = function(e) NULL)
  is.finite(b[["FP_between"]])
}

permute_and_refit_car <- function(data, model_formula, adjMatrix, phases,
                                  rectangle_col = "stat_rec",
                                  fp_col = "FP_between") {
  data_perm <- permute_fp_between_data(data, rectangle_col, fp_col)
  fit <- fit_car_on_data(data_perm, model_formula, adjMatrix)
  if (!is_car_fit_ok(fit)) {
    return(list(
      coefs = na_coef_template(phases),
      failed = TRUE,
      fail_reason = if (inherits(fit, "car_fit_fail")) fit$error else "non-finite FP_between"
    ))
  }
  list(
    coefs = extract_fp_between_coefs_spamm(fit, phases),
    failed = FALSE,
    fail_reason = NA_character_
  )
}

empirical_p <- function(null_vals, obs) {
  x <- null_vals[is.finite(null_vals)]
  if (!length(x) || !is.finite(obs)) return(NA_real_)
  (sum(abs(x) >= abs(obs)) + 1) / (length(x) + 1)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return(sprintf("%.4g (< 0.001)", p))
  sprintf("%.4g", p)
}

null_summary_stats <- function(null_vals) {
  x <- null_vals[is.finite(null_vals)]
  c(
    mean = mean(x),
    sd = stats::sd(x),
    iqr = as.numeric(stats::IQR(x)),
    q025 = as.numeric(stats::quantile(x, 0.025)),
    q975 = as.numeric(stats::quantile(x, 0.975))
  )
}

interpret_h2_spatial <- function(obs, q025, q975, null_mean) {
  inside <- is.finite(obs) && obs >= q025 && obs <= q975
  direction <- if (!is.finite(obs) || !is.finite(null_mean)) {
    "indeterminate"
  } else if (abs(obs) > abs(null_mean)) {
    "permutation null closer to zero than observed (spatial arrangement strengthens |coef|)"
  } else if (abs(obs) < abs(null_mean)) {
    "permutation null farther from zero than observed (spatial arrangement weakens |coef|)"
  } else {
    "similar magnitude"
  }
  statement <- if (inside) {
    paste0(
      "H2 effect does not depend on the spatial arrangement of FP_between. ",
      "No evidence the CAR spatial intercept and FP_between are substituting for ",
      "one another. Reported effect is not an artifact of spatial confounding."
    )
  } else {
    paste0(
      "The CAR spatial intercept and FP_between are competing for the same spatial ",
      "signal when FP_between's real spatial arrangement is present. The reported ",
      "H2 effect is at least partly attributable to this overlap, not purely to ",
      "the fishing-pressure effect itself."
    )
  }
  list(inside_95 = inside, direction = direction, statement = statement)
}

load_adjMatrix_primary_car <- function(path_round2, path_reporting) {
  if (!file.exists(path_round2)) {
    stop(
      "Primary CAR adjacency source missing: ", path_round2,
      " (same file used by run_h2h3_phase_v2_reporting.R)"
    )
  }
  round2 <- readRDS(path_round2)
  if (is.null(round2$adjMatrix)) {
    stop("Round 2 RDS missing adjMatrix: ", path_round2)
  }
  adj <- round2$adjMatrix
  confirm_note <- "loaded from h2h3_feasibility_round2_model_objects.rds$adjMatrix"
  if (file.exists(path_reporting)) {
    rep <- readRDS(path_reporting)
    if (!is.null(rep$adjMatrix)) {
      same <- identical(adj, rep$adjMatrix)
      if (!same) {
        stop(
          "adjMatrix in round2 RDS differs from phase_v2_reporting_model_objects.rds; ",
          "refusing to proceed with mismatched adjacency."
        )
      }
      confirm_note <- paste0(
        confirm_note,
        "; confirmed identical to phase_v2_reporting_model_objects.rds$adjMatrix"
      )
    }
  }
  list(adjMatrix = adj, source = path_round2, confirm_note = confirm_note)
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# FP_between spatial confounding bootstrap (CAR) — run log")
logmsg("")
logmsg(
  "Rectangle-level permutation of FP_between under the primary CAR within-between ",
  "model. Destroys the spatial arrangement of the H2 covariate while leaving ",
  "residual, FP_within, phase_v2, rectangle IDs, and the queen-contiguity ",
  "adjacency matrix untouched."
)
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf("spaMM %s", as.character(utils::packageVersion("spaMM"))))
logmsg(sprintf("N_BOOT = %d; SEED = %d; N_CORES = %d", N_BOOT, SEED, N_CORES))

# ---------------------------------------------------------------------------
# Load data + adjacency + observed CAR fit
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
stopifnot(file.exists(path_models_v2))
v2 <- readRDS(path_models_v2)
if (is.null(v2$data)) stop("primary_model_v2.rds missing data.")

phases <- if (!is.null(v2$phase_v2_labels)) {
  as.character(v2$phase_v2_labels)
} else {
  PHASE_V2_LEVELS
}
phase_col <- "phase_v2"
slope_names <- unname(vapply(phases, slope_col, character(1)))
coef_names <- c("FP_between", slope_names)

adj_pack <- load_adjMatrix_primary_car(path_round2, path_reporting)
adjMatrix <- adj_pack$adjMatrix
logmsg("Adjacency: ", adj_pack$confirm_note)
logmsg(sprintf(
  "adjMatrix dim = %d x %d; class = %s",
  nrow(adjMatrix), ncol(adjMatrix), paste(class(adjMatrix), collapse = "/")
))

formula_car <- residual ~ FP_between * phase_v2 + FP_within * phase_v2 +
  adjacency(1 | stat_rec)

dat <- v2$data
dat$stat_rec <- factor(as.character(dat$stat_rec), levels = rownames(adjMatrix))
if (anyNA(dat$stat_rec)) {
  stop("Some haul stat_rec values are missing from adjMatrix rownames.")
}
if (!phase_col %in% names(dat)) stop("Analysis data missing phase_v2.")
dat[[phase_col]] <- factor(dat[[phase_col]], levels = phases)
if (anyNA(dat[[phase_col]])) stop("NA in phase_v2 after relevel.")

n_hauls <- nrow(dat)
n_rect <- nlevels(dat$stat_rec)
fp_map <- unique(dat[, c("stat_rec", "FP_between")])
stopifnot(nrow(fp_map) == n_rect)

logmsg("Loaded data from: ", path_models_v2)
logmsg(sprintf("Phase scheme: phase_v2 (column `%s`)", phase_col))
logmsg(sprintf("Phase levels: %s", paste(phases, collapse = " | ")))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  n_hauls, n_rect, min(dat$year), max(dat$year)
))
logmsg("Formula: ", paste(deparse(formula_car), collapse = " "))
logmsg(sprintf(
  "FP_between unique per rectangle: %d values; range [%.3f, %.3f]",
  nrow(fp_map), min(fp_map$FP_between), max(fp_map$FP_between)
))

# Observed CAR fit: prefer saved primary CAR; else refit once
fit_obs <- NULL
obs_source <- NA_character_
if (!is.null(v2$fit_wb_car_v2)) {
  fit_obs <- v2$fit_wb_car_v2
  obs_source <- "primary_model_v2.rds$fit_wb_car_v2"
} else if (file.exists(path_reporting)) {
  rep <- readRDS(path_reporting)
  if (!is.null(rep$fit_wb_car_v2)) {
    fit_obs <- rep$fit_wb_car_v2
    obs_source <- "phase_v2_reporting_model_objects.rds$fit_wb_car_v2"
  }
}
if (is.null(fit_obs)) {
  logmsg("Saved CAR fit missing — refitting observed CAR once for baseline coefficients.")
  fit_obs <- spaMM::fitme(
    formula_car, data = dat, adjMatrix = adjMatrix, method = "REML"
  )
  obs_source <- "refit in this script (same formula + adjMatrix)"
}
logmsg("Observed CAR source: ", obs_source)

# Sanity: permuting FP_between must not alter FP_within / phase / residual / levels
set.seed(SEED)
dat_check <- permute_fp_between_data(dat)
stopifnot(identical(dat_check$FP_within, dat$FP_within))
stopifnot(identical(dat_check$residual, dat$residual))
stopifnot(identical(as.character(dat_check$stat_rec), as.character(dat$stat_rec)))
stopifnot(identical(as.character(dat_check[[phase_col]]), as.character(dat[[phase_col]])))
stopifnot(identical(levels(dat_check$stat_rec), levels(dat$stat_rec)))
stopifnot(!isTRUE(all.equal(dat_check$FP_between, dat$FP_between)))
logmsg(
  "Permutation sanity check passed: only FP_between changes; ",
  "FP_within / residual / stat_rec / ", phase_col, " / adjMatrix untouched."
)

# ---------------------------------------------------------------------------
# Step 1 — observed H2 coefficients (CAR)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1 — Observed H2 coefficients (CAR baseline)")
obs_coefs <- extract_fp_between_coefs_spamm(fit_obs, phases)
obs_main <- unname(obs_coefs[["FP_between"]])

obs_slopes_tab <- extract_wb_phase_slopes_spamm(
  fit_obs,
  term_name = "FP_between",
  model_id = "wb_car_v2",
  hypothesis_group = "H2_spatial_between",
  phases = phases
)
obs_fe <- obs_slopes_tab[1, ]

logmsg(sprintf(
  "Observed fixef FP_between (reference phase %s): %+0.6f",
  phases[[1]], obs_main
))
logmsg(sprintf(
  "  SE = %.6f; 95%% CI [%.6f, %.6f]; z = %.3f; p = %.4g",
  obs_fe$fp_slope_se, obs_fe$fp_slope_lo, obs_fe$fp_slope_hi,
  obs_fe$statistic, obs_fe$p_value
))
for (i in seq_len(nrow(obs_slopes_tab))) {
  r <- obs_slopes_tab[i, ]
  logmsg(sprintf(
    "  Phase %s slope = %+0.6f (SE %.6f; 95%% CI [%.6f, %.6f]; p = %.4g)",
    r$phase, r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi, r$p_value
  ))
}

# ---------------------------------------------------------------------------
# Steps 2–3 — permutation loop
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Steps 2–3 — Rectangle-level FP_between permutation + CAR refit")
logmsg(
  "Each replicate: sample(FP_between) across the 158 rectangles once, reassign ",
  "to all hauls in that rectangle, refit CAR with spaMM::fitme(..., method = REML)."
)
logmsg("")
logmsg("### Matched permutation procedure (explicit)")
logmsg(
  "For each replicate: (1) one global shuffle of the rectangle-level FP_between ",
  "values; (2) one refit of the full CAR formula with fixed adjMatrix; (3) one ",
  "fixef() extract from that single fit, recording the reference-phase ",
  "FP_between main effect and the FP_between × phase_v2 interaction increments ",
  "(combined to phase slopes as main + interaction). Not four independent ",
  "shuffle-and-refit procedures. The resulting null distributions are columns ",
  "of one n_boot × K matrix — matched randomizations."
)

t_start <- proc.time()[[3]]
set.seed(SEED)

fail_reasons <- rep(NA_character_, N_BOOT)

if (N_CORES > 1L && .Platform$OS.type == "unix") {
  logmsg(sprintf(
    "Running with parallel::mclapply, mc.cores = %d (spaMM; watch for fork issues)",
    N_CORES
  ))
  seeds <- SEED + seq_len(N_BOOT)
  null_list <- parallel::mclapply(seq_len(N_BOOT), function(i) {
    set.seed(seeds[[i]])
    permute_and_refit_car(dat, formula_car, adjMatrix, phases)
  }, mc.cores = N_CORES)
  null_mat <- do.call(rbind, lapply(null_list, function(x) x$coefs))
  fail_reasons <- vapply(null_list, function(x) {
    if (isTRUE(x$failed)) as.character(x$fail_reason) else NA_character_
  }, character(1))
} else {
  logmsg("Running sequentially (replicate)")
  null_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = length(coef_names))
  colnames(null_mat) <- coef_names
  for (i in seq_len(N_BOOT)) {
    res <- permute_and_refit_car(dat, formula_car, adjMatrix, phases)
    null_mat[i, ] <- res$coefs
    if (isTRUE(res$failed)) fail_reasons[[i]] <- res$fail_reason
    if (i %% 50L == 0L || i == N_BOOT) {
      cat(sprintf("  ... completed %d / %d\n", i, N_BOOT))
    }
  }
}
if (is.null(dim(null_mat))) {
  null_mat <- matrix(null_mat, nrow = N_BOOT, byrow = TRUE)
}
colnames(null_mat) <- coef_names

t_end <- proc.time()[[3]]
runtime_sec <- t_end - t_start

failed_idx <- which(!is.finite(null_mat[, "FP_between"]))
n_failed <- length(failed_idx)
n_ok <- N_BOOT - n_failed
logmsg(sprintf(
  "Finished: n_boot = %d; n_failed = %d; n_ok = %d; runtime = %.1f sec (%.2f sec/rep)",
  N_BOOT, n_failed, n_ok, runtime_sec, runtime_sec / N_BOOT
))
if (n_failed > 0L) {
  logmsg(
    "FAILED iterations (NOT silently dropped from n_boot; recorded as NA in null matrix):"
  )
  for (i in failed_idx) {
    logmsg(sprintf(
      "  perm_id=%d reason=%s",
      i, ifelse(is.na(fail_reasons[[i]]), "non-finite coef", fail_reasons[[i]])
    ))
  }
} else {
  logmsg("All iterations produced finite FP_between coefficients.")
}
if (runtime_sec > 1800) {
  logmsg(
    "FLAG: runtime > 30 min. Consider reducing N_BOOT or increasing N_CORES ",
    "(env vars N_BOOT / N_CORES)."
  )
} else {
  logmsg(sprintf(
    paste0(
      "Runtime practicality: %.1f min for %d CAR refits (%.2f sec/rep) — ",
      "1000 iterations is practical at this wall-clock cost."
    ),
    runtime_sec / 60, N_BOOT, runtime_sec / N_BOOT
  ))
}

# ---------------------------------------------------------------------------
# Step 4 — compare observed to null
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 4 — Observed vs permutation null (CAR)")

targets <- c(
  FP_between = sprintf(
    "FP_between (reference-phase main effect, %s)", phases[[1]]
  ),
  stats::setNames(
    paste0("Phase ", phases, " FP_between slope"),
    slope_names
  )
)

summary_rows <- lapply(names(targets), function(nm) {
  obs <- unname(obs_coefs[[nm]])
  null_x <- null_mat[, nm]
  null_clean <- null_x[is.finite(null_x)]
  ns <- null_summary_stats(null_clean)
  p_emp <- empirical_p(null_clean, obs)
  inter <- interpret_h2_spatial(obs, ns[["q025"]], ns[["q975"]], ns[["mean"]])
  logmsg(sprintf(
    "%s: obs=%+.6f; null mean=%+.6f sd=%.6f iqr=%.6f; 95%% null [%.6f, %.6f]; p_emp=%s; inside_95=%s",
    targets[[nm]], obs, ns[["mean"]], ns[["sd"]], ns[["iqr"]],
    ns[["q025"]], ns[["q975"]],
    fmt_p(p_emp), inter$inside_95
  ))
  logmsg(sprintf("  Direction: %s", inter$direction))
  logmsg(sprintf("  Interpretation: %s", inter$statement))
  data.frame(
    target = nm,
    label = targets[[nm]],
    phase = NA_character_,
    obs_coef = obs,
    null_mean = unname(ns[["mean"]]),
    null_sd = unname(ns[["sd"]]),
    null_iqr = unname(ns[["iqr"]]),
    null_q025 = unname(ns[["q025"]]),
    null_q975 = unname(ns[["q975"]]),
    p_empirical = p_emp,
    inside_null_95 = inter$inside_95,
    spatially_confounded = !inter$inside_95,
    direction = inter$direction,
    interpretation = inter$statement,
    stringsAsFactors = FALSE
  )
})
summary_df <- dplyr::bind_rows(summary_rows)

# Cleaner phase labels for slope rows
summary_df$phase <- dplyr::case_when(
  summary_df$target == "FP_between" ~ phases[[1]],
  summary_df$target %in% slope_names ~ phases[match(summary_df$target, slope_names)],
  TRUE ~ summary_df$phase
)

any_outside <- any(summary_df$target %in% slope_names & !summary_df$inside_null_95)

# Presentation table (phase slopes only)
phase_table <- summary_df %>%
  dplyr::filter(target %in% slope_names) %>%
  dplyr::mutate(
    phase_label = phase,
    confounded_lab = ifelse(spatially_confounded, "Yes", "No")
  )

# ---------------------------------------------------------------------------
# RE vs CAR null-width comparison
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## RE vs CAR null-width comparison")

compare_df <- NULL
if (file.exists(path_re_objects)) {
  re_obj <- readRDS(path_re_objects)
  re_null <- re_obj$null_mat
  re_summary <- re_obj$summary
  compare_rows <- lapply(slope_names, function(nm) {
    ph <- phases[match(nm, slope_names)]
    car_row <- summary_df[summary_df$target == nm, , drop = FALSE]
    re_x <- re_null[, nm]
    re_x <- re_x[is.finite(re_x)]
    re_ns <- null_summary_stats(re_x)
    re_inside <- if (!is.null(re_summary) && nm %in% re_summary$target) {
      isTRUE(re_summary$inside_null_95[re_summary$target == nm][[1]])
    } else {
      NA
    }
    data.frame(
      phase = ph,
      car_obs = car_row$obs_coef[[1]],
      car_null_sd = car_row$null_sd[[1]],
      car_null_iqr = car_row$null_iqr[[1]],
      car_null_q025 = car_row$null_q025[[1]],
      car_null_q975 = car_row$null_q975[[1]],
      car_inside_95 = car_row$inside_null_95[[1]],
      car_spatially_confounded = car_row$spatially_confounded[[1]],
      re_obs = if (!is.null(re_obj$obs_coefs)) unname(re_obj$obs_coefs[[nm]]) else NA_real_,
      re_null_sd = unname(re_ns[["sd"]]),
      re_null_iqr = unname(re_ns[["iqr"]]),
      re_null_q025 = unname(re_ns[["q025"]]),
      re_null_q975 = unname(re_ns[["q975"]]),
      re_inside_95 = re_inside,
      re_spatially_confounded = if (is.na(re_inside)) NA else !re_inside,
      sd_ratio_car_over_re = car_row$null_sd[[1]] / unname(re_ns[["sd"]]),
      iqr_ratio_car_over_re = car_row$null_iqr[[1]] / unname(re_ns[["iqr"]]),
      stringsAsFactors = FALSE
    )
  })
  compare_df <- dplyr::bind_rows(compare_rows)
  readr::write_csv(compare_df, path_out_compare)
  logmsg("Wrote: ", path_out_compare)
  for (i in seq_len(nrow(compare_df))) {
    r <- compare_df[i, ]
    logmsg(sprintf(
      "  %s: RE null SD=%.6f IQR=%.6f | CAR null SD=%.6f IQR=%.6f | SD ratio CAR/RE=%.3f",
      r$phase, r$re_null_sd, r$re_null_iqr, r$car_null_sd, r$car_null_iqr,
      r$sd_ratio_car_over_re
    ))
  }
} else {
  logmsg(
    "RE archive objects not found at ", path_re_objects,
    " — skipping side-by-side null-width table."
  )
}

# ---------------------------------------------------------------------------
# Deliverables
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Deliverables")

results_df <- data.frame(
  record_type = "permutation",
  perm_id = seq_len(N_BOOT),
  converged = is.finite(null_mat[, "FP_between"]),
  fail_reason = fail_reasons,
  stringsAsFactors = FALSE
)
for (nm in coef_names) {
  results_df[[nm]] <- as.numeric(null_mat[, nm])
}

meta_row <- data.frame(
  record_type = "metadata",
  perm_id = NA_integer_,
  converged = NA,
  fail_reason = NA_character_,
  phase_scheme = "phase_v2",
  phase_col = phase_col,
  phases = paste(phases, collapse = "|"),
  n_boot = N_BOOT,
  n_failed = n_failed,
  n_ok = n_ok,
  seed = SEED,
  runtime_sec = runtime_sec,
  n_cores = N_CORES,
  adjMatrix_source = adj_pack$source,
  stringsAsFactors = FALSE
)
for (nm in coef_names) meta_row[[nm]] <- NA_real_

# Keep ALL iterations (failed = converged FALSE with NA coefs) — do not drop
results_out <- dplyr::bind_rows(meta_row, results_df)
readr::write_csv(results_out, path_out_results)
logmsg("Wrote: ", path_out_results)
logmsg(sprintf(
  "Results CSV retains all %d iterations (converged + failed); failed not dropped.",
  N_BOOT
))

label_map <- c(
  FP_between = paste0("Main effect\n(ref. ", phases[[1]], ")"),
  stats::setNames(paste0("Phase\n", phases), slope_names)
)
results_ok <- results_df %>% dplyr::filter(converged)
plot_df <- tidyr::pivot_longer(
  results_ok,
  cols = tidyselect::all_of(coef_names),
  names_to = "target",
  values_to = "coef"
) %>%
  dplyr::mutate(
    target_label = factor(
      unname(label_map[target]),
      levels = unname(label_map[coef_names])
    )
  )

obs_plot <- data.frame(
  target = names(obs_coefs),
  obs = as.numeric(obs_coefs),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    target_label = factor(
      unname(label_map[target]),
      levels = levels(plot_df$target_label)
    )
  )

p <- ggplot(plot_df, aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey75", colour = "white") +
  geom_density(colour = "grey25", linewidth = 0.6) +
  geom_vline(
    data = obs_plot, aes(xintercept = obs),
    colour = "#B2182B", linewidth = 0.9, linetype = "solid"
  ) +
  facet_wrap(~target_label, scales = "free", nrow = 2) +
  labs(
    title = "FP_between spatial permutation null (CAR primary model)",
    subtitle = sprintf(
      "Red line = observed coefficient; n_boot = %d, n_failed = %d, seed = %d; phases = %s",
      N_BOOT, n_failed, SEED, paste(phases, collapse = " / ")
    ),
    x = "Coefficient under rectangle-level FP_between permutation",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    panel.grid.minor = element_blank()
  )
ggplot2::ggsave(path_out_fig, p, width = 11, height = 7, dpi = 150)
logmsg("Wrote: ", path_out_fig)

form_txt <- paste(deparse(formula_car), collapse = " ")
summary_lines <- c(
  "# FP_between spatial confounding bootstrap (CAR) — summary",
  "",
  "Model-appropriate spatial confounding check for H2 under the **CAR** primary",
  "specification used for reported H2 spatial contrasts.",
  "",
  sprintf("- Formula: `%s` (REML, spaMM::fitme)", form_txt),
  "- Phase scheme: **phase_v2** (`phase_v2`)",
  sprintf("- Phase levels: %s", paste(phases, collapse = " / ")),
  sprintf("- Adjacency: %s", adj_pack$confirm_note),
  "",
  "The superseded RE `(1 | stat_rec)` bootstrap is archived at",
  "`exploratory/pipeline/permutation_bootstrap_FP_between_RE.R` /",
  "`exploratory/outputs/permutation_bootstrap_FP_between_RE_*`.",
  "",
  "## Settings",
  "",
  sprintf("- Seed: `%d`", SEED),
  sprintf("- `n_boot`: %d", N_BOOT),
  sprintf("- `n_failed` (error / non-finite; **not dropped**): %d", n_failed),
  sprintf("- `n_ok`: %d", n_ok),
  sprintf("- Cores: %d", N_CORES),
  sprintf("- Runtime: %.1f sec (%.2f sec/replicate; %.2f min total)",
          runtime_sec, runtime_sec / N_BOOT, runtime_sec / 60),
  sprintf("- Rectangles: %d; hauls: %d", n_rect, n_hauls),
  sprintf("- Observed CAR source: `%s`", obs_source),
  sprintf("- adjMatrix source: `%s`", adj_pack$source),
  "",
  "### Matched permutation procedure",
  "",
  "For each of the replicates: **one global shuffle** of rectangle-level",
  "`FP_between` → **one CAR refit** → **one `fixef()` extract** recording the",
  "reference-phase main effect and the three `FP_between:phase_v2` interaction",
  "increments (phase slopes = main + interaction). Adjacency matrix fixed.",
  "",
  "## Observed H2 coefficients (CAR)",
  "",
  sprintf(
    "- Reference-phase main effect `FP_between` (%s): **%+.6f** (SE %.6f; 95%% CI [%.6f, %.6f]; p = %.4g)",
    phases[[1]], obs_fe$fp_slope, obs_fe$fp_slope_se, obs_fe$fp_slope_lo,
    obs_fe$fp_slope_hi, obs_fe$p_value
  ),
  "",
  "Phase-specific `FP_between` slopes (presented H2 effects under CAR):",
  ""
)
for (i in seq_len(nrow(obs_slopes_tab))) {
  r <- obs_slopes_tab[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "- %s: **%+.6f** (SE %.6f; 95%% CI [%.6f, %.6f]; p = %.4g)",
      r$phase, r$fp_slope, r$fp_slope_se, r$fp_slope_lo, r$fp_slope_hi, r$p_value
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "## Presentation table (CAR)",
  "",
  "| Phase | H2 (`FP_between`) coefficient — CAR | Spatially confounded? |",
  "|---|---|---|"
)
for (i in seq_len(nrow(phase_table))) {
  r <- phase_table[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.3f | %s |",
      r$phase_label, r$obs_coef, r$confounded_lab
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "**Spatially confounded** = Yes if the observed coefficient falls outside the",
  "null distribution's 95% interval.",
  "",
  "## Null distribution vs observed (CAR)",
  "",
  "| Target | Observed | Null mean | Null SD | Null IQR | Null 2.5% | Null 97.5% | Empirical p | Inside null 95%? |",
  "|--------|----------|-----------|---------|----------|-----------|------------|-------------|------------------|"
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.6f | %+.6f | %.6f | %.6f | %.6f | %.6f | %s | %s |",
      r$label, r$obs_coef, r$null_mean, r$null_sd, r$null_iqr, r$null_q025, r$null_q975,
      fmt_p(r$p_empirical), ifelse(r$inside_null_95, "yes", "no")
    )
  )
}

if (!is.null(compare_df)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## RE vs CAR null-width comparison (phase slopes)",
    "",
    "Archived RE bootstrap (`exploratory/outputs/permutation_bootstrap_FP_between_RE_*`)",
    "alongside this CAR run — same four phases, same reference phase (1985–1991).",
    "",
    "| Phase | RE obs | CAR obs | RE null SD | CAR null SD | RE null IQR | CAR null IQR | SD ratio (CAR/RE) | RE confounded? | CAR confounded? |",
    "|-------|--------|---------|------------|-------------|-------------|--------------|-------------------|----------------|-----------------|"
  )
  for (i in seq_len(nrow(compare_df))) {
    r <- compare_df[i, ]
    summary_lines <- c(
      summary_lines,
      sprintf(
        "| %s | %+.3f | %+.3f | %.6f | %.6f | %.6f | %.6f | %.3f | %s | %s |",
        r$phase, r$re_obs, r$car_obs, r$re_null_sd, r$car_null_sd,
        r$re_null_iqr, r$car_null_iqr, r$sd_ratio_car_over_re,
        ifelse(isTRUE(r$re_spatially_confounded), "Yes", "No"),
        ifelse(isTRUE(r$car_spatially_confounded), "Yes", "No")
      )
    )
  }
}

summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation",
  "",
  "Rule:",
  "",
  "- Observed coefficient **inside** the permutation null 95% interval → H2 effect does",
  "  not depend on the spatial arrangement of `FP_between`; no evidence that the",
  "  CAR spatial term and `FP_between` are substituting for one another.",
  "- Observed coefficient **outside** the null 95% interval → CAR spatial term and",
  "  `FP_between` compete for the same spatial signal under the real map; reported H2",
  "  is at least partly attributable to that overlap.",
  "",
  "Per-target statements:",
  ""
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf("### %s", r$label),
    "",
    sprintf("- Direction: %s", r$direction),
    sprintf("- Empirical p = %s; inside null 95%% = %s", fmt_p(r$p_empirical), r$inside_null_95),
    sprintf("- %s", r$interpretation),
    ""
  )
}

headline <- if (any_outside) {
  paste0(
    "**Headline (phase-specific H2 slopes, phase_v2, CAR):** at least one ",
    "presented phase-specific `FP_between` slope falls outside its spatial-permutation ",
    "null 95% interval under the CAR model, indicating that the reported H2 effect ",
    "depends on the real spatial arrangement of fishing pressure and is at least ",
    "partly entangled with the CAR spatial term's signal."
  )
} else {
  paste0(
    "**Headline (phase-specific H2 slopes, phase_v2, CAR):** all presented ",
    "phase-specific `FP_between` slopes fall inside their spatial-permutation null ",
    "95% intervals under the CAR model. No evidence that the reported H2 effects ",
    "are an artifact of spatial confounding between `FP_between` and ",
    "`adjacency(1 | stat_rec)`."
  )
}
summary_lines <- c(summary_lines, "## Headline", "", headline, "")
summary_lines <- c(
  summary_lines,
  "## Outputs",
  "",
  sprintf("- Results CSV: `%s`", path_out_results),
  sprintf("- Null figure: `%s`", path_out_fig),
  sprintf("- Run log: `%s`", path_out_run_log),
  sprintf("- Session info: `%s`", path_out_session),
  sprintf("- RE vs CAR null-width CSV: `%s`", path_out_compare),
  sprintf("- Archived RE bootstrap: `%s`", path_re_summary),
  ""
)
writeLines(summary_lines, path_out_summary)
logmsg("Wrote: ", path_out_summary)

saveRDS(
  list(
    seed = SEED,
    n_boot = N_BOOT,
    n_failed = n_failed,
    n_ok = n_ok,
    n_cores = N_CORES,
    runtime_sec = runtime_sec,
    fail_reasons = fail_reasons,
    failed_perm_ids = failed_idx,
    adjMatrix_source = adj_pack$source,
    adjMatrix_confirm = adj_pack$confirm_note,
    obs_source = obs_source,
    scheme = "phase_v2",
    phase_col = phase_col,
    phases = phases,
    formula = formula_car,
    obs_coefs = obs_coefs,
    obs_slopes = obs_slopes_tab,
    null_mat = null_mat,
    summary = summary_df,
    phase_table = phase_table,
    compare_re_car = compare_df,
    results_all = results_df
  ),
  path_out_rds
)
logmsg("Wrote: ", path_out_rds)

logmsg("")
logmsg("## Headline")
logmsg(headline)

writeLines(run_log, path_out_run_log)
cat("Run log written to: ", path_out_run_log, "\n", sep = "")
