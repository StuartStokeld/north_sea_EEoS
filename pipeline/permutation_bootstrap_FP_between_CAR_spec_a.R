# Spec A on CAR — Task 2 permutation bootstrap
#
# Base model (from Task 1, gate CLEAR):
#   residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2
#              + FP_within * phase_v2 + adjacency(1 | stat_rec)
#
# Permutation logic (must recompute lag — not a formula swap on CAR-alone):
#   1. Shuffle rectangle-level FP_between
#   2. Recompute FP_between_lag from shuffled map via fixed KNN(k=4)
#   3. Refit CAR+lag; extract phase-specific FP_between slopes
#
# Run: Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R
# Optional: N_BOOT=1000 SEED=42 N_CORES=4

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
  stop("Run from pipeline/ or Rscript pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir

local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))

if (!requireNamespace("spaMM", quietly = TRUE)) {
  stop("Package 'spaMM' required (check .R_libs/spaMM).")
}
suppressPackageStartupMessages(library(spaMM))

PHASE_V2_LEVELS <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

path_car_lag <- file.path(project_root, "outputs", "primary_model_v2_spec_a_car.rds")
path_gate <- file.path(project_root, "outputs", "spec_a_car_identifiability.md")
path_knn <- file.path(project_root, "outputs", "knn_listw_k4.rds")
path_car_alone_summary <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_summary.md"
)

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_results <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_spec_a_results.csv"
)
path_out_summary <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_spec_a_summary.md"
)
path_out_fig <- file.path(
  fig_dir, "permutation_bootstrap_FP_between_CAR_spec_a_null.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_spec_a_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_spec_a_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_spec_a_objects.rds"
)
path_out_compare <- file.path(
  project_root, "outputs", "permutation_bootstrap_FP_between_CAR_vs_CAR_spec_a.csv"
)

N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
N_CORES <- as.integer(Sys.getenv("N_CORES", unset = "4"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

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
    stop("FP_between is not unique per rectangle.")
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
      model_formula, data = data, adjMatrix = adjMatrix, method = "REML"
    ),
    error = function(e) structure(list(error = conditionMessage(e)), class = "car_fit_fail")
  )
}

is_car_fit_ok <- function(fit) {
  if (inherits(fit, "car_fit_fail") || is.null(fit)) return(FALSE)
  b <- tryCatch(spaMM::fixef(fit), error = function(e) NULL)
  is.finite(b[["FP_between"]])
}

# Manual lag check against helper
manual_lag_from_map <- function(fp_map, nb, region_id) {
  fp_map <- fp_map[match(region_id, normalize_stat_rec(fp_map$stat_rec)), , drop = FALSE]
  compute_spatial_lag(fp_map$FP_between, nb)
}

permute_and_refit_car_spec_a <- function(data, model_formula, adjMatrix, phases, nb) {
  data_perm <- permute_fp_between_data(data)
  data_perm <- recompute_fp_between_lag_on_data(data_perm, nb)
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
    q025 = as.numeric(stats::quantile(x, 0.025)),
    q975 = as.numeric(stats::quantile(x, 0.975))
  )
}

# ---------------------------------------------------------------------------
logmsg("# FP_between spatial confounding bootstrap — CAR Spec A — run log")
logmsg("")
logmsg(
  "Task 2: permute FP_between, recompute FP_between_lag from KNN(k=4), ",
  "refit CAR+lag. Gate CLEAR required from Task 1."
)

if (!file.exists(path_car_lag)) {
  stop("Task 1 artifact missing: ", path_car_lag)
}
gate_txt <- if (file.exists(path_gate)) paste(readLines(path_gate, warn = FALSE), collapse = "\n") else ""
if (!grepl("GATE: CLEAR", gate_txt, fixed = TRUE)) {
  stop("Task 1 gate is not CLEAR — do not run Task 2. See ", path_gate)
}

sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg(sprintf("spaMM %s; N_BOOT=%d; SEED=%d; N_CORES=%d",
               as.character(utils::packageVersion("spaMM")), N_BOOT, SEED, N_CORES))

obj <- readRDS(path_car_lag)
fit_obs <- obj$fit_car_lag
dat <- obj$data
adjMatrix <- obj$adjMatrix
form <- obj$formula
phases <- PHASE_V2_LEVELS
knn <- readRDS(path_knn)
nb <- knn$nb
region_id <- normalize_stat_rec(attr(nb, "region.id"))

dat$stat_rec <- factor(as.character(dat$stat_rec), levels = rownames(adjMatrix))
dat$phase_v2 <- factor(dat$phase_v2, levels = phases)

slope_names <- unname(vapply(phases, slope_col, character(1)))
coef_names <- c("FP_between", slope_names)

logmsg("")
logmsg("## Inputs")
logmsg("CAR+lag fit: ", path_car_lag)
logmsg("k-NN: ", path_knn)
logmsg("Formula: ", paste(deparse(form), collapse = " "))
logmsg(sprintf("n=%d hauls; %d rectangles", nrow(dat), nlevels(dat$stat_rec)))

# ---------------------------------------------------------------------------
# Spot-check lag recomputation
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Lag recomputation spot-check")
set.seed(SEED)
dat_check <- permute_fp_between_data(dat)
dat_check <- recompute_fp_between_lag_on_data(dat_check, nb)
fp_map <- unique(as.data.frame(dat_check)[, c("stat_rec", "FP_between")])
fp_map$stat_rec <- normalize_stat_rec(fp_map$stat_rec)
manual <- manual_lag_from_map(fp_map, nb, region_id)
lag_from_data <- dat_check$FP_between_lag[
  match(region_id, normalize_stat_rec(dat_check$stat_rec))
]
# unique match: first haul per rectangle
lag_from_data <- vapply(region_id, function(id) {
  dat_check$FP_between_lag[normalize_stat_rec(dat_check$stat_rec) == id][1]
}, numeric(1))
max_abs_diff <- max(abs(manual - lag_from_data), na.rm = TRUE)
logmsg(sprintf("max |helper lag − manual lag| after one shuffle = %.3g", max_abs_diff))
if (max_abs_diff > 1e-10) stop("Lag recomputation spot-check FAILED")
stopifnot(!isTRUE(all.equal(dat_check$FP_between, dat$FP_between)))
stopifnot(!isTRUE(all.equal(dat_check$FP_between_lag, dat$FP_between_lag)))
stopifnot(identical(dat_check$FP_within, dat$FP_within))
logmsg("Spot-check passed: lag recomputed from permuted FP_between; FP_within unchanged.")

# ---------------------------------------------------------------------------
# Observed coefficients
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Observed coefficients (CAR+lag)")
obs_coefs <- extract_fp_between_coefs_spamm(fit_obs, phases)
for (ph in phases) {
  logmsg(sprintf("  %s FP_between slope = %+.6f", ph, obs_coefs[[slope_col(ph)]]))
}

# ---------------------------------------------------------------------------
# Runtime probe (5 reps) before full run
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Runtime probe (5 replicates)")
t_probe <- system.time({
  set.seed(SEED)
  for (i in 1:5) {
    invisible(permute_and_refit_car_spec_a(dat, form, adjMatrix, phases, nb))
  }
})
sec_per_rep <- t_probe[["elapsed"]] / 5
proj_sec <- sec_per_rep * N_BOOT / max(1, N_CORES)
logmsg(sprintf(
  "Probe: %.2f sec / 5 = %.2f sec/rep; projected full run ≈ %.1f sec (%.1f min) at N_CORES=%d",
  t_probe[["elapsed"]], sec_per_rep, proj_sec, proj_sec / 60, N_CORES
))
# CAR-alone baseline ~58s / 4 cores. Flag if >> 5× that wall time.
if (proj_sec > 58 * 5) {
  logmsg(sprintf(
    "FLAG: projected runtime %.1f sec well beyond CAR-alone ~58s/4-core baseline. Continuing; report this."
  , proj_sec))
}

# ---------------------------------------------------------------------------
# Full permutation loop
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Permutation loop")
t_start <- proc.time()[[3]]
set.seed(SEED)

fail_reasons <- character(N_BOOT)
if (N_CORES > 1L && .Platform$OS.type == "unix") {
  logmsg(sprintf("Running parallel::mclapply, mc.cores = %d", N_CORES))
  seeds <- SEED + seq_len(N_BOOT)
  null_list <- parallel::mclapply(seq_len(N_BOOT), function(i) {
    set.seed(seeds[[i]])
    permute_and_refit_car_spec_a(dat, form, adjMatrix, phases, nb)
  }, mc.cores = N_CORES)
} else {
  logmsg("Running sequentially")
  null_list <- vector("list", N_BOOT)
  progress_every <- max(1L, floor(N_BOOT / 10L))
  for (i in seq_len(N_BOOT)) {
    null_list[[i]] <- permute_and_refit_car_spec_a(dat, form, adjMatrix, phases, nb)
    if (i == 1L || i %% progress_every == 0L || i == N_BOOT) {
      logmsg(sprintf("  progress %d/%d (%.1f sec)", i, N_BOOT, proc.time()[[3]] - t_start))
    }
  }
}

null_mat <- do.call(rbind, lapply(null_list, function(x) x$coefs))
colnames(null_mat) <- coef_names
failed <- vapply(null_list, function(x) isTRUE(x$failed), logical(1))
fail_reasons <- vapply(null_list, function(x) {
  if (isTRUE(x$failed)) as.character(x$fail_reason) else NA_character_
}, character(1))

runtime_sec <- proc.time()[[3]] - t_start
n_failed <- sum(failed)
n_ok <- N_BOOT - n_failed
logmsg(sprintf(
  "Finished: n_boot=%d; n_failed=%d; n_ok=%d; runtime=%.1f sec (%.2f sec/rep wall)",
  N_BOOT, n_failed, n_ok, runtime_sec, runtime_sec / N_BOOT
))
if (n_failed > 0L) {
  tab_fail <- sort(table(fail_reasons[failed]), decreasing = TRUE)
  logmsg("Failure reasons:")
  for (nm in names(tab_fail)) {
    logmsg(sprintf("  %s: %d", nm, tab_fail[[nm]]))
  }
}

# ---------------------------------------------------------------------------
# Compare to null + CAR-alone
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Observed vs null (CAR+lag)")

summary_rows <- lapply(coef_names, function(nm) {
  obs <- unname(obs_coefs[[nm]])
  null_x <- null_mat[, nm]
  null_clean <- null_x[is.finite(null_x)]
  ns <- null_summary_stats(null_clean)
  p_emp <- empirical_p(null_clean, obs)
  inside <- is.finite(obs) && obs >= ns[["q025"]] && obs <= ns[["q975"]]
  label <- if (nm == "FP_between") {
    sprintf("FP_between (ref %s)", phases[[1]])
  } else {
    paste0("Phase ", gsub("slope_", "", nm), " FP_between slope")
  }
  logmsg(sprintf(
    "%s: obs=%+.6f; null 95%% [%.6f, %.6f]; p=%s; inside=%s",
    label, obs, ns[["q025"]], ns[["q975"]], fmt_p(p_emp), inside
  ))
  data.frame(
    target = nm, label = label, obs_coef = obs,
    null_mean = unname(ns[["mean"]]), null_sd = unname(ns[["sd"]]),
    null_q025 = unname(ns[["q025"]]), null_q975 = unname(ns[["q975"]]),
    p_empirical = p_emp, inside_null_95 = inside,
    stringsAsFactors = FALSE
  )
})
summary_df <- dplyr::bind_rows(summary_rows)

# CAR-alone comparison from archived summary table
car_alone_flag <- c(
  "1985-1991" = "Yes", "1992-2001" = "Yes",
  "2002-2007" = "Yes", "2008-2015" = "No"
)
compare_rows <- lapply(phases, function(ph) {
  sc <- slope_col(ph)
  r <- summary_df[summary_df$target == sc, ]
  data.frame(
    phase = ph,
    obs_car_alone = c(-0.062797, -0.024352, 0.084112, 0.013256)[match(ph, phases)],
    confounded_car_alone = car_alone_flag[[ph]],
    obs_car_spec_a = r$obs_coef,
    null_q025 = r$null_q025,
    null_q975 = r$null_q975,
    p_empirical = r$p_empirical,
    confounded_car_spec_a = ifelse(r$inside_null_95, "No", "Yes"),
    stringsAsFactors = FALSE
  )
})
compare_df <- dplyr::bind_rows(compare_rows)
write_csv(compare_df, path_out_compare)

logmsg("")
logmsg("## Comparison vs CAR-alone confounded phases")
logmsg("| Phase | CAR-alone confounded? | CAR+lag confounded? | CAR+lag obs |")
logmsg("|---|---|---|---|")
for (i in seq_len(nrow(compare_df))) {
  r <- compare_df[i, ]
  logmsg(sprintf(
    "| %s | %s | %s | %+.4f |",
    r$phase, r$confounded_car_alone, r$confounded_car_spec_a, r$obs_car_spec_a
  ))
}

# ---------------------------------------------------------------------------
# Deliverables
# ---------------------------------------------------------------------------
results_df <- data.frame(
  perm_id = seq_len(N_BOOT),
  converged = !failed,
  fail_reason = fail_reasons,
  stringsAsFactors = FALSE
)
for (nm in coef_names) results_df[[nm]] <- as.numeric(null_mat[, nm])
write_csv(results_df, path_out_results)

results_ok <- results_df %>% dplyr::filter(converged)
label_map <- c(
  FP_between = paste0("Main effect\n(ref. ", phases[[1]], ")"),
  stats::setNames(paste0("Phase\n", phases), slope_names)
)
if (nrow(results_ok) > 0L) {
  plot_df <- tidyr::pivot_longer(
    results_ok, cols = tidyselect::all_of(coef_names),
    names_to = "target", values_to = "coef"
  ) %>%
    dplyr::mutate(
      target_label = factor(unname(label_map[target]), levels = unname(label_map[coef_names]))
    )
  obs_plot <- data.frame(
    target = names(obs_coefs), obs = as.numeric(obs_coefs), stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      target_label = factor(unname(label_map[target]), levels = levels(plot_df$target_label))
    )
  p <- ggplot(plot_df, aes(x = coef)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey75", colour = "white") +
    geom_density(colour = "grey25", linewidth = 0.6) +
    geom_vline(data = obs_plot, aes(xintercept = obs), colour = "#B2182B", linewidth = 0.9) +
    facet_wrap(~target_label, scales = "free", nrow = 2) +
    labs(
      title = "FP_between spatial permutation null (CAR + Spec A lag)",
      subtitle = sprintf("Lag recomputed per shuffle; n_boot=%d; n_failed=%d; seed=%d",
                         N_BOOT, n_failed, SEED),
      x = "Coefficient under rectangle-level FP_between permutation",
      y = "Density"
    ) +
    theme_bw(base_size = 11)
  ggplot2::ggsave(path_out_fig, p, width = 11, height = 7, dpi = 150)
}

any_outside <- any(compare_df$confounded_car_spec_a == "Yes")
# Does lag clear the CAR-alone flagged phases?
car_flagged <- phases[car_alone_flag[phases] == "Yes"]
cleared <- all(compare_df$confounded_car_spec_a[compare_df$phase %in% car_flagged] == "No")
headline <- if (cleared) {
  paste0(
    "**Headline (CAR Spec A):** lag term clears all CAR-alone confounded phases ",
    "(1985–1991, 1992–2001, 2002–2007)."
  )
} else if (!any_outside) {
  paste0(
    "**Headline (CAR Spec A):** all phase-specific FP_between slopes inside null 95% ",
    "under CAR+lag (including phases that were not flagged under CAR-alone)."
  )
} else {
  still <- compare_df$phase[compare_df$confounded_car_spec_a == "Yes"]
  paste0(
    "**Headline (CAR Spec A):** lag does **not** clear CAR's confounded phases. ",
    "Still outside null: ", paste(still, collapse = ", "), "."
  )
}

summary_lines <- c(
  "# FP_between spatial confounding bootstrap — CAR Spec A summary",
  "",
  sprintf("Formula: `%s`", paste(deparse(form), collapse = " ")),
  "",
  "Each replicate: shuffle rectangle-level `FP_between`, **recompute ",
  "`FP_between_lag` from the shuffled map** (fixed KNN k=4), refit CAR+lag.",
  "",
  sprintf("- Seed: `%d`; n_boot: %d; n_failed: %d; n_ok: %d; cores: %d; runtime: %.1f sec",
          SEED, N_BOOT, n_failed, n_ok, N_CORES, runtime_sec),
  sprintf("- Probe sec/rep: %.2f; projected vs CAR-alone ~58s/4-core baseline",
          sec_per_rep),
  "",
  "## Comparison vs CAR-alone",
  "",
  "| Phase | CAR-alone obs | CAR-alone confounded? | CAR+lag obs | Null 2.5% | Null 97.5% | CAR+lag confounded? |",
  "|-------|---------------|------------------------|-------------|-----------|------------|---------------------|"
)
for (i in seq_len(nrow(compare_df))) {
  r <- compare_df[i, ]
  summary_lines <- c(summary_lines, sprintf(
    "| %s | %+.4f | %s | %+.4f | %.4f | %.4f | %s |",
    r$phase, r$obs_car_alone, r$confounded_car_alone,
    r$obs_car_spec_a, r$null_q025, r$null_q975, r$confounded_car_spec_a
  ))
}
summary_lines <- c(
  summary_lines, "",
  "## Headline", "", headline, "",
  "## Routing hint (Task 3)", "",
  if (cleared) {
    "Task 3 branch A candidate: promote Spec A (CAR-scoped) as mechanism check."
  } else {
    "Task 3 branch B candidate: archive Spec A — does not explain CAR's confounded phases."
  },
  "",
  "## Outputs",
  "",
  sprintf("- `%s`", path_out_results),
  sprintf("- `%s`", path_out_summary),
  sprintf("- `%s`", path_out_compare),
  ""
)
writeLines(summary_lines, path_out_summary)

saveRDS(
  list(
    seed = SEED, n_boot = N_BOOT, n_failed = n_failed, n_ok = n_ok,
    n_cores = N_CORES, runtime_sec = runtime_sec, sec_per_rep_probe = sec_per_rep,
    formula = form, obs_coefs = obs_coefs, null_mat = null_mat,
    summary = summary_df, compare = compare_df, failed = failed,
    fail_reasons = fail_reasons, headline = headline, cleared_car_phases = cleared
  ),
  path_out_rds
)

logmsg("")
logmsg(headline)
writeLines(run_log, path_out_run_log)
cat("Run log:", path_out_run_log, "\n")
