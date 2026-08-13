# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Spatial confounding bootstrap on FP_between — Spec A (HISTORICAL — RE).
#
# SUPERSEDED as the presented Spec A claim. Current Spec A bootstrap is:
#   pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R
# See outputs/spec_a_car_routing_note.md.
#
# Historical formula (same protocol as archived RE bootstrap):
#   residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2
#            + FP_within * phase_v2 + (1 | stat_rec)  [REML]
#
# Primary CAR confounding check (no lag):
#   pipeline/permutation_bootstrap_FP_between_CAR.R
#
# After each FP_between shuffle, FP_between_lag is recomputed from the shuffled
# map using fixed k-NN neighbour structure (k=4).
#
# Run: Rscript --vanilla pipeline/run_h2h3_fp_between_confounding_bootstrap_spec_a.R
#
# Optional env: N_BOOT=1000  SEED=42  N_CORES=1

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
  stop("Run from pipeline/ or Rscript pipeline/run_h2h3_fp_between_confounding_bootstrap_spec_a.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2_common.R"))
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))
source(file.path(script_dir, "R", "h2h3_knn_spatial_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: ",
    "Rscript --vanilla pipeline/run_h2h3_fp_between_confounding_bootstrap_spec_a.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

path_spec_a <- file.path(project_root, "outputs", "primary_model_v2_spec_a.rds")
path_knn <- file.path(project_root, "outputs", "knn_listw_k4.rds")
path_orig_summary <- file.path(
  project_root, "exploratory", "outputs",
  "permutation_bootstrap_FP_between_RE_summary.md"
)

fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_results <- file.path(
  project_root, "outputs", "fp_between_confounding_bootstrap_spec_a_results.csv"
)
path_out_summary <- file.path(
  project_root, "outputs", "fp_between_confounding_bootstrap_spec_a_summary.md"
)
path_out_fig <- file.path(
  fig_dir, "fp_between_confounding_bootstrap_spec_a_null.png"
)
path_out_run_log <- file.path(
  project_root, "outputs", "fp_between_confounding_bootstrap_spec_a_run_log.md"
)
path_out_session <- file.path(
  project_root, "outputs", "fp_between_confounding_bootstrap_spec_a_sessionInfo.txt"
)
path_out_rds <- file.path(
  project_root, "outputs", "fp_between_confounding_bootstrap_spec_a_objects.rds"
)

N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
N_CORES <- as.integer(Sys.getenv("N_CORES", unset = "1"))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

PHASE_V2_LEVELS <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")

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

extract_fp_between_coefs <- function(fit, phases) {
  out <- na_coef_template(phases)
  b <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
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

fit_primary_on_data <- function(data, model_formula) {
  tryCatch(
    glmmTMB::glmmTMB(
      model_formula,
      data = data,
      family = stats::gaussian(),
      REML = TRUE
    ),
    error = function(e) NULL
  )
}

permute_and_refit_spec_a <- function(data, model_formula, phases, nb_knn) {
  data_perm <- permute_fp_between_data(data)
  data_perm <- recompute_fp_between_lag_on_data(data_perm, nb_knn)
  fit <- fit_primary_on_data(data_perm, model_formula)
  if (is.null(fit)) return(na_coef_template(phases))
  extract_fp_between_coefs(fit, phases)
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
      "No evidence the rectangle intercept and FP_between are substituting for ",
      "one another. Reported effect is not an artifact of spatial confounding."
    )
  } else {
    paste0(
      "The rectangle intercept and FP_between are competing for the same spatial ",
      "signal when FP_between's real spatial arrangement is present. The reported ",
      "H2 effect is at least partly attributable to this overlap, not purely to ",
      "the fishing-pressure effect itself."
    )
  }
  list(inside_95 = inside, direction = direction, statement = statement)
}

load_spec_a_model <- function(path) {
  obj <- readRDS(path)
  if (is.null(obj$primary_model_spec_a) || is.null(obj$data) || is.null(obj$formula_spec_a)) {
    stop("primary_model_v2_spec_a.rds missing primary_model_spec_a, data, or formula_spec_a.")
  }
  phases <- if (!is.null(obj$phase_v2_labels)) {
    as.character(obj$phase_v2_labels)
  } else {
    PHASE_V2_LEVELS
  }
  list(
    source = path,
    scheme = "phase_v2_spec_a",
    fit = obj$primary_model_spec_a,
    data = obj$data,
    formula = obj$formula_spec_a,
    phase_col = "phase_v2",
    phases = phases,
    model_id = "wb_primary_phase_v2_spec_a"
  )
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# FP_between spatial confounding bootstrap — Spec A — run log")
logmsg("")
logmsg(
  "Rectangle-level permutation of FP_between under the Spec A model. ",
  "FP_between_lag recomputed from shuffled FP_between after each permutation; ",
  "k-NN neighbour structure fixed. Targets: same FP_between × phase_v2 ",
  "coefficients as the original bootstrap for comparison."
)

logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo: ", path_out_session)
logmsg(sprintf("glmmTMB %s", as.character(utils::packageVersion("glmmTMB"))))
logmsg(sprintf("N_BOOT = %d; SEED = %d; N_CORES = %d", N_BOOT, SEED, N_CORES))

stopifnot(file.exists(path_spec_a), file.exists(path_knn))

primary <- load_spec_a_model(path_spec_a)
knn_pack <- readRDS(path_knn)
nb_knn <- knn_pack$nb

fit_obs <- primary$fit
dat <- primary$data
form <- primary$formula
phase_col <- primary$phase_col
phases <- unname(as.character(primary$phases))
slope_names <- unname(vapply(phases, slope_col, character(1)))
coef_names <- c("FP_between", slope_names)

dat$stat_rec <- factor(dat$stat_rec)
dat[[phase_col]] <- factor(dat[[phase_col]], levels = phases)

n_hauls <- nrow(dat)
n_rect <- length(unique(dat$stat_rec))

logmsg("")
logmsg("## Inputs")
logmsg("Loaded Spec A: ", path_spec_a)
logmsg("k-NN weights: ", path_knn)
logmsg(sprintf("Phase levels: %s", paste(phases, collapse = " | ")))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  n_hauls, n_rect, min(dat$year), max(dat$year)
))
logmsg("Formula: ", paste(deparse(form), collapse = " "))

# Sanity check
set.seed(SEED)
dat_check <- permute_fp_between_data(dat)
dat_check <- recompute_fp_between_lag_on_data(dat_check, nb_knn)
stopifnot(identical(dat_check$FP_within, dat$FP_within))
stopifnot(identical(dat_check$residual, dat$residual))
stopifnot(!isTRUE(all.equal(dat_check$FP_between, dat$FP_between)))
stopifnot(!isTRUE(all.equal(dat_check$FP_between_lag, dat$FP_between_lag)))
logmsg(
  "Permutation sanity check passed: FP_between and FP_between_lag change; ",
  "FP_within / residual / stat_rec / phase_v2 unchanged."
)

# Observed coefficients
logmsg("")
logmsg("## Step 1 — Observed H2 coefficients (Spec A baseline)")
obs_coefs <- extract_fp_between_coefs(fit_obs, phases)
obs_slopes_tab <- extract_wb_phase_slopes(
  fit_obs,
  term_name = "FP_between",
  model_id = primary$model_id,
  hypothesis_group = "H2_spatial_between",
  phases = phases
)
obs_fe <- obs_slopes_tab[1, ]

# Permutation loop
logmsg("")
logmsg("## Steps 2–3 — Permutation + refit (lag recomputed each replicate)")
t_start <- proc.time()[[3]]
set.seed(SEED)

if (N_CORES > 1L && .Platform$OS.type == "unix") {
  logmsg(sprintf("Running with parallel::mclapply, mc.cores = %d", N_CORES))
  seeds <- SEED + seq_len(N_BOOT)
  null_list <- parallel::mclapply(seq_len(N_BOOT), function(i) {
    set.seed(seeds[[i]])
    permute_and_refit_spec_a(dat, form, phases, nb_knn)
  }, mc.cores = N_CORES)
  null_mat <- do.call(rbind, null_list)
} else {
  logmsg("Running sequentially")
  null_mat <- matrix(NA_real_, nrow = N_BOOT, ncol = length(coef_names))
  colnames(null_mat) <- coef_names
  progress_every <- max(1L, floor(N_BOOT / 10L))
  for (i in seq_len(N_BOOT)) {
    null_mat[i, ] <- permute_and_refit_spec_a(dat, form, phases, nb_knn)
    if (i == 1L || i %% progress_every == 0L || i == N_BOOT) {
      elapsed_i <- proc.time()[[3]] - t_start
      logmsg(sprintf(
        "  progress %d/%d (%.1f sec elapsed, ~%.2f sec/rep)",
        i, N_BOOT, elapsed_i, elapsed_i / i
      ))
    }
  }
}
if (is.null(dim(null_mat))) {
  null_mat <- matrix(null_mat, nrow = N_BOOT, byrow = TRUE)
}
colnames(null_mat) <- coef_names

t_end <- proc.time()[[3]]
runtime_sec <- t_end - t_start
n_failed <- sum(!is.finite(null_mat[, "FP_between"]))
n_ok <- N_BOOT - n_failed
logmsg(sprintf(
  "Finished: n_boot = %d; n_failed = %d; n_ok = %d; runtime = %.1f sec",
  N_BOOT, n_failed, n_ok, runtime_sec
))

# Compare to null
logmsg("")
logmsg("## Step 4 — Observed vs permutation null (Spec A)")

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
    "%s: obs=%+.6f; null 95%% [%.6f, %.6f]; p_emp=%s; inside_95=%s",
    targets[[nm]], obs, ns[["q025"]], ns[["q975"]], fmt_p(p_emp), inter$inside_95
  ))
  data.frame(
    target = nm,
    label = targets[[nm]],
    obs_coef = obs,
    null_mean = unname(ns[["mean"]]),
    null_sd = unname(ns[["sd"]]),
    null_q025 = unname(ns[["q025"]]),
    null_q975 = unname(ns[["q975"]]),
    p_empirical = p_emp,
    inside_null_95 = inter$inside_95,
    direction = inter$direction,
    interpretation = inter$statement,
    stringsAsFactors = FALSE
  )
})
summary_df <- dplyr::bind_rows(summary_rows)
any_outside <- any(summary_df$target %in% slope_names & !summary_df$inside_null_95)

# Original bootstrap comparison (phase slopes from archived summary table)
if (file.exists(path_orig_summary)) {
  orig_txt <- readLines(path_orig_summary, warn = FALSE)
  table_rows <- grep("^\\| Phase .* FP_between slope \\|", orig_txt, value = TRUE)
  logmsg("")
  logmsg("## Comparison vs original primary bootstrap")
  for (ph in phases) {
    sc <- slope_col(ph)
    spec_inside <- summary_df$inside_null_95[summary_df$target == sc]
    row <- grep(ph, table_rows, fixed = TRUE, value = TRUE)
    orig_flag <- if (length(row) == 1L) {
      # last non-empty cell is yes/no for inside_null_95
      cells <- trimws(strsplit(row, "\\|", fixed = FALSE)[[1]])
      cells <- cells[nzchar(cells)]
      cells[[length(cells)]]
    } else {
      "NA"
    }
    logmsg(sprintf(
      "Phase %s: original inside_95=%s; Spec A inside_95=%s",
      ph,
      orig_flag,
      if (length(spec_inside) && isTRUE(spec_inside)) "yes" else "no"
    ))
  }
}

# Deliverables
results_df <- data.frame(
  perm_id = seq_len(N_BOOT),
  converged = is.finite(null_mat[, "FP_between"]),
  stringsAsFactors = FALSE
)
for (nm in coef_names) {
  results_df[[nm]] <- as.numeric(null_mat[, nm])
}
results_ok <- results_df %>% dplyr::filter(converged)

readr::write_csv(
  dplyr::bind_rows(
    data.frame(record_type = "metadata", perm_id = NA_integer_, n_boot = N_BOOT,
               seed = SEED, runtime_sec = runtime_sec, stringsAsFactors = FALSE),
    results_ok %>% dplyr::mutate(record_type = "permutation")
  ),
  path_out_results
)

label_map <- c(
  FP_between = paste0("Main effect\n(ref. ", phases[[1]], ")"),
  stats::setNames(paste0("Phase\n", phases), slope_names)
)
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
    colour = "#B2182B", linewidth = 0.9
  ) +
  facet_wrap(~target_label, scales = "free", nrow = 2) +
  labs(
    title = "FP_between spatial permutation null (Spec A: + FP_between_lag * phase_v2)",
    subtitle = sprintf(
      "Lag recomputed per shuffle; n_boot = %d; seed = %d",
      N_BOOT, SEED
    ),
    x = "Coefficient under rectangle-level FP_between permutation",
    y = "Density"
  ) +
  theme_bw(base_size = 11)
ggplot2::ggsave(path_out_fig, p, width = 11, height = 7, dpi = 150)

form_txt <- paste(deparse(form), collapse = " ")
headline <- if (any_outside) {
  paste0(
    "**Headline (Spec A):** at least one phase-specific `FP_between` slope falls ",
    "outside its spatial-permutation null 95% interval under the lag-augmented model."
  )
} else {
  paste0(
    "**Headline (Spec A):** all phase-specific `FP_between` slopes fall inside ",
    "their spatial-permutation null 95% intervals under the lag-augmented model."
  )
}

summary_lines <- c(
  "# FP_between spatial confounding bootstrap — Spec A summary",
  "",
  sprintf("Formula: `%s` (REML)", form_txt),
  "",
  "Each replicate: shuffle rectangle-level `FP_between` once, **recompute ",
  "`FP_between_lag` from the shuffled map** (fixed k-NN structure), refit Spec A.",
  "",
  sprintf("- Seed: `%d`; n_boot: %d; n_failed: %d; runtime: %.1f sec", SEED, N_BOOT, n_failed, runtime_sec),
  "",
  "## Null distribution vs observed",
  "",
  "| Target | Observed | Null 2.5% | Null 97.5% | Empirical p | Inside null 95%? |",
  "|--------|----------|-----------|------------|-------------|------------------|"
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.6f | %.6f | %.6f | %s | %s |",
      r$label, r$obs_coef, r$null_q025, r$null_q975,
      fmt_p(r$p_empirical), ifelse(r$inside_null_95, "yes", "no")
    )
  )
}
summary_lines <- c(summary_lines, "", "## Headline", "", headline, "")
writeLines(summary_lines, path_out_summary)

saveRDS(
  list(
    seed = SEED, n_boot = N_BOOT, n_failed = n_failed, runtime_sec = runtime_sec,
    formula = form, obs_coefs = obs_coefs, null_mat = null_mat, summary = summary_df
  ),
  path_out_rds
)

logmsg("")
logmsg("## Deliverables")
logmsg("Wrote: ", path_out_results)
logmsg("Wrote: ", path_out_summary)
logmsg("Wrote: ", path_out_fig)
logmsg("Wrote: ", path_out_rds)
logmsg("")
logmsg(headline)

writeLines(run_log, path_out_run_log)
cat("Run log:", path_out_run_log, "\n")
