# ARCHIVED (exploratory): this was the original spatial-confounding permutation
# test, run against the RE (1 | stat_rec) specification before CAR was adopted
# as the primary H2 model; see pipeline/permutation_bootstrap_FP_between_CAR.R
# for the current version.
#
# Spatial confounding bootstrap on FP_between (H2) — RE specification
#
# PURPOSE: test whether the reported H2 effect depends on the spatial arrangement
# of rectangle-level FP_between, relative to the shared (1 | stat_rec) random
# intercept. Shuffles FP_between once per rectangle (preserving the set of
# values), leaves residual / FP_within / phase factor / rectangle IDs untouched,
# and refits the within-between glmmTMB model.
#
# Model (policy-anchored phase_v2; superseded for H2 spatial contrasts by CAR):
#   residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)  [REML]
#   phases: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015
#
# Prefers outputs/primary_model_v2.rds; falls back to h2h3_wb_model_objects.rds
# (legacy data-driven phase) if v2 is absent.
#
# Run: Rscript --vanilla exploratory/pipeline/permutation_bootstrap_FP_between_RE.R
# Writes under exploratory/outputs/.
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
  stop("Run: Rscript --vanilla exploratory/pipeline/permutation_bootstrap_FP_between_RE.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2h3_within_between_helpers.R"))

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop(
    "Package 'glmmTMB' required. Run with: ",
    "Rscript --vanilla exploratory/pipeline/permutation_bootstrap_FP_between_RE.R"
  )
}
suppressPackageStartupMessages(library(glmmTMB))

# ---------------------------------------------------------------------------
# Paths / settings (archived under exploratory/; live primary model stays top-level)
# ---------------------------------------------------------------------------
path_models_v2 <- file.path(project_root, "outputs", "primary_model_v2.rds")
path_models_legacy <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")

out_root <- file.path(project_root, "exploratory", "outputs")
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_results <- file.path(
  out_root, "permutation_bootstrap_FP_between_RE_results.csv"
)
path_out_summary <- file.path(
  out_root, "permutation_bootstrap_FP_between_RE_summary.md"
)
path_out_fig <- file.path(
  fig_dir, "permutation_bootstrap_FP_between_RE_null.png"
)
path_out_run_log <- file.path(
  out_root, "permutation_bootstrap_FP_between_RE_run_log.md"
)
path_out_session <- file.path(
  out_root, "permutation_bootstrap_FP_between_RE_sessionInfo.txt"
)
path_out_rds <- file.path(
  out_root, "permutation_bootstrap_FP_between_RE_objects.rds"
)

N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
# Default sequential: TMB/glmmTMB is unsafe under forked mclapply.
n_cores_env <- Sys.getenv("N_CORES", unset = "1")
N_CORES <- as.integer(n_cores_env)
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- 1L

PHASE_V2_LEVELS <- c("1985-1991", "1992-2001", "2002-2007", "2008-2015")
PHASE_LEGACY_LEVELS <- c("1985-1988", "1989-2000", "2001-2007", "2008-2015")

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

permute_and_refit <- function(data, model_formula, phases,
                              rectangle_col = "stat_rec",
                              fp_col = "FP_between") {
  data_perm <- permute_fp_between_data(data, rectangle_col, fp_col)
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

load_primary_model <- function(path_v2, path_legacy) {
  if (file.exists(path_v2)) {
    obj <- readRDS(path_v2)
    if (is.null(obj$primary_model_v2) || is.null(obj$data) || is.null(obj$formula_v2)) {
      stop("primary_model_v2.rds missing primary_model_v2, data, or formula_v2.")
    }
    phases <- if (!is.null(obj$phase_v2_labels)) {
      as.character(obj$phase_v2_labels)
    } else {
      PHASE_V2_LEVELS
    }
    list(
      source = path_v2,
      scheme = "phase_v2",
      fit = obj$primary_model_v2,
      data = obj$data,
      formula = obj$formula_v2,
      phase_col = "phase_v2",
      phases = phases,
      model_id = "wb_primary_phase_v2"
    )
  } else if (file.exists(path_legacy)) {
    obj <- readRDS(path_legacy)
    if (is.null(obj$fit_wb) || is.null(obj$data) || is.null(obj$formula_wb)) {
      stop("h2h3_wb_model_objects.rds missing fit_wb, data, or formula_wb.")
    }
    list(
      source = path_legacy,
      scheme = "phase_legacy",
      fit = obj$fit_wb,
      data = obj$data,
      formula = obj$formula_wb,
      phase_col = "phase",
      phases = PHASE_LEGACY_LEVELS,
      model_id = "wb_primary"
    )
  } else {
    stop("No primary model found. Expected ", path_v2, " or ", path_legacy)
  }
}

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
logmsg("# FP_between spatial confounding bootstrap — run log")
logmsg("")
logmsg(
  "Rectangle-level permutation of FP_between under the current primary ",
  "within-between glmmTMB model. Destroys the spatial arrangement of the H2 ",
  "covariate while leaving residual, FP_within, the phase factor, and ",
  "(1 | stat_rec) untouched."
)
logmsg("")
logmsg("## Session")
sink(path_out_session)
print(utils::sessionInfo())
sink()
logmsg("sessionInfo written to: ", path_out_session)
logmsg(sprintf("glmmTMB %s", as.character(utils::packageVersion("glmmTMB"))))
logmsg(sprintf("N_BOOT = %d; SEED = %d; N_CORES = %d", N_BOOT, SEED, N_CORES))

# ---------------------------------------------------------------------------
# Load saved primary fit + analysis data
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Inputs")
primary <- load_primary_model(path_models_v2, path_models_legacy)
fit_obs <- primary$fit
dat <- primary$data
form <- primary$formula
phase_col <- primary$phase_col
phases <- unname(as.character(primary$phases))
slope_names <- unname(vapply(phases, slope_col, character(1)))
coef_names <- c("FP_between", slope_names)

dat$stat_rec <- factor(dat$stat_rec)
if (!phase_col %in% names(dat)) {
  stop("Analysis data missing phase column: ", phase_col)
}
dat[[phase_col]] <- factor(dat[[phase_col]], levels = phases)

n_hauls <- nrow(dat)
n_rect <- length(unique(dat$stat_rec))
fp_map <- unique(dat[, c("stat_rec", "FP_between")])
stopifnot(nrow(fp_map) == n_rect)

logmsg("Loaded: ", primary$source)
logmsg(sprintf("Phase scheme: %s (column `%s`)", primary$scheme, phase_col))
logmsg(sprintf("Phase levels: %s", paste(phases, collapse = " | ")))
logmsg(sprintf(
  "Analysis data: %d hauls, %d rectangles, years %d–%d",
  n_hauls, n_rect, min(dat$year), max(dat$year)
))
logmsg("Formula: ", paste(deparse(form), collapse = " "))
logmsg(sprintf(
  "FP_between unique per rectangle: %d values; range [%.3f, %.3f]",
  nrow(fp_map), min(fp_map$FP_between), max(fp_map$FP_between)
))

# Sanity: permuting FP_between must not alter FP_within / phase / residual
set.seed(SEED)
dat_check <- permute_fp_between_data(dat)
stopifnot(identical(dat_check$FP_within, dat$FP_within))
stopifnot(identical(dat_check$residual, dat$residual))
stopifnot(identical(as.character(dat_check$stat_rec), as.character(dat$stat_rec)))
stopifnot(identical(as.character(dat_check[[phase_col]]), as.character(dat[[phase_col]])))
stopifnot(!isTRUE(all.equal(dat_check$FP_between, dat$FP_between)))
logmsg(
  "Permutation sanity check passed: only FP_between changes; ",
  "FP_within / residual / stat_rec / ", phase_col, " unchanged."
)

# ---------------------------------------------------------------------------
# Step 1 — observed H2 coefficients
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 1 — Observed H2 coefficients (baseline)")
obs_coefs <- extract_fp_between_coefs(fit_obs, phases)
obs_main <- unname(obs_coefs[["FP_between"]])

obs_slopes_tab <- extract_wb_phase_slopes(
  fit_obs,
  term_name = "FP_between",
  model_id = primary$model_id,
  hypothesis_group = "H2_spatial_between",
  phases = phases
)

# Reference-phase SE/CI from slope table row 1
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
logmsg("## Steps 2–3 — Rectangle-level FP_between permutation + refit")
logmsg(
  "Each replicate: sample(FP_between) across the 158 rectangles once, reassign ",
  "to all hauls in that rectangle, refit primary model with REML = TRUE."
)
logmsg("")
logmsg("### Matched permutation procedure (explicit)")
logmsg(
  "For each replicate: (1) one global shuffle of the rectangle-level FP_between ",
  "values; (2) one refit of the full primary formula; (3) one fixef() extract ",
  "from that single fit, recording the reference-phase FP_between main effect ",
  "and the FP_between × phase interaction increments (combined to phase slopes ",
  "as main + interaction). Not four independent shuffle-and-refit procedures. ",
  "The resulting null distributions are columns of one n_boot × K matrix — ",
  "matched randomizations — so the empirical p-values come from a coherent ",
  "matched permutation procedure."
)

t_start <- proc.time()[[3]]
set.seed(SEED)

if (N_CORES > 1L && .Platform$OS.type == "unix") {
  logmsg(sprintf("Running with parallel::mclapply, mc.cores = %d", N_CORES))
  seeds <- SEED + seq_len(N_BOOT)
  null_list <- parallel::mclapply(seq_len(N_BOOT), function(i) {
    set.seed(seeds[[i]])
    permute_and_refit(dat, form, phases)
  }, mc.cores = N_CORES)
  null_mat <- do.call(rbind, null_list)
} else {
  logmsg("Running sequentially (replicate)")
  null_mat <- t(replicate(
    N_BOOT,
    permute_and_refit(dat, form, phases),
    simplify = TRUE
  ))
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
  "Finished: n_boot = %d; n_failed = %d; n_ok = %d; runtime = %.1f sec (%.2f sec/rep)",
  N_BOOT, n_failed, n_ok, runtime_sec, runtime_sec / N_BOOT
))
if (n_failed > 0L) {
  logmsg(
    "NOTE: failed permutations recorded as NA (not silently dropped from n_boot). ",
    "A high failure rate can indicate estimation problems under spatially shuffled FP_between."
  )
}
if (runtime_sec > 1800) {
  logmsg(
    "FLAG: runtime > 30 min. Consider reducing N_BOOT or increasing N_CORES ",
    "(env vars N_BOOT / N_CORES)."
  )
}

# ---------------------------------------------------------------------------
# Step 4 — compare observed to null
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Step 4 — Observed vs permutation null")

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
    "%s: obs=%+.6f; null mean=%+.6f sd=%.6f; 95%% null [%.6f, %.6f]; p_emp=%s; inside_95=%s",
    targets[[nm]], obs, ns[["mean"]], ns[["sd"]], ns[["q025"]], ns[["q975"]],
    fmt_p(p_emp), inter$inside_95
  ))
  logmsg(sprintf("  Direction: %s", inter$direction))
  logmsg(sprintf("  Interpretation: %s", inter$statement))
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

# ---------------------------------------------------------------------------
# Deliverables
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Deliverables")

results_df <- data.frame(
  perm_id = seq_len(N_BOOT),
  converged = is.finite(null_mat[, "FP_between"]),
  stringsAsFactors = FALSE
)
for (nm in coef_names) {
  results_df[[nm]] <- as.numeric(null_mat[, nm])
}
results_ok <- results_df %>% dplyr::filter(converged)

meta_row <- data.frame(
  record_type = "metadata",
  perm_id = NA_integer_,
  converged = NA,
  phase_scheme = primary$scheme,
  phase_col = phase_col,
  phases = paste(phases, collapse = "|"),
  n_boot = N_BOOT,
  n_failed = n_failed,
  n_ok = n_ok,
  seed = SEED,
  runtime_sec = runtime_sec,
  n_cores = N_CORES,
  stringsAsFactors = FALSE
)
# Align coefficient columns as NA in metadata row
for (nm in coef_names) meta_row[[nm]] <- NA_real_

results_out <- dplyr::bind_rows(
  meta_row,
  results_ok %>%
    dplyr::mutate(
      record_type = "permutation",
      phase_scheme = NA_character_,
      phase_col = NA_character_,
      phases = NA_character_,
      n_boot = NA_integer_,
      n_failed = NA_integer_,
      n_ok = NA_integer_,
      seed = NA_integer_,
      runtime_sec = NA_real_,
      n_cores = NA_integer_
    )
)
readr::write_csv(results_out, path_out_results)
logmsg("Wrote: ", path_out_results)

# Plot labels
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
    colour = "#B2182B", linewidth = 0.9, linetype = "solid"
  ) +
  facet_wrap(~target_label, scales = "free", nrow = 2) +
  labs(
    title = sprintf(
      "FP_between spatial permutation null (%s primary model)",
      primary$scheme
    ),
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

form_txt <- paste(deparse(form), collapse = " ")
summary_lines <- c(
  "# FP_between spatial confounding bootstrap — summary",
  "",
  "Model-appropriate spatial confounding check for H2: randomly reassigns which",
  "rectangle receives which `FP_between` value, then refits the current primary model.",
  "",
  sprintf("- Formula: `%s` (REML)", form_txt),
  sprintf("- Phase scheme: **%s** (`%s`)", primary$scheme, phase_col),
  sprintf("- Phase levels: %s", paste(phases, collapse = " / ")),
  "",
  "## Settings",
  "",
  sprintf("- Seed: `%d`", SEED),
  sprintf("- `n_boot`: %d", N_BOOT),
  sprintf("- `n_failed` (non-converged / NA): %d", n_failed),
  sprintf("- `n_ok`: %d", n_ok),
  sprintf("- Cores: %d", N_CORES),
  sprintf("- Runtime: %.1f sec (%.2f sec/replicate)", runtime_sec, runtime_sec / N_BOOT),
  sprintf("- Rectangles: %d; hauls: %d", n_rect, n_hauls),
  sprintf("- Source fit: `%s`", primary$source),
  "",
  "## Observed H2 coefficients",
  "",
  sprintf(
    "- Reference-phase main effect `FP_between` (%s): **%+.6f** (SE %.6f; 95%% CI [%.6f, %.6f]; p = %.4g)",
    phases[[1]], obs_fe$fp_slope, obs_fe$fp_slope_se, obs_fe$fp_slope_lo,
    obs_fe$fp_slope_hi, obs_fe$p_value
  ),
  "",
  "Phase-specific `FP_between` slopes (presented H2 effects):",
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
  "## Null distribution vs observed",
  "",
  "| Target | Observed | Null mean | Null SD | Null 2.5% | Null 97.5% | Empirical p | Inside null 95%? |",
  "|--------|----------|-----------|---------|-----------|------------|-------------|------------------|"
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  summary_lines <- c(
    summary_lines,
    sprintf(
      "| %s | %+.6f | %+.6f | %.6f | %.6f | %.6f | %s | %s |",
      r$label, r$obs_coef, r$null_mean, r$null_sd, r$null_q025, r$null_q975,
      fmt_p(r$p_empirical), ifelse(r$inside_null_95, "yes", "no")
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "## Interpretation (Step 5)",
  "",
  "Rule:",
  "",
  "- Observed coefficient **inside** the permutation null 95% interval → H2 effect does",
  "  not depend on the spatial arrangement of `FP_between`; no evidence that the",
  "  rectangle intercept and `FP_between` are substituting for one another.",
  "- Observed coefficient **outside** the null 95% interval → rectangle intercept and",
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
    "**Headline (phase-specific H2 slopes, ", primary$scheme, "):** at least one ",
    "presented phase-specific `FP_between` slope falls outside its spatial-permutation ",
    "null 95% interval, indicating that the reported H2 effect depends on the real ",
    "spatial arrangement of fishing pressure and is at least partly entangled with ",
    "the rectangle intercept's spatial signal."
  )
} else {
  paste0(
    "**Headline (phase-specific H2 slopes, ", primary$scheme, "):** all presented ",
    "phase-specific `FP_between` slopes fall inside their spatial-permutation null ",
    "95% intervals. No evidence that the reported H2 effects are an artifact of ",
    "spatial confounding between `FP_between` and `(1 | stat_rec)`."
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
    source = primary$source,
    scheme = primary$scheme,
    phase_col = phase_col,
    phases = phases,
    formula = form,
    obs_coefs = obs_coefs,
    obs_slopes = obs_slopes_tab,
    null_mat = null_mat,
    summary = summary_df,
    results_ok = results_ok
  ),
  path_out_rds
)
logmsg("Wrote: ", path_out_rds)

logmsg("")
logmsg("## Headline")
logmsg(headline)

writeLines(run_log, path_out_run_log)
cat("Run log written to: ", path_out_run_log, "\n", sep = "")
