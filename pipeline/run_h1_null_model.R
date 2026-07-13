# H1 null model: randomise B_obs while holding B_pred (from observed S, N, E) fixed.
# Primary null scheme is chosen from log(B_obs) distribution shape.

# H1 B-only null model: randomise B_obs while B_pred (from S, N, E) stays fixed.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h1_null_model.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "h1_null_helpers.R"))

path_preds <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_decision <- file.path(project_root, "outputs", "null_sampling_decision.rds")
path_dist_summary <- file.path(project_root, "outputs", "null_distribution_summary.csv")
path_perms <- file.path(project_root, "outputs", "null_permutations.rds")
path_summary <- file.path(project_root, "outputs", "null_summary.rds")
path_fig_dir <- file.path(project_root, "outputs", "figures")

dir.create(path_fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(path_preds))

N_PERM <- 999L
SEED <- 42L

haul <- readRDS(path_preds)
B_obs <- haul$B_obs
B_pred <- haul$B_pred
log_obs <- log(B_obs)
log_pred <- log(B_pred)

T_obs <- log_r2(log_obs, log_pred)
median_abs_res_obs <- median(abs(log_obs - log_pred), na.rm = TRUE)

decision <- assess_b_null_sampling(B_obs)
write.csv(decision$metrics, path_dist_summary, row.names = FALSE)

decision_record <- list(
  assessment = decision,
  selected_at = Sys.time(),
  n_hauls = nrow(haul),
  T_obs = T_obs,
  median_abs_residual_obs = median_abs_res_obs
)
saveRDS(decision_record, path_decision)

cat("=== B null sampling decision ===\n")
cat(decision$rationale, "\n")
cat("Primary method:", decision$primary_method, "\n")
cat("Robustness method:", decision$robustness_method, "\n\n")

run_one_method <- function(method, label) {
  T_null <- run_null_simulation(
    B_obs_g = B_obs,
    B_pred_g = B_pred,
    method = method,
    n_perm = N_PERM,
    log_bounds_95 = decision$log_bounds_95,
    seed = SEED
  )

  abs_res_null <- vapply(
    seq_len(N_PERM),
    function(i) {
      if (method == "b_shuffle") {
        B_null <- sample(B_obs, replace = FALSE)
      } else {
        B_null <- exp(runif(
          length(B_obs),
          min = decision$log_bounds_95[1L],
          max = decision$log_bounds_95[2L]
        ))
      }
      median(abs(log(B_null) - log_pred), na.rm = TRUE)
    },
    numeric(1)
  )

  sum_tab <- summarise_null_test(T_obs, T_null, label)
  sum_tab$median_abs_residual_obs <- median_abs_res_obs
  sum_tab$median_abs_residual_null <- median(abs_res_null)

  fig_path <- file.path(
    path_fig_dir,
    paste0("null_r2_", gsub("_", "-", method), ".png")
  )
  plot_null_r2_distribution(T_obs, T_null, label, fig_path)

  list(
    method = method,
    label = label,
    T_null = T_null,
    summary = sum_tab,
    figure = fig_path
  )
}

primary_res <- run_one_method(
  decision$primary_method,
  paste0("primary: ", decision$primary_method)
)
robust_res <- run_one_method(
  decision$robustness_method,
  paste0("robustness: ", decision$robustness_method)
)

perm_record <- list(
  n_perm = N_PERM,
  seed = SEED,
  T_obs = T_obs,
  primary = primary_res,
  robustness = robust_res,
  decision = decision
)
saveRDS(perm_record, path_perms)

summary_record <- list(
  decision = decision,
  T_obs = T_obs,
  median_abs_residual_obs = median_abs_res_obs,
  primary = primary_res$summary,
  robustness = robust_res$summary,
  figures = c(primary_res$figure, robust_res$figure)
)
saveRDS(summary_record, path_summary)

cat("=== Observed H1 performance ===\n")
cat("R² (log scale, B_obs vs B_pred):", round(T_obs, 4), "\n")
cat("Median |residual|:", round(median_abs_res_obs, 4), "\n\n")

print_summary <- function(s, title) {
  cat("=== ", title, " ===\n", sep = "")
  cat("Null R² median:", round(s$null_median, 4), "\n")
  cat("Null R² 2.5%–97.5%:", round(s$null_q025, 4), "–", round(s$null_q975, 4), "\n")
  cat("p-value (one-sided, P(null R² >= observed)):", round(s$p_value_one_sided, 4), "\n")
  cat("p-value (two-sided):", round(s$p_value_two_sided, 4), "\n")
  cat("Median |residual| under null:", round(s$median_abs_residual_null, 4), "\n\n")
}

print_summary(primary_res$summary, "Primary null")
print_summary(robust_res$summary, "Robustness null")

cat("Saved:\n")
cat(" ", path_decision, "\n")
cat(" ", path_dist_summary, "\n")
cat(" ", path_perms, "\n")
cat(" ", path_summary, "\n")
cat(" ", primary_res$figure, "\n")
cat(" ", robust_res$figure, "\n")
