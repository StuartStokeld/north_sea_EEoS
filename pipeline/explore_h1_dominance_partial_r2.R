# H1 dominance / size-homogeneity: partial-R^2 follow-up
# (follow-up to cursor_briefing_dominance_confound.md)
#
# Quick statistical follow-up — does not touch the primary H1 pipeline. Uses
# the already-built per-haul dominance table (outputs/h1_dominance_haul_table.csv);
# no new data or recomputation of N/E/B_pred.
#
# Question: how much residual variance do D and size_CV jointly explain,
# beyond log(B_obs) alone? Replaces "two separate binned gradients" with one
# defensible partial-R^2 statistic.

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/explore_h1_dominance_partial_r2.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
project_root <- get_project_root_from(script_dir)

path_haul_tbl <- file.path(project_root, "outputs", "h1_dominance_haul_table.csv")
path_out_coef <- file.path(project_root, "outputs", "h1_dominance_partial_r2_coefficients.csv")
path_out_summary <- file.path(project_root, "outputs", "h1_dominance_partial_r2_summary.csv")
path_out_anova <- file.path(project_root, "outputs", "h1_dominance_partial_r2_anova.csv")

stopifnot(file.exists(path_haul_tbl))

haul <- read.csv(path_haul_tbl, stringsAsFactors = FALSE)
n_total <- nrow(haul)

# ---------------------------------------------------------------------------
# Missingness / usability check
# ---------------------------------------------------------------------------
required_cols <- c("ln_ratio", "D", "size_CV", "B_obs")
usable <- haul %>%
  filter(
    if_all(all_of(required_cols), ~ is.finite(.x)),
    B_obs > 0
  )
n_used <- nrow(usable)
n_dropped <- n_total - n_used

cat("Loaded", path_haul_tbl, "\n")
cat("n rows in table:", n_total, "\n")
cat("n usable for regression (finite ln_ratio/D/size_CV, B_obs > 0):", n_used, "\n")
cat("n dropped for missingness/non-finite values:", n_dropped, "\n\n")

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
m0 <- lm(ln_ratio ~ log(B_obs), data = usable)
m1 <- lm(ln_ratio ~ log(B_obs) + D + size_CV, data = usable)
m2 <- lm(ln_ratio ~ log(B_obs) + D * size_CV, data = usable)

r2_m0 <- summary(m0)$r.squared
r2_m1 <- summary(m1)$r.squared
r2_m2 <- summary(m2)$r.squared

partial_r2_m1_over_m0 <- (r2_m1 - r2_m0) / (1 - r2_m0)

anova_m1_m2 <- anova(m1, m2)

# ---------------------------------------------------------------------------
# Output tables
# ---------------------------------------------------------------------------
model_summary <- tibble(
  model = c("m0: ln_ratio ~ log(B_obs)",
            "m1: ln_ratio ~ log(B_obs) + D + size_CV",
            "m2: ln_ratio ~ log(B_obs) + D * size_CV"),
  r_squared = c(r2_m0, r2_m1, r2_m2),
  adj_r_squared = c(summary(m0)$adj.r.squared, summary(m1)$adj.r.squared, summary(m2)$adj.r.squared),
  n = n_used
)

coef_m1 <- summary(m1)$coefficients
coef_m2 <- summary(m2)$coefficients

coef_tbl <- bind_rows(
  tibble(
    model = "m1",
    term = rownames(coef_m1),
    estimate = coef_m1[, "Estimate"],
    std_error = coef_m1[, "Std. Error"],
    t_value = coef_m1[, "t value"],
    p_value = coef_m1[, "Pr(>|t|)"]
  ),
  tibble(
    model = "m2",
    term = rownames(coef_m2),
    estimate = coef_m2[, "Estimate"],
    std_error = coef_m2[, "Std. Error"],
    t_value = coef_m2[, "t value"],
    p_value = coef_m2[, "Pr(>|t|)"]
  )
)

anova_tbl <- tibble(
  model = c("m1", "m2"),
  res_df = anova_m1_m2[["Res.Df"]],
  rss = anova_m1_m2[["RSS"]],
  df = anova_m1_m2[["Df"]],
  sum_of_sq = anova_m1_m2[["Sum of Sq"]],
  f_value = anova_m1_m2[["F"]],
  p_value = anova_m1_m2[["Pr(>F)"]]
)

write.csv(model_summary, path_out_summary, row.names = FALSE)
write.csv(coef_tbl, path_out_coef, row.names = FALSE)
write.csv(anova_tbl, path_out_anova, row.names = FALSE)

# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------
cat("=== R^2 by model ===\n")
print(as.data.frame(model_summary))

cat("\nHeadline: partial R^2 of {D, size_CV} beyond m0 =",
    sprintf("%.4f", partial_r2_m1_over_m0), "\n")
cat("  (= (R2_m1 - R2_m0) / (1 - R2_m0) = (", sprintf("%.4f", r2_m1), "-", sprintf("%.4f", r2_m0),
    ") / (1 -", sprintf("%.4f", r2_m0), "))\n\n")

cat("=== m1 coefficients (D, size_CV) ===\n")
print(coef_m1[c("D", "size_CV"), ])

cat("\n=== m2 interaction term (D:size_CV) ===\n")
print(coef_m2["D:size_CV", ])

cat("\n=== anova(m1, m2) ===\n")
print(anova_m1_m2)

cat("\nSaved:\n")
cat(" ", path_out_summary, "\n")
cat(" ", path_out_coef, "\n")
cat(" ", path_out_anova, "\n")
