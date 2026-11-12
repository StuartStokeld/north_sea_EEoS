# Diagnostics for the primary CAR model fitted to the absolute residual |ln B_obs - ln B_pred|
# (the corrected specification, F-A5) alongside the signed residual (current code).
# Read-only. Writes audit/outputs/abs_residual/ and audit/figures/abs_residual_diagnostics.png.
#
# Run from the repo root: Rscript audit/scripts/abs_residual_diagnostics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("audit/scripts/audit_refit_helpers.R")
out_dir <- "audit/outputs/abs_residual"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("audit/figures", showWarnings = FALSE)

haul <- readRDS("outputs/haul_eeos_predictions.rds")
dat <- audit_build_h2h3_data(audit_recompute_residuals(haul))

moments <- function(x) {
  m <- mean(x); s <- sd(x)
  c(n = length(x), mean = m, sd = s, median = median(x), min = min(x), max = max(x),
    skewness = mean((x - m)^3) / s^3, excess_kurtosis = mean((x - m)^4) / s^4 - 3)
}
set.seed(42)
rows <- list()
fits <- list()
for (rsp in c("residual", "abs_residual")) {
  s <- audit_fit_primary_car(dat, rsp)
  fit <- attr(s, "fit")
  fits[[rsp]] <- fit
  e <- as.numeric(residuals(fit))
  fv <- as.numeric(fitted(fit))
  sw <- shapiro.test(sample(e, 5000))
  qq <- qqnorm(e, plot.it = FALSE)
  rows[[rsp]] <- tibble(
    response = rsp,
    stat = c(names(moments(dat[[rsp]])), "share_response_lt_0",
             paste0("model_resid_", names(moments(e))),
             "shapiro_W_5000", "shapiro_p_5000", "qq_correlation",
             "spearman_absresid_vs_fitted", "spearman_absresid_vs_fitted_p",
             "share_fitted_lt_0", "rho_CAR", "phi_residual_var", "lambda_rect_var"),
    value = c(moments(dat[[rsp]]), mean(dat[[rsp]] < 0), moments(e),
              sw$statistic, sw$p.value, cor(qq$x, qq$y),
              cor(abs(e), fv, method = "spearman"),
              cor.test(abs(e), fv, method = "spearman", exact = FALSE)$p.value,
              mean(fv < 0), attr(s, "rho"),
              spaMM::residVar(fit, which = "fit")[1],
              as.numeric(spaMM::VarCorr(fit)[1, "Variance"]))
  )
}
tab <- bind_rows(rows) %>% tidyr::pivot_wider(names_from = response, values_from = value)
print(as.data.frame(tab), digits = 4)
write_csv(tab, file.path(out_dir, "diagnostics.csv"))

png("audit/figures/abs_residual_diagnostics.png", width = 1800, height = 1200, res = 150)
op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1))
for (rsp in c("residual", "abs_residual")) {
  fit <- fits[[rsp]]
  e <- as.numeric(residuals(fit))
  fv <- as.numeric(fitted(fit))
  lab <- if (rsp == "residual") "signed residual (current code)" else "absolute residual (methods)"
  hist(dat[[rsp]], breaks = 80, main = paste("Response:", lab), xlab = rsp, col = "grey80", border = "white")
  qqnorm(e, main = paste("Normal QQ, model residuals:", lab), pch = 16, cex = 0.3)
  qqline(e, col = "red")
  plot(fv, e, pch = 16, cex = 0.3, col = grDevices::adjustcolor("black", 0.3),
       xlab = "fitted", ylab = "model residual", main = paste("Residual vs fitted:", lab))
  abline(h = 0, col = "red")
}
par(op)
dev.off()
cat("Saved audit/figures/abs_residual_diagnostics.png\n")
