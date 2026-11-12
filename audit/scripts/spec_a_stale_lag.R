# F-B15: the live Spec A (spatial lag) model joins FP_between_lag from outputs/fp_between_lag_rectangle.rds,
# written only by an exploratory script before the DataType fix. Recompute the k = 4 lag from the
# current FP_between (same neighbour lists) and refit Spec A to quantify the effect. Read-only.
suppressPackageStartupMessages(library(dplyr))
source("audit/scripts/audit_refit_helpers.R")
source("pipeline/R/h2h3_knn_spatial_helpers.R")
v2 <- readRDS("outputs/primary_model_v2.rds")
dat <- as.data.frame(v2$data)
lag_file <- readRDS("outputs/fp_between_lag_rectangle.rds")
cur <- dat %>% distinct(stat_rec = as.character(stat_rec), FP_between)
nb <- strsplit(lag_file$neighbour_stat_rec, "|", fixed = TRUE)
lag_now <- vapply(nb, function(n) mean(cur$FP_between[match(n, cur$stat_rec)]), 0)
cmp <- tibble(stat_rec = lag_file$stat_rec, lag_file = lag_file$FP_between_lag, lag_now = lag_now)
cat(sprintf("FP_between_lag stale vs recomputed: max |d| = %.3f, median |d| = %.3f, r = %.4f\n",
            max(abs(cmp$lag_file - cmp$lag_now)), median(abs(cmp$lag_file - cmp$lag_now)),
            cor(cmp$lag_file, cmp$lag_now)))
adj <- readRDS("outputs/h2h3_feasibility_round2_model_objects.rds")$adjMatrix
script <- "pipeline/run_h2h3_spec_a_car_identifiability.R"
env <- new.env(); exprs <- as.list(parse(script, keep.source = FALSE))
f <- eval(Filter(function(e) is.call(e) && length(e) >= 2 && identical(deparse(e[[2]]), "formula_car_lag"), exprs)[[1]][[3]])
cat("Spec A formula:", deparse(f), "\n")
fit_one <- function(lagvals, response) {
  d <- dat
  d$FP_between_lag <- lagvals[match(as.character(d$stat_rec), cmp$stat_rec)]
  d$stat_rec <- factor(as.character(d$stat_rec), levels = rownames(adj))
  d$phase_v2 <- factor(as.character(d$phase_v2), levels = AUDIT_PHASE_V2)
  d$y <- if (response == "abs") abs(d$residual) else d$residual
  ff <- update(f, y ~ .)
  fit <- spaMM::fitme(ff, data = d, adjMatrix = adj, method = "REML")
  bind_rows(
    extract_wb_phase_slopes_spamm(fit, "FP_between", "specA", "H2", phases = AUDIT_PHASE_V2),
    extract_wb_phase_slopes_spamm(fit, "FP_between_lag", "specA", "lag", phases = AUDIT_PHASE_V2)
  ) %>% transmute(component, phase, slope = fp_slope, se = fp_slope_se, p = p_value)
}
out <- bind_rows(
  fit_one(cmp$lag_file, "signed") %>% mutate(lag = "stale (file)", response = "signed"),
  fit_one(cmp$lag_now, "signed") %>% mutate(lag = "recomputed", response = "signed"),
  fit_one(cmp$lag_file, "abs") %>% mutate(lag = "stale (file)", response = "absolute"),
  fit_one(cmp$lag_now, "abs") %>% mutate(lag = "recomputed", response = "absolute")
)
dir.create("audit/outputs/spec_a", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out, "audit/outputs/spec_a/spec_a_stale_vs_recomputed_lag.csv")
print(as.data.frame(out %>% mutate(across(c(slope, se), ~ round(.x, 4)), p = signif(p, 3))))
