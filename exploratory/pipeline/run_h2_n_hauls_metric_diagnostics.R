# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Haul-count vs residual-metric diagnostics for H2 (base R; no tidyverse).
#
# Supervisor check: relationship between n_hauls per rectangle and metric
# output. Partial pooling handles unequal-n for inference; this script screens
# whether metric *values* (central tendency) are artefacts of sample size, and
# documents the expected wedge/funnel in dispersion as averages converge.
#
# Decision rule:
#   - Wedge in SD/SE vs n_hauls: expected → document, no action
#   - Linear/monotone trend in mean metric vs log(n_hauls): escalate only if
#     |Spearman| >= 0.25 AND loess slope is visually clear (flagged in run log)
#
# Run: Rscript --vanilla pipeline/run_h2_n_hauls_metric_diagnostics.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/run_h2_n_hauls_metric_diagnostics.R")
}

project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

SPEARMAN_TREND_THRESHOLD <- 0.25
PHASE_LEVELS <- c("1985-1988", "1989-2000", "2001-2007", "2008-2015")
H2_YEAR_MIN <- 1985L
H2_YEAR_MAX <- 2015L

normalize_stat_rec <- function(x) {
  x <- as.character(x)
  x <- gsub('"', "", x, fixed = TRUE)
  trimws(x)
}

build_phase_factor <- function(year) {
  breaks <- c(-Inf, 1989, 2001, 2008, Inf)
  labels <- PHASE_LEVELS
  cut(year, breaks = breaks, labels = labels, right = FALSE)
}

path_panel_csv <- file.path(project_root, "outputs", "h2_rectangle_panel.csv")
path_panel_rds <- file.path(project_root, "outputs", "h2_rectangle_panel.rds")
path_haul <- file.path(project_root, "outputs", "haul_eeos_predictions.rds")
path_couce <- file.path(project_root, "outputs", "h2_couce_year_effort.rds")
path_pooling <- file.path(project_root, "outputs", "h2h3_wb_partial_pooling.csv")

path_wb_models <- file.path(project_root, "outputs", "h2h3_wb_model_objects.rds")
path_diag <- file.path(project_root, "outputs", "h2_n_hauls_metric_diagnostics.csv")
path_bins <- file.path(project_root, "outputs", "h2_n_hauls_metric_bins.csv")
path_sens <- file.path(project_root, "outputs", "h2_n_hauls_metric_sensitivity.csv")
path_run_log <- file.path(project_root, "outputs", "h2_n_hauls_metric_run_log.md")
fig_dir <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
path_fig_metric <- file.path(fig_dir, "h2_n_hauls_vs_metric.png")
path_fig_wedge <- file.path(fig_dir, "h2_n_hauls_wedge_sd.png")
path_fig_phase <- file.path(fig_dir, "h2_n_hauls_vs_metric_by_phase.png")

stopifnot(
  file.exists(path_haul),
  file.exists(path_couce),
  file.exists(path_pooling),
  file.exists(path_wb_models),
  file.exists(path_panel_csv) || file.exists(path_panel_rds)
)

ols_fp_row <- function(y, fp, weights = NULL, label, phase, n_rect, response) {
  ok <- is.finite(y) & is.finite(fp)
  if (!is.null(weights)) ok <- ok & is.finite(weights) & weights > 0
  y <- y[ok]
  fp <- fp[ok]
  w <- if (is.null(weights)) NULL else weights[ok]
  if (length(y) < 5L) {
    return(data.frame(
      sensitivity = label, phase = phase, response = response,
      n_rectangles = length(y), slope = NA_real_, se = NA_real_,
      p_value = NA_real_, r2 = NA_real_, stringsAsFactors = FALSE
    ))
  }
  fit <- if (is.null(w)) lm(y ~ fp) else lm(y ~ fp, weights = w)
  sm <- summary(fit)
  data.frame(
    sensitivity = label,
    phase = phase,
    response = response,
    n_rectangles = length(y),
    slope = unname(coef(fit)[["fp"]]),
    se = unname(sm$coefficients["fp", "Std. Error"]),
    p_value = unname(sm$coefficients["fp", "Pr(>|t|)"]),
    r2 = sm$r.squared,
    stringsAsFactors = FALSE
  )
}

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

trend_stats <- function(x, y, metric_name, layer, phase = "all") {
  ok <- is.finite(x) & is.finite(y) & x > 0
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  empty <- data.frame(
    layer = layer,
    phase = phase,
    metric = metric_name,
    n = n,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    ols_slope_log_n = NA_real_,
    ols_se = NA_real_,
    ols_p = NA_real_,
    ols_r2 = NA_real_,
    mean_trend_flag = FALSE,
    stringsAsFactors = FALSE
  )
  if (n < 5L) return(empty)

  ct <- suppressWarnings(cor.test(log(x), y, method = "spearman", exact = FALSE))
  fit <- lm(y ~ log(x))
  sm <- summary(fit)
  slope <- unname(coef(fit)[["log(x)"]])
  se <- unname(sm$coefficients["log(x)", "Std. Error"])
  p_ols <- unname(sm$coefficients["log(x)", "Pr(>|t|)"])
  rho <- unname(ct$estimate)
  flag <- is.finite(rho) && abs(rho) >= SPEARMAN_TREND_THRESHOLD
  data.frame(
    layer = layer,
    phase = phase,
    metric = metric_name,
    n = n,
    spearman_rho = rho,
    spearman_p = ct$p.value,
    ols_slope_log_n = slope,
    ols_se = se,
    ols_p = p_ols,
    ols_r2 = sm$r.squared,
    mean_trend_flag = flag,
    stringsAsFactors = FALSE
  )
}

ntile5 <- function(x) {
  r <- rank(x, ties.method = "first")
  as.integer(cut(r, breaks = quantile(r, probs = seq(0, 1, 0.2)),
                 include.lowest = TRUE, labels = FALSE))
}

iqr_val <- function(x) {
  as.numeric(diff(quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)))
}

bin_summary <- function(df, layer, phase = "all") {
  df$n_haul_bin <- ntile5(df$n_hauls)
  do.call(rbind, lapply(sort(unique(df$n_haul_bin)), function(b) {
    sub <- df[df$n_haul_bin == b, , drop = FALSE]
    data.frame(
      layer = layer,
      phase = phase,
      n_haul_bin = b,
      n_rectangles = nrow(sub),
      n_hauls_min = min(sub$n_hauls),
      n_hauls_max = max(sub$n_hauls),
      n_hauls_median = median(sub$n_hauls),
      mean_abs_residual_mean = mean(sub$mean_abs_residual, na.rm = TRUE),
      mean_abs_residual_iqr = iqr_val(sub$mean_abs_residual),
      mean_residual_mean = mean(sub$mean_residual, na.rm = TRUE),
      mean_residual_iqr = iqr_val(sub$mean_residual),
      sd_abs_residual_mean = mean(sub$sd_abs_residual, na.rm = TRUE),
      se_abs_residual_mean = mean(sub$se_abs_residual, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

aggregate_rect <- function(dat, by_phase = FALSE) {
  if (by_phase) {
    keys <- interaction(dat$stat_rec, dat$phase, drop = TRUE)
    split_idx <- split(seq_len(nrow(dat)), keys)
  } else {
    split_idx <- split(seq_len(nrow(dat)), dat$stat_rec)
  }
  rows <- lapply(split_idx, function(idx) {
    sub <- dat[idx, , drop = FALSE]
    n <- nrow(sub)
    sd_abs <- if (n > 1L) sd(sub$abs_residual) else NA_real_
    sd_res <- if (n > 1L) sd(sub$residual) else NA_real_
    out <- data.frame(
      stat_rec = sub$stat_rec[1L],
      n_hauls = n,
      mean_abs_residual = mean(sub$abs_residual, na.rm = TRUE),
      mean_residual = mean(sub$residual, na.rm = TRUE),
      sd_abs_residual = sd_abs,
      sd_residual = sd_res,
      se_abs_residual = sd_abs / sqrt(n),
      stringsAsFactors = FALSE
    )
    if (by_phase) out$phase <- as.character(sub$phase[1L])
    out
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Layer 1 — rectangle panel
# ---------------------------------------------------------------------------
if (file.exists(path_panel_rds)) {
  panel <- readRDS(path_panel_rds)
} else {
  panel <- read.csv(path_panel_csv, stringsAsFactors = FALSE)
}
panel$stat_rec <- normalize_stat_rec(panel$stat_rec)
panel$se_abs_residual <- panel$sd_abs_residual / sqrt(panel$n_hauls)

logmsg("## H2 haul-count vs metric diagnostics")
logmsg("")
logmsg(sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
logmsg("")
logmsg(sprintf(
  "Rectangle panel: %d rectangles, n_hauls range %d–%d (median %s).",
  nrow(panel),
  min(panel$n_hauls),
  max(panel$n_hauls),
  format(median(panel$n_hauls), digits = 4)
))
logmsg(sprintf(
  "Mean-trend escalate threshold: |Spearman rho| >= %.2f (central tendency only).",
  SPEARMAN_TREND_THRESHOLD
))
logmsg("")

diag_rows <- rbind(
  trend_stats(panel$n_hauls, panel$mean_abs_residual, "mean_abs_residual", "rectangle_panel"),
  trend_stats(panel$n_hauls, panel$mean_residual, "mean_residual", "rectangle_panel"),
  trend_stats(panel$n_hauls, panel$sd_abs_residual, "sd_abs_residual", "rectangle_panel"),
  trend_stats(panel$n_hauls, panel$se_abs_residual, "se_abs_residual", "rectangle_panel")
)
bins_out <- bin_summary(panel, "rectangle_panel")

# ---------------------------------------------------------------------------
# Layer 2 — live WB universe + phase
# ---------------------------------------------------------------------------
haul <- readRDS(path_haul)
couce_year <- readRDS(path_couce)
haul$stat_rec <- normalize_stat_rec(haul$stat_rec)
couce_year$stat_rec <- normalize_stat_rec(couce_year$stat_rec)

haul <- haul[
  haul$stat_rec %in% panel$stat_rec &
    is.finite(haul$residual) &
    is.finite(haul$abs_residual) &
    haul$year >= H2_YEAR_MIN &
    haul$year <= H2_YEAR_MAX,
  ,
  drop = FALSE
]

# Inner join Couce year effort
couce_key <- paste(couce_year$stat_rec, couce_year$year, sep = "\r")
haul_key <- paste(haul$stat_rec, haul$year, sep = "\r")
keep <- haul_key %in% couce_key
dat <- haul[keep, , drop = FALSE]
dat$phase <- build_phase_factor(dat$year)

logmsg(sprintf(
  "Live WB universe: %d hauls in %d rectangles.",
  nrow(dat),
  length(unique(dat$stat_rec))
))

wb_rect <- aggregate_rect(dat, by_phase = FALSE)
diag_rows <- rbind(
  diag_rows,
  trend_stats(wb_rect$n_hauls, wb_rect$mean_abs_residual, "mean_abs_residual", "wb_universe"),
  trend_stats(wb_rect$n_hauls, wb_rect$mean_residual, "mean_residual", "wb_universe"),
  trend_stats(wb_rect$n_hauls, wb_rect$sd_abs_residual, "sd_abs_residual", "wb_universe"),
  trend_stats(wb_rect$n_hauls, wb_rect$sd_residual, "sd_residual", "wb_universe"),
  trend_stats(wb_rect$n_hauls, wb_rect$se_abs_residual, "se_abs_residual", "wb_universe")
)
bins_out <- rbind(bins_out, bin_summary(wb_rect, "wb_universe"))

wb_phase_rect <- aggregate_rect(dat, by_phase = TRUE)
for (ph in PHASE_LEVELS) {
  sub <- wb_phase_rect[wb_phase_rect$phase == ph, , drop = FALSE]
  diag_rows <- rbind(
    diag_rows,
    trend_stats(sub$n_hauls, sub$mean_abs_residual, "mean_abs_residual", "wb_phase", ph),
    trend_stats(sub$n_hauls, sub$mean_residual, "mean_residual", "wb_phase", ph),
    trend_stats(sub$n_hauls, sub$sd_abs_residual, "sd_abs_residual", "wb_phase", ph),
    trend_stats(sub$n_hauls, sub$se_abs_residual, "se_abs_residual", "wb_phase", ph)
  )
}

# ---------------------------------------------------------------------------
# Layer 3 — partial pooling
# ---------------------------------------------------------------------------
pooling <- read.csv(path_pooling, stringsAsFactors = FALSE)
if ("model_id" %in% names(pooling)) {
  primary <- grepl("primary|wb_primary", pooling$model_id)
  if (any(primary)) pooling <- pooling[primary, , drop = FALSE]
}
pool_ct <- suppressWarnings(
  cor.test(log(pooling$n_hauls), pooling$shrinkage_ratio, method = "spearman", exact = FALSE)
)
logmsg("")
logmsg("### Partial pooling (existing artefact)")
logmsg(sprintf(
  "Spearman(log n_hauls, shrinkage_ratio) = %.3f (p = %.3g, n = %d).",
  unname(pool_ct$estimate),
  pool_ct$p.value,
  nrow(pooling)
))
logmsg("Shrinkage rising with n is the model answer to unequal precision;")
logmsg("it is separate from whether the DV mean trends with n.")
logmsg("")

diag_rows <- rbind(
  diag_rows,
  data.frame(
    layer = "partial_pooling",
    phase = "all",
    metric = "shrinkage_ratio",
    n = nrow(pooling),
    spearman_rho = unname(pool_ct$estimate),
    spearman_p = pool_ct$p.value,
    ols_slope_log_n = NA_real_,
    ols_se = NA_real_,
    ols_p = NA_real_,
    ols_r2 = NA_real_,
    mean_trend_flag = FALSE,
    stringsAsFactors = FALSE
  )
)

# ---------------------------------------------------------------------------
# Decision rule + sensitivity (only if mean trend flagged)
# ---------------------------------------------------------------------------
mean_metrics <- c("mean_abs_residual", "mean_residual")
mean_flags <- diag_rows[
  diag_rows$metric %in% mean_metrics &
    diag_rows$layer %in% c("rectangle_panel", "wb_universe", "wb_phase"),
  ,
  drop = FALSE
]
any_mean_trend <- any(mean_flags$mean_trend_flag, na.rm = TRUE)

wedge_rows <- diag_rows[
  diag_rows$metric %in% c("se_abs_residual", "sd_abs_residual", "sd_residual") &
    diag_rows$layer %in% c("rectangle_panel", "wb_universe"),
  ,
  drop = FALSE
]
wedge_ok <- any(wedge_rows$spearman_rho < 0, na.rm = TRUE)

logmsg("### Decision rule")
logmsg(sprintf(
  "Wedge / convergence (SE or SD vs n declining): %s",
  if (wedge_ok) "YES — expected pattern present" else "NOT CLEAR — still document averaging expectation"
))
logmsg(sprintf(
  "Linear/monotone mean-metric trend (|rho| >= %.2f): %s",
  SPEARMAN_TREND_THRESHOLD,
  if (any_mean_trend) "YES — escalate (sensitivity)" else "NO — stop; no further adjustment"
))
logmsg("")

sens_rows <- NULL
if (any_mean_trend) {
  flagged <- mean_flags[mean_flags$mean_trend_flag, , drop = FALSE]
  logmsg("Flagged mean-trend rows:")
  for (i in seq_len(nrow(flagged))) {
    logmsg(sprintf(
      "  - %s / %s / %s: Spearman rho = %.3f (p = %.3g)",
      flagged$layer[i], flagged$phase[i], flagged$metric[i],
      flagged$spearman_rho[i], flagged$spearman_p[i]
    ))
  }
  logmsg("")
  logmsg("### Sensitivity (escalation)")
  logmsg(
    "Rectangle-level OLS of residual metric ~ FP_between: unweighted, ",
    "n-weighted, inverse-SE^2 weighted, and lowest-n-decile dropped. ",
    "Descriptive robustness only — not a replacement for the haul-level WB primary."
  )
  logmsg("")

  wb_mod <- readRDS(path_wb_models)
  wb_dat <- as.data.frame(wb_mod$data)
  wb_dat$stat_rec <- normalize_stat_rec(wb_dat$stat_rec)
  wb_dat$phase <- as.character(wb_dat$phase)

  # FP_between is rectangle-constant; take one value per rectangle
  fp_map <- wb_dat[!duplicated(wb_dat$stat_rec), c("stat_rec", "FP_between")]
  panel_fp <- merge(panel, fp_map, by = "stat_rec", all.x = TRUE)

  # Legacy panel: overall weighted OLS (always useful context when escalating)
  sens_rows <- rbind(
    ols_fp_row(panel_fp$mean_abs_residual, panel_fp$FP_between, NULL,
               "panel_unweighted", "all", nrow(panel_fp), "mean_abs_residual"),
    ols_fp_row(panel_fp$mean_abs_residual, panel_fp$FP_between, panel_fp$n_hauls,
               "panel_weight_n", "all", nrow(panel_fp), "mean_abs_residual"),
    ols_fp_row(panel_fp$mean_abs_residual, panel_fp$FP_between,
               1 / (panel_fp$se_abs_residual^2),
               "panel_weight_inv_se2", "all", nrow(panel_fp), "mean_abs_residual"),
    ols_fp_row(panel_fp$mean_residual, panel_fp$FP_between, NULL,
               "panel_unweighted", "all", nrow(panel_fp), "mean_residual"),
    ols_fp_row(panel_fp$mean_residual, panel_fp$FP_between, panel_fp$n_hauls,
               "panel_weight_n", "all", nrow(panel_fp), "mean_residual"),
    ols_fp_row(panel_fp$mean_residual, panel_fp$FP_between,
               1 / (panel_fp$se_abs_residual^2),
               "panel_weight_inv_se2", "all", nrow(panel_fp), "mean_residual")
  )

  flagged_phases <- unique(as.character(flagged$phase[flagged$layer == "wb_phase"]))
  # Also escalate overall layers if flagged
  if (any(flagged$layer %in% c("rectangle_panel", "wb_universe"))) {
    flagged_phases <- unique(c(flagged_phases, PHASE_LEVELS))
  }

  for (ph in flagged_phases) {
    sub_h <- wb_dat[wb_dat$phase == ph, , drop = FALSE]
    # rectangle means within phase
    keys <- unique(sub_h$stat_rec)
    rect <- do.call(rbind, lapply(keys, function(sr) {
      h <- sub_h[sub_h$stat_rec == sr, , drop = FALSE]
      n <- nrow(h)
      sd_abs <- if (n > 1L) sd(h$abs_residual) else NA_real_
      data.frame(
        stat_rec = sr,
        n_hauls = n,
        mean_abs_residual = mean(h$abs_residual, na.rm = TRUE),
        mean_residual = mean(h$residual, na.rm = TRUE),
        se_abs_residual = sd_abs / sqrt(n),
        FP_between = h$FP_between[1L],
        stringsAsFactors = FALSE
      )
    }))
    n_cut <- as.numeric(quantile(rect$n_hauls, 0.1, names = FALSE, type = 7))
    keep <- rect$n_hauls > n_cut
    rect_drop <- rect[keep, , drop = FALSE]

    for (resp in c("mean_abs_residual", "mean_residual")) {
      sens_rows <- rbind(
        sens_rows,
        ols_fp_row(rect[[resp]], rect$FP_between, NULL,
                   "phase_unweighted", ph, nrow(rect), resp),
        ols_fp_row(rect[[resp]], rect$FP_between, rect$n_hauls,
                   "phase_weight_n", ph, nrow(rect), resp),
        ols_fp_row(rect[[resp]], rect$FP_between, 1 / (rect$se_abs_residual^2),
                   "phase_weight_inv_se2", ph, nrow(rect), resp),
        ols_fp_row(rect_drop[[resp]], rect_drop$FP_between, NULL,
                   "phase_drop_lowest_n_decile", ph, nrow(rect_drop), resp)
      )
    }

    # Re-check mean vs n after dropping lowest decile
    post_abs <- trend_stats(
      rect_drop$n_hauls, rect_drop$mean_abs_residual,
      "mean_abs_residual", "wb_phase_drop_low_n", ph
    )
    post_signed <- trend_stats(
      rect_drop$n_hauls, rect_drop$mean_residual,
      "mean_residual", "wb_phase_drop_low_n", ph
    )
    diag_rows <- rbind(diag_rows, post_abs, post_signed)

    logmsg(sprintf("Phase %s: n_hauls decile cut = %.1f (drop n <= cut); retained %d / %d rectangles.",
                   ph, n_cut, nrow(rect_drop), nrow(rect)))
    logmsg(sprintf(
      "  After drop — mean_abs vs log n: Spearman rho = %.3f (flag = %s)",
      post_abs$spearman_rho, post_abs$mean_trend_flag
    ))
    logmsg(sprintf(
      "  After drop — mean_residual vs log n: Spearman rho = %.3f (flag = %s)",
      post_signed$spearman_rho, post_signed$mean_trend_flag
    ))
  }

  logmsg("")
  logmsg("Sensitivity slopes (residual metric ~ FP_between):")
  for (i in seq_len(nrow(sens_rows))) {
    logmsg(sprintf(
      "  %s | %s | %s: slope = %.4f (SE = %.4f, p = %.3g, n = %d)",
      sens_rows$sensitivity[i], sens_rows$phase[i], sens_rows$response[i],
      sens_rows$slope[i], sens_rows$se[i], sens_rows$p_value[i],
      sens_rows$n_rectangles[i]
    ))
  }

  # Robustness verdict: compare unweighted vs weighted / drop for flagged phases
  robust_ok <- TRUE
  for (ph in flagged_phases) {
    for (resp in c("mean_abs_residual", "mean_residual")) {
      base <- sens_rows$slope[
        sens_rows$sensitivity == "phase_unweighted" &
          sens_rows$phase == ph & sens_rows$response == resp
      ]
      alts <- sens_rows$slope[
        sens_rows$sensitivity %in% c("phase_weight_n", "phase_weight_inv_se2",
                                     "phase_drop_lowest_n_decile") &
          sens_rows$phase == ph & sens_rows$response == resp
      ]
      if (length(base) == 1L && length(alts) > 0L && is.finite(base)) {
        # Material change: sign flip or |alt - base| > |base|
        if (any(sign(alts) != sign(base) & abs(alts) > 1e-8, na.rm = TRUE) ||
            any(abs(alts - base) > abs(base) + 1e-8, na.rm = TRUE)) {
          robust_ok <- FALSE
        }
      }
    }
  }
  logmsg("")
  if (robust_ok) {
    logmsg(
      "Robustness: weighted / drop-low-n slopes keep the same sign and similar ",
      "magnitude as unweighted phase OLS — H2 FP_between association for flagged ",
      "phase(s) is not an artefact of sparse rectangles. No change to primary WB model."
    )
  } else {
    logmsg(
      "Robustness: at least one weighted / drop-low-n slope changed sign or more ",
      "than doubled |magnitude| relative to unweighted — flag for discussion; ",
      "consider reporting phase-specific min-n sensitivity alongside primary results."
    )
  }
} else {
  logmsg("No central-tendency trend exceeds the pre-registered threshold.")
  logmsg("Action: document expected wedge; cite partial pooling for unequal-n inference;")
  logmsg("no min-n / weighted sensitivity required for this check.")
}
logmsg("")

# ---------------------------------------------------------------------------
# Figures (base graphics)
# ---------------------------------------------------------------------------
add_loess <- function(x, y, col = "#c0392b", lwd = 2) {
  ok <- is.finite(x) & is.finite(y) & x > 0
  if (sum(ok) < 10L) return(invisible(NULL))
  xok <- x[ok]
  yok <- y[ok]
  lx <- log10(xok)
  fit <- tryCatch(loess(yok ~ lx, span = 0.75), error = function(e) NULL)
  if (is.null(fit)) return(invisible(NULL))
  xx <- seq(min(xok), max(xok), length.out = 100)
  yy <- predict(fit, newdata = data.frame(lx = log10(xx)))
  lines(xx, yy, col = col, lwd = lwd)
}

png(path_fig_metric, width = 1000, height = 450, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1.2), oma = c(3.5, 0, 2, 0))
plot(panel$n_hauls, panel$mean_abs_residual,
     log = "x", pch = 16, col = adjustcolor("#2c3e50", 0.45), cex = 0.7,
     xlab = "N hauls (log)", ylab = "mean_abs_residual",
     main = "mean |log residual|")
add_loess(panel$n_hauls, panel$mean_abs_residual)
plot(panel$n_hauls, panel$mean_residual,
     log = "x", pch = 16, col = adjustcolor("#2c3e50", 0.45), cex = 0.7,
     xlab = "N hauls (log)", ylab = "mean_residual",
     main = "mean signed log residual")
add_loess(panel$n_hauls, panel$mean_residual)
mtext("Rectangle residual metrics vs haul count (H2 panel)",
      outer = TRUE, font = 2, cex = 1.05)
mtext(
  paste0(
    "Loess on log n. Escalate only if mean trends with n (|Spearman| >= ",
    SPEARMAN_TREND_THRESHOLD, "). Partial pooling handles unequal-n inference separately."
  ),
  side = 1, outer = TRUE, cex = 0.7, line = 1.8, adj = 0
)
par(op)
dev.off()

png(path_fig_wedge, width = 1000, height = 450, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1.2), oma = c(3.5, 0, 2, 0))
plot(panel$n_hauls, panel$sd_abs_residual,
     log = "x", pch = 16, col = adjustcolor("#2c3e50", 0.45), cex = 0.7,
     xlab = "N hauls (log)", ylab = "sd_abs_residual",
     main = "Within-rectangle SD")
add_loess(panel$n_hauls, panel$sd_abs_residual, col = "#2980b9")
plot(panel$n_hauls, panel$se_abs_residual,
     log = "x", pch = 16, col = adjustcolor("#2c3e50", 0.45), cex = 0.7,
     xlab = "N hauls (log)", ylab = "se_abs_residual",
     main = "SE = SD / sqrt(n)")
add_loess(panel$n_hauls, panel$se_abs_residual, col = "#2980b9")
mtext("Wedge / convergence: residual dispersion vs haul count",
      outer = TRUE, font = 2, cex = 1.05)
mtext(
  "Narrowing SE with n is expected as averages converge. Wedge alone does not trigger sensitivity analysis.",
  side = 1, outer = TRUE, cex = 0.7, line = 1.8, adj = 0
)
par(op)
dev.off()

png(path_fig_phase, width = 1100, height = 700, res = 120)
op <- par(mfrow = c(2, 4), mar = c(3.5, 3.5, 2.5, 0.8), oma = c(2, 2, 3, 0))
for (metric in c("mean_abs_residual", "mean_residual")) {
  for (ph in PHASE_LEVELS) {
    sub <- wb_phase_rect[wb_phase_rect$phase == ph, , drop = FALSE]
    y <- sub[[metric]]
    plot(sub$n_hauls, y,
         log = "x", pch = 16, col = adjustcolor("#2c3e50", 0.35), cex = 0.55,
         xlab = if (metric == "mean_residual") "N hauls" else "",
         ylab = if (ph == PHASE_LEVELS[1L]) metric else "",
         main = if (metric == "mean_abs_residual") ph else "")
    add_loess(sub$n_hauls, y, lwd = 1.6)
  }
}
mtext("Phase-stratified residual metrics vs haul count (live WB universe)",
      outer = TRUE, font = 2, cex = 1.05, line = 1)
par(op)
dev.off()

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
write.csv(diag_rows, path_diag, row.names = FALSE)
write.csv(bins_out, path_bins, row.names = FALSE)
if (!is.null(sens_rows)) {
  write.csv(sens_rows, path_sens, row.names = FALSE)
}

logmsg("### Key Spearman results (central tendency)")
for (i in seq_len(nrow(mean_flags))) {
  logmsg(sprintf(
    "  %s | %s | %s: rho = %.3f (p = %.3g), OLS slope = %.4f, flag = %s",
    mean_flags$layer[i], mean_flags$phase[i], mean_flags$metric[i],
    mean_flags$spearman_rho[i], mean_flags$spearman_p[i],
    mean_flags$ols_slope_log_n[i], mean_flags$mean_trend_flag[i]
  ))
}
logmsg("")
logmsg("### Key Spearman results (dispersion / wedge)")
for (i in seq_len(nrow(wedge_rows))) {
  logmsg(sprintf(
    "  %s | %s: rho = %.3f (p = %.3g)",
    wedge_rows$layer[i], wedge_rows$metric[i],
    wedge_rows$spearman_rho[i], wedge_rows$spearman_p[i]
  ))
}
logmsg("")
logmsg("### Outputs")
logmsg(sprintf("- %s", path_diag))
logmsg(sprintf("- %s", path_bins))
if (!is.null(sens_rows)) logmsg(sprintf("- %s", path_sens))
logmsg(sprintf("- %s", path_fig_metric))
logmsg(sprintf("- %s", path_fig_wedge))
logmsg(sprintf("- %s", path_fig_phase))
logmsg("")
logmsg("### Verdict summary")
overall_null <- !any(
  mean_flags$mean_trend_flag[
    mean_flags$layer %in% c("rectangle_panel", "wb_universe")
  ],
  na.rm = TRUE
)
if (!any_mean_trend) {
  logmsg(
    "Unequal haul counts are handled by partial pooling for inference. ",
    "Rectangle-level residual metrics show the expected wedge in dispersion ",
    "with increasing n. No systematic linear association between haul count ",
    "and residual metric central tendency exceeded |Spearman| >= 0.25; ",
    "no further adjustment was applied."
  )
} else if (overall_null) {
  logmsg(
    "Unequal haul counts are handled by partial pooling for inference. ",
    "Overall rectangle metrics show no mean trend with n_hauls and the expected ",
    "SE wedge. A phase-specific mean trend (see flagged rows) triggered ",
    "weighted / drop-low-n sensitivity of residual ~ FP_between; see sensitivity ",
    "table and robustness note above. Primary haul-level WB model unchanged."
  )
} else {
  logmsg(
    "Unequal haul counts are handled by partial pooling for inference. ",
    "A mean-metric trend with n_hauls exceeded the pre-registered threshold; ",
    "see flagged rows and sensitivity table."
  )
}

writeLines(run_log, path_run_log)
message("Wrote: ", path_run_log)
message("Wrote: ", path_diag)
if (!is.null(sens_rows)) message("Wrote: ", path_sens)
message("Wrote: ", path_fig_metric)
message("Wrote: ", path_fig_wedge)
message("Wrote: ", path_fig_phase)

dd_fig <- file.path(project_root, "display_discussion", "figures")
dir.create(dd_fig, recursive = TRUE, showWarnings = FALSE)
for (fp in c(path_fig_metric, path_fig_wedge, path_fig_phase)) {
  file.copy(fp, file.path(dd_fig, basename(fp)), overwrite = TRUE)
}
message("Copied figures to display_discussion/figures/")
