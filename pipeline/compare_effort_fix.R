# Compare live outputs across the DATRAS haul-effort fix: old (pre-fix), primary
# (exactly 30-min hauls) and sensitivity (26-34 min window).
#
# Runs (output folders; each holds the same file layout as outputs/):
#   old      outputs/archive_pre_effort_fix/      pipeline before the fix (gitignored copy)
#   primary  outputs/                             effort fix, HaulDur == 30 (primary)
#   sens     outputs/sensitivity_window_26_34/    effort fix, 26 <= HaulDur <= 34 (gitignored)
#
# Run from repo root:
#   Rscript --vanilla pipeline/compare_effort_fix.R [--timings=<timings.tsv>]
#   Rscript --vanilla pipeline/compare_effort_fix.R --old=<dir> --primary=<dir> --sens=<dir> --out=<md>
#
# Writes outputs/effort_fix_comparison.md (default) and its tables to
# outputs/effort_fix_comparison_tables/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

cli_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cli_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("pipeline", winslash = "/", mustWork = TRUE)
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
# Same project-local library the H2/H3 scripts prepend (spaMM lives there)
local_lib <- file.path(project_root, ".R_libs")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(glmmTMB)
  library(spaMM)
})
source(file.path(r_dir, "h2h3_results_helpers.R"))

arg_val <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
}
dirs <- c(
  old = arg_val("old", file.path(project_root, "outputs", "archive_pre_effort_fix")),
  primary = arg_val("primary", file.path(project_root, "outputs")),
  sens = arg_val("sens", file.path(project_root, "outputs", "sensitivity_window_26_34"))
)
run_ids <- names(dirs)
run_short <- c(old = "old", primary = "30 min", sens = "26–34")
run_long <- c(old = "old", primary = "30 min (primary)", sens = "26–34 (sensitivity)")
path_md <- arg_val("out", file.path(project_root, "outputs", "effort_fix_comparison.md"))
path_timings <- arg_val("timings", NULL)
dir_tables <- file.path(dirname(path_md), sub("\\.md$", "_tables", basename(path_md)))
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(dir.exists(dirs)))
pth <- function(run, f) file.path(dirs[[run]], f)
rel <- function(d) sub(paste0(project_root, "/"), "", d)

# --- helpers ---------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "–", formatC(x, digits = d, format = "f"))
fmt_p <- function(p) ifelse(is.na(p), "–", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))
fmt_n <- function(n) format(n, big.mark = ",")
md_table <- function(df) {
  df <- as.data.frame(df)
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, rows)
}
haul_country <- function(haul_id) vapply(strsplit(haul_id, " "), `[`, "", 4L)
haul_key6 <- function(haul_id) vapply(strsplit(haul_id, " "), function(p) paste(p[1], p[2], p[3], p[4], p[5], p[8], sep = "_"), "")
read_rds_quiet <- function(path) suppressWarnings(readRDS(path))
write_tab <- function(df, name) write.csv(df, file.path(dir_tables, paste0(name, ".csv")), row.names = FALSE)
per_run <- function(f) setNames(lapply(run_ids, f), run_ids)
rename_runs <- function(df, labels = run_long) rename(df, any_of(setNames(run_ids, labels[run_ids])))
na0 <- function(x) ifelse(is.na(x), 0L, x)

hh <- readRDS(file.path(project_root, "outputs", "datras_hh_raw.rds")) %>%
  filter(hh_key_n == 1L) %>%
  select(haul_key, DataType, HaulDur)

md <- c(
  "# DATRAS haul-effort fix: old vs exactly 30 min (primary) vs 26–34 min (sensitivity)",
  "",
  sprintf("Generated %s by `pipeline/compare_effort_fix.R` on branch `fix/haul-effort-datatype`.", Sys.Date()),
  "",
  sprintf("- **old**: `%s`. Pipeline before the fix (DataType C per-hour counts used as raw counts; all hauls, any duration).", rel(dirs[["old"]])),
  sprintf("- **30 min (primary)**: `%s`. Effort fix: valid hauls (HaulVal V), DataType C/R/P/S, **HaulDur == 30 exactly**; C counts × HaulDur/60.", rel(dirs[["primary"]])),
  sprintf("- **26–34 (sensitivity)**: `%s`. Same fix with the window 26 ≤ HaulDur ≤ 34 (run first as the primary, then moved there unchanged).", rel(dirs[["sens"]])),
  "- Model formulas, phases, priors and fishing-pressure construction are unchanged. `B_obs` (FishGlob raw haul weight) is unchanged.",
  ""
)

# =================================================================================
# H1
# =================================================================================
pred <- per_run(function(r) readRDS(pth(r, "haul_eeos_predictions.rds")))
bm <- per_run(function(r) readRDS(pth(r, "haul_h1_benchmarks.rds")))

same_ids <- Reduce(intersect, lapply(bm, `[[`, "haul_id"))
n_set <- vapply(bm, nrow, 1L)
n_primary_not_in <- vapply(run_ids, function(r) length(setdiff(bm$primary$haul_id, bm[[r]]$haul_id)), 1L)

pattern_offset <- function(lo, lp) {
  ss_tot <- sum((lo - mean(lo))^2)
  bias <- mean(lp - lo)
  tibble(
    n = length(lo),
    log_r2_1to1 = 1 - sum((lo - lp)^2) / ss_tot,
    pearson_r = cor(lo, lp),
    ols_slope = unname(coef(lm(lo ~ lp))[2]),
    median_ratio = median(exp(lp - lo)),
    mean_log_bias = bias,
    log_r2_bias_removed = 1 - sum((lo - (lp - bias))^2) / ss_tot
  )
}
# Benchmarks as defined in run_h1_harte_baseline.R / run_h1_lne_reference.R
bench_defs <- list(
  eeos_biomass = list(label = "EEoS (S, N, E), unfitted", obs = function(d) d$ln_B_obs, pred = function(d) d$ln_B_pred),
  productivity_1to1 = list(label = "Productivity 1:1, log(E × m_min), headline baseline", obs = function(d) d$ln_B_obs, pred = function(d) d$ln_E_calibrated),
  productivity_1to1_uncalibrated = list(label = "Productivity 1:1 uncalibrated, log E_raw", obs = function(d) d$ln_B_obs, pred = function(d) d$ln_E_raw),
  ln_e_ols = list(label = "ln(E) OLS, fitted", obs = function(d) d$ln_B_obs, pred = function(d) d$ln_B_pred_E),
  productivity_ratio = list(label = "Productivity ratio E/B^(3/4) (Harte Fig 2; different target)", obs = function(d) log(d$ratio_obs), pred = function(d) log(d$ratio_pred))
)
score_runs <- function(sets) {
  bind_rows(lapply(names(bench_defs), function(id) {
    def <- bench_defs[[id]]
    bind_rows(lapply(run_ids, function(r) {
      d <- sets[[r]]
      pattern_offset(def$obs(d), def$pred(d)) %>% mutate(run = r)
    })) %>% mutate(model_id = id, model = def$label, .before = 1)
  }))
}
h1_own <- score_runs(bm)
bm_same <- lapply(bm, function(b) b[match(same_ids, b$haul_id), ])
stopifnot(all(abs(bm_same$old$B_obs - bm_same$primary$B_obs) < 1e-6),
          all(abs(bm_same$sens$B_obs - bm_same$primary$B_obs) < 1e-6))
h1_same <- score_runs(bm_same)
write_tab(h1_own, "h1_pattern_vs_offset_own_sets")
write_tab(h1_same, "h1_pattern_vs_offset_same_hauls")

# Sanity: full-set log_r2 reproduces each run's h1_model_comparison.csv
for (r in run_ids) {
  mc <- read.csv(pth(r, "h1_model_comparison.csv"), stringsAsFactors = FALSE)
  for (id in c("eeos_biomass", "productivity_1to1", "productivity_1to1_uncalibrated", "ln_e_ols")) {
    d <- bench_defs[[id]]
    stopifnot(abs(pattern_offset(d$obs(bm[[r]]), d$pred(bm[[r]]))$log_r2_1to1 - mc$log_r2[mc$model_id == id]) < 1e-8)
  }
}

bench_md <- function(tab) {
  md_table(tab %>% transmute(
    model, run = run_short[run], n = fmt_n(n),
    `1:1 log_r2` = fmt(log_r2_1to1), `Pearson r` = fmt(pearson_r), `OLS slope` = fmt(ols_slope),
    `median ratio` = fmt(median_ratio), `mean log bias` = fmt(mean_log_bias),
    `log_r2 bias-removed` = fmt(log_r2_bias_removed)
  ))
}

h1_metrics <- function(p) {
  tibble(
    n_hauls = nrow(p),
    median_B_pred_over_B_obs = median(p$B_pred / p$B_obs),
    log_r2 = log_r2(log(p$B_obs), log(p$B_pred)),
    mean_residual = mean(p$residual),
    residual_sd = sd(p$residual)
  )
}
h1_full <- bind_rows(
  bind_rows(lapply(run_ids, function(r) h1_metrics(pred[[r]]) %>% mutate(run = r, hauls = "own set"))),
  bind_rows(lapply(run_ids, function(r) h1_metrics(pred[[r]] %>% filter(haul_id %in% same_ids)) %>% mutate(run = r, hauls = "same (30-min) hauls")))
) %>% relocate(hauls, run)
write_tab(h1_full, "h1_eeos_residual_summary")

with_meta <- function(p) {
  p %>%
    mutate(haul_key = haul_key6(haul_id), Country = haul_country(haul_id)) %>%
    left_join(hh, by = "haul_key")
}
rm_ <- lapply(pred, with_meta)
dt_counts <- bind_rows(lapply(run_ids, function(r) rm_[[r]] %>% count(DataType) %>% mutate(run = r))) %>%
  mutate(DataType = coalesce(DataType, "(no unique HH)")) %>%
  pivot_wider(names_from = run, values_from = n, values_fill = 0L)
write_tab(dt_counts, "h1_hauls_by_datatype")
dur_counts <- bind_rows(lapply(run_ids, function(r) rm_[[r]] %>%
  mutate(dur = case_when(is.na(HaulDur) ~ "(no unique HH)", HaulDur == 30 ~ "30", HaulDur >= 26 & HaulDur <= 34 ~ "26–34, not 30",
                         HaulDur >= 55 ~ "≥ 55", TRUE ~ "other")) %>%
  count(dur) %>% mutate(run = r))) %>%
  pivot_wider(names_from = run, values_from = n, values_fill = 0L)
write_tab(dur_counts, "h1_hauls_by_duration")

# Same hauls: input changes by DataType (old -> 30 min)
common <- inner_join(
  pred$old %>% select(haul_id, B_pred_old = B_pred, residual_old = residual, N_old = N, E_old = E, S_old = S),
  pred$primary %>% select(haul_id, B_pred_new = B_pred, residual_new = residual, N_new = N, E_new = E, S_new = S),
  by = "haul_id"
) %>%
  mutate(haul_key = haul_key6(haul_id)) %>%
  left_join(hh, by = "haul_key")
common_tab <- common %>%
  group_by(DataType) %>%
  summarise(
    n_hauls = n(),
    median_N_new_over_old = median(N_new / N_old),
    median_E_new_over_old = median(E_new / E_old),
    share_S_unchanged = mean(S_new == S_old),
    median_Bpred_new_over_old = median(B_pred_new / B_pred_old),
    mean_residual_old = mean(residual_old),
    mean_residual_new = mean(residual_new),
    .groups = "drop"
  )
write_tab(common_tab, "h1_common_hauls_by_datatype")
# 30 min vs 26-34 on the shared 30-min hauls: inputs should be identical
ps <- inner_join(pred$primary %>% select(haul_id, S, N, E, B_pred),
                 pred$sens %>% select(haul_id, S2 = S, N2 = N, E2 = E, B_pred2 = B_pred), by = "haul_id")
ps_max_diff <- max(abs(ps$S - ps$S2), abs(ps$N / ps$N2 - 1), abs(ps$E / ps$E2 - 1), abs(ps$B_pred / ps$B_pred2 - 1))

null_tabs <- bind_rows(lapply(run_ids, function(r) {
  f <- pth(r, "null_distribution_summary.csv")
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) %>% mutate(run = r, .before = 1) else NULL
}))

md <- c(md,
  "## H1 — haul-level EEoS",
  "",
  "### Haul sets",
  "",
  sprintf("H1 hauls: old %s; 30 min %s; 26–34 %s. The 30-min set nests in both others (%d of its hauls missing from old, %d from 26–34).",
          fmt_n(n_set[["old"]]), fmt_n(n_set[["primary"]]), fmt_n(n_set[["sens"]]),
          n_primary_not_in[["old"]], n_primary_not_in[["sens"]]),
  "",
  md_table(dt_counts %>% rename_runs(run_short)),
  "",
  "By recorded HaulDur (minutes):",
  "",
  md_table(dur_counts %>% rename(HaulDur = dur) %>% rename_runs(run_short)),
  "",
  "Definitions (log scale; obs = ln B_obs, pred = ln B_pred):",
  "",
  "- **1:1 log_r2** = 1 − Σ(obs − pred)² / Σ(obs − mean obs)². Penalises pattern *and* offset; can be negative.",
  "- **Pearson r**: pattern only (invariant to offset and scale).",
  "- **OLS slope** of obs on pred: 1 means predicted spread matches observed spread; < 1 means predictions too spread out.",
  "- **median ratio** = median(B_pred / B_obs); **mean log bias** = mean(pred − obs) (> 0 = over-prediction; exp(bias) = geometric-mean ratio).",
  "- **log_r2 bias-removed** = 1:1 log_r2 after subtracting the mean log bias from pred: what 1:1 log_r2 would be with the offset gone.",
  "",
  "### Benchmarks, each run on its own haul set (as each run would be reported)",
  "",
  bench_md(h1_own),
  "",
  sprintf("### Benchmarks on the same %s hauls (the 30-min set) in every run", fmt_n(length(same_ids))),
  "",
  "Separates the effect of the count conversion from the change of haul set. On these hauls the 30-min and 26–34",
  sprintf("inputs are identical (max relative difference in S, N, E, B_pred: %s), so the 26–34 rows differ from the 30-min rows",
          formatC(ps_max_diff, format = "e", digits = 1)),
  "only for ln(E) OLS, which is fitted to each run's own full set.",
  "",
  bench_md(h1_same),
  "",
  "Notes: the productivity ratio (Harte Fig 2) compares E/B^(3/4) predicted vs observed, a different target from ln B.",
  "ln(E) OLS is fitted (to each run's full H1 set), so its mean bias is ~0 on its own set by construction.",
  "Full-set log_r2 values for EEoS, productivity 1:1 (calibrated and uncalibrated) and ln(E) OLS were checked to reproduce",
  "each run's `h1_model_comparison.csv`.",
  "",
  "### EEoS residual summary",
  "",
  md_table(h1_full %>% transmute(hauls, run = run_short[run], n_hauls = fmt_n(n_hauls), `median B_pred/B_obs` = fmt(median_B_pred_over_B_obs),
                                 log_r2 = fmt(log_r2), `mean residual (= −bias)` = fmt(mean_residual), `residual SD` = fmt(residual_sd))),
  "",
  "### Same hauls, old → 30 min: what changed in the inputs (by DataType)",
  "",
  md_table(common_tab %>% transmute(
    DataType, n_hauls,
    `median N new/old` = fmt(median_N_new_over_old), `median E new/old` = fmt(median_E_new_over_old),
    `S unchanged` = paste0(fmt(100 * share_S_unchanged, 1), "%"),
    `median B_pred new/old` = fmt(median_Bpred_new_over_old),
    `mean resid old` = fmt(mean_residual_old), `mean resid 30 min` = fmt(mean_residual_new)
  )),
  ""
)
if (nrow(null_tabs)) {
  write_tab(null_tabs, "h1_null_model")
  md <- c(md,
    "### Uniform-draw null model (`null_distribution_summary.csv`, each run's own H1 set)",
    "",
    md_table(null_tabs %>% mutate(run = run_short[run]) %>% mutate(across(where(is.numeric) & !any_of("n"), ~ fmt(.x, 3)))),
    ""
  )
}

# =================================================================================
# Residual means by country x period (C -> R switch steps)
# =================================================================================
dt_year <- readRDS(file.path(project_root, "outputs", "datras_hh_raw.rds")) %>%
  filter(HaulVal == "V") %>%
  group_by(Country, Year) %>%
  summarise(share_C = mean(DataType == "C", na.rm = TRUE), n = n(), .groups = "drop")
write_tab(dt_year, "datatype_share_C_by_country_year")

switches <- tibble(
  Country = c("DK", "NL", "GB-SCT"),
  switch_year = c(2000L, 2000L, 2014L)
) %>%
  rowwise() %>%
  mutate(
    share_C_before = mean(dt_year$share_C[dt_year$Country == Country & dt_year$Year %in% (switch_year - 4L):(switch_year - 1L)]),
    share_C_after = mean(dt_year$share_C[dt_year$Country == Country & dt_year$Year %in% switch_year:(switch_year + 3L)])
  ) %>%
  ungroup()

step_table <- function(r, run) {
  bind_rows(lapply(seq_len(nrow(switches)), function(i) {
    cc <- switches$Country[[i]]; y <- switches$switch_year[[i]]
    yrs_b <- (y - 4L):(y - 1L); yrs_a <- y:min(y + 3L, 2015L)
    own_b <- r$residual[r$Country == cc & r$year %in% yrs_b]
    own_a <- r$residual[r$Country == cc & r$year %in% yrs_a]
    oth_b <- r$residual[r$Country != cc & r$year %in% yrs_b]
    oth_a <- r$residual[r$Country != cc & r$year %in% yrs_a]
    tibble(
      run = run, Country = cc, switch_year = y,
      years_before = sprintf("%d-%d", min(yrs_b), max(yrs_b)),
      years_after = sprintf("%d-%d", min(yrs_a), max(yrs_a)),
      n_before = length(own_b), n_after = length(own_a),
      mean_before = mean(own_b), mean_after = mean(own_a),
      step = mean(own_a) - mean(own_b),
      step_other_countries = mean(oth_a) - mean(oth_b),
      step_net = (mean(own_a) - mean(own_b)) - (mean(oth_a) - mean(oth_b))
    )
  }))
}
steps <- bind_rows(lapply(run_ids, function(r) step_table(rm_[[r]], r)))
write_tab(steps, "switch_steps")
steps_wide <- steps %>%
  select(Country, switch_year, run, n_before, n_after, step, step_net) %>%
  pivot_wider(names_from = run, values_from = c(n_before, n_after, step, step_net))

phase_v2_of <- function(y) cut(y, c(1984, 1991, 2001, 2007, 2015), labels = c("1985-1991", "1992-2001", "2002-2007", "2008-2015"))
cp <- bind_rows(lapply(run_ids, function(r) rm_[[r]] %>% mutate(run = r))) %>%
  mutate(period = phase_v2_of(year)) %>%
  group_by(run, Country, period) %>%
  summarise(n = n(), mean_residual = mean(residual), share_C = mean(DataType == "C", na.rm = TRUE), .groups = "drop")
write_tab(cp, "residual_by_country_phase_v2")
cp_wide <- cp %>%
  select(run, Country, period, n, mean_residual) %>%
  pivot_wider(names_from = run, values_from = c(n, mean_residual)) %>%
  left_join(cp %>% filter(run == "old") %>% select(Country, period, share_C_old = share_C), by = c("Country", "period")) %>%
  arrange(Country, period)

cy <- bind_rows(lapply(run_ids, function(r) rm_[[r]] %>% mutate(run = r))) %>%
  group_by(run, Country, year) %>%
  summarise(n = n(), mean_residual = mean(residual), .groups = "drop")
write_tab(cy, "residual_by_country_year")
cy_switch <- cy %>%
  inner_join(switches %>% select(Country, switch_year), by = "Country") %>%
  filter(year >= switch_year - 4L, year <= switch_year + 3L) %>%
  select(run, Country, year, n, mean_residual) %>%
  pivot_wider(names_from = run, values_from = c(n, mean_residual)) %>%
  left_join(dt_year %>% select(Country, year = Year, share_C), by = c("Country", "year")) %>%
  arrange(Country, year)

md <- c(md,
  "## Residual means by country × period (DataType C → R switches)",
  "",
  "Before the fix, DataType C hauls had ~2× N and E, so B_pred was too high and residuals too low while a",
  "country reported C. A C → R switch then shows as an upward step in that country's mean residual.",
  "`step` = mean residual in the 4 years from the switch minus the 4 years before (GB-SCT: 2014–2015 after);",
  "`net` subtracts the same-years change in all other countries (removes North-Sea-wide year effects).",
  "Each run on its own H1 haul set.",
  "",
  "Switch years confirmed from HH (share of valid hauls with DataType C, 4 years before → after):",
  "",
  md_table(switches %>% transmute(Country, switch_year, `share C before` = fmt(share_C_before, 2), `share C after` = fmt(share_C_after, 2))),
  "",
  md_table(steps_wide %>% transmute(
    Country, `switch year` = switch_year,
    `step old` = fmt(step_old), `step 30 min` = fmt(step_primary), `step 26–34` = fmt(step_sens),
    `net old` = fmt(step_net_old), `net 30 min` = fmt(step_net_primary), `net 26–34` = fmt(step_net_sens),
    `n before/after old` = sprintf("%d / %d", n_before_old, n_after_old),
    `n before/after 30 min` = sprintf("%d / %d", n_before_primary, n_after_primary),
    `n before/after 26–34` = sprintf("%d / %d", n_before_sens, n_after_sens)
  )),
  "",
  "### Year by year around each switch (H1 hauls)",
  "",
  md_table(cy_switch %>% transmute(
    Country, year, `share C` = fmt(share_C, 2),
    `n old` = na0(n_old), `n 30 min` = na0(n_primary), `n 26–34` = na0(n_sens),
    `resid old` = fmt(mean_residual_old), `resid 30 min` = fmt(mean_residual_primary), `resid 26–34` = fmt(mean_residual_sens)
  )),
  "",
  "### All countries × `phase_v2` period",
  "",
  md_table(cp_wide %>% transmute(
    Country, period, `share C (old hauls)` = fmt(share_C_old, 2),
    `n old` = na0(n_old), `n 30 min` = na0(n_primary), `n 26–34` = na0(n_sens),
    `resid old` = fmt(mean_residual_old), `resid 30 min` = fmt(mean_residual_primary), `resid 26–34` = fmt(mean_residual_sens)
  )),
  ""
)

# =================================================================================
# H2 / H3
# =================================================================================
pm <- per_run(function(r) read_rds_quiet(pth(r, "primary_model_v2.rds")))

panel_n <- function(p) {
  tibble(n_hauls = nrow(p$data), n_rectangles = n_distinct(p$data$stat_rec)) %>%
    bind_cols(as_tibble(as.list(setNames(p$phase_sizes$n_obs, paste0("n_", p$phase_sizes$phase_v2)))))
}
pn <- bind_rows(lapply(run_ids, function(r) panel_n(pm[[r]]) %>% mutate(run = r, .before = 1)))
write_tab(pn, "h2h3_panel_sizes")
rects <- lapply(pm, function(p) sort(unique(as.character(p$data$stat_rec))))
rects_same <- identical(rects$old, rects$primary) && identical(rects$old, rects$sens)

r2_row <- function(p) {
  car <- tryCatch(nakagawa_r2_spamm(p$fit_wb_car_v2, p$fit_wb_car_v2$data), error = function(e) NULL)
  re <- tryCatch({
    # nakagawa_r2_glmmtmb() without newdata needs the fit's original data object in scope
    fit <- p$primary_model_v2
    fe <- as.numeric(stats::predict(fit, newdata = p$data, re.form = NA))
    var_f <- stats::var(fe)
    var_r <- as.numeric(attr(glmmTMB::VarCorr(fit)$cond$stat_rec, "stddev")[1])^2
    var_e <- as.numeric(stats::sigma(fit))^2
    list(r2_marginal = var_f / (var_f + var_r + var_e), r2_conditional = (var_f + var_r) / (var_f + var_r + var_e))
  }, error = function(e) NULL)
  tibble(
    car_r2_marginal = if (is.null(car)) NA_real_ else car$r2_marginal,
    car_r2_conditional = if (is.null(car)) NA_real_ else car$r2_conditional,
    re_r2_marginal = if (is.null(re)) NA_real_ else re$r2_marginal,
    re_r2_conditional = if (is.null(re)) NA_real_ else re$r2_conditional
  )
}
r2 <- bind_rows(lapply(run_ids, function(r) r2_row(pm[[r]]) %>% mutate(run = r, .before = 1)))
write_tab(r2, "h2h3_r2")

fe_all <- function(file, model_id) {
  bind_rows(lapply(run_ids, function(r) {
    read.csv(pth(r, file), stringsAsFactors = FALSE) %>%
      transmute(model_id = model_id, run = r, term, estimate, se = std_error, p = p_value)
  }))
}
fe <- bind_rows(fe_all("phase_v2_car_fixed_effects.csv", "wb_car_v2"), fe_all("phase_v2_primary_fixed_effects.csv", "wb_primary_v2"))
write_tab(fe, "h2h3_fixed_effects")

slopes <- bind_rows(lapply(run_ids, function(r) {
  read.csv(pth(r, "phase_v2_fp_slopes_by_phase.csv"), stringsAsFactors = FALSE) %>%
    transmute(model_id, component, phase, run = r, slope = fp_slope, se = fp_slope_se, lo = fp_slope_lo, hi = fp_slope_hi,
              p = p_value, p_bonf = pmin(1, 4 * p_value))
}))
write_tab(slopes, "h2h3_slopes_by_phase")

slope_md <- function(mid) {
  w <- slopes %>%
    filter(model_id == mid) %>%
    mutate(cell = sprintf("%s (%s; %s)%s", fmt(slope, 4), fmt_p(p), fmt_p(p_bonf), ifelse(p_bonf < 0.05, " †", ""))) %>%
    select(component, phase, run, cell) %>%
    pivot_wider(names_from = run, values_from = cell)
  md_table(w %>%
    mutate(component = ifelse(component == "FP_between", "H2 FP_between", "H3 FP_within")) %>%
    rename(hypothesis = component) %>%
    rename_runs())
}
fe_md <- function(mid) {
  w <- fe %>%
    filter(model_id == mid, term != "(Intercept)") %>%
    mutate(cell = sprintf("%s (%s; p %s)", fmt(estimate, 4), fmt(se, 4), fmt_p(p))) %>%
    select(term, run, cell) %>%
    pivot_wider(names_from = run, values_from = cell)
  md_table(w %>% rename_runs())
}
sig_summary <- slopes %>%
  filter(model_id == "wb_car_v2") %>%
  group_by(component, phase) %>%
  summarise(
    sign = paste(ifelse(slope[match(run_ids, run)] > 0, "+", "−"), collapse = " / "),
    p_lt_05 = paste(ifelse(p[match(run_ids, run)] < 0.05, "yes", "no"), collapse = " / "),
    bonf_lt_05 = paste(ifelse(p_bonf[match(run_ids, run)] < 0.05, "yes", "no"), collapse = " / "),
    .groups = "drop"
  )

# Permutation null (CAR, FP_between): empirical two-sided p recomputed from the null draws
perm_p <- function(dir) {
  path <- file.path(dir, "permutation_bootstrap_FP_between_CAR_results.csv")
  if (!file.exists(path)) return(NULL)
  res <- read.csv(path, stringsAsFactors = FALSE) %>% filter(record_type == "permutation", converged %in% c(TRUE, "TRUE"))
  obs <- read.csv(file.path(dir, "phase_v2_fp_slopes_by_phase.csv"), stringsAsFactors = FALSE) %>%
    filter(model_id == "wb_car_v2", component == "FP_between")
  cols <- c("1985-1991" = "slope_1985_1991", "1992-2001" = "slope_1992_2001", "2002-2007" = "slope_2002_2007", "2008-2015" = "slope_2008_2015")
  bind_rows(lapply(names(cols), function(ph) {
    null <- res[[cols[[ph]]]]; o <- obs$fp_slope[obs$phase == ph]
    tibble(phase = ph, observed = o, null_lo = quantile(null, 0.025), null_hi = quantile(null, 0.975),
           p_perm = (1 + sum(abs(null) >= abs(o))) / (length(null) + 1), n_perm = length(null))
  }))
}
perm <- bind_rows(lapply(run_ids, function(r) { x <- perm_p(dirs[[r]]); if (is.null(x)) NULL else mutate(x, run = r, .before = 1) }))
if (nrow(perm)) write_tab(perm, "h2_permutation_car")

subs <- function(dir, h) {
  f <- file.path(dir, sprintf("h2h3_car_subsampling_summary_%s.csv", h))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE) %>% filter(abs(retention - 0.8) < 1e-9) %>%
    transmute(hypothesis = h, phase, pct_sign_match, pct_sign_and_sig_match)
}
sub_tab <- bind_rows(lapply(run_ids, function(r) {
  x <- bind_rows(subs(dirs[[r]], "H2"), subs(dirs[[r]], "H3"))
  if (nrow(x)) mutate(x, run = r) else NULL
}))
if (nrow(sub_tab)) write_tab(sub_tab, "h2h3_car_subsampling_r08")

speca <- function(dir) {
  f <- file.path(dir, "spec_a_car_identifiability.csv")
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) %>% select(any_of(c("gate", "rho_primary", "max_se_ratio_vs_fp_between", "n_warnings"))) else NULL
}
sa <- bind_rows(lapply(run_ids, function(r) { x <- speca(dirs[[r]]); if (is.null(x)) NULL else mutate(x, run = r, .before = 1) }))

md <- c(md,
  "## H2 / H3 — primary within-between model (`phase_v2`)",
  "",
  "`residual ~ FP_between * phase_v2 + FP_within * phase_v2`; CAR `adjacency(1|stat_rec)` (`wb_car_v2`, primary for H2 and H3)",
  "and companion RE `(1|stat_rec)` (`wb_primary_v2`). Each run's own haul panel.",
  "",
  md_table(pn %>% mutate(run = run_short[run]) %>% mutate(across(everything(), as.character))),
  "",
  sprintf("Same rectangle set in all three runs: %s (%s rectangles).",
          if (rects_same) "yes" else "**no**", paste(unique(lengths(rects)), collapse = " / ")),
  "",
  md_table(r2 %>% transmute(run = run_short[run], `CAR R²m` = fmt(car_r2_marginal), `CAR R²c` = fmt(car_r2_conditional),
                            `RE R²m` = fmt(re_r2_marginal), `RE R²c` = fmt(re_r2_conditional))),
  "",
  "Nakagawa R² as in `nakagawa_r2_spamm()` / `nakagawa_r2_glmmtmb()` (`pipeline/R/h2h3_results_helpers.R`).",
  "",
  "### Phase-specific slopes — CAR `wb_car_v2` (presented H2 and H3 effects)",
  "",
  "Each cell: slope (p; Bonferroni p). Bonferroni p = min(1, 4 × p) within each m = 4 family (H2, H3 separately);",
  "† = Bonferroni p < 0.05. SEs and 95% CIs are in `effort_fix_comparison_tables/h2h3_slopes_by_phase.csv`.",
  "",
  slope_md("wb_car_v2"),
  "",
  "Sign and significance, old / 30 min / 26–34:",
  "",
  md_table(sig_summary %>% transmute(hypothesis = ifelse(component == "FP_between", "H2", "H3"), phase, sign,
                                     `p < 0.05` = p_lt_05, `Bonferroni p < 0.05` = bonf_lt_05)),
  "",
  "### Fixed effects — CAR `wb_car_v2`",
  "",
  "Each cell: estimate (SE; p).",
  "",
  fe_md("wb_car_v2"),
  "",
  "### Phase-specific slopes — companion RE `wb_primary_v2`",
  "",
  slope_md("wb_primary_v2"),
  "",
  "### Fixed effects — companion RE `wb_primary_v2`",
  "",
  fe_md("wb_primary_v2"),
  ""
)
if (nrow(perm)) {
  md <- c(md,
    "### Robustness — FP_between spatial permutation (CAR, 1,000 shuffles)",
    "",
    "Each cell: two-sided p_perm, recomputed from the null draws.",
    "",
    md_table(perm %>% mutate(v = fmt_p(p_perm)) %>% select(phase, run, v) %>%
               pivot_wider(names_from = run, values_from = v) %>% rename_runs()),
    ""
  )
}
if (nrow(sub_tab)) {
  md <- c(md,
    "### Robustness — CAR rectangle subsampling (80% retention, 1,000 draws)",
    "",
    "Each cell: % of draws matching the full-data sign and significance (% matching sign only).",
    "",
    md_table(sub_tab %>%
      mutate(v = sprintf("%s (%s)", fmt(pct_sign_and_sig_match, 1), fmt(pct_sign_match, 1))) %>%
      select(hypothesis, phase, run, v) %>%
      pivot_wider(names_from = run, values_from = v) %>%
      rename_runs()),
    ""
  )
}
if (nrow(sa)) {
  md <- c(md,
    "### Robustness — Spec A (CAR + k-NN spatial lag) identifiability gate",
    "",
    md_table(sa %>% mutate(run = run_short[run]) %>% mutate(across(where(is.numeric), ~ fmt(.x, 3)))),
    ""
  )
}

# =================================================================================
# Run notes
# =================================================================================
md <- c(md, "## Run notes", "")
if (!is.null(path_timings) && file.exists(path_timings)) {
  tm <- read.delim(path_timings, stringsAsFactors = FALSE)
  md <- c(md,
    sprintf("Exactly-30-min live pipeline run (CLAUDE.md order): %d scripts, total %.1f min; exit codes: %s.",
            nrow(tm), sum(tm$seconds) / 60, paste(unique(tm$exit), collapse = ", ")),
    "",
    md_table(tm %>% transmute(script, exit, seconds)),
    ""
  )
}
md <- c(md,
  "- `run_h2h3_within_between.R` did not add the project `.R_libs` to its library path (unlike the other H2/H3 scripts), so",
  "  spaMM was not found and its original-phase CAR sensitivity was skipped in the pipeline runs above. The path was added and",
  "  the script re-run on its own at exactly 30 min, regenerating `outputs/h2h3_wb_car_fixed_effects.csv`; its log is otherwise",
  "  identical to the pipeline run (same 8,670 hauls; same RE fit read by `run_h2h3_phase_v2_refit.R`), so nothing downstream",
  "  changes. The primary `wb_car_v2` fit (in `run_h2h3_phase_v2_reporting.R`) was never affected. The 26–34 sensitivity",
  "  folder does not contain this original-phase CAR table.",
  ""
)

si_version <- function(dir, pkg) {
  files <- list.files(dir, pattern = "sessionInfo\\.txt$", full.names = TRUE)
  hits <- unlist(lapply(files, function(f) regmatches(readLines(f, warn = FALSE), gregexpr(paste0("\\b", pkg, "_[0-9][0-9.-]*"), readLines(f, warn = FALSE)))))
  v <- unique(sub(paste0("^", pkg, "_"), "", hits))
  if (length(v)) paste(sort(v), collapse = ", ") else "–"
}
si_r <- function(dir) {
  files <- list.files(dir, pattern = "sessionInfo\\.txt$", full.names = TRUE)
  v <- unique(unlist(lapply(files, function(f) grep("^R version", readLines(f, warn = FALSE), value = TRUE))))
  if (length(v)) paste(sub(" \\(.*", "", v), collapse = "; ") else "–"
}
lock_txt <- paste(readLines(file.path(project_root, "renv.lock"), warn = FALSE), collapse = "\n")
lock_version <- function(pkg) {
  m <- regmatches(lock_txt, regexpr(paste0('"Package": "', pkg, '",\\s*"Version": "[^"]+"'), lock_txt))
  if (length(m)) sub('.*"Version": "([^"]+)"', "\\1", m) else "not in lock"
}
lock_r <- sub('.*"R": \\{\\s*"Version": "([^"]+)".*', "\\1", lock_txt)
pkgs <- c("spaMM", "glmmTMB", "dplyr", "rlang", "Matrix", "TMB")
pk <- tibble(package = pkgs, renv_lock = vapply(pkgs, lock_version, ""))
for (r in run_ids) pk[[run_short[[r]]]] <- vapply(pkgs, function(p) si_version(dirs[[r]], p), "")
md <- c(md,
  sprintf("**Package versions differ from `renv.lock`.** `renv.lock` records R %s; old run: %s; 30 min: %s; 26–34: %s.",
          lock_r, si_r(dirs[["old"]]), si_r(dirs[["primary"]]), si_r(dirs[["sens"]])),
  "`renv::restore()` fails on R 4.6 (pinned versions such as rlang 1.1.6 do not compile), so the installed",
  "current versions were used and `renv.lock` was not modified. Versions below are read from the `*_sessionInfo.txt`",
  "files the H2/H3 scripts write (H1 scripts do not write one); several values mean different scripts saw different libraries.",
  "",
  md_table(pk),
  ""
)

writeLines(md, path_md)
cat("Saved:", path_md, "\nTables:", dir_tables, "\n")
