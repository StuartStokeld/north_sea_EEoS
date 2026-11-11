# GATE 1 validation of the DATRAS haul-effort alignment (read-only w.r.t. the live outputs).
#
# Rebuilds the effort-aligned DATRAS state IN MEMORY and checks it against
# expectations, without overwriting any live pipeline output. Writes only to
# outputs/datras_effort_validation/.
#
# Run from repo root (after pipeline/import_datras_hh.R):
#   Rscript pipeline/validate_haul_effort_alignment.R
#
# Reports:
#   1. HL -> HH match rate
#   2. hauls kept / dropped by reason x DataType
#   3. kept-haul duration range; GB-SCT 1985-1998 long tows kept (expect 0)
#   4. median DATRAS-implied biomass / FishGlob B_obs by DataType x duration band,
#      before and after converting C counts per hour -> per haul
#   5. share of non-integer raw counts after conversion
#   6. haul counts in H1 terms (new in-memory panel; previous H1 / H2-H3 haul sets kept)
#   7. HL-import fix: CSV HL block (new reader) vs tracked datras_hl_raw.rds

suppressPackageStartupMessages({
  library(dplyr)
})

cli_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cli_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run from pipeline/ or Rscript pipeline/validate_haul_effort_alignment.R")
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "h1_common.R"))
project_root <- get_project_root_from(script_dir)
source(file.path(r_dir, "datras_constants.R"))
source(file.path(r_dir, "datras_hl_helpers.R"))
source(file.path(r_dir, "datras_state_helpers.R"))
source(file.path(r_dir, "datras_csv_import.R"))
source(file.path(r_dir, "datras_hh_helpers.R"))
source(file.path(r_dir, "h1_join_helpers.R"))

out_dir <- file.path(project_root, "outputs", "datras_effort_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_hl <- file.path(project_root, "outputs", "datras_hl_raw.rds")
path_hh <- file.path(project_root, "outputs", "datras_hh_raw.rds")
path_lw <- file.path(project_root, "outputs", "fishbase_lw_lookup_v2.csv")
path_fishglob <- file.path(project_root, "FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
# "Previous" = pre-fix run, archived in outputs/archive_pre_effort_fix/ (gitignored);
# falls back to the live outputs if the archive is absent.
old_dir <- file.path(project_root, "outputs", "archive_pre_effort_fix")
if (!dir.exists(old_dir)) old_dir <- file.path(project_root, "outputs")
path_old_preds <- file.path(old_dir, "haul_eeos_predictions.rds")
path_old_primary <- file.path(old_dir, "primary_model_v2.rds")
stopifnot(file.exists(path_hl), file.exists(path_hh), file.exists(path_lw), file.exists(path_fishglob))

hl <- readRDS(path_hl)
hh <- readRDS(path_hh)
lw <- read.csv(path_lw, stringsAsFactors = FALSE)
window <- haul_dur_window()
# Band label for the active window: "30 min" (primary) or e.g. "26-34 min" (sensitivity)
window_label <- if (window[["min"]] == window[["max"]]) sprintf("%g min", window[["min"]]) else sprintf("%g-%g min", window[["min"]], window[["max"]])

# --- 1-3, 5: primary alignment ---------------------------------------------------
aligned <- align_hl_to_hh(hl, hh, dur_min = window[["min"]], dur_max = window[["max"]])
kept <- aligned$hauls %>% filter(effort_reason == "kept")

dur_check <- tibble(
  check = c("kept_HaulDur_min", "kept_HaulDur_max", "kept_outside_window",
            "kept_GB_SCT_1985_1998_total", "kept_GB_SCT_1985_1998_ge_55min",
            "kept_DataTypes"),
  value = c(
    min(kept$HaulDur), max(kept$HaulDur),
    sum(kept$HaulDur < window[["min"]] | kept$HaulDur > window[["max"]]),
    sum(kept$Country == "GB-SCT" & kept$Year <= 1998L),
    sum(kept$Country == "GB-SCT" & kept$Year <= 1998L & kept$HaulDur >= 55),
    paste(sort(unique(kept$DataType)), collapse = "/")
  )
)
sct_long <- aligned$hauls %>%
  filter(Country == "GB-SCT", Year <= 1998L) %>%
  mutate(dur_band = if_else(HaulDur >= 55, ">=55 min",
    if_else(HaulDur >= window[["min"]] & HaulDur <= window[["max"]], window_label, "other"))) %>%
  count(dur_band, effort_reason, name = "n_hauls")

# --- 4: biomass ratio by DataType x duration band ---------------------------------
fishglob <- build_fishglob_haul_table(path_fishglob)

# All valid hauls with any DataType/duration, conversion applied, to compare bands
wide <- align_hl_to_hh(hl, hh, dur_min = 0, dur_max = Inf, datatypes = c("C", "R", "P", "S"))
bins_wide <- add_lw_mass_to_bins(clean_hl_raw(wide$hl), lw)
b_wide <- bins_wide %>%
  group_by(haul_key) %>%
  summarise(B_datras_raw = sum(NumberAtLength * mass_g), .groups = "drop") %>%
  inner_join(wide$hauls %>% select(haul_key, Country, Year, DataType, HaulDur, effort_f), by = "haul_key") %>%
  mutate(B_datras_reported = B_datras_raw / effort_f) %>%
  inner_join(fishglob %>% select(join_key, B_obs), by = c("haul_key" = "join_key")) %>%
  filter(B_obs > 0) %>%
  mutate(
    ratio_reported = B_datras_reported / B_obs,
    ratio_converted = B_datras_raw / B_obs,
    dur_band = case_when(
      HaulDur >= window[["min"]] & HaulDur <= window[["max"]] ~ window_label,
      HaulDur >= 55 & HaulDur <= 65 ~ "55-65 min",
      TRUE ~ "other"
    ),
    in_window = HaulDur >= window[["min"]] & HaulDur <= window[["max"]]
  )

ratio_table <- bind_rows(
  b_wide %>% group_by(DataType, dur_band),
  b_wide %>% filter(in_window) %>% mutate(dur_band = "ALL KEPT (in window)") %>% group_by(DataType, dur_band)
) %>%
  summarise(
    n_hauls = n(),
    median_ratio_before_conversion = round(median(ratio_reported), 3),
    median_ratio_after_conversion = round(median(ratio_converted), 3),
    iqr_after_lo = round(quantile(ratio_converted, 0.25), 3),
    iqr_after_hi = round(quantile(ratio_converted, 0.75), 3),
    .groups = "drop"
  ) %>%
  arrange(dur_band != "ALL KEPT (in window)", DataType, dur_band)

# --- 6: haul counts in H1 terms ---------------------------------------------------
built <- build_state_from_hl(aligned$hl, lw)
state <- built$haul_state %>%
  left_join(aligned$hauls %>% select(haul_key, DataType, HaulDur), by = "haul_key")
state_valid <- normalize_haul_E(filter_valid_haul_state(state))
joined <- join_fishglob_datras(fishglob, prepare_datras_for_join(state_valid))
h1_new <- filter_eeos_inputs(joined)

fg_key <- function(id) vapply(strsplit(id, " "), function(p) paste(p[1], p[2], p[3], p[4], p[5], p[8], sep = "_"), "")
old_h1 <- if (file.exists(path_old_preds)) readRDS(path_old_preds) %>% mutate(haul_key = fg_key(haul_id)) else NULL
old_h23 <- if (file.exists(path_old_primary)) suppressWarnings(readRDS(path_old_primary))$data %>% mutate(haul_key = fg_key(haul_id)) else NULL

reason_of <- function(keys) {
  tibble(haul_key = keys) %>%
    left_join(aligned$hauls %>% select(haul_key, effort_reason, DataType), by = "haul_key") %>%
    mutate(effort_reason = coalesce(effort_reason, "not_in_HL"), DataType = coalesce(DataType, "(none)"))
}

counts <- tibble(
  quantity = c(
    "HL hauls (Q1 1985-2015)",
    "effort-kept HL hauls",
    "DATRAS state hauls (valid LW/E), new",
    "FishGlob-joined hauls, new",
    "H1 hauls passing EEoS input filters, new (pre-EEoS solve)",
    "previous H1 hauls (haul_eeos_predictions.rds)",
    "previous H1 hauls kept by effort alignment",
    "previous H2/H3 hauls (primary_model_v2.rds$data)",
    "previous H2/H3 hauls kept by effort alignment"
  ),
  n = c(
    nrow(aligned$hauls), nrow(kept), nrow(state_valid), nrow(joined), nrow(h1_new),
    if (is.null(old_h1)) NA else nrow(old_h1),
    if (is.null(old_h1)) NA else sum(old_h1$haul_key %in% kept$haul_key),
    if (is.null(old_h23)) NA else nrow(old_h23),
    if (is.null(old_h23)) NA else sum(old_h23$haul_key %in% kept$haul_key)
  ),
  expected = c(NA, NA, NA, NA, NA, 12069, 10941, 10464, 9402)
) %>%
  mutate(pct_of_previous = case_when(
    grepl("previous H1 hauls kept", quantity) ~ round(100 * n / n[quantity == "previous H1 hauls (haul_eeos_predictions.rds)"], 1),
    grepl("previous H2/H3 hauls kept", quantity) ~ round(100 * n / n[quantity == "previous H2/H3 hauls (primary_model_v2.rds$data)"], 1),
    grepl("H1 hauls passing", quantity) ~ round(100 * n / n[quantity == "previous H1 hauls (haul_eeos_predictions.rds)"], 1),
    TRUE ~ NA_real_
  ))

old_h1_reasons <- if (is.null(old_h1)) NULL else reason_of(old_h1$haul_key) %>% count(effort_reason, DataType, name = "n_hauls")
old_h23_reasons <- if (is.null(old_h23)) NULL else reason_of(old_h23$haul_key) %>% count(effort_reason, DataType, name = "n_hauls")
new_h1_by_type <- h1_new %>% count(DataType, name = "n_hauls")

# Hauls newly lost at the LW/EEoS-filter stage because counts shrank (C x 0.5)
lost_after_kept <- if (is.null(old_h1)) NULL else {
  old_kept <- old_h1$haul_key[old_h1$haul_key %in% kept$haul_key]
  lost <- setdiff(old_kept, h1_new$join_key)
  joined_reason <- joined %>%
    filter(join_key %in% lost) %>%
    mutate(exclusion_reason = classify_eeos_exclusion(.)) %>%
    count(exclusion_reason, name = "n_hauls")
  bind_rows(
    joined_reason,
    tibble(exclusion_reason = "not joined / invalid state", n_hauls = length(lost) - sum(joined_reason$n_hauls))
  )
}

# --- 7: HL-import fix check -------------------------------------------------------
csv_path <- find_ices_hl_csv(project_root)
hl_import_check <- if (is.null(csv_path)) {
  tibble(check = "csv_present", value = "no ICES CSV on disk — skipped")
} else {
  hl_csv <- import_hl_from_ices_csv(csv_path)
  key6 <- function(d) make_haul_key(d$Survey, d$Year, d$Quarter, d$Country, d$Ship, d$HaulNo)
  k_csv <- unique(key6(hl_csv)); k_rds <- unique(key6(hl))
  row_id <- function(d) paste(key6(d), d$Valid_Aphia, d$LngtCode, d$LngtClass, d$HLNoAtLngt, d$SubFactor)
  r_csv <- row_id(hl_csv); r_rds <- row_id(hl)
  tibble(
    check = c("rows_csv", "rows_rds", "hauls_csv", "hauls_rds", "hauls_in_both",
              "rds_rows_found_in_csv_pct", "sum_raised_counts_csv", "sum_raised_counts_rds"),
    value = as.character(c(
      nrow(hl_csv), nrow(hl), length(k_csv), length(k_rds), length(intersect(k_csv, k_rds)),
      round(100 * mean(r_rds %in% r_csv), 3),
      sum(hl_csv$HLNoAtLngt * hl_csv$SubFactor), sum(hl$HLNoAtLngt * hl$SubFactor)
    ))
  )
}

# --- Write ------------------------------------------------------------------------
write.csv(aligned$summary, file.path(out_dir, "alignment_summary.csv"), row.names = FALSE)
write.csv(aligned$raw_count_check, file.path(out_dir, "raw_count_check.csv"), row.names = FALSE)
write.csv(dur_check, file.path(out_dir, "duration_checks.csv"), row.names = FALSE)
write.csv(sct_long, file.path(out_dir, "gb_sct_1985_1998.csv"), row.names = FALSE)
write.csv(ratio_table, file.path(out_dir, "biomass_ratio_by_datatype.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "haul_counts.csv"), row.names = FALSE)
if (!is.null(old_h1_reasons)) write.csv(old_h1_reasons, file.path(out_dir, "previous_h1_hauls_by_reason.csv"), row.names = FALSE)
if (!is.null(old_h23_reasons)) write.csv(old_h23_reasons, file.path(out_dir, "previous_h2h3_hauls_by_reason.csv"), row.names = FALSE)
if (!is.null(lost_after_kept)) write.csv(lost_after_kept, file.path(out_dir, "previous_h1_kept_but_now_filtered.csv"), row.names = FALSE)
write.csv(hl_import_check, file.path(out_dir, "hl_import_fix_check.csv"), row.names = FALSE)

show <- function(title, df) {
  cat("\n===", title, "===\n")
  print(as.data.frame(df), row.names = FALSE)
}
cat(sprintf("\nWindow: %g-%g min | HL->HH match rate: %.3f%%\n", window[["min"]], window[["max"]], 100 * aligned$match_rate))
show("Hauls by effort reason x DataType (all HL hauls)", aligned$summary)
show("Duration checks (kept hauls)", dur_check)
show("GB-SCT 1985-1998 by duration band", sct_long)
show("Median DATRAS-implied biomass / B_obs", ratio_table)
show("Raw-count integer check (kept hauls)", aligned$raw_count_check)
show("Haul counts", counts)
if (!is.null(old_h1_reasons)) show("Previous H1 hauls by effort reason", old_h1_reasons)
if (!is.null(old_h23_reasons)) show("Previous H2/H3 hauls by effort reason", old_h23_reasons)
if (!is.null(lost_after_kept)) show("Previous H1 hauls effort-kept but now failing state/EEoS filters", lost_after_kept)
show("New H1 panel by DataType", new_h1_by_type)
show("HL-import fix: CSV (new block reader) vs tracked datras_hl_raw.rds", hl_import_check)
cat("\nWritten to:", out_dir, "\n")
