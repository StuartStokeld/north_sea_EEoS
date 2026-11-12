# Diagnosis: H1 hauls whose FishGlob B_obs omits species with NA weight.
#
# build_fishglob_haul_table() computes B_obs = sum(wgt, na.rm = TRUE) * 1000, so a species
# whose FishGlob `wgt` is NA contributes 0 g. This script characterises why `wgt` is NA,
# how much biomass is missing, whether it relates to the residual or to fishing pressure,
# and how the H1 metrics and primary CAR slopes change under four treatments.
# Read-only with respect to the pipeline: writes to audit/outputs/ only.
#
# Run from the repo root: Rscript audit/scripts/diagnose_missing_weights.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})
source("audit/scripts/audit_refit_helpers.R")
for (f in c("datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
            "datras_csv_import.R", "datras_hh_helpers.R")) {
  source(file.path("pipeline", "R", f))
}
out_dir <- "audit/outputs/missing_weights"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
haul <- readRDS("outputs/haul_eeos_predictions.rds")
haul$haul_key <- audit_haul_key(haul$haul_id)
stopifnot(!anyDuplicated(haul$haul_key))

env <- new.env()
load("FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData", envir = env)
fg <- env$data %>% filter(quarter == 1, year >= 1985, year <= 2015)
fg$haul_key <- audit_haul_key(fg$haul_id)
fg <- fg %>% filter(haul_key %in% haul$haul_key)

lw_fg <- read_csv("FishGlob_data/length_weight/length.weight_DATRAS_3August2023.csv",
                  show_col_types = FALSE) %>%
  distinct(taxa, .keep_all = TRUE) %>%
  select(taxa, a_fg = a, b_fg = b)
lw_ours <- read_csv("outputs/fishbase_lw_lookup_v2.csv", show_col_types = FALSE) %>%
  filter(is.finite(a_final), is.finite(b_final)) %>%
  select(aphia_id, a_ours = a_final, b_ours = b_final)

al <- load_aligned_hl(".")
hl_raw <- al$hl
haul_dt <- hl_raw %>% distinct(haul_key, DataType)
stopifnot(!anyDuplicated(haul_dt$haul_key))

# FishGlob rows are one per haul x accepted taxon (several DATRAS AphiaIDs can map to one taxon)
# FishGlob joins its LW table on taxa + family + genus + rank (get_datras.R l.773-812);
# flag rows where that full-key join finds no a/b although the taxon name is in the table.
lw_fg_full <- read_csv("FishGlob_data/length_weight/length.weight_DATRAS_3August2023.csv",
                       show_col_types = FALSE) %>%
  select(taxa, family, genus, rank, a, b)
fg_key <- fg %>%
  distinct(accepted_name, family, genus, rank) %>%
  left_join(lw_fg_full %>% distinct(taxa, family, genus, rank, .keep_all = TRUE),
            by = c("accepted_name" = "taxa", "family", "genus", "rank")) %>%
  transmute(accepted_name, family, genus, rank, lw_fullkey_ok = is.finite(a) & is.finite(b),
            lw_name_in_table = accepted_name %in% lw_fg_full$taxa)
say("accepted names with >1 family/genus/rank combination in FishGlob: %d",
    sum(duplicated(fg_key$accepted_name)))

fg_sp <- fg %>%
  left_join(fg_key, by = c("accepted_name", "family", "genus", "rank")) %>%
  group_by(haul_key, accepted_name) %>%
  summarise(
    wgt = sum(wgt), num = sum(num),
    lw_fullkey_ok = all(lw_fullkey_ok), lw_name_in_table = all(lw_name_in_table),
    aphia_ids = list(unique(verbatim_aphia_id)),
    country = first(country), year = first(year),
    .groups = "drop"
  ) %>%
  left_join(lw_fg, by = c("accepted_name" = "taxa"))

aphia_map <- fg %>% distinct(verbatim_aphia_id, accepted_name)
stopifnot(!anyDuplicated(aphia_map$verbatim_aphia_id))

# ---------------------------------------------------------------------------
# 1. Mechanism
# ---------------------------------------------------------------------------
# Our HL (April 2026 snapshot), effort-aligned, all rows (including NA lengths/counts)
hl_sp_raw <- hl_raw %>%
  filter(haul_key %in% haul$haul_key) %>%
  inner_join(aphia_map, by = c("Valid_Aphia" = "verbatim_aphia_id")) %>%
  group_by(haul_key, accepted_name) %>%
  summarise(
    hl_rows = n(),
    hl_rows_na_count = sum(is.na(HLNoAtLngt) | HLNoAtLngt <= 0),
    hl_rows_na_length = sum(is.na(LngtClass) | is.na(LngtCode)),
    .groups = "drop"
  )

na_rows <- fg_sp %>%
  filter(is.na(wgt)) %>%
  left_join(hl_sp_raw, by = c("haul_key", "accepted_name")) %>%
  left_join(haul_dt, by = "haul_key") %>%
  mutate(
    mechanism = case_when(
      !lw_name_in_table ~ "M1b taxon name absent from FishGlob LW table",
      !lw_fullkey_ok ~ "M1a FishGlob LW join fails on family/genus/rank (name present)",
      is.na(num) ~ "M2 no length counts in FishGlob (num NA)",
      TRUE ~ "M3 counts present but >= 1 length row NA (sum() without na.rm)"
    ),
    our_hl = case_when(
      is.na(hl_rows) ~ "species absent from our HL",
      hl_rows_na_count == hl_rows ~ "our HL: all counts NA",
      hl_rows_na_length > 0 ~ "our HL: some lengths NA",
      hl_rows_na_count > 0 ~ "our HL: some counts NA",
      TRUE ~ "our HL: complete lengths"
    )
  )

say("FishGlob haul x taxon rows in H1 hauls: %d; NA wgt: %d (%.2f%%); hauls affected: %d of %d",
    nrow(fg_sp), nrow(na_rows), 100 * nrow(na_rows) / nrow(fg_sp),
    n_distinct(na_rows$haul_key), nrow(haul))
mech <- na_rows %>% count(mechanism, our_hl, name = "haul_taxon_rows") %>% arrange(mechanism)
print(as.data.frame(mech))
write_csv(mech, file.path(out_dir, "mechanism.csv"))

by_species <- na_rows %>%
  group_by(accepted_name, mechanism) %>%
  summarise(haul_taxon_rows = n(), hauls = n_distinct(haul_key), .groups = "drop") %>%
  left_join(fg_sp %>% count(accepted_name, name = "hauls_present"), by = "accepted_name") %>%
  mutate(pct_of_hauls_present = round(100 * hauls / hauls_present, 1)) %>%
  arrange(desc(hauls))
write_csv(by_species, file.path(out_dir, "by_species.csv"))

haul_meta <- haul %>%
  select(haul_key, year) %>%
  mutate(country = sapply(strsplit(haul_key, "_"), `[`, 4)) %>%
  left_join(haul_dt, by = "haul_key") %>%
  mutate(affected = haul_key %in% na_rows$haul_key)
by_country <- haul_meta %>% group_by(country) %>%
  summarise(hauls = n(), pct = round(100 * mean(affected), 1), affected = sum(affected), .groups = "drop") %>%
  arrange(desc(affected))
by_year <- haul_meta %>% group_by(year) %>%
  summarise(hauls = n(), pct = round(100 * mean(affected), 1), affected = sum(affected), .groups = "drop")
by_datatype <- haul_meta %>% group_by(DataType) %>%
  summarise(hauls = n(), pct = round(100 * mean(affected), 1), affected = sum(affected), .groups = "drop")
by_country_year <- haul_meta %>% filter(affected) %>% count(country, year, name = "affected_hauls") %>%
  arrange(desc(affected_hauls))
write_csv(by_country, file.path(out_dir, "by_country.csv"))
write_csv(by_year, file.path(out_dir, "by_year.csv"))
write_csv(by_datatype, file.path(out_dir, "by_datatype.csv"))
write_csv(by_country_year, file.path(out_dir, "by_country_year.csv"))

# ---------------------------------------------------------------------------
# 2. Magnitude: estimated weight of every haul x taxon from our HL lengths
# ---------------------------------------------------------------------------
# FishGlob convention: Length = LngtClass/10 for "." and "0" (mm), LngtClass for "1"/"5" (cm),
# lower bound of the class (no midpoint), FishGlob LW table. Fallback: our lookup with the
# pipeline's length midpoints (lngt_midpoint_cm), for taxa FishGlob has no LW for.
bins <- clean_hl_raw(hl_raw, non_fish = integer(0)) %>%
  filter(haul_key %in% haul$haul_key) %>%
  inner_join(aphia_map, by = c("Valid_Aphia" = "verbatim_aphia_id")) %>%
  left_join(lw_fg, by = c("accepted_name" = "taxa")) %>%
  left_join(lw_ours, by = c("AphiaID" = "aphia_id")) %>%
  mutate(
    L_fg = ifelse(lngt_code %in% c(".", "0"), LngtClass * 0.1, LngtClass),
    w_fg_kg = NumberAtLength * a_fg * L_fg^b_fg / 1000,
    w_ours_kg = NumberAtLength * a_ours * length_cm^b_ours / 1000
  )
est_sp <- bins %>%
  group_by(haul_key, accepted_name) %>%
  summarise(
    n_ind = sum(NumberAtLength),
    w_fg_kg = if (all(is.finite(w_fg_kg))) sum(w_fg_kg) else NA_real_,
    w_ours_kg = if (all(is.finite(w_ours_kg))) sum(w_ours_kg) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    w_est_kg = coalesce(w_fg_kg, w_ours_kg),
    w_est_source = case_when(
      !is.na(w_fg_kg) ~ "FishGlob LW + FishGlob length convention",
      !is.na(w_ours_kg) ~ "pipeline LW lookup + pipeline midpoints",
      TRUE ~ "unresolved"
    )
  )

# Fallback for taxa without length rows: reported species-category weight from the raw
# DATRAS CSV (9 Nov 2026 snapshot; SpeciesCategoryWeight = CatCatchWgt, g; per hour for C).
csv_path <- list.files(".", pattern = "^Unaggregated trawl and biological information.*\\.csv$",
                       recursive = TRUE, full.names = TRUE)
csv_path <- csv_path[!grepl("^\\./(\\.venv|renv)/", csv_path)][1]
ccw <- NULL
if (!is.na(csv_path)) {
  blocks <- read_ices_csv_blocks(csv_path)
  hlc <- blocks$HL
  hhc <- blocks$HH
  hhc$haul_key <- paste("NS-IBTS", hhc$Year, hhc$Quarter, hhc$Country, hhc$Platform, hhc$HaulNumber, sep = "_")
  hlc$haul_key <- paste("NS-IBTS", hlc$Year, hlc$Quarter, hlc$Country, hlc$Platform, hlc$HaulNumber, sep = "_")
  hh_dur <- hhc %>% distinct(haul_key, DataType_hh = DataType, HaulDuration) %>%
    group_by(haul_key) %>% filter(n() == 1L) %>% ungroup()
  ccw <- hlc %>%
    filter(haul_key %in% haul$haul_key) %>%
    distinct(haul_key, AphiaID, SpeciesSex, SpeciesCategory, SpeciesCategoryWeight) %>%
    mutate(SpeciesCategoryWeight = ifelse(SpeciesCategoryWeight < 0, NA, SpeciesCategoryWeight)) %>%
    inner_join(aphia_map, by = c("AphiaID" = "verbatim_aphia_id")) %>%
    left_join(hh_dur, by = "haul_key") %>%
    mutate(w = ifelse(DataType_hh == "C", SpeciesCategoryWeight * HaulDuration / 60, SpeciesCategoryWeight)) %>%
    group_by(haul_key, accepted_name) %>%
    summarise(w_catch_kg = if (all(is.finite(w))) sum(w) / 1000 else NA_real_, .groups = "drop")
  say("Raw CSV used for CatCatchWgt fallback: %s", csv_path)
}

sp_all <- fg_sp %>%
  select(haul_key, accepted_name, wgt) %>%
  left_join(est_sp, by = c("haul_key", "accepted_name")) %>%
  left_join(if (is.null(ccw)) tibble(haul_key = character(), accepted_name = character(), w_catch_kg = numeric()) else ccw,
            by = c("haul_key", "accepted_name")) %>%
  mutate(
    missing = is.na(wgt),
    w_est_source = case_when(
      !is.na(w_est_kg) ~ w_est_source,
      !is.na(w_catch_kg) ~ "reported catch weight (raw CSV, Nov 2026)",
      !is.na(w_est_source) ~ "unresolved (lengths present, no LW)",
      TRUE ~ "unresolved (no lengths, no catch weight)"
    ),
    w_est_kg = coalesce(w_est_kg, w_catch_kg)
  )

fill_src <- sp_all %>% filter(missing) %>%
  count(w_est_source, name = "haul_taxon_rows")
print(as.data.frame(fill_src))
write_csv(fill_src, file.path(out_dir, "fill_source.csv"))

mag <- sp_all %>%
  group_by(haul_key) %>%
  summarise(
    n_taxa = n(),
    n_missing_taxa = sum(missing),
    n_ind_total = sum(n_ind, na.rm = TRUE),
    n_ind_missing = sum(n_ind[missing], na.rm = TRUE),
    B_obs_kg_reported = sum(wgt, na.rm = TRUE),
    w_missing_kg = sum(w_est_kg[missing], na.rm = TRUE),
    n_missing_unresolved = sum(missing & is.na(w_est_kg)),
    .groups = "drop"
  ) %>%
  filter(n_missing_taxa > 0) %>%
  mutate(
    share_ind_missing = n_ind_missing / n_ind_total,
    share_biomass_missing = w_missing_kg / (B_obs_kg_reported + w_missing_kg)
  )
stopifnot(nrow(mag) == n_distinct(na_rows$haul_key))
write_csv(mag, file.path(out_dir, "per_haul_magnitude.csv"))
qs <- c(0, .1, .25, .5, .75, .9, .95, .99, 1)
dist_tab <- tibble(
  quantile = paste0(100 * qs, "%"),
  share_individuals_pct = 100 * quantile(mag$share_ind_missing, qs, na.rm = TRUE),
  share_biomass_pct = 100 * quantile(mag$share_biomass_missing, qs, na.rm = TRUE)
)
print(as.data.frame(dist_tab), digits = 3)
write_csv(dist_tab, file.path(out_dir, "magnitude_quantiles.csv"))
bands <- mag %>%
  mutate(band = cut(100 * share_biomass_missing, c(-Inf, 0.1, 1, 5, 10, 25, 50, Inf),
                    labels = c("<=0.1%", "0.1-1%", "1-5%", "5-10%", "10-25%", "25-50%", ">50%"))) %>%
  count(band, name = "hauls")
print(as.data.frame(bands))
write_csv(bands, file.path(out_dir, "magnitude_bands.csv"))
say("hauls with >=1 missing taxon unresolved (no weight estimate): %d", sum(mag$n_missing_unresolved > 0))

# ---------------------------------------------------------------------------
# 3. How it shows up in the current results
# ---------------------------------------------------------------------------
haul2 <- haul %>%
  left_join(mag %>% select(haul_key, share_biomass_missing, share_ind_missing, w_missing_kg,
                           n_missing_unresolved), by = "haul_key") %>%
  mutate(affected = !is.na(share_biomass_missing))
cmp_h1 <- haul2 %>% group_by(affected) %>%
  summarise(
    hauls = n(),
    mean_signed = mean(residual), median_signed = median(residual),
    mean_abs = mean(abs_residual), median_abs = median(abs_residual),
    median_ln_B_obs = median(ln_B_obs), .groups = "drop"
  )
print(as.data.frame(cmp_h1), digits = 4)
write_csv(cmp_h1, file.path(out_dir, "affected_vs_unaffected_h1.csv"))
wt_signed <- wilcox.test(residual ~ affected, data = haul2)
wt_abs <- wilcox.test(abs_residual ~ affected, data = haul2)
say("Wilcoxon affected vs unaffected: signed p = %.3g; absolute p = %.3g", wt_signed$p.value, wt_abs$p.value)
aff <- haul2 %>% filter(affected)
ct_s <- cor.test(aff$share_biomass_missing, aff$residual, method = "spearman", exact = FALSE)
ct_a <- cor.test(aff$share_biomass_missing, aff$abs_residual, method = "spearman", exact = FALSE)
say("Among affected hauls: Spearman(missing biomass share, signed residual) = %.3f (p = %.3g); with absolute residual = %.3f (p = %.3g)",
    ct_s$estimate, ct_s$p.value, ct_a$estimate, ct_a$p.value)
assoc <- tibble(
  test = c("Wilcoxon signed residual, affected vs not", "Wilcoxon absolute residual, affected vs not",
           "Spearman share vs signed residual (affected)", "Spearman share vs absolute residual (affected)"),
  estimate = c(NA, NA, unname(ct_s$estimate), unname(ct_a$estimate)),
  p_value = c(wt_signed$p.value, wt_abs$p.value, ct_s$p.value, ct_a$p.value)
)
write_csv(assoc, file.path(out_dir, "residual_association.csv"))

dat_a <- audit_build_h2h3_data(audit_recompute_residuals(haul2))
dat_a$affected <- dat_a$haul_key %in% aff$haul_key
by_phase <- dat_a %>% group_by(phase_v2) %>%
  summarise(hauls = n(), pct = round(100 * mean(affected), 1), affected = sum(affected), .groups = "drop")
print(as.data.frame(by_phase))
write_csv(by_phase, file.path(out_dir, "h2h3_by_phase.csv"))
rect <- dat_a %>% group_by(stat_rec) %>%
  summarise(hauls = n(), pct = 100 * mean(affected), affected = sum(affected),
            FP_between = first(FP_between), .groups = "drop")
say("H2/H3: %d of %d hauls affected; %d of %d rectangles have >=1 affected haul; max rectangle share %.1f%%",
    sum(dat_a$affected), nrow(dat_a), sum(rect$affected > 0), nrow(rect), max(rect$pct))
ct_rect <- cor.test(rect$FP_between, rect$pct, method = "spearman", exact = FALSE)
say("Rectangle level: Spearman(FP_between, %% hauls affected) = %.3f (p = %.3g)", ct_rect$estimate, ct_rect$p.value)
write_csv(rect, file.path(out_dir, "h2h3_by_rectangle.csv"))
fp_cmp <- dat_a %>% group_by(affected) %>%
  summarise(hauls = n(), mean_FP_between = mean(FP_between), mean_FP_within = mean(FP_within),
            sd_FP_within = sd(FP_within), .groups = "drop")
print(as.data.frame(fp_cmp), digits = 4)
write_csv(fp_cmp, file.path(out_dir, "h2h3_fp_affected.csv"))
glm_aff <- glm(affected ~ FP_between * phase_v2 + FP_within * phase_v2, family = binomial, data = dat_a)
glm_tab <- as.data.frame(summary(glm_aff)$coefficients)
glm_tab$term <- rownames(glm_tab)
print(glm_tab[, c("term", "Estimate", "Std. Error", "Pr(>|z|)")], digits = 3)
write_csv(glm_tab, file.path(out_dir, "h2h3_affected_glm.csv"))

# ---------------------------------------------------------------------------
# 4. Treatments
# ---------------------------------------------------------------------------
share <- setNames(haul2$share_biomass_missing, haul2$haul_key)
unres <- haul2$haul_key[!is.na(haul2$n_missing_unresolved) & haul2$n_missing_unresolved > 0]
# Thresholds use the estimated share (a lower bound where a missing taxon is unresolved).
excl_above <- function(thr) {
  haul2 %>% filter(!(affected & share_biomass_missing > thr))
}
treat <- list(
  "a current (NA weight counted as 0)" = haul2,
  "b1 exclude hauls missing >1% biomass" = excl_above(0.01),
  "b5 exclude hauls missing >5% biomass" = excl_above(0.05),
  "b10 exclude hauls missing >10% biomass" = excl_above(0.10),
  "c exclude all affected hauls" = haul2 %>% filter(!affected),
  "d fill missing taxa from DATRAS lengths x LW" = haul2 %>%
    mutate(B_obs = B_obs + 1000 * coalesce(w_missing_kg, 0))
)
res <- lapply(names(treat), function(nm) audit_run_treatment(treat[[nm]], nm))
h1_tab <- bind_rows(lapply(res, `[[`, "h1"))
slopes <- bind_rows(lapply(res, `[[`, "slopes"))
print(as.data.frame(h1_tab), digits = 4)
write_csv(h1_tab, file.path(out_dir, "treatments_h1.csv"))
write_csv(slopes, file.path(out_dir, "treatments_car_slopes.csv"))
say("Hauls whose missing taxa include an unresolved weight (share is a lower bound): %d", length(unres))

writeLines(c(
  "## H1 metrics by treatment", "", audit_md_table(h1_tab, 4), "",
  "## CAR slopes by treatment (slope (SE); * raw p<0.05; ** Bonferroni p<0.05)", "",
  audit_md_table(
    audit_format_slopes(slopes) %>%
      select(treatment, response, component, phase, cell) %>%
      pivot_wider(names_from = phase, values_from = cell) %>%
      left_join(slopes %>% distinct(treatment, response, n_hauls_h2h3, n_rect), by = c("treatment", "response"))
  )
), file.path(out_dir, "treatments_tables.md"))
say("Done. Outputs in %s", out_dir)
