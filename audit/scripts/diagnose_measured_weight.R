# Diagnosis: DATRAS measured catch weight (SpeciesCategoryWeight = CatCatchWgt, g) as an alternative B_obs.
# Source: the raw ICES CSV (HH + HL blocks, downloaded 9 Nov 2026; local only, not in git).
# Read-only; writes audit/outputs/measured_weight/. Run from the repo root.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})
source("audit/scripts/audit_refit_helpers.R")
for (f in c("datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
            "datras_csv_import.R", "datras_hh_helpers.R")) source(file.path("pipeline", "R", f))
out_dir <- "audit/outputs/measured_weight"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

haul <- readRDS("outputs/haul_eeos_predictions.rds")
haul$haul_key <- audit_haul_key(haul$haul_id)
lw <- read_csv("outputs/fishbase_lw_lookup_v2.csv", show_col_types = FALSE)
pool_ids <- lw$aphia_id[is.finite(lw$a_final) & is.finite(lw$b_final)]

csv <- list.files(".", pattern = "^Unaggregated trawl and biological information.*\\.csv$",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(csv))
b <- read_ices_csv_blocks(csv)
key <- function(d) paste("NS-IBTS", d$Year, d$Quarter, d$Country, d$Platform, d$HaulNumber, sep = "_")
hh <- b$HH %>% mutate(haul_key = key(.)) %>% filter(haul_key %in% haul$haul_key) %>%
  distinct(haul_key, DataType, HaulDuration, Country, Year) %>%
  group_by(haul_key) %>% filter(n() == 1L) %>% ungroup()
hl <- b$HL %>% mutate(haul_key = key(.)) %>% filter(haul_key %in% hh$haul_key)

# 1. How weights are recorded: one SpeciesCategoryWeight per haul x species x sex x category
#    (repeated on every length row). Check it is constant within that key.
cat_key <- c("haul_key", "AphiaID", "SpeciesSex", "SpeciesCategory")
w_var <- hl %>% group_by(across(all_of(cat_key))) %>%
  summarise(n_w = n_distinct(SpeciesCategoryWeight), .groups = "drop")
say("haul x species x sex x category keys: %d; with >1 distinct weight: %d", nrow(w_var), sum(w_var$n_w > 1))
cats <- hl %>% distinct(across(all_of(cat_key)), SpeciesCategoryWeight) %>%
  group_by(across(all_of(cat_key))) %>% slice(1) %>% ungroup() %>%
  left_join(hh, by = "haul_key") %>%
  mutate(valid = !is.na(SpeciesCategoryWeight) & SpeciesCategoryWeight >= 0,
         w_g = ifelse(valid, ifelse(DataType == "C", SpeciesCategoryWeight * HaulDuration / 60,
                                    SpeciesCategoryWeight), NA_real_),
         in_pool = AphiaID %in% pool_ids,
         non_fish = AphiaID %in% NON_FISH_APHIAIDS)
say("category records: %d; valid weight: %.1f%%; zero weight: %d", nrow(cats), 100 * mean(cats$valid),
    sum(cats$valid & cats$SpeciesCategoryWeight == 0))
p_sub <- cats %>% filter(DataType == "P") %>% count(haul_key, AphiaID, SpeciesSex) %>% filter(n > 1)
say("P hauls: species x sex with >1 category (subcategories, summed): %d", nrow(p_sub))

sp <- cats %>% group_by(haul_key, AphiaID, in_pool, non_fish) %>%
  summarise(valid = all(valid), w_g = if (all(valid)) sum(w_g) else NA_real_, .groups = "drop")
say("species x haul records: %d; valid: %.1f%%; in S/N/E pool: %d (valid %.1f%%)", nrow(sp),
    100 * mean(sp$valid), sum(sp$in_pool), 100 * mean(sp$valid[sp$in_pool]))

hm <- sp %>% group_by(haul_key) %>%
  summarise(n_pool = sum(in_pool), n_pool_valid = sum(in_pool & valid),
            complete_pool = all(valid[in_pool]),
            B_meas_pool_g = sum(w_g[in_pool], na.rm = TRUE),
            w_excluded_g = sum(w_g[!in_pool], na.rm = TRUE),
            w_total_g = sum(w_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(share_excluded = w_excluded_g / w_total_g) %>%
  left_join(hh, by = "haul_key")
haul_m <- haul %>% select(haul_key, B_obs, B_pred, year, stat_rec, S) %>% left_join(hm, by = "haul_key")
haul_m$complete_pool[is.na(haul_m$complete_pool)] <- FALSE
say("H1 hauls: %d; in Nov-2026 CSV: %d; complete pool weights: %d (%.1f%%)", nrow(haul),
    sum(!is.na(haul_m$n_pool)), sum(haul_m$complete_pool), 100 * mean(haul_m$complete_pool))

cov_dt <- haul_m %>% group_by(DataType) %>%
  summarise(hauls = n(), pct_complete = round(100 * mean(complete_pool), 1),
            pct_species_records_valid = round(100 * sum(n_pool_valid, na.rm = TRUE) / sum(n_pool, na.rm = TRUE), 1),
            .groups = "drop")
cov_cy <- haul_m %>% group_by(Country, Year) %>%
  summarise(hauls = n(), pct_complete = round(100 * mean(complete_pool), 1), .groups = "drop")
cov_c <- haul_m %>% group_by(Country) %>%
  summarise(hauls = n(), pct_complete = round(100 * mean(complete_pool), 1), .groups = "drop")
cov_y <- haul_m %>% group_by(Year) %>%
  summarise(hauls = n(), pct_complete = round(100 * mean(complete_pool), 1), .groups = "drop")
print(as.data.frame(cov_dt)); print(as.data.frame(cov_c))
write_csv(cov_dt, file.path(out_dir, "coverage_by_datatype.csv"))
write_csv(cov_cy, file.path(out_dir, "coverage_by_country_year.csv"))
write_csv(cov_c, file.path(out_dir, "coverage_by_country.csv"))
write_csv(cov_y, file.path(out_dir, "coverage_by_year.csv"))
say("excluded-species weight share (complete hauls): median %.2f%%, 90%% %.2f%%, max %.1f%%",
    100 * median(haul_m$share_excluded[haul_m$complete_pool]),
    100 * quantile(haul_m$share_excluded[haul_m$complete_pool], 0.9),
    100 * max(haul_m$share_excluded[haul_m$complete_pool]))

# 3. Agreement with FishGlob B_obs
cmp <- haul_m %>% filter(complete_pool, B_meas_pool_g > 0) %>% mutate(ratio = B_meas_pool_g / B_obs)
q <- quantile(cmp$ratio, c(.01, .05, .25, .5, .75, .95, .99))
print(round(q, 3))
ag <- tibble(stat = names(q), ratio_measured_over_fishglob = unname(q))
write_csv(ag, file.path(out_dir, "agreement_quantiles.csv"))
out <- cmp %>% mutate(outlier = abs(log(ratio)) > log(2))
say("complete hauls compared: %d; within 10%%: %.1f%%; outside x2: %d", nrow(cmp),
    100 * mean(abs(cmp$ratio - 1) < 0.1), sum(out$outlier))
ol <- out %>% group_by(Country) %>% summarise(hauls = n(), outliers = sum(outlier),
                                              median_ratio = median(ratio), .groups = "drop")
ol_y <- out %>% group_by(Country, Year) %>% summarise(hauls = n(), outliers = sum(outlier),
                                                      median_ratio = median(ratio), .groups = "drop") %>%
  filter(outliers > 0) %>% arrange(desc(outliers))
print(as.data.frame(ol)); print(head(as.data.frame(ol_y), 12))
write_csv(ol, file.path(out_dir, "agreement_by_country.csv"))
write_csv(ol_y, file.path(out_dir, "outliers_by_country_year.csv"))
write_csv(out %>% select(haul_key, Country, Year, DataType, B_obs, B_meas_pool_g, ratio, share_excluded),
          file.path(out_dir, "per_haul_comparison.csv"))

# 4. H1 + primary CAR on complete hauls: FishGlob B_obs vs measured B_obs (same hauls)
sub <- haul %>% filter(haul_key %in% cmp$haul_key)
meas <- sub %>% mutate(B_obs = cmp$B_meas_pool_g[match(haul_key, cmp$haul_key)])
ra <- audit_run_treatment(sub, "FishGlob B_obs, complete-weight hauls")
rb <- audit_run_treatment(meas, "measured B_obs, complete-weight hauls")
h1 <- bind_rows(audit_h1_metrics(audit_recompute_residuals(haul)) %>% mutate(treatment = "FishGlob B_obs, all H1 hauls", .before = 1),
                ra$h1, rb$h1)
sl <- bind_rows(ra$slopes, rb$slopes)
print(as.data.frame(h1), digits = 4)
write_csv(h1, file.path(out_dir, "h1_metrics.csv"))
write_csv(sl, file.path(out_dir, "car_slopes.csv"))
writeLines(c("## H1", "", audit_md_table(h1, 4), "", "## CAR slopes", "",
             audit_md_table(audit_format_slopes(sl) %>% select(treatment, response, component, phase, cell, p_raw, p_bonf, n_hauls_h2h3, n_rect))),
           file.path(out_dir, "tables.md"))

# 5. The factor-100 cluster: classify country-year blocks, test a block-level x100 unit correction.
cy <- cmp %>% mutate(low = ratio < 0.05) %>% group_by(Country, Year) %>%
  summarise(hauls = n(), share_low = mean(low), median_ratio = median(ratio), .groups = "drop")
say("country-years compared: %d; all-low: %d; all-normal: %d; mixed: %d", nrow(cy),
    sum(cy$share_low == 1), sum(cy$share_low == 0), sum(cy$share_low > 0 & cy$share_low < 1))
print(as.data.frame(cy %>% filter(share_low > 0 & share_low < 1)))
write_csv(cy, file.path(out_dir, "factor100_by_country_year.csv"))
blocks <- cy %>% filter(share_low > 0.5) %>% transmute(Country, Year, x100 = TRUE)
print(as.data.frame(blocks %>% group_by(Country) %>% summarise(years = paste(range(Year), collapse = "-"), n_years = n())))
cmp2 <- cmp %>% left_join(blocks, by = c("Country", "Year")) %>%
  mutate(B_meas_corr = ifelse(coalesce(x100, FALSE), 100 * B_meas_pool_g, B_meas_pool_g),
         ratio_corr = B_meas_corr / B_obs)
q2 <- quantile(cmp2$ratio_corr, c(.01, .05, .25, .5, .75, .95, .99)); print(round(q2, 3))
say("after block x100: within 10%%: %.1f%%; within 25%%: %.1f%%; outside x2: %d", 100 * mean(abs(cmp2$ratio_corr - 1) < 0.1),
    100 * mean(abs(cmp2$ratio_corr - 1) < 0.25), sum(abs(log(cmp2$ratio_corr)) > log(2)))
write_csv(tibble(stat = names(q2), ratio = unname(q2)), file.path(out_dir, "agreement_quantiles_x100.csv"))
write_csv(cmp2 %>% group_by(Country) %>% summarise(hauls = n(), median_ratio = median(ratio_corr),
          outside_x2 = sum(abs(log(ratio_corr)) > log(2)), .groups = "drop"), file.path(out_dir, "agreement_by_country_x100.csv"))
meas2 <- sub %>% mutate(B_obs = cmp2$B_meas_corr[match(haul_key, cmp2$haul_key)])
rc <- audit_run_treatment(meas2, "measured B_obs x100 blocks, complete-weight hauls")
h1b <- bind_rows(h1, rc$h1); slb <- bind_rows(sl, rc$slopes)
print(as.data.frame(h1b), digits = 4)
write_csv(h1b, file.path(out_dir, "h1_metrics.csv")); write_csv(slb, file.path(out_dir, "car_slopes.csv"))
writeLines(c("## H1", "", audit_md_table(h1b, 4), "", "## CAR slopes", "",
             audit_md_table(audit_format_slopes(slb) %>% select(treatment, response, component, phase, cell, p_raw, p_bonf, n_hauls_h2h3, n_rect))),
           file.path(out_dir, "tables.md"))
