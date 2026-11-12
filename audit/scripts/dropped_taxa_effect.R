# F-B8: taxa in the pipeline's HL rows that the LW lookup drops. Lists them, compares the rebuilt
# lookup (audit/scripts/rebuild_lw_lookup.R, FishBase 26.06) with the current one, and estimates the
# effect on S, N, E and H1. Read-only; writes audit/outputs/dropped_taxa/.
# Run from the repo root after rebuild_lw_lookup.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("audit/scripts/audit_refit_helpers.R")
for (f in c("datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
            "datras_csv_import.R", "datras_hh_helpers.R")) source(file.path("pipeline", "R", f))
out_dir <- "audit/outputs/dropped_taxa"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

cur <- read_csv("outputs/fishbase_lw_lookup_v2.csv", show_col_types = FALSE)
new <- read_csv("audit/outputs/lw_rebuild/fishbase_lw_lookup_rebuilt_26.06.csv", show_col_types = FALSE)
haul <- readRDS("outputs/haul_eeos_predictions.rds")
haul$haul_key <- audit_haul_key(haul$haul_id)
hl <- load_aligned_hl(".")$hl
bins <- clean_hl_raw(hl) %>% filter(haul_key %in% haul$haul_key)

have <- cur$aphia_id[is.finite(cur$a_final) & is.finite(cur$b_final)]
dropped <- bins %>% filter(!AphiaID %in% have) %>%
  group_by(aphia_id = AphiaID) %>%
  summarise(individuals = sum(NumberAtLength), hauls = n_distinct(haul_key), .groups = "drop") %>%
  mutate(in_lookup_without_coef = aphia_id %in% cur$aphia_id,
         pct_of_N = 100 * individuals / sum(bins$NumberAtLength)) %>%
  left_join(new %>% select(aphia_id, worms_name, worms_status, worms_valid_name, worms_rank, fb_name,
                           fallback_level, a_final, b_final), by = "aphia_id") %>%
  arrange(desc(individuals))
# share of S: per haul, dropped species / (kept + dropped species)
s_all <- bins %>% group_by(haul_key) %>% summarise(S_all = n_distinct(AphiaID), .groups = "drop")
s_drop <- bins %>% filter(!AphiaID %in% have) %>% group_by(haul_key) %>%
  summarise(S_drop = n_distinct(AphiaID), .groups = "drop")
sh <- s_all %>% left_join(s_drop, by = "haul_key") %>% mutate(S_drop = coalesce(S_drop, 0L))
say("Dropped taxa in H1 hauls: %d (%d absent from lookup, %d in lookup without a/b); %.3f%% of N",
    nrow(dropped), sum(!dropped$in_lookup_without_coef), sum(dropped$in_lookup_without_coef), sum(dropped$pct_of_N))
say("H1 hauls losing >= 1 species: %d of %d; among them median share of S lost %.1f%%, max %.1f%%",
    sum(sh$S_drop > 0), nrow(sh), 100 * median((sh$S_drop / sh$S_all)[sh$S_drop > 0]),
    100 * max(sh$S_drop / sh$S_all))
say("Resolved by the rebuilt lookup: %d of %d", sum(is.finite(dropped$a_final)), nrow(dropped))
write_csv(dropped, file.path(out_dir, "dropped_taxa.csv"))

# Parameter comparison for taxa with coefficients in the current lookup
cmp <- cur %>% filter(is.finite(a_final)) %>%
  select(aphia_id, accepted_name, a_cur = a_final, b_cur = b_final, level_cur = fallback_level) %>%
  inner_join(new %>% select(aphia_id, a_new = a_final, b_new = b_final, level_new = fallback_level, fb_name),
             by = "aphia_id") %>%
  mutate(w30_ratio = (a_new * 30^b_new) / (a_cur * 30^b_cur),
         identical = abs(a_new - a_cur) < 1e-12 & abs(b_new - b_cur) < 1e-12)
say("Current taxa with coefficients: %d; present in rebuild: %d; identical a,b: %d; |W(30 cm) ratio - 1| > 10%%: %d",
    sum(is.finite(cur$a_final)), nrow(cmp), sum(cmp$identical, na.rm = TRUE),
    sum(abs(cmp$w30_ratio - 1) > 0.1, na.rm = TRUE))
write_csv(cmp, file.path(out_dir, "lookup_parameter_comparison.csv"))

# Effect on S, N, E, H1: (i) current lookup + resolved dropped taxa; (ii) fully rebuilt lookup
lk_i <- bind_rows(cur %>% filter(is.finite(a_final)),
                  dropped %>% filter(is.finite(a_final)) %>% select(aphia_id, a_final, b_final))
lk_ii <- new %>% filter(is.finite(a_final))
st0 <- build_state_from_hl(hl, cur)$haul_state %>% filter(haul_key %in% haul$haul_key)
run_variant <- function(lk, label) {
  st <- build_state_from_hl(hl, lk)$haul_state %>% filter(haul_key %in% haul$haul_key) %>%
    mutate(E = E_raw / m_min^0.75)
  m <- st0 %>% mutate(E0 = E_raw / m_min^0.75) %>% select(haul_key, S0 = S, N0 = N, E0, m0 = m_min) %>%
    inner_join(st %>% select(haul_key, S, N, E, m_min), by = "haul_key")
  ch <- m %>% filter(S != S0 | abs(N / N0 - 1) > 1e-12 | abs(E / E0 - 1) > 1e-12)
  say("[%s] hauls changed: %d; S up in %d (median +%g); median N ratio %.4f; median E ratio %.4f",
      label, nrow(ch), sum(ch$S > ch$S0), median((ch$S - ch$S0)[ch$S > ch$S0]), median(ch$N / ch$N0), median(ch$E / ch$E0))
  pr <- audit_eeos_predict(ch$S, ch$N, ch$E)
  ch <- bind_cols(ch, pr)
  h2 <- haul %>% left_join(ch %>% select(haul_key, Bn = B_pred_norm, mm = m_min), by = "haul_key") %>%
    mutate(B_pred = ifelse(is.na(Bn), B_pred, Bn * mm))
  write_csv(ch, file.path(out_dir, sprintf("changed_hauls_%s.csv", gsub("[^a-z]", "", tolower(label)))))
  audit_h1_metrics(audit_recompute_residuals(h2)) %>%
    mutate(variant = label, hauls_changed = nrow(ch),
           failed_solves = sum(ch$max_abs_constraint_dev > 0.05, na.rm = TRUE), .before = 1)
}
h1 <- bind_rows(
  audit_h1_metrics(haul) %>% mutate(variant = "current", hauls_changed = 0L, failed_solves = NA, .before = 1),
  run_variant(lk_i, "current + resolved dropped taxa"),
  run_variant(lk_ii, "fully rebuilt lookup (FishBase 26.06)")
)
print(as.data.frame(h1), digits = 4)
write_csv(h1, file.path(out_dir, "h1_effect.csv"))
