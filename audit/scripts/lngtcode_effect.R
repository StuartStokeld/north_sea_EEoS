# Effect of the LngtCode "0" / "5" midpoint errors (findings F-A1, F-A2) on H1 and the primary CAR.
#
# Rebuilds haul state with the ICES-correct midpoints ("." 1 mm in mm; "0" 5 mm in mm; "1" 1 cm in cm;
# "5" 5 cm in cm), re-runs EEoS for every haul whose S, N, E or m_min changed, and refits the primary
# CAR on the signed and absolute residual. The live pipeline is not modified: the corrected midpoint
# function is swapped in only inside this script. EEoS results are cached in audit/outputs/lngtcode/.
#
# Run from the repo root: Rscript audit/scripts/lngtcode_effect.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("audit/scripts/audit_refit_helpers.R")
for (f in c("datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
            "datras_csv_import.R", "datras_hh_helpers.R")) {
  source(file.path("pipeline", "R", f))
}
out_dir <- "audit/outputs/lngtcode"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")

midpoint_ices <- function(code, cls) {
  code <- as.character(code)
  cls <- as.numeric(cls)
  out <- rep(NA_real_, length(cls))
  out[code == "."] <- (cls[code == "."] + 0.5) / 10
  out[code == "0"] <- (cls[code == "0"] + 2.5) / 10
  out[code == "1"] <- cls[code == "1"] + 0.5
  out[code == "5"] <- cls[code == "5"] + 2.5
  out
}

haul <- readRDS("outputs/haul_eeos_predictions.rds")
haul$haul_key <- audit_haul_key(haul$haul_id)
lw <- read_csv("outputs/fishbase_lw_lookup_v2.csv", show_col_types = FALSE)
hl <- load_aligned_hl(".")$hl

live_mid <- lngt_midpoint_cm
cur <- build_state_from_hl(hl, lw)$haul_state
lngt_midpoint_cm <- midpoint_ices
fix <- build_state_from_hl(hl, lw)$haul_state
lngt_midpoint_cm <- live_mid

st <- cur %>% select(haul_key, S, N, E_raw, m_min) %>%
  inner_join(fix %>% select(haul_key, S2 = S, N2 = N, E_raw2 = E_raw, m_min2 = m_min), by = "haul_key") %>%
  mutate(E = E_raw / m_min^0.75, E2 = E_raw2 / m_min2^0.75) %>%
  filter(haul_key %in% haul$haul_key)
stopifnot(nrow(st) == nrow(haul))
chk <- haul %>% select(haul_key, E_live = E) %>% inner_join(st, by = "haul_key")
stopifnot(max(abs(chk$E / chk$E_live - 1)) < 1e-9)

changed <- st %>% filter(abs(E2 / E - 1) > 1e-12 | S2 != S | abs(N2 - N) > 1e-9)
say("H1 hauls whose state changes: %d of %d", nrow(changed), nrow(st))

cache <- file.path(out_dir, "eeos_corrected_midpoints.csv")
done <- if (file.exists(cache)) read_csv(cache, show_col_types = FALSE) else tibble(haul_key = character())
todo <- changed %>% filter(!haul_key %in% done$haul_key)
if (nrow(todo) > 0L) {
  say("Running EEoS for %d hauls", nrow(todo))
  t0 <- Sys.time()
  pr <- audit_eeos_predict(todo$S2, todo$N2, todo$E2)
  say("EEoS time: %.1f min", as.numeric(Sys.time() - t0, units = "mins"))
  done <- bind_rows(done, todo %>% select(haul_key, S2, N2, E2, m_min2) %>% bind_cols(pr))
  write_csv(done, cache)
}

fixed_haul <- haul %>%
  left_join(done %>% select(haul_key, m_min2, B_pred_norm2 = B_pred_norm, dev2 = max_abs_constraint_dev),
            by = "haul_key") %>%
  mutate(
    B_pred = ifelse(is.na(B_pred_norm2), B_pred, B_pred_norm2 * m_min2),
    m_min = coalesce(m_min2, m_min)
  )
say("Corrected hauls with NA prediction: %d; with |constraint dev| > 5%%: %d",
    sum(!is.na(fixed_haul$m_min2) & is.na(fixed_haul$B_pred_norm2)),
    sum(fixed_haul$dev2 > 0.05, na.rm = TRUE))

ra <- audit_run_treatment(haul, "current midpoints")
rb <- audit_run_treatment(fixed_haul, "ICES midpoints")
h1 <- bind_rows(ra$h1, rb$h1)
sl <- bind_rows(ra$slopes, rb$slopes)
print(as.data.frame(h1), digits = 4)
write_csv(h1, file.path(out_dir, "h1_metrics.csv"))
write_csv(sl, file.path(out_dir, "car_slopes.csv"))
ratio <- fixed_haul$B_pred / haul$B_pred
say("B_pred ratio corrected/current: median %.3f; hauls changed >10%%: %d", median(ratio), sum(abs(ratio - 1) > 0.1))
writeLines(c(
  "## H1", "", audit_md_table(h1, 4), "",
  "## Primary CAR phase slopes", "",
  audit_md_table(audit_format_slopes(sl) %>% select(treatment, response, component, phase, cell, p_raw, p_bonf), 3)
), file.path(out_dir, "tables.md"))
