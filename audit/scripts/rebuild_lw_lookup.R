# DRAFT (audit): re-runnable rebuild of the FishBase length-weight lookup.
#
# Replaces the non-runnable supplementary/fishbase_lookup.Rmd (undefined `all_spp`, hard-coded
# paths) with a script that:
#   1. builds the taxon list from the pipeline's own HL rows (effort-aligned, after clean_hl_raw()),
#      keyed by the DATRAS Valid_Aphia the pipeline joins on (not FishGlob names);
#   2. resolves each AphiaID's name with WoRMS (cached CSV, so re-runs are offline);
#   3. matches WoRMS valid / verbatim names to FishBase valid names (FishBase synonym step: follow-up);
#   4. selects one LW record per species with the notebook's rule (TL > FL > SL > NA type > other,
#      then largest Number);
#   5. falls back to genus mean, then family mean, of a and b over TL / NA-type records
#      (notebook rule: pool = records of the species in the list), then to the same means over all
#      FishBase species of that genus / family (extended pool, flagged separately).
#
# Not wired into the pipeline. Outputs go to audit/outputs/lw_rebuild/.
#
# Run from the repo root:
#   Rscript audit/scripts/rebuild_lw_lookup.R                 # FishBase release 25.04 and latest
#   EEOS_FB_VERSIONS=26.06 Rscript audit/scripts/rebuild_lw_lookup.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rfishbase)
})
for (f in c("h1_common.R", "datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
            "datras_csv_import.R", "datras_hh_helpers.R")) {
  source(file.path("pipeline", "R", f))
}
out_dir <- "audit/outputs/lw_rebuild"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")
versions <- strsplit(Sys.getenv("EEOS_FB_VERSIONS", "25.04,26.06"), ",")[[1]]

# ---------------------------------------------------------------------------
# 1. Taxon list from the pipeline's HL rows
# ---------------------------------------------------------------------------
al <- load_aligned_hl(".")
bins <- clean_hl_raw(al$hl)
taxa_list <- bins %>%
  group_by(aphia_id = AphiaID) %>%
  summarise(total_N = sum(NumberAtLength), n_rows = n(), n_hauls = n_distinct(haul_key),
            .groups = "drop") %>%
  arrange(desc(total_N))
say("Fish taxa in pipeline HL rows: %d", nrow(taxa_list))

# ---------------------------------------------------------------------------
# 2. WoRMS records (cached)
# ---------------------------------------------------------------------------
worms_cache <- file.path(out_dir, "worms_records.csv")
`%||%` <- function(a, b) if (is.null(a)) b else a
worms <- if (file.exists(worms_cache)) read_csv(worms_cache, show_col_types = FALSE) else NULL
todo <- setdiff(taxa_list$aphia_id, worms$AphiaID)
if (length(todo) > 0L) {
  say("Querying WoRMS for %d AphiaIDs (only the IDs are sent)", length(todo))
  fetch <- function(id) {
    url <- sprintf("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/%d", id)
    for (attempt in 1:3) {
      r <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
      if (!is.null(r)) break
      Sys.sleep(2 * attempt)
    }
    Sys.sleep(0.2)
    if (is.null(r)) return(NULL) # not cached: retried on the next run
    tibble(
      AphiaID = id, scientificname = r$scientificname %||% NA, rank = r$rank %||% NA,
      status = r$status %||% NA, valid_AphiaID = r$valid_AphiaID %||% NA,
      valid_name = r$valid_name %||% NA, genus = r$genus %||% NA, family = r$family %||% NA,
      class = r$class %||% NA, source = "WoRMS REST"
    )
  }
  probe <- tryCatch(jsonlite::fromJSON("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/126436"),
                    error = function(e) NULL)
  got <- if (is.null(probe)) {
    say("WoRMS REST is unavailable; skipping (re-run later to fill the cache)")
    tibble()
  } else {
    bind_rows(lapply(todo, fetch))
  }
  if (nrow(got) > 0L) {
    worms <- bind_rows(worms, got)
    write_csv(worms, worms_cache)
  }
  say("WoRMS records obtained: %d of %d", nrow(got), length(todo))
}
# Fallback for IDs WoRMS did not return: DATRAS species name from the raw ICES CSV (if present)
still <- setdiff(taxa_list$aphia_id, worms$AphiaID)
if (length(still) > 0L) {
  csv <- list.files(".", pattern = "^Unaggregated trawl and biological information.*\\.csv$",
                    recursive = TRUE, full.names = TRUE)
  if (length(csv) > 0L) {
    hl_names <- read_ices_csv_blocks(csv[1])$HL %>% distinct(AphiaID, Species) %>%
      filter(AphiaID %in% still) %>% group_by(AphiaID) %>% slice(1) %>% ungroup()
    worms <- bind_rows(worms, hl_names %>% transmute(
      AphiaID, scientificname = Species, valid_name = Species,
      rank = ifelse(grepl(" ", Species), "Species", NA_character_),
      genus = sub(" .*", "", Species), source = "DATRAS CSV name (WoRMS unavailable)"
    ))
    say("WoRMS unavailable for %d IDs: used DATRAS species names", length(still))
  }
}
taxa_list <- taxa_list %>% left_join(worms, by = c("aphia_id" = "AphiaID"))

# ---------------------------------------------------------------------------
# 3-5. Per FishBase release
# ---------------------------------------------------------------------------
pick_best <- function(lw) {
  lw %>%
    mutate(
      type_priority = case_when(
        Type == "TL" ~ 1, Type == "FL" ~ 2, Type == "SL" ~ 3, is.na(Type) ~ 4,
        Type %in% c("NG", "OT", "WD") ~ 5, TRUE ~ 6
      ),
      sample_size = suppressWarnings(as.numeric(Number))
    ) %>%
    filter(is.finite(a), is.finite(b)) %>%
    group_by(Species) %>%
    arrange(type_priority, desc(sample_size), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(lw_type_flag = case_when(
      Type == "TL" ~ "TL_direct", Type == "FL" ~ "FL_needs_conversion",
      Type == "SL" ~ "SL_needs_conversion", is.na(Type) ~ "type_unknown", TRUE ~ "other_review"
    ))
}

build_one <- function(version) {
  say("--- FishBase release %s ---", version)
  fb_taxa <- as.data.frame(load_taxa(version = version)) %>% select(Species, Genus, Family)
  fb_species <- unique(fb_taxa$Species)

  # Name resolution: WoRMS valid name, then WoRMS verbatim name, then FishBase synonyms
  sp <- taxa_list %>% filter(rank == "Species")
  cand <- sp %>%
    mutate(
      fb_name = case_when(
        valid_name %in% fb_species ~ valid_name,
        scientificname %in% fb_species ~ scientificname,
        TRUE ~ NA_character_
      ),
      fb_match = case_when(
        valid_name %in% fb_species ~ "WoRMS valid name",
        scientificname %in% fb_species ~ "WoRMS verbatim name",
        TRUE ~ NA_character_
      )
    )
  # FishBase synonym resolution: deferred follow-up (scope decision 12 Nov 2026)
  lw_raw <- as.data.frame(length_weight(unique(na.omit(cand$fb_name)), version = version))
  best <- pick_best(lw_raw)

  direct <- cand %>%
    left_join(best %>% select(Species, a, b, Type, lw_type_flag, Number),
              by = c("fb_name" = "Species"))

  # Fallback pools (TL or NA-type records)
  pool_list <- lw_raw %>% filter(Type == "TL" | is.na(Type), is.finite(a), is.finite(b)) %>%
    left_join(fb_taxa, by = "Species")
  genus_list <- pool_list %>% group_by(Genus) %>%
    summarise(a_g = mean(a), b_g = mean(b), n_g = n(), .groups = "drop")
  family_list <- pool_list %>% group_by(Family) %>%
    summarise(a_f = mean(a), b_f = mean(b), n_f = n(), .groups = "drop")

  others <- taxa_list %>% filter(is.na(rank) | rank != "Species") %>%
    mutate(fb_name = NA_character_, fb_match = NA_character_)
  all_t <- bind_rows(direct, others) %>%
    mutate(
      genus_key = ifelse(!is.na(rank) & rank == "Genus", scientificname, genus),
      family_key = ifelse(!is.na(rank) & rank == "Family", scientificname, family)
    )
  # FishBase family for the genus where WoRMS and FishBase disagree
  fb_gen_fam <- fb_taxa %>% distinct(Genus, Family) %>% group_by(Genus) %>% filter(n() == 1L) %>% ungroup()
  all_t <- all_t %>% left_join(fb_gen_fam %>% rename(fb_family = Family), by = c("genus_key" = "Genus")) %>%
    mutate(family_key = coalesce(fb_family, family_key))

  # Extended pool: all FishBase species in the genus / family (only fetched for unresolved taxa)
  unresolved <- all_t %>% filter(is.na(a))
  in_list_g <- unresolved$genus_key %in% genus_list$Genus
  in_list_f <- unresolved$family_key %in% family_list$Family
  ext_g <- unique(na.omit(unresolved$genus_key[!in_list_g]))
  ext_f <- unique(na.omit(unresolved$family_key[!in_list_g & !in_list_f]))
  ext_species <- fb_taxa %>% filter(Genus %in% ext_g | Family %in% ext_f) %>% pull(Species) %>% unique()
  genus_ext <- family_ext <- NULL
  if (length(ext_species) > 0L) {
    pool_ext <- as.data.frame(length_weight(ext_species, version = version)) %>%
      filter(Type == "TL" | is.na(Type), is.finite(a), is.finite(b)) %>%
      left_join(fb_taxa, by = "Species")
    genus_ext <- pool_ext %>% group_by(Genus) %>%
      summarise(a_ge = mean(a), b_ge = mean(b), n_ge = n(), .groups = "drop")
    family_ext <- pool_ext %>% group_by(Family) %>%
      summarise(a_fe = mean(a), b_fe = mean(b), n_fe = n(), .groups = "drop")
  }
  all_t <- all_t %>%
    left_join(genus_list, by = c("genus_key" = "Genus")) %>%
    left_join(family_list, by = c("family_key" = "Family"))
  if (!is.null(genus_ext)) {
    all_t <- all_t %>% left_join(genus_ext, by = c("genus_key" = "Genus")) %>%
      left_join(family_ext, by = c("family_key" = "Family"))
  } else {
    all_t <- all_t %>% mutate(a_ge = NA_real_, b_ge = NA_real_, a_fe = NA_real_, b_fe = NA_real_)
  }
  all_t %>%
    mutate(
      fallback_level = case_when(
        !is.na(a) ~ lw_type_flag,
        !is.na(a_g) ~ "genus_mean",
        !is.na(a_f) ~ "family_mean",
        !is.na(a_ge) ~ "genus_mean_fishbase",
        !is.na(a_fe) ~ "family_mean_fishbase",
        TRUE ~ "MISSING"
      ),
      a_final = coalesce(a, a_g, a_f, a_ge, a_fe),
      b_final = coalesce(b, b_g, b_f, b_ge, b_fe),
      accepted_name = coalesce(fb_name, valid_name, scientificname),
      fishbase_version = version,
      rfishbase_version = as.character(utils::packageVersion("rfishbase"))
    ) %>%
    select(accepted_name, aphia_id, total_N, n_hauls, a_final, b_final, Type, fallback_level,
           worms_name = scientificname, worms_status = status, worms_valid_name = valid_name,
           worms_rank = rank, fb_name, fb_match, genus_key, family_key, Number,
           fishbase_version, rfishbase_version) %>%
    arrange(desc(total_N))
}

for (v in versions) {
  lk <- build_one(v)
  write_csv(lk, file.path(out_dir, sprintf("fishbase_lw_lookup_rebuilt_%s.csv", v)))
  say("release %s: %d taxa; by fallback level:", v, nrow(lk))
  print(as.data.frame(lk %>% count(fallback_level)))
}
say("Done. Outputs in %s", out_dir)
