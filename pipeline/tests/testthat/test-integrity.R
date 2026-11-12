# Layer 4.5: pipeline integrity - accounting, join cardinality, provenance.

suppressPackageStartupMessages(library(dplyr))

h1_rebuild <- function() {
  if (!is.null(.audit_cache$h1)) return(.audit_cache$h1)
  root <- audit_repo_root()
  for (f in c("h1_common.R", "datras_constants.R", "datras_hl_helpers.R", "datras_state_helpers.R",
              "datras_csv_import.R", "datras_hh_helpers.R", "h1_join_helpers.R")) {
    source_pipeline(f, env = globalenv())
  }
  fg_path <- skip_if_no_file("FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
  al <- suppressMessages(load_aligned_hl(root))
  lw <- read.csv(file.path(root, "outputs", "fishbase_lw_lookup_v2.csv"))
  st <- build_state_from_hl(al$hl, lw)$haul_state %>% filter_valid_haul_state() %>% normalize_haul_E()
  ready <- prepare_datras_for_join(st)
  fg <- build_fishglob_haul_table(fg_path)
  joined <- join_fishglob_datras(fg, ready)
  eeos <- filter_eeos_inputs(joined)
  .audit_cache$h1 <- list(al = al, st = st, ready = ready, fg = fg, joined = joined, eeos = eeos)
  .audit_cache$h1
}

test_that("accounting: every HL haul ends in exactly one effort bucket; counts reconcile", {
  h <- h1_rebuild()
  hl_hauls <- cached_rds("outputs", "datras_hl_raw.rds") %>%
    distinct(Survey, Year, Quarter, Country, Ship, HaulNo) %>% nrow()
  expect_equal(nrow(h$al$hauls), hl_hauls) # one row per HL haul
  expect_false(anyDuplicated(h$al$hauls$haul_key) > 0)
  expect_true(all(h$al$hauls$effort_reason %in% c("kept", "no_hh", "ambiguous_hh_key", "invalid", "datatype", "duration")))
  tab <- table(h$al$hauls$effort_reason)
  # documented counts (CLAUDE.md / outputs/datras_effort_validation): 12,389 = 10,406 + 9 + 33 + 1,941
  expect_equal(hl_hauls, 12389L)
  expect_equal(as.integer(tab[["kept"]]), 10406L)
  expect_equal(sum(tab), hl_hauls)
  expect_setequal(unique(h$al$hl$haul_key), h$al$hauls$haul_key[h$al$hauls$effort_reason == "kept"])
})

test_that("accounting: H1 stage counts 10,406 -> 10,172 -> 10,128 and the outputs match the rebuild", {
  h <- h1_rebuild()
  expect_equal(nrow(h$st), 10406L)
  expect_equal(nrow(h$joined), 10172L)
  expect_equal(nrow(h$eeos), 10128L)
  pred <- cached_rds("outputs", "haul_eeos_predictions.rds")
  expect_setequal(h$eeos$haul_id, pred$haul_id)
  m <- inner_join(h$eeos %>% select(haul_id, S, N, E, B_obs, m_min),
                  pred %>% select(haul_id, S, N, E, B_obs, m_min), by = "haul_id")
  expect_equal(m$S.x, m$S.y)
  expect_equal(m$N.x, m$N.y)
  expect_equal(m$E.x, m$E.y)
  expect_equal(m$B_obs.x, m$B_obs.y)
})

test_that("joins neither duplicate nor lose rows", {
  h <- h1_rebuild()
  expect_false(anyDuplicated(h$fg$join_key) > 0)
  expect_false(anyDuplicated(h$ready$haul_key) > 0)
  expect_equal(nrow(h$joined), length(intersect(h$fg$join_key, h$ready$haul_key)))
  expect_false(anyDuplicated(h$joined$join_key) > 0)
})

test_that("accounting: H2/H3 panel = 8,670 hauls in 158 rectangles, each with >= 10 H1 hauls", {
  v2 <- cached_rds("outputs", "primary_model_v2.rds")
  d <- v2$data
  expect_equal(nrow(d), 8670L)
  expect_equal(length(unique(d$stat_rec)), 158L)
  pred <- cached_rds("outputs", "haul_eeos_predictions.rds")
  n_h1 <- table(pred$stat_rec)[unique(as.character(d$stat_rec))]
  expect_true(all(n_h1 >= 10))
  expect_true(all(d$haul_id %in% pred$haul_id))
})

test_that("[F-B13] provenance: every input of a live script is produced by a live script or is a documented raw input", {
  g <- dg_build(audit_repo_root())
  li <- dg_live_inputs(g, DG_RAW_INPUTS)
  flagged <- li[grepl("FLAG", li$status), ]
  expect_equal(nrow(flagged), 0L,
               info = paste(unique(paste(flagged$script, "<-", flagged$path)), collapse = "\n"))
})

test_that("[F-B13] provenance: no live script reads from exploratory/", {
  g <- dg_build(audit_repo_root())
  bad <- g[g$script_class == "live" & g$direction == "read" & grepl("^exploratory/", g$path), ]
  expect_equal(nrow(bad), 0L, info = paste(bad$script, bad$path))
})

test_that("provenance: no hard-coded absolute paths in pipeline code", {
  files <- c(list.files(audit_pipeline_dir(), pattern = "\\.R$", full.names = TRUE),
             list.files(file.path(audit_pipeline_dir(), "R"), pattern = "\\.R$", full.names = TRUE))
  hits <- unlist(lapply(files, function(f) {
    x <- readLines(f, warn = FALSE)
    if (any(grepl("/Users/|C:\\\\", x))) basename(f) else NULL
  }))
  expect_length(hits, 0L)
})

# Determinism checks (seeds, two-run identity) deferred by scope decision (12 Nov 2026).
