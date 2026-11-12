# Layer 4.2: input data contracts for FishGlob, the length-weight lookup, Couce effort,
# the CAR adjacency and the EEoS input table.

suppressPackageStartupMessages(library(dplyr))

fishglob_q1 <- function() {
  p <- skip_if_no_file("FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
  if (is.null(.audit_cache$fg)) {
    e <- new.env()
    load(p, envir = e)
    .audit_cache$fg <- e$data %>% filter(quarter == 1, year >= 1985, year <= 2015)
  }
  .audit_cache$fg
}

# --- FishGlob ----------------------------------------------------------------

test_that("PROVENANCE (change detector): FishGlob RData is the traced FishGlob commit acde3dc / v2.1", {
  p <- skip_if_no_file("FishGlob_data", "outputs", "Cleaned_data", "NS-IBTS_clean.RData")
  expect_identical(digest::digest(file = p, algo = "sha256"),
                   "bd98d2eae367eba9aa82934996478c765a4eb974d95e8ec0edb0a9ec6b7dad13")
})

test_that("FishGlob schema and haul_id format (8 space-separated fields, unique 6-field key)", {
  fg <- fishglob_q1()
  expect_true(all(c("haul_id", "year", "quarter", "stat_rec", "wgt", "num", "accepted_name",
                    "verbatim_aphia_id") %in% names(fg)))
  ids <- unique(fg$haul_id)
  expect_true(all(lengths(strsplit(ids, " ")) == 8L))
  expect_false(anyDuplicated(fishglob_key(ids)) > 0)
  expect_true(all(grepl("^[0-9]{2}[A-M][0-9]$", unique(fg$stat_rec))))
})

test_that("FishGlob wgt is non-negative kg", {
  fg <- fishglob_q1()
  expect_true(all(fg$wgt >= 0, na.rm = TRUE))
  per_haul_kg <- tapply(fg$wgt, fg$haul_id, sum, na.rm = TRUE)
  # a 30-60 min GOV haul: median catch of order 10-1000 kg (g or t would be off by 1e3)
  expect_gt(median(per_haul_kg), 10)
  expect_lt(median(per_haul_kg), 1000)
})

test_that("[F-B9] every FishGlob taxon weight in the H1 hauls is present (no NA counted as 0 in B_obs)", {
  fg <- fishglob_q1()
  pred <- cached_rds("outputs", "haul_eeos_predictions.rds")
  fg1 <- fg[fg$haul_id %in% pred$haul_id, ]
  expect_equal(sum(is.na(fg1$wgt)), 0L)
})

# --- Length-weight lookup ----------------------------------------------------

lw_lookup <- function() {
  p <- skip_if_no_file("outputs", "fishbase_lw_lookup_v2.csv")
  read.csv(p, stringsAsFactors = FALSE)
}

test_that("LW lookup schema: one row per AphiaID, a > 0, b plausible", {
  lw <- lw_lookup()
  expect_true(all(c("accepted_name", "aphia_id", "a_final", "b_final", "Type", "fallback_level") %in% names(lw)))
  expect_false(anyDuplicated(lw$aphia_id) > 0)
  ok <- is.finite(lw$a_final) & is.finite(lw$b_final)
  expect_true(all(lw$a_final[ok] > 0 & lw$a_final[ok] < 1))
  expect_true(all(lw$b_final[ok] >= 2 & lw$b_final[ok] <= 4))
})

test_that("LW b outside 2.5-3.5 only for the documented records (audit F-A11)", {
  lw <- lw_lookup()
  ok <- is.finite(lw$b_final)
  out <- sort(lw$aphia_id[ok & (lw$b_final < 2.5 | lw$b_final > 3.5)])
  documented <- sort(c(126868L, 105814L, 127312L, 154675L, 126904L, 127072L, 127193L))
  expect_identical(as.integer(out), documented)
})

test_that("[F-A11] LW records used without conversion are total-length (TL) relationships", {
  lw <- lw_lookup()
  used <- lw[is.finite(lw$a_final), ]
  allowed <- c("TL_direct", "genus_mean", "family_mean")
  expect_true(all(used$fallback_level %in% allowed),
              info = paste(unique(used$fallback_level[!used$fallback_level %in% allowed]), collapse = ", "))
})

test_that("[F-B8] every fish taxon in the pipeline's HL rows has LW coefficients", {
  hl <- cached_rds("outputs", "datras_hl_raw.rds")
  source_pipeline("datras_constants.R")
  lw <- lw_lookup()
  fish <- unique(hl$Valid_Aphia[!hl$Valid_Aphia %in% NON_FISH_APHIAIDS])
  have <- lw$aphia_id[is.finite(lw$a_final) & is.finite(lw$b_final)]
  missing <- setdiff(fish, have)
  expect_length(missing, 0L)
})

test_that("[F-B12] the non-fish exclusion list has no duplicates and covers bivalves/cephalopods in HL", {
  source_pipeline("datras_constants.R")
  expect_false(anyDuplicated(NON_FISH_APHIAIDS) > 0)
  # Ostrea edulis (140658) and Todarodes (138281) occur in HL and are not fish
  expect_true(all(c(140658L, 138281L) %in% NON_FISH_APHIAIDS))
})

# --- Couce fishing effort ----------------------------------------------------

test_that("Couce rectangle-year effort: hours >= 0, otter + beam = total, unique keys, in scope", {
  cy <- cached_rds("outputs", "h2_couce_year_effort.rds")
  expect_true(all(c("stat_rec", "year", "hours_otter", "hours_beam", "hours_total") %in% names(cy)))
  expect_false(anyDuplicated(cy[, c("stat_rec", "year")]) > 0)
  expect_true(all(cy$year >= 1985L & cy$year <= 2015L))
  expect_true(all(cy$hours_otter >= 0 & cy$hours_beam >= 0, na.rm = TRUE))
  expect_false(anyNA(cy$hours_total))
  expect_equal(cy$hours_total, cy$hours_otter + cy$hours_beam)
  expect_true(all(grepl("^[0-9]{2}[A-M][0-9]$", cy$stat_rec)))
})

test_that("Couce rectangle summary agrees with the rectangle-year table", {
  cy <- cached_rds("outputs", "h2_couce_year_effort.rds")
  cr <- cached_rds("outputs", "h2_couce_rectangle_effort.rds")
  agg <- cy %>% group_by(stat_rec) %>% summarise(tot = sum(hours_total), n = n(), .groups = "drop")
  m <- inner_join(cr, agg, by = "stat_rec")
  expect_equal(nrow(m), nrow(cr))
  expect_equal(m$total_hours, m$tot)
})

# --- CAR adjacency -----------------------------------------------------------

test_that("CAR adjacency: binary, symmetric, zero diagonal, no islands, one component", {
  adj <- cached_rds("outputs", "h2h3_feasibility_round2_model_objects.rds")$adjMatrix
  A <- as.matrix(adj)
  attributes(A) <- list(dim = dim(A), dimnames = dimnames(A))
  expect_true(all(A %in% c(0, 1)))
  expect_true(isSymmetric(unname(A)))
  expect_true(all(diag(A) == 0))
  expect_true(all(rowSums(A) >= 1))
  # connectivity by breadth-first search
  seen <- 1L
  frontier <- 1L
  while (length(frontier)) {
    nb <- setdiff(which(colSums(A[frontier, , drop = FALSE]) > 0), seen)
    seen <- c(seen, nb)
    frontier <- nb
  }
  expect_equal(length(seen), nrow(A))
})

test_that("CAR adjacency equals queen contiguity derived from the rectangle codes", {
  adj <- cached_rds("outputs", "h2h3_feasibility_round2_model_objects.rds")$adjMatrix
  A <- as.matrix(adj)
  rc <- ices_rect_rowcol(rownames(A))
  Q <- outer(seq_len(nrow(rc)), seq_len(nrow(rc)), function(i, j) {
    as.numeric(abs(rc[i, "row"] - rc[j, "row"]) <= 1 & abs(rc[i, "col"] - rc[j, "col"]) <= 1 & i != j)
  })
  expect_equal(dim(A), dim(Q))
  expect_equal(as.vector(A), as.vector(Q))
})

test_that("CAR adjacency rectangles are exactly the H2/H3 analysis rectangles", {
  adj <- cached_rds("outputs", "h2h3_feasibility_round2_model_objects.rds")$adjMatrix
  panel <- cached_rds("outputs", "h2_rectangle_panel.rds")
  expect_setequal(rownames(adj), panel$stat_rec)
})

# --- EEoS inputs / predictions table ----------------------------------------

test_that("EEoS haul table satisfies the solver's domain and the documented definitions", {
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  expect_false(anyDuplicated(h$haul_id) > 0)
  expect_true(all(h$S >= 2 & h$N >= h$S & h$E > h$N))
  expect_true(all(h$m_min > 0 & h$B_obs > 0 & is.finite(h$B_pred) & h$B_pred > 0))
  expect_equal(h$E, h$E_raw / h$m_min^0.75)
  expect_equal(h$B_pred, h$B_pred_norm * h$m_min)
  expect_equal(h$residual, log(h$B_obs) - log(h$B_pred))
  expect_true(all(h$year >= 1985L & h$year <= 2015L))
})

test_that("[F-B14] EEoS haul table carries the haul_key used for joins", {
  h <- cached_rds("outputs", "haul_eeos_predictions.rds")
  expect_true("haul_key" %in% names(h))
})
