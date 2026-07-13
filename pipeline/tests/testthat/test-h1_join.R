test_that("prepare_datras_for_join removes NA E", {
  skip_if_not_installed("dplyr")
  source(file.path("..", "..", "R", "h1_join_helpers.R"))

  datras <- tibble(
    haul_key = c("a", "b", "c"),
    S = c(3L, 2L, 4L),
    N = c(10L, 5L, 20L),
    E_raw = c(100, NA_real_, 200),
    E = c(50, NA_real_, 80),
    m_min = c(0.1, 0.2, NA_real_)
  )

  ready <- prepare_datras_for_join(datras)
  expect_equal(nrow(ready), 1L)
  expect_false(any(is.na(ready$E)))
})

test_that("join_fishglob_datras fails on null E", {
  skip_if_not_installed("dplyr")
  source(file.path("..", "..", "R", "h1_join_helpers.R"))

  fg <- tibble(
    haul_id = "h1",
    join_key = "k1",
    B_obs = 1000,
    stat_rec = "44G1",
    year = 2000L
  )
  datras <- tibble(
    haul_key = "k1",
    haul_id = "h1",
    S = 3L,
    N = 10L,
    E_raw = 100,
    E = NA_real_,
    m_min = 0.1,
    m_min_epsilon = 0.1,
    min_epsilon = 1,
    n_species_with_lw = 3L,
    m_min_method = "min_length_bin"
  )

  expect_error(join_fishglob_datras(fg, datras), "NA E")
})

test_that("classify_eeos_exclusion labels zero biomass", {
  skip_if_not_installed("dplyr")
  source(file.path("..", "..", "R", "h1_join_helpers.R"))

  haul <- tibble(
    E = c(20, 20),
    E_raw = c(100, 100),
    m_min = c(0.1, 0.1),
    S = c(3L, 3L),
    N = c(10L, 10L),
    B_obs = c(100, 0)
  )
  reasons <- classify_eeos_exclusion(haul)
  expect_equal(reasons, c("passed", "bad_B_obs"))
})
