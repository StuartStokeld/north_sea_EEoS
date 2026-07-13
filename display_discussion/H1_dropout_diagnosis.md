# H1 dropout diagnosis

Generated: 2026-07-13 from `pipeline/run_h1_dropout_diagnosis.R`.

## Root cause of historical 1998 / 2013–14 spikes

Earlier drafts showed large year-level dropout in 1998 and 2013–14. That was **not** missing DATRAS HL years.
Bins with missing FishBase LW parameters previously propagated `NA` through `sum()` / `min()`, producing hauls with `NA E_raw` / `NA m_min` that survived the join and appeared as mass dropout.

**Fix (implemented):** `add_lw_mass_to_bins()` drops bad LW bins; `filter_valid_haul_state()` excludes incomplete hauls before `datras_haul_state.rds` is saved; `prepare_datras_for_join()` enforces valid normalised E before the FishGlob join (`stopifnot(null E == 0)`).

## Spike-year check (1998, 2013, 2014)

| Year | HL hauls | DATRAS state | FishGlob joined | EEoS predictions | HL→DATRAS drop | DATRAS→join drop | Join→pred drop | Status |
|------|----------|--------------|-----------------|------------------|----------------|------------------|----------------|--------|
| 1998 | 404 | 404 | 404 | 404 | 0 | 0 | 0 | no_dropout_spike |
| 2013 | 394 | 394 | 393 | 393 | 0 | 1 | 0 | DATRAS hauls without FishGlob B_obs |
| 2014 | 342 | 342 | 340 | 340 | 0 | 2 | 0 | DATRAS hauls without FishGlob B_obs |

## Current largest filter-year drop

Year **1991**: 424 joined hauls → 384 predictions (40 excluded at EEoS filters).

The dominant exclusion reason is **`bad_B_obs`** (FishGlob `sum(wgt) == 0`), concentrated in **1991** (40 hauls). This is intentional: `B_obs > 0` is required for log-scale comparison.

## Filter exclusions (all years)

| Reason | Hauls |
|--------|------:|
| bad_B_obs | 44 |
| S_lt_2 | 4 |

## Audit verdict

| Check | Status | Detail |
|-------|--------|--------|
| spike_years_hl_retention | PASS | 1998, 2013, 2014 |
| spike_years_join_retention | PASS | max datras->join drop = 2 (2014) |
| spike_years_prediction_retention | PASS | max join->pred drop = 0 (1998, 2013, 2014) |

**Panel:** 12,069 hauls with EEoS predictions (12389 DATRAS state hauls; 12117 FishGlob-matched).

*Outputs: `outputs/h1_dropout_funnel_by_year.csv`, `outputs/h1_filter_exclusions.csv`, `outputs/h1_join_gaps.csv`.*
