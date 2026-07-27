# Pre-H3 exploratory visualisation & variance decomposition — run log

Descriptive/exploratory only. No H3 model fit, no H3 conclusions drawn.

## Provisional constants
- SPARSE_HAUL_THRESHOLD = 5 (hauls/rectangle/year; below this flagged 'sparse', not dropped, in count visuals)
- MIN_YEARS_PER_RECT = 10 (years at >=SPARSE_HAUL_THRESHOLD hauls needed for a rectangle to be 'usable for within-rectangle temporal analysis')
- DECADE_BINS = 1985-1994, 1995-2004, 2005-2015 (map faceting only, not modelling; provisional — supervisor may prefer different cut points)

## Sign-convention flag: signed residual
This script defines `resid_signed = log(B_pred) - log(B_obs)` per the briefing. This is the *negative* of the H1/H2 pipeline's primary `residual` column (residual = log(B_obs) - log(B_pred); see pipeline/README.md 'Key conventions'). `resid_magnitude` (= existing `abs_residual`) is identical under either convention. Flagged here, not resolved — downstream readers should check sign before comparing to H2 tables.

## Section A — temporal distribution of hauls
Loaded 12117 Q1 hauls (1985-2015) from NS-IBTS_clean.RData via build_fishglob_haul_table() — this is the full haul set for Section A (independent of DATRAS join / EEoS filter success used later in B/C).
Distinct ICES rectangles with >=1 haul in the study period: 197
A.3: 1 of 197 rectangles have >=10 years with >=5 hauls/year ('usable for within-rectangle temporal analysis'). THIS IS THE DENOMINATOR for any H3 within-rectangle design.
Context (not an exclusion, just distribution): hauls/rectangle/year across all 5392 rectangle-year cells with >=1 haul — median=2, 90th pct=3, 95th pct=4, max=12. The NS-IBTS Q1 design samples most rectangles once or twice a year, which is why SPARSE_HAUL_THRESHOLD=5 (a within-YEAR count) is a demanding bar; flagged here as context for A.3, not as a recommendation to change the provisional constants.

## Section F — structural exclusions
Rectangles with no Couce et al. (2020) fishing-hours coverage at all (Skagerrak/Kattegat, eastern English Channel per the known structural exclusion): 32 of 197. Excluded from fishing-pressure-dependent visuals/analysis (D, E) only; included in A-C (haul/biomass/residual visuals do not depend on Couce coverage).
Rectangles usable_for_fishing_analysis (usable_temporal AND has_couce_coverage), i.e. the D.3/E analysis set: 0 of 197.
Reconciliation note: this 32-rectangle count is taken over the full 197-rectangle universe used in this script's Section A (any haul, no pre-filter). The briefing's stated '26 rectangles' figure (see outputs/step0_exclusion_counts.csv, outputs/step0_robustness_check.md) is computed over a smaller 187-rectangle universe that Step 0 first restricted to >=5 total hauls pooled across all years (a different, coarser threshold to this script's per-year SPARSE_HAUL_THRESHOLD). The two counts are consistent once that universe difference is accounted for; not re-derived further here.
Saved rectangle usability flags: /Users/stuartstokeld/north_sea_eeos/outputs/h3_pre_rectangle_usability_flags.csv

## Section B — predicted vs observed biomass over time
12069 of 12117 Q1 hauls (99.6%) have a successful EEoS prediction (existing H1 pipeline output, outputs/haul_eeos_predictions.rds); 48 hauls lack B_pred and are excluded from the B/C biomass and residual visuals below. Itemised H1 exclusion reasons are not re-derived here — see outputs/h1_dropout_summary.csv and outputs/h1_filter_exclusions.csv for the existing audit trail.
B.2 (facet by region/subarea) SKIPPED: no usable region/subarea field exists in NS-IBTS_clean.RData. `sub_area` is present but 100% NA (0 of 364196 species-level rows have a value); `continent` is constant ('europe'); `stratum`/`season` are 100% NA; `survey_unit` only distinguishes Q1 vs Q3 survey (constant within this Q1-only analysis). Per the briefing, no new regional grouping was constructed for this task.

## Section C — EEoS residual over time and space
C.2: 0 haul-level rows fall outside DECADE_BINS and are excluded from the decade maps (none expected given DECADE_BINS spans 1985-2015).
C.2: 0 rectangles with residual data could not be matched to the ICES rectangle shapefile geometry and are excluded from the map (retained in the CSV-based analyses).
C.3: rectangle-year cells restricted to n_hauls >= SPARSE_HAUL_THRESHOLD (5) before fitting slopes, i.e. the same 'qualifying year' definition used in A.3 — a documented methodological choice, not silent filtering.
C.3: fitted signed-residual slope for 1 of 1 usable_temporal rectangles (remainder had <2 qualifying years).

## Section D — fishing pressure over time and space
D uses outputs/h2_couce_year_effort.rds directly (6353 rectangle-year rows, 215 rectangles, 1985-2015); rectangles with no Couce record are absent from this table by construction (the F exclusion).
D.2: 0 Couce-covered rectangles could not be matched to shapefile geometry and are excluded from the map.
D.3: fitted fishing-pressure slope for 0 of 0 usable_for_fishing_analysis rectangles (same slope method as C.3).

## Section E — variance decomposition
E restricted to the 0 rectangles flagged usable_for_fishing_analysis (usable_temporal AND has_couce_coverage) — the same rectangle set is used for fishing pressure AND both residual decompositions below, so the three ICC values are directly comparable.
E.1 fishing-pressure input: 0 rectangle-year rows across 0 rectangles (all available Couce years used, not restricted to qualifying survey years, since Couce coverage is a modelled reconstruction independent of haul sparsity).
E.2 residual input: 0 rectangle-year rows across 0 rectangles, restricted to qualifying years (n_hauls >= SPARSE_HAUL_THRESHOLD), consistent with A.3/C.3.
Saved variance decomposition summary (E.1-E.3 combined): /Users/stuartstokeld/north_sea_eeos/outputs/h3_pre_variance_decomposition.csv
ICC values (uninterpreted): 
  - fishing_pressure_hours_total: ICC = NA, n_rectangles = 0, sd_within/sd_between = NA
  - residual_magnitude: ICC = NA, n_rectangles = 0, sd_within/sd_between = NA
  - residual_signed: ICC = NA, n_rectangles = 0, sd_within/sd_between = NA

## Headline numbers (reported, not interpreted)
- Rectangles usable_temporal (A.3, the core feasibility denominator): 1 of 197
- Rectangles usable_temporal AND has_couce_coverage (usable_for_fishing_analysis, D.3/E denominator): 0 of 197
- usable_temporal rectangle(s): 43G1
- ICC (fishing_pressure_hours_total): NA (see note on n_rectangles above)
- ICC (residual_magnitude): NA
- ICC (residual_signed): NA

## Figure index
- **A.1** `h3_pre_A1_haul_count_per_year.png` — Haul count per year, 1985-2015, all Q1 NS-IBTS hauls (bar chart).
- **A.2** `h3_pre_A2_haul_year_heatmap.png` — Rectangle x year haul-count heatmap (primary visual for spatial-temporal coverage gaps); cell = haul count, 0-filled for years with no haul.
- **B.1** `h3_pre_B1_biomass_pred_vs_obs_timeseries.png` — Median log(B_pred) and log(B_obs) per year, with IQR ribbons.
- **C.1** `h3_pre_C1_residual_timeseries.png` — Mean signed residual (top) and mean residual magnitude (bottom) per year, with 95% CI; two panels.
- **C.2** `h3_pre_C2_residual_maps_by_decade.png` — Mean signed residual by ICES rectangle, one map per decade bin (1985-94 / 1995-2004 / 2005-15), common colour scale across all three.
- **C.3** `h3_pre_C3_residual_slope_map.png` — Rectangle-level OLS slope of mean signed residual vs year, usable_temporal rectangles only; direction + magnitude only, no p-values.
- **D.1** `h3_pre_D1_fishing_pressure_timeseries.png` — Mean Couce fishing hours per year across all rectangles with coverage, with IQR ribbon.
- **D.2** `h3_pre_D2_fishing_pressure_maps_by_decade.png` — Mean Couce fishing hours by ICES rectangle, one map per decade bin, common log1p colour scale across all three (same faceting logic as C.2).
- **D.3** `h3_pre_D3_fishing_pressure_slope_map.png` — Rectangle-level OLS slope of annual fishing hours vs year, usable_for_fishing_analysis rectangles only; direction + magnitude only, no p-values.
- **B.2** SKIPPED — no figure produced (see Section B note above).

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_pre_rectangle_usability_flags.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_pre_variance_decomposition.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_pre_exploration_run_log.md (this file)
- 9 figures in /Users/stuartstokeld/north_sea_eeos/outputs/figures (h3_pre_*.png)
