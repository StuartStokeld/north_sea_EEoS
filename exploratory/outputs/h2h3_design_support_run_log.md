# Quantitative support for the H2/H3 design decision — run log

Reports numbers only. No temporal phase boundary, spatial scheme, or resolution is recommended anywhere in this script or its outputs — that decision is made separately.

## Deviation from the brief's output format
The brief requests 'one CSV per section (A, B, C, D)'. Sections A, B and D each specify multiple genuinely different table shapes (e.g. A: a year x metric series, a phase-boundary-change table, a within-phase-trend table). Rather than force heterogeneous tables into one wide/long CSV (which would be harder to read and use than the source tables), each sub-table is saved as its own CSV, prefixed by section (A1/A2/A3, B1-B4), plus a single CSV for C (which genuinely is one table) and one for D. This mirrors the existing repo convention (see outputs/h3_pre_*.csv, outputs/h3_policy_*.csv — always one file per logical table, never force-merged). Flagged here as a deviation from the literal instruction, not silently done.

## Rectangle universes (reconciliation)
THREE distinct rectangle universes appear across this task's tables — none is more 'correct' than another, they answer different questions:
  (1) **215-rectangle full Couce universe** (`outputs/h2_couce_year_effort.rds`'s own distinct stat_rec set) — every rectangle Couce et al. (2020) provide a fishing-hours reconstruction for, REGARDLESS of whether NS-IBTS Q1 ever hauled there. This is the SAME universe used by the existing h3_pre_D1/D2 figures ('across all rectangles with coverage' — see outputs/exploratory_review/index.md Task 1 D.1/D.2 rows) and is the universe used in **Section A (A1-A3) and Section B (B1-B4)** below, since those sections describe the Couce fishing-pressure dataset on its own terms, extending those two existing figures directly.
  (2) **197-rectangle NS-IBTS haul-bearing universe** (`outputs/h3_pre_rectangle_usability_flags.csv`, all rows) — every rectangle with >=1 Q1 haul, 1985-2015, regardless of Couce coverage. Used to build Scheme A (2x2/3x3 block-merge) spatial units in **Section C and Section D** (Scheme A is blind to fishing pressure by construction, per run_h3_policy_zones.R, so it is built on the full haul-bearing universe, not just the Couce-covered subset).
  (3) **165-rectangle Couce-covered subset of the 197** (`has_couce_coverage` flag in the same CSV) — rectangles with BOTH >=1 Q1 haul AND >=1 Couce fishing-hours record. Used to build Scheme B (3-/4-tier pressure zones) in **Section C and Section D** (Scheme B needs a fishing-pressure value per rectangle to assign tiers, per run_h3_policy_zones.R), AND as the relaxed-rule individual-rectangle universe for **Section C item 1** (see below — the relaxed rule turns out to select exactly this same 165-rectangle set).
50 rectangles are in the full Couce universe (1) but NOT in the 197-rectangle haul-bearing universe (2) at all, i.e. Couce provides a fishing-hours reconstruction for water NS-IBTS Q1 never sampled in this period — consistent with the existing h3_pre run log's note that D uses outputs/h2_couce_year_effort.rds directly, unrestricted by the haul universe.

## Section A — whole-study-area temporal trend in fishing pressure
Universe: full Couce dataset (215 rectangles), matching h3_pre_D1 exactly. Phase boundaries (1985-1990, 1990-2000, 2000-2010, 2010-2015) are the visual/approximate boundaries named in the briefing (from h3_pre_D1/D3), NOT statistically derived — treated as provisional throughout this section. Boundary years are shared between adjacent phases (e.g. 1990 is the last year of phase 1 AND the first year of phase 2), a provisional modelling choice, not an error.
A1 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA1_year_fishing_summary.csv (31 year rows, 1985-2015; n_rect_contributing ranges 198-210 rectangles/year).
A2 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA2_phase_boundary_changes.csv (5 rows: 4 phase-boundary transitions + full 1985 vs 2015).
  - 1985-1990: mean hours 11926.87 -> 17160.84 (abs change +5233.97, +43.9%)
  - 1990-2000: mean hours 17160.84 -> 14968.94 (abs change -2191.90, -12.8%)
  - 2000-2010: mean hours 14968.94 -> 8714.51 (abs change -6254.42, -41.8%)
  - 2010-2015: mean hours 8714.51 -> 7206.00 (abs change -1508.51, -17.3%)
  - 1985-2015 (full period): mean hours 11926.87 -> 7206.00 (abs change -4720.87, -39.6%)
A3 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA3_phase_trends.csv (within-phase linear trends, one row per phase; NOT one 30-year trend).
  - 1985-1990 (n=6 yrs): slope=917.885 hours/yr (increasing, 6.56%/yr of phase mean), R^2=0.879
  - 1990-2000 (n=11 yrs): slope=-295.308 hours/yr (decreasing, -1.80%/yr of phase mean), R^2=0.849
  - 2000-2010 (n=11 yrs): slope=-756.547 hours/yr (decreasing, -6.79%/yr of phase mean), R^2=0.933
  - 2010-2015 (n=6 yrs): slope=-218.654 hours/yr (decreasing, -2.88%/yr of phase mean), R^2=0.493
Follow-up flag (NOT run here, per the briefing): a formal structural-break detection (e.g. `strucchange::breakpoints()` on the annual mean-hours series) could sharpen these four visual phase boundaries into statistically located break years, rather than the eyeballed 1990/2000/2010 boundaries used above. Flagged as a possible follow-up only.

## Section B — spatial heterogeneity of fishing pressure by decade
Universe: full Couce dataset (215 rectangles), matching h3_pre_D2 exactly. DECADE_BINS reused as-is from run_h3_pre_exploration.R: 1985-1994, 1995-2004, 2005-2015.
B1 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB1_rect_decade_wide.csv (215 rectangles x 3 decade columns [1985-1994, 1995-2004, 2005-2015]; NA = no Couce record for that rectangle in that decade).
B2 comparison 'decade2_minus_decade1': n_rectangles_used = 210 (both decades non-NA).
B2 comparison 'decade3_minus_decade2': n_rectangles_used = 212 (both decades non-NA).
B2 comparison 'decade3_minus_decade1': n_rectangles_used = 211 (both decades non-NA).
B2 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB2_rect_decade_change.csv (633 rows across 3 decade-pair comparisons, absolute + percent change).
B3 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB3_change_classification.csv — tercile classification of |decade3-decade1 absolute change| via dplyr::ntile(., 3) on 211 rectangles with both decades present. PROVISIONAL exact boundary values: 'stable' <= 3318.77 hours; 'large_change' >= 10261.00 hours (bottom/top tercile of |absolute change|, NOT a formal cutoff).
  table(change_class): stable=71, moderate=70, large_change=70
B4 saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB4_decade_persistence_correlation.csv — spatial persistence (Pearson + Spearman correlation of rectangle-level fishing pressure across decade pairs):
  - decade1_vs_decade2: n=210, Pearson r=0.716, Spearman rho=0.866
  - decade2_vs_decade3: n=212, Pearson r=0.876, Spearman rho=0.936
  - decade1_vs_decade3: n=211, Pearson r=0.470, Spearman rho=0.801

## Section C — ICC by spatial-unit definition (5 resolutions)
Relaxed inclusion rule for individual-rectangle resolution (#1): has_couce_coverage AND n_years_present >= 2 (n_years_present = years with >=1 haul, ANY count — NOT the stricter usable_for_fishing_analysis flag, which required >=10 years with >=5 hauls/year and returned 0 rectangles; see h3_pre_exploration_run_log.md). 165 of 197 rectangles qualify under this relaxed rule.
This relaxed rule turns out to select EXACTLY the same set as has_couce_coverage alone (165 of 197) — i.e. essentially every Couce-covered rectangle already has >=2 years of Q1 haul data (median n_years_present across all 197 rectangles = 31 of a possible 31), so the >=2-year bar adds no further restriction beyond Couce coverage itself. Reported plainly, not treated as a coincidence requiring adjustment.
Input data for ALL FIVE resolutions below is the SAME: 5112 rectangle-year fishing-pressure rows and 4799 rectangle-year residual-mean rows (165 rectangles), both restricted to the 165-rectangle relaxed universe. NOTE: unlike the pre-H3 script's Section E, residual rectangle-year means here use ALL years with >=1 haul (no SPARSE_HAUL_THRESHOLD >=5-hauls/year gate) — a deliberate further relaxation, consistent with the brief's relaxed-rule instruction. Only the GROUPING (unit_id) changes across the five rows below — this isolates the effect of spatial coarsening from any change in the underlying data.
Scheme A spatial units are constructed over the full 197-rectangle haul-bearing universe (blind to fishing pressure, matching run_h3_policy_zones.R); Scheme B over the 165-rectangle Couce-covered subset (needs a pressure value to assign tiers). Both then joined against the SAME 165-rectangle relaxed input data above — a block/zone that happens to span rectangles outside the relaxed universe simply contributes fewer observations, not zero (reported via n_units_with_*_data).

C saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designC_icc_by_spatial_unit.csv — ordered finest to coarsest:
  - 1. Individual ICES rectangle (relaxed rule): n_units_constructed=165, n_units_with_data(fish/resid)=165/165 | ICC fishing=0.5893, resid_mag=0.2274, resid_signed=0.2372
  - 2. Scheme A, 2x2 block: n_units_constructed=62, n_units_with_data(fish/resid)=50/50 | ICC fishing=0.4527, resid_mag=0.1708, resid_signed=0.1746
  - 3. Scheme A, 3x3 block: n_units_constructed=34, n_units_with_data(fish/resid)=26/26 | ICC fishing=0.3772, resid_mag=0.1365, resid_signed=0.1420
  - 4. Scheme B, 3-tier zone: n_units_constructed=22, n_units_with_data(fish/resid)=22/22 | ICC fishing=0.5304, resid_mag=0.0979, resid_signed=0.1039
  - 5. Scheme B, 4-tier zone: n_units_constructed=49, n_units_with_data(fish/resid)=49/49 | ICC fishing=0.5241, resid_mag=0.1359, resid_signed=0.1425

## Section D — haul-effort unevenness vs fishing pressure
Reuses the SAME Scheme A (2x2/3x3, built over the 197-rectangle haul-bearing universe) and Scheme B (3-/4-tier, built over the 165-rectangle Couce-covered universe) spatial units as Section C above (identical construction to run_h3_policy_zones.R), with the same POLICY_BREAK_YEAR = 2003 pre/post split reused as-is from that script.
D saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designD_haul_unevenness.csv:
  - Scheme A (block merge) / 2x2 / pre: n_units=62, CV=61.0%, Gini=0.342, max/min ratio (strict/excl-zero)=Inf/92.7, cor(fishing,haulcount) Pearson=-0.042 Spearman=0.076 (n=50)
  - Scheme A (block merge) / 2x2 / post: n_units=62, CV=50.2%, Gini=0.279, max/min ratio (strict/excl-zero)=Inf/22.4, cor(fishing,haulcount) Pearson=-0.362 Spearman=-0.142 (n=50)
  - Scheme A (block merge) / 3x3 / pre: n_units=34, CV=74.9%, Gini=0.420, max/min ratio (strict/excl-zero)=Inf/503.0, cor(fishing,haulcount) Pearson=0.112 Spearman=0.191 (n=26)
  - Scheme A (block merge) / 3x3 / post: n_units=34, CV=63.3%, Gini=0.356, max/min ratio (strict/excl-zero)=160.5/160.5, cor(fishing,haulcount) Pearson=-0.109 Spearman=0.005 (n=26)
  - Scheme B (pressure zones) / 3-tier / pre: n_units=22, CV=207.4%, Gini=0.760, max/min ratio (strict/excl-zero)=533.2/533.2, cor(fishing,haulcount) Pearson=0.056 Spearman=0.316 (n=22)
  - Scheme B (pressure zones) / 3-tier / post: n_units=22, CV=202.3%, Gini=0.745, max/min ratio (strict/excl-zero)=Inf/70.3, cor(fishing,haulcount) Pearson=0.102 Spearman=-0.124 (n=22)
  - Scheme B (pressure zones) / 4-tier / pre: n_units=49, CV=183.4%, Gini=0.648, max/min ratio (strict/excl-zero)=Inf/294.0, cor(fishing,haulcount) Pearson=0.038 Spearman=0.175 (n=49)
  - Scheme B (pressure zones) / 4-tier / post: n_units=49, CV=166.2%, Gini=0.613, max/min ratio (strict/excl-zero)=Inf/374.5, cor(fishing,haulcount) Pearson=0.143 Spearman=0.067 (n=49)
D4 (optional byproduct) figures saved: h2h3_designD4_A_2x2_haulcount_vs_pressure.png, h2h3_designD4_A_3x3_haulcount_vs_pressure.png, h2h3_designD4_B_3tier_haulcount_vs_pressure.png, h2h3_designD4_B_4tier_haulcount_vs_pressure.png

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA1_year_fishing_summary.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA2_phase_boundary_changes.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA3_phase_trends.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB1_rect_decade_wide.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB2_rect_decade_change.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB3_change_classification.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designB4_decade_persistence_correlation.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designC_icc_by_spatial_unit.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designD_haul_unevenness.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designD4_A_2x2_haulcount_vs_pressure.png
  - /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designD4_A_3x3_haulcount_vs_pressure.png
  - /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designD4_B_3tier_haulcount_vs_pressure.png
  - /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designD4_B_4tier_haulcount_vs_pressure.png (optional D4 byproduct figures)
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_design_support_run_log.md (this file)
