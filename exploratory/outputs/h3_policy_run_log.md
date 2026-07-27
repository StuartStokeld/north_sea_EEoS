# H3 strategy feasibility (policy-period split x coarse spatial zones) — run log

Feasibility/visualisation check only. No H3 model fit, no scheme/parameter selected as final, no significance testing.

## Provisional constants
- POLICY_BREAK_YEAR = 2003 (pre = 1985-2002 [18 yrs]; post = 2003-2015 [13 yrs]; midpoint of the 2002-2004 cod-crisis reform window; PROVISIONAL pending exact date confirmation with Jake/HH; periods are unbalanced by construction — reported, not corrected for.)
- BLOCK_SIZES = 2, 3 (Scheme A block_size x block_size rectangle merge)
- PRESSURE_TIERS = 3, 4 (Scheme B quantile-based tier counts, via dplyr::ntile on mean_annual_hours_total pooled over the full period)
- MIN_HAULS_PER_CELL = 5 (reused from SPARSE_HAUL_THRESHOLD in the pre-H3 feasibility task; EXPLORATORY REFERENCE POINT ONLY — retained in the feasibility summary for continuity with the pre-H3 task, NOT the primary analysis filter; see 'Sample-size adequacy' section below for the statistically defensible alternative)
- FIXED_CANDIDATE_THRESHOLDS = 1, 2, 3, 5, 10, 15, 20, 30, 40, 50 (sensitivity grid, 'explore beyond >=5 hauls')
- RELIABILITY_TARGETS = 0.7, 0.8, 0.9 (Spearman-Brown reliability targets for the cell mean; conventional psychometric benchmarks, not a new invention: ~0.7 'acceptable', 0.8 'good', 0.9 'excellent' per Nunnally 1978)

## Universe of rectangles
This task uses the SAME 197-rectangle haul-bearing universe as Section A of the pre-H3 feasibility script (outputs/h3_pre_rectangle_usability_flags.csv), NOT Step 0's 187-rectangle universe (which pre-filtered on total pooled hauls >= 5) — see that script's run log for the reconciliation between the two. Scheme A (blind to fishing pressure) uses the full 197-rectangle universe. Scheme B (fishing-pressure tiers) additionally requires Couce coverage, restricting to 165 of 197 rectangles (the same has_couce_coverage flag from the pre-H3 task).
Rectangles without matching ICES shapefile geometry (excluded from block/zone construction): 0 of 197 (Scheme A universe); 0 of 165 (Scheme B universe).

## Scheme A — geographic block merge
Building Scheme A blocks at 2x2...
  197 rectangles merged into 62 blocks (sizes: 1-4 rectangles/block).
A_2x2: 55 of 62 units have data in both periods for mean residual; 50 of 62 for mean fishing hours (units missing a period, e.g. zero contributing rectangles that period, excluded from the scatter, not imputed).
  Sample-size adequacy (Scheme A (block merge) / 2x2): n_cells=117, n_hauls=12069, ICC=0.1774 (haul-level residual variance between-cell vs within-cell); required n/cell for reliability 0.7/0.8/0.9 -> 11/19/42.
Building Scheme A blocks at 3x3...
  197 rectangles merged into 34 blocks (sizes: 1-9 rectangles/block).
A_3x3: 30 of 34 units have data in both periods for mean residual; 26 of 34 for mean fishing hours (units missing a period, e.g. zero contributing rectangles that period, excluded from the scatter, not imputed).
  Sample-size adequacy (Scheme A (block merge) / 3x3): n_cells=64, n_hauls=12069, ICC=0.1513 (haul-level residual variance between-cell vs within-cell); required n/cell for reliability 0.7/0.8/0.9 -> 14/23/51.

## Scheme B — contiguous fishing-pressure zones
Method: (1) classify each of the 165 Couce-covered rectangles into a quantile-based pressure tier (dplyr::ntile on mean_annual_hours_total, pooled 1985-2015, from outputs/h2_couce_rectangle_effort.rds); (2) build a rook (edge-sharing, not diagonal) adjacency list directly from ICES-grid row/col indices (derived from the shapefile's own SOUTH/WEST corner fields — no spdep/igraph dependency needed, since grid adjacency is directly computable from integer row/col indices); (3) find connected components via a plain breadth-first search restricted to same-tier neighbours (R/h3_policy_zones_helpers.R::build_contiguous_zones()). This was straightforward to implement with the grid-index approach (no igraph/spdep required).
Building Scheme B zones at 3 tiers...
  3-tier split: table(tier) = tier1=55, tier2=55, tier3=55; resulting in 22 contiguous zones (14 singleton [n=1 rectangle] zones — same-tier rectangles with no same-tier rook-neighbour; NOT pooled with other same-tier-but-non-adjacent rectangles). Largest zone: 52 rectangles.
B_3tier: 21 of 22 units have data in both periods for mean residual; 22 of 22 for mean fishing hours (units missing a period, e.g. zero contributing rectangles that period, excluded from the scatter, not imputed).
  Sample-size adequacy (Scheme B (pressure zones) / 3-tier): n_cells=43, n_hauls=10498, ICC=0.0975 (haul-level residual variance between-cell vs within-cell); required n/cell for reliability 0.7/0.8/0.9 -> 22/38/84.
Building Scheme B zones at 4 tiers...
  4-tier split: table(tier) = tier1=42, tier2=41, tier3=41, tier4=41; resulting in 49 contiguous zones (31 singleton [n=1 rectangle] zones — same-tier rectangles with no same-tier rook-neighbour; NOT pooled with other same-tier-but-non-adjacent rectangles). Largest zone: 29 rectangles.
B_4tier: 45 of 49 units have data in both periods for mean residual; 49 of 49 for mean fishing hours (units missing a period, e.g. zero contributing rectangles that period, excluded from the scatter, not imputed).
  Sample-size adequacy (Scheme B (pressure zones) / 4-tier): n_cells=94, n_hauls=10498, ICC=0.1250 (haul-level residual variance between-cell vs within-cell); required n/cell for reliability 0.7/0.8/0.9 -> 17/28/63.

## Feasibility summary
Saved feasibility summary: /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_feasibility_summary.csv
  - Scheme A (block merge) / 2x2: n_units=62, n_units_usable_both_periods=55 (>= 5 hauls in both periods); pre mean/median/min = 114.9/141.5/0; post mean/median/min = 80.5/94.5/0
  - Scheme A (block merge) / 3x3: n_units=34, n_units_usable_both_periods=29 (>= 5 hauls in both periods); pre mean/median/min = 209.5/182.0/0; post mean/median/min = 146.9/151.0/2
  - Scheme B (pressure zones) / 3-tier: n_units=22, n_units_usable_both_periods=21 (>= 5 hauls in both periods); pre mean/median/min = 287.6/41.0/4; post mean/median/min = 191.6/30.5/0
  - Scheme B (pressure zones) / 4-tier: n_units=49, n_units_usable_both_periods=45 (>= 5 hauls in both periods); pre mean/median/min = 129.1/40.0/0; post mean/median/min = 86.0/27.0/0

Best-performing combination by n_units_usable_both_periods (at the exploratory MIN_HAULS_PER_CELL=5 reference point): Scheme A (block merge) / 2x2 (55 of 62 units, 88.7%). Reported, not a recommendation — comparison material for the supervisor discussion.

## Sample-size adequacy for MIN_HAULS_PER_CELL (statistically defensible thresholds)
The exploratory MIN_HAULS_PER_CELL = 5 gate above is a reference point only, not the primary feasibility filter. Method: for each scheme/parameter combination, label every haul with its cell (spatial unit x pre/post period), then decompose the haul-level residual's variance into between-cell and within-cell (haul-to-haul) components (one-way random-effects ANOVA, method-of-moments estimator for unbalanced groups — R/h3_pre_exploration_helpers.R::variance_components_anova(), same estimator used in the pre-H3 task's Section E). The resulting ICC (intraclass correlation = between-cell variance / total variance) is the fraction of haul-to-haul residual variation that reflects genuine between-cell signal rather than sampling noise. Inverting the Spearman-Brown formula (reliability of an n-haul cell mean = n*ICC / (1+(n-1)*ICC)) gives the number of hauls a cell needs for its mean to reach a target reliability — this is the same logic used to decide how many test items or raters are 'enough' in measurement theory, applied here to deciding how many hauls are 'enough' for one cell mean. If ICC <= 0 (no detectable between-cell signal at all), required n is reported as Inf — no amount of extra replication manufactures a signal that isn't there.
Saved sample-size adequacy table: /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_sample_size_adequacy.csv
  - Scheme A (block merge) / 2x2: n_cells=117, n_hauls=12069, ICC=0.1774 (var_between=0.0424, var_within=0.1968); R0.7->n=11, R0.8->n=19, R0.9->n=42
  - Scheme A (block merge) / 3x3: n_cells=64, n_hauls=12069, ICC=0.1513 (var_between=0.0362, var_within=0.2034); R0.7->n=14, R0.8->n=23, R0.9->n=51
  - Scheme B (pressure zones) / 3-tier: n_cells=43, n_hauls=10498, ICC=0.0975 (var_between=0.0236, var_within=0.2181); R0.7->n=22, R0.8->n=38, R0.9->n=84
  - Scheme B (pressure zones) / 4-tier: n_cells=94, n_hauls=10498, ICC=0.1250 (var_between=0.0300, var_within=0.2102); R0.7->n=17, R0.8->n=28, R0.9->n=63
Saved feasibility-by-threshold sensitivity table: /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_feasibility_by_threshold.csv
This table repeats the n_units_usable_both_periods calculation across the fixed candidate grid (1, 2, 3, 5, 10, 15, 20, 30, 40, 50) AND each combo's own reliability-derived threshold(s) — the full 'explore beyond >=5 hauls' picture, not a single number.

## Secondary/robustness variant (symmetric trimmed window)
Primary split shows promise (best combination has 88.7% of units usable in both periods, >= 50% gate) — building the symmetric trimmed-window variant for the best-performing combination only.
Symmetric window: 13 years either side of 2003 (pre = 1990-2002, post = 2003-2015) — NOT built as a full figure set here; flagged as the next step if the supervisor wants it pursued (out of scope to build a second full figure set speculatively beyond confirming the primary split's promise).

## Figure index
- **A_2x2_units_map** `h3_policy_A_2x2_units_map.png` — Scheme A (block merge) / 2x2 spatial units (n=62), for visual sanity-check before any statistic is computed on them.
- **A_2x2_haulcount** `h3_policy_A_2x2_haulcount_prepost.png` — Haul count per 2x2 spatial unit, pre-2003 vs post-2003, side by side.
- **A_2x2_prepost** `h3_policy_A_2x2_prepost_comparison.png` — Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 2x2 spatial unit; descriptive only.
- **A_3x3_units_map** `h3_policy_A_3x3_units_map.png` — Scheme A (block merge) / 3x3 spatial units (n=34), for visual sanity-check before any statistic is computed on them.
- **A_3x3_haulcount** `h3_policy_A_3x3_haulcount_prepost.png` — Haul count per 3x3 spatial unit, pre-2003 vs post-2003, side by side.
- **A_3x3_prepost** `h3_policy_A_3x3_prepost_comparison.png` — Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 3x3 spatial unit; descriptive only.
- **B_3tier_units_map** `h3_policy_B_3tier_units_map.png` — Scheme B (pressure zones) / 3-tier spatial units (n=22), for visual sanity-check before any statistic is computed on them.
- **B_3tier_haulcount** `h3_policy_B_3tier_haulcount_prepost.png` — Haul count per 3-tier spatial unit, pre-2003 vs post-2003, side by side.
- **B_3tier_prepost** `h3_policy_B_3tier_prepost_comparison.png` — Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 3-tier spatial unit; descriptive only.
- **B_4tier_units_map** `h3_policy_B_4tier_units_map.png` — Scheme B (pressure zones) / 4-tier spatial units (n=49), for visual sanity-check before any statistic is computed on them.
- **B_4tier_haulcount** `h3_policy_B_4tier_haulcount_prepost.png` — Haul count per 4-tier spatial unit, pre-2003 vs post-2003, side by side.
- **B_4tier_prepost** `h3_policy_B_4tier_prepost_comparison.png` — Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 4-tier spatial unit; descriptive only.
- **sample_size_sensitivity** `h3_policy_sample_size_sensitivity.png` — Percent of spatial units usable in both periods vs haul-count threshold, one line per scheme/parameter, with each combo's reliability-derived thresholds marked.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_feasibility_summary.csv (single-threshold reference table, MIN_HAULS_PER_CELL=5 only — exploratory)
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_sample_size_adequacy.csv (variance decomposition + reliability-derived required-n per combo — statistically grounded)
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_feasibility_by_threshold.csv (feasibility sensitivity across the full threshold grid — the primary deliverable for this question)
- /Users/stuartstokeld/north_sea_eeos/outputs/h3_policy_run_log.md (this file)
- 13 figures in /Users/stuartstokeld/north_sea_eeos/outputs/figures (h3_policy_*.png)
