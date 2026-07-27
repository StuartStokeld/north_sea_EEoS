# Exploratory outputs — master index

Collation only. No synthesis, comparison, or conclusions across tasks are drawn here —
that happens separately, in the review conversation. Every caption/description below is
quoted or paraphrased directly from the source task's own code/docs/run log; nothing here
is a new interpretation. Files are **copied** (not moved) from their original locations in
`outputs/` / `outputs/figures/` / `display_discussion/`, so nothing elsewhere is broken by
this collation.

See `../../COMPLETION_NOTE.md` (i.e. `outputs/exploratory_review_completion_note.md`) for
counts and gaps.

---

## Task 1 — Pre-H3 exploratory visualisation & variance decomposition

**What it covered:** Exploratory/descriptive check of the spatial and temporal
distribution of hauls, predicted biomass, EEoS residuals, and Couce fishing pressure,
to assess whether fishing pressure varies enough *within* ICES rectangles across years to
give H3 statistical traction (H2 already used between-rectangle variance over the full
30-year period). No H3 model fit, no H3 conclusions drawn.
**Run log:** `run_logs/h3_pre_exploration_run_log.md`

## Task 2 — H1 thesis figures

**What it covered:** The core H1 haul-level test of whether the Ecological Equation of
State (EEoS) predicts community biomass from state variables (S, N, E) — the two direct
reproductions of Harte et al. (2022)'s validation figures (Test 1, Test 2), the extended
fitted ln(E) comparison, the null-permutation check, and two parallel diagnostic tracks
(catchability scaling by biomass quartile; haul-level dominance/size-homogeneity as a
second axis of prediction failure).
**Run log / companion docs:** no single run log; `run_logs/h1_catchability_scaling_exploration.md`
and `run_logs/h1_dominance_size_homogeneity_exploration.md` are auto-generated companion
logs for the two diagnostic tracks. Headline H1 methodology/results are narrated in
`display_discussion/H1_methods_draft.md` and `display_discussion/H1_results_summary.Rmd`
(not copied here — source docs, not run logs; referenced for caption provenance only).

## Task 3 — Core visual suite

**What it covered (per the task brief):** haul map, temporal coverage, composition/
big-shoal figures, H1 residual map, reframed H2/H3 strategy exploration.
**Status: two of the five items were generated in a follow-up task** ("Missing Visuals —
Haul Map & Composition Figures") after this collation first flagged them as gaps: a
standalone geographic haul map (`haul_map_1`), and the plain S/D/size_CV composition
distributions (`composition_1`/`2`/`3`) — as opposed to the residual-ratio overlays
already indexed under Task 2 (`h1_dominance_ratio_vs_D.png`,
`h1_dominance_ratio_vs_sizeCV.png`), which remain distinct and are not duplicated below.
The remaining three items are still gaps: a standalone "H1 residual map" and "temporal
coverage" figure (distinct from Task 1's A.1/A.2, already indexed there) could not be
found anywhere on disk; "reframed H2/H3 strategy exploration" most plausibly refers to
Task 1 and/or Task 4 themselves (already indexed separately), not a distinct 5th figure
set. The H2 residual map was explicitly deferred (out of scope for the follow-up task).
**Run log:** `run_logs/task3_missing_visuals_run_log.md`

## Task 4 — Policy-period × spatial-zone feasibility check (with sample-size-adequacy extension)

**What it covered:** Whether coarsening both temporal (pre/post 2003 policy-break split)
and spatial (2×2/3×3 rectangle-block merge, or 3-/4-tier contiguous fishing-pressure
zones) resolution gives H3 more statistical traction than the annual within-rectangle
design tested in Task 1. Extended with a reliability-based sample-size adequacy analysis
(ICC / Spearman-Brown-derived required-n) that supersedes the original exploratory
"≥5 hauls in both periods" filter as the primary feasibility filter.
**Run log:** `run_logs/h3_policy_run_log.md`

---

## Index table

| Task | Figure/Stats ID | File path (within `exploratory_review/`) | Caption / description | Notes |
|---|---|---|---|---|
| 1 | A.1 | `figures/h3_pre_A1_haul_count_by_year.png` | Haul count per year, 1985-2015, all Q1 NS-IBTS hauls (bar chart). | universe = 197 rectangles, 12,117 Q1 hauls |
| 1 | A.2 | `figures/h3_pre_A2_rect_year_heatmap.png` | Rectangle x year haul-count heatmap (primary visual for spatial-temporal coverage gaps); cell = haul count, 0-filled for years with no haul. | median hauls/rectangle/year = 2, 90th pct = 3, 95th pct = 4, max = 12 |
| 1 | B.1 | `figures/h3_pre_B1_biomass_pred_vs_obs_timeseries.png` | Median log(B_pred) and log(B_obs) per year, with IQR ribbons. | 12,069 of 12,117 hauls (99.6%) have a successful EEoS prediction |
| 1 | B.2 | — not produced (no file) | N/A | B.2 region facet skipped, sub_area 100% NA; continent/stratum/season also unusable for grouping |
| 1 | C.1 | `figures/h3_pre_C1_residual_timeseries.png` | Mean signed residual (top) and mean residual magnitude (bottom) per year, with 95% CI; two panels. | uses `resid_signed = log(B_pred) - log(B_obs)`, the *negative* of the H1/H2 pipeline's primary `residual` convention — flagged, not resolved |
| 1 | C.2 | `figures/h3_pre_C2_residual_maps_by_decade.png` | Mean signed residual by ICES rectangle, one map per decade bin (1985-94 / 1995-2004 / 2005-15), common colour scale across all three. | 0 rows excluded for falling outside decade bins; 0 rectangles unmatched to shapefile geometry |
| 1 | C.3 | `figures/h3_pre_C3_residual_slope_map.png` | Rectangle-level OLS slope of mean signed residual vs year, usable_temporal rectangles only; direction + magnitude only, no p-values. | slope fitted for 1 of 1 usable_temporal rectangles (43G1) |
| 1 | D.1 | `figures/h3_pre_D1_fishing_pressure_timeseries.png` | Mean Couce fishing hours per year across all rectangles with coverage, with IQR ribbon. | 6,353 rectangle-year rows, 215 rectangles, 1985-2015 |
| 1 | D.2 | `figures/h3_pre_D2_fishing_pressure_maps_by_decade.png` | Mean Couce fishing hours by ICES rectangle, one map per decade bin, common log1p colour scale across all three (same faceting logic as C.2). | 0 Couce-covered rectangles unmatched to shapefile geometry |
| 1 | D.3 | `figures/h3_pre_D3_fishing_pressure_slope_map.png` | Rectangle-level OLS slope of annual fishing hours vs year, usable_for_fishing_analysis rectangles only; direction + magnitude only, no p-values. | fitted for 0 of 0 usable_for_fishing_analysis rectangles — empty map |
| 1 | rectangle usability flags | `stats/h3_pre_rectangle_usability_flags.csv` | Per-rectangle haul counts, qualifying-year counts, `usable_temporal` and `has_couce_coverage`/`usable_for_fishing_analysis` flags (A.3/F). | feasibility: 1/197 rectangles qualify as usable_temporal; 0/197 usable_for_fishing_analysis; 32/197 lack Couce coverage |
| 1 | variance decomposition | `stats/h3_pre_variance_decomposition.csv` | E.1-E.3 combined: between/within-rectangle variance components and ICC for fishing pressure, residual magnitude, and signed residual, restricted to usable_for_fishing_analysis rectangles. | ICC = NA for all three metrics (0 usable_for_fishing_analysis rectangles — denominator is empty) |
| 1 | run log | `run_logs/h3_pre_exploration_run_log.md` | Full run log: provisional constants, sign-convention flag, exclusion counts, headline numbers, figure index. | universe = 197 rectangles |
| 2 | Test 1 | `figures/h1_test1_harte_fig1_logB_pred_vs_obs.png` | "Harte Fig 1: EEoS predicted vs observed biomass" (plot title, `h1_harte_baseline_helpers.R::save_harte_baseline_figures()`). | log_r2 = −0.223 (Rmd-reported headline number); unfitted, parameter-free reproduction of Harte et al. (2022) Fig 1 |
| 2 | Test 1 baseline | `figures/h1_productivity_1to1_calibrated.png` | "Unfitted productivity 1:1: log(E × m_min) vs log(B_obs)" (plot title, same helper). | headline Test-1 baseline; log_r2 = 0.736 (Rmd-reported) |
| 2 | Test 1 diagnostic | `figures/h1_productivity_1to1_uncalibrated_diagnostic.png` | "Diagnostic (uncalibrated): log(E_raw) vs log(B_obs)" (plot title, same helper). | diagnostic only, retained as `model_id = productivity_1to1_uncalibrated`; NOT the Test-1 headline baseline |
| 2 | Test 2 | `figures/h1_test2_harte_fig2_productivity_ratio.png` | "Harte Fig 2: productivity ratio (sqrt display)" (plot title, `h1_harte_baseline_helpers.R::plot_harte_fig2()`). | unfitted reproduction of Harte et al. (2022) Fig 2 (E/B^0.75 ratio); Harte's own reference value = 0.600 |
| 2 | extended comparison | `figures/h1_extended_lnE_ols_reference.png` | "ln(E) OLS (fitted correlative): ln(B_obs) ~ ln(E)" (plot title, `run_h1_lne_reference.R`). | extended comparison only, reported after the two unfitted tests; fitted model, not on the same footing as Test 1/2 |
| 2 | null model (primary scheme) | `figures/h1_null_r2_uniform_log_b95.png` | "Null distribution of R² (log scale): primary: uniform_log_b_95" (plot title, `h1_null_helpers.R::plot_null_r2_distribution()`); corresponds to the "Uniform on 95% log(B_obs)" scheme in `H1_results_summary.Rmd`'s null table. | confirmed primary via `outputs/null_sampling_decision.rds` (`primary_method = "uniform_log_b_95"`, chosen because log(B_obs) distribution is peaked); null median log_r2 = −2.02 per Rmd table |
| 2 | null model (robustness scheme) | `figures/h1_null_r2_b_shuffle.png` | "Null distribution of R² (log scale): robustness: b_shuffle" (same generator); corresponds to "B shuffle" scheme in the Rmd's null table. | confirmed robustness scheme via `outputs/null_sampling_decision.rds` (`robustness_method = "b_shuffle"`); null median log_r2 = −2.60 per Rmd table; both null p < 0.001 vs observed log_r2 |
| 2 | catchability exploration | `figures/h1_catchability_ratio_vs_bobs.png` | "Ratio vs log B_obs" (image alt text, `H1_catchability_scaling_exploration.md`). | exploratory only; no correction applied to the pipeline |
| 2 | dominance figure (D) | `figures/h1_dominance_ratio_vs_D.png` | "Figure 4a legend: log(B_pred/B_obs) vs Berger-Parker dominance D; loess smooth; N = 12,069 hauls." (`display_discussion/H1_H2_methods_results_draft.md`, L110). | diagnostic only, parallel to the B_obs-quartile magnitude-bias result; no correction applied |
| 2 | dominance figure (size_CV) | `figures/h1_dominance_ratio_vs_sizeCV.png` | "Figure 4b legend: log(B_pred/B_obs) vs dominant-species size_CV; loess smooth." (`display_discussion/H1_H2_methods_results_draft.md`, L114). | diagnostic only, parallel to the B_obs-quartile magnitude-bias result; no correction applied |
| 2 | biomass-quartile gradient | `stats/h1_catchability_by_quartile.csv` | Median/mean B_pred/B_obs ratio by observed-biomass quartile (columns: `b_quartile, n, ln_B_obs_min, ln_B_obs_max, median_ratio, mean_ln_ratio`). | exists only as a table/CSV, not a standalone figure — see gap note in completion note |
| 2 | Test 1/2/extended metrics | `stats/h1_harte_baseline_metrics.csv` | Per-model metrics (log_r2, cor2, log_rmse, raw_rmse, median_ratio, Fig-2 Pearson R² all/trimmed) for the Harte-baseline models. | — |
| 2 | full model comparison | `stats/h1_model_comparison.csv` | All H1 model comparisons (EEoS, productivity 1:1 calibrated/uncalibrated, productivity ratio, ln(E) OLS): fit stats, SS ratios, Harte-criterion flags. | — |
| 2 | ln(E) OLS coefficients | `stats/h1_lne_coefficients.csv` | Fitted `lm(log B_obs ~ log E)` coefficient table (term, Estimate, Std. Error, t, p). | extended-comparison model only |
| 2 | Fig 2 leverage diagnostic | `stats/h1_fig2_leverage_hauls.csv` | High-leverage hauls excluded for the trimmed Test-2 Pearson R² (haul_id, S, N, E, B_obs, B_pred, predicted/observed ratio, leverage rank). | scoped to the Test-2 ratio-space regression only; not dropped from the main H1 dataset |
| 2 | dropout summary | `stats/h1_dropout_summary.csv` | Haul counts by pipeline stage (columns: `stage, n_hauls`). | data description |
| 2 | dropout by reason | `stats/h1_dropout_by_reason.csv` | Haul counts by exclusion reason (columns: `exclusion_reason, n_hauls`). | data description |
| 2 | dropout by rectangle | `stats/h1_dropout_by_stat_rec.csv` | Join/prediction/dropout counts per ICES rectangle (columns: `stat_rec, n_joined, n_predicted, n_dropped, pct_dropped`). | data description |
| 2 | dropout by year | `stats/h1_dropout_by_year.csv` | Join/prediction/dropout counts per year (columns: `year, n_joined, n_predicted, n_dropped, pct_dropped`). | data description |
| 2 | dropout funnel by year | `stats/h1_dropout_funnel_by_year.csv` | Full pipeline funnel per year: HL hauls -> DATRAS state -> FishGlob join -> EEoS predictions, with drop counts at each stage. | data description |
| 2 | dropout spike years | `stats/h1_dropout_spike_years.csv` | Same funnel columns as above, restricted to years flagged as dropout spikes, plus `spike_status`. | data description |
| 2 | filter exclusions | `stats/h1_filter_exclusions.csv` | Per-haul exclusion audit trail (haul_id, join_key, year, stat_rec, S, N, E, B_obs, exclusion_reason). | data description |
| 2 | join gaps | `stats/h1_join_gaps.csv` | Per-record join-gap audit trail (gap_side, join_key, year, B_obs, stat_rec, S, N, E, E_raw, m_min). | data description |
| 2 | catchability scaling summary | `stats/h1_catchability_scaling_summary.csv` | Key-value summary (n_hauls, median ratio, mean log-ratio, power-law slope/intercept, log_r2 unadjusted/after-division, interpretation flag). | exploratory only; no correction applied |
| 2 | dominance by D (bins) | `stats/h1_dominance_by_D_bins.csv` | Quartile/decile binned median B_pred/B_obs ratio by Berger-Parker dominance D. | — |
| 2 | dominance by D x B_obs quartile | `stats/h1_dominance_by_D_x_bobs_quartile.csv` | Same D-binned ratio table, nested within each B_obs quartile. | — |
| 2 | dominance by size_CV (bins) | `stats/h1_dominance_by_sizeCV_bins.csv` | Quartile/decile binned median ratio by dominant-species size-homogeneity (size_CV). | — |
| 2 | dominance by size_CV x B_obs quartile | `stats/h1_dominance_by_sizeCV_x_bobs_quartile.csv` | Same size_CV-binned ratio table, nested within each B_obs quartile. | — |
| 2 | dominance conditional correlation | `stats/h1_dominance_conditional_correlation.csv` | Pearson/Spearman r of D and size_CV with ln_ratio, computed within each B_obs quartile. | — |
| 2 | dominance correlation matrix | `stats/h1_dominance_correlation_matrix.csv` | Pairwise Pearson/Spearman r and |r|>0.4 flag across D, size_CV, N_haul, B_obs, n_bins_dominant_species. | only 2 of 10 pairs exceed |r|>0.4, neither involving D or size_CV vs B_obs/N |
| 2 | dominance data quality by rectangle | `stats/h1_dominance_dataquality_by_stat_rec.csv` | Low-size_CV-decile share per rectangle vs baseline share, with `flagged` indicator. | — |
| 2 | dominance data quality by year | `stats/h1_dominance_dataquality_by_year.csv` | Same cross-tab by year. | — |
| 2 | dominance data quality summary | `stats/h1_dominance_dataquality_summary.csv` | Rolled-up flagged/unflagged share of the low-size_CV decile, by grouping. | flagged years/rectangles show ~no enrichment (+0.3pp / −1.7pp) in the low-size_CV decile |
| 2 | dominance haul table | `stats/h1_dominance_haul_table.csv` | Full per-haul table: haul_id, D, dominant species, size_CV, n_bins_dominant_species, N, B_obs, B_pred, ratio, residual, B_obs_quartile. | 12,069 hauls, matched back to haul_key with zero loss |
| 2 | dominance partial R2 (ANOVA) | `stats/h1_dominance_partial_r2_anova.csv` | Nested-model ANOVA table (res_df, RSS, df, sum_of_sq, F, p) for the dominance partial-R² check. | — |
| 2 | dominance partial R2 coefficients | `stats/h1_dominance_partial_r2_coefficients.csv` | Coefficient table for the dominance partial-R² models. | — |
| 2 | dominance partial R2 summary | `stats/h1_dominance_partial_r2_summary.csv` | R²/adj-R²/n per model in the dominance partial-R² comparison. | — |
| 2 | dominance taxonomic breakdown | `stats/h1_dominance_taxonomic_breakdown.csv` | Dominant-species identity and combined share, by metric extreme decile (descriptive only). | descriptive only, not used to construct any species-based flag/filter |
| 2 | null distribution summary | `stats/h1_null_distribution_summary.csv` | Distribution-shape diagnostics on log(B_obs) used to choose the primary null scheme (median, IQR, range, skewness, kurtosis, histogram peak ratio, 95% log bounds). | — |
| 2 | pipeline audit results | `stats/h1_pipeline_audit_results.csv` | Named data-integrity checks with pass/fail status and notes (check_id, description, status, notes). | data description |
| 2 | DATRAS build diagnostics | `stats/h1_datras_build_diagnostics.csv` | Key-value diagnostics from the DATRAS HL/state build step. | data description |
| 2 | DATRAS length-code audit | `stats/h1_datras_lngt_code_audit.csv` | Length-code (LngtCode) audit: row counts and example bin midpoints per code, flagged as known/unknown. | data description |
| 2 | run log (catchability) | `run_logs/h1_catchability_scaling_exploration.md` | Auto-generated exploration doc/log for the catchability-scaling diagnostic (exploratory only, no correction applied). | — |
| 2 | run log (dominance) | `run_logs/h1_dominance_size_homogeneity_exploration.md` | Auto-generated exploration doc/log for the dominance/size-homogeneity diagnostic. | — |
| 3 | haul_map_1 | `figures/task3_haul_map_1_total_hauls_per_rectangle.png` | "haul_map_1: Total haul count per ICES rectangle, 1985-2015" (plot title, `run_missing_visuals.R`); choropleth of total Q1 NS-IBTS haul count per rectangle, full period; the 32 rectangles without Couce et al. (2020) fishing-pressure coverage are outlined in red. | shapefile: `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp` via `load_ices_rectangles_sf()`; all 197 of 197 usability-flags-table rectangles matched to shapefile geometry |
| 3 | composition_1 | `figures/task3_composition_1_S_D_sizeCV_histograms.png` | "composition_1: Distribution of S, D, and size_CV across hauls" (plot title, `run_missing_visuals.R`); histograms of species richness (S), Berger-Parker dominance (D), and dominant-species size_CV. | N = 12,069 hauls; plain distributions, not the ratio overlays already indexed under Task 2 |
| 3 | composition_2 | `figures/task3_composition_2_D_vs_sizeCV_scatter.png` | "composition_2: D vs size_CV, plain distribution" (plot title, `run_missing_visuals.R`); plain scatter of D vs size_CV with density contours, no ratio/residual colouring. | N = 12,069 hauls; distinct from Task 2's `h1_dominance_ratio_vs_D.png`/`h1_dominance_ratio_vs_sizeCV.png`, which plot log(B_pred/B_obs) against these same axes |
| 3 | composition_3 | `figures/task3_composition_3_big_shoal_reference_share.png` | "composition_3: Share of hauls above the D top-decile reference line" (plot title, `run_missing_visuals.R`); D vs size_CV scatter with hauls at/above the D top-decile boundary highlighted. | no formal D x size_CV "big shoal" cutoff was found anywhere in the H1 dominance analysis (checked `explore_h1_haul_dominance.R` and all `display_discussion/*.md`); reused the pre-existing top-decile boundary from `stats/h1_dominance_by_D_bins.csv` (D >= 0.9006, n=1,206, 10.0% of hauls) as a labelled reference line only, per the task brief's fallback instruction — not a formal classification |
| 4 | Scheme A 2x2 units map | `figures/policy_zone_schemeA_2x2_units_map.png` | "Scheme A (block merge) / 2x2 spatial units (n=62), for visual sanity-check before any statistic is computed on them." | universe = 197 rectangles; 62 blocks (1-4 rectangles/block) |
| 4 | Scheme A 2x2 haul count pre/post | `figures/policy_zone_schemeA_2x2_haulcount_prepost.png` | "Haul count per 2x2 spatial unit, pre-2003 vs post-2003, side by side." | — |
| 4 | Scheme A 2x2 pre/post comparison | `figures/policy_zone_schemeA_2x2_prepost_comparison.png` | "Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 2x2 spatial unit; descriptive only." | 55 of 62 units have data in both periods for mean residual; 50 of 62 for mean fishing hours |
| 4 | Scheme A 3x3 units map | `figures/policy_zone_schemeA_3x3_units_map.png` | "Scheme A (block merge) / 3x3 spatial units (n=34), for visual sanity-check before any statistic is computed on them." | 34 blocks (1-9 rectangles/block) |
| 4 | Scheme A 3x3 haul count pre/post | `figures/policy_zone_schemeA_3x3_haulcount_prepost.png` | "Haul count per 3x3 spatial unit, pre-2003 vs post-2003, side by side." | — |
| 4 | Scheme A 3x3 pre/post comparison | `figures/policy_zone_schemeA_3x3_prepost_comparison.png` | "Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 3x3 spatial unit; descriptive only." | 30 of 34 units have data in both periods for mean residual; 26 of 34 for mean fishing hours |
| 4 | Scheme B 3-tier units map | `figures/policy_zone_schemeB_3tier_units_map.png` | "Scheme B (pressure zones) / 3-tier spatial units (n=22), for visual sanity-check before any statistic is computed on them." | universe = 165 Couce-covered rectangles; 22 contiguous zones (14 singleton zones) |
| 4 | Scheme B 3-tier haul count pre/post | `figures/policy_zone_schemeB_3tier_haulcount_prepost.png` | "Haul count per 3-tier spatial unit, pre-2003 vs post-2003, side by side." | — |
| 4 | Scheme B 3-tier pre/post comparison | `figures/policy_zone_schemeB_3tier_prepost_comparison.png` | "Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 3-tier spatial unit; descriptive only." | 21 of 22 units have data in both periods for mean residual; 22 of 22 for mean fishing hours |
| 4 | Scheme B 4-tier units map | `figures/policy_zone_schemeB_4tier_units_map.png` | "Scheme B (pressure zones) / 4-tier spatial units (n=49), for visual sanity-check before any statistic is computed on them." | 49 contiguous zones (31 singleton zones) |
| 4 | Scheme B 4-tier haul count pre/post | `figures/policy_zone_schemeB_4tier_haulcount_prepost.png` | "Haul count per 4-tier spatial unit, pre-2003 vs post-2003, side by side." | — |
| 4 | Scheme B 4-tier pre/post comparison | `figures/policy_zone_schemeB_4tier_prepost_comparison.png` | "Pre-2003 vs post-2003 mean residual and mean fishing hours, one point per 4-tier spatial unit; descriptive only." | 45 of 49 units have data in both periods for mean residual; 49 of 49 for mean fishing hours |
| 4 | sample-size sensitivity | `figures/policy_zone_sample_size_sensitivity.png` | "Percent of spatial units usable in both periods vs haul-count threshold, one line per scheme/parameter, with each combo's reliability-derived thresholds marked." | ICC 0.10-0.18 across the 4 combos; required-n 11-84 hauls/cell depending on combo and reliability target (0.7/0.8/0.9) |
| 4 | feasibility summary (single threshold) | `stats/policy_zone_feasibility_summary.csv` | n_units_usable_both_periods and haul-count distribution per scheme/parameter combo, at MIN_HAULS_PER_CELL=5 only. | ≥5-haul filter is exploratory-only, superseded by reliability-derived thresholds — not the primary analysis filter |
| 4 | sample-size adequacy | `stats/policy_zone_sample_size_adequacy.csv` | Variance decomposition (between/within-cell) and reliability-derived required-n per scheme/parameter combo, from haul-level residual ICC. | ICC 0.0975-0.1774; required-n for 0.8 reliability = 19-38 hauls/cell depending on combo |
| 4 | feasibility by threshold | `stats/policy_zone_feasibility_by_threshold.csv` | Feasibility (n_units_usable_both_periods, pct) swept across a fixed threshold grid (1-50) plus each combo's own reliability-derived thresholds. | the primary sensitivity deliverable for the ≥5-haul-filter follow-up |
| 4 | run log | `run_logs/h3_policy_run_log.md` | Full run log: provisional constants, universe of rectangles, Scheme A/B construction method, feasibility summary, sample-size-adequacy method and results, figure index. | universe = 197 rectangles (Scheme A), 165 Couce-covered rectangles (Scheme B) |
