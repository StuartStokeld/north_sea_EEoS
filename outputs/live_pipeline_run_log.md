# Live pipeline — record for supervisor review

This file states which scripts currently produce the **presented** results.
Superseded / exploratory work lives under [`exploratory/`](../exploratory/) and is
not part of this list.

Generated: 2026-07-27 (repo reorganization; no new analysis).

---

## Live pipeline

### H1 — haul-level EEoS

| Script | Role |
|--------|------|
| `pipeline/build_datras_state_variables.R` | State variables (S, N, E, m_min) |
| `pipeline/build_eeos_predictions.R` | FishGlob join + EEoS B_pred + residuals |
| `pipeline/run_h1_harte_baseline.R` | Harte unfitted baseline (Tests 1–2) |
| `pipeline/run_h1_lne_reference.R` | ln(E) OLS + model comparison |
| `pipeline/run_h1_null_model.R` | B-only null permutations |
| `pipeline/run_pipeline_diagnostics.R` | Pipeline audit |
| `pipeline/run_h1_dropout_diagnosis.R` | Dropout funnel / spike years |

Supporting H1 diagnostics (still under `pipeline/`, not superseded):
`explore_h1_catchability_scaling.R`, `explore_h1_haul_dominance.R`,
`explore_h1_dominance_partial_r2.R`.

### H2 inputs (rectangle panel + fishing effort)

| Script | Role |
|--------|------|
| `pipeline/import_couce_fishing_effort.R` | Couce fishing hours |
| `pipeline/build_h2_rectangle_panel.R` | Rectangle-level residual panel |

### Structural-break check (phase boundaries for H2/H3)

| Script | Role |
|--------|------|
| `pipeline/run_h2h3_structbreak_check.R` | BIC structural breaks → phases 1989 / 2001 / 2008 |
| `pipeline/R/h2h3_structbreak_helpers.R` | Helpers |

Outputs: `outputs/h2h3_designA4_structbreak_*.csv`,
`outputs/figures/h2h3_designA4_structbreak_series.png`,
`outputs/h2h3_designA4_structbreak_run_log.md`.

### H2/H3 — final within-between decomposed model (biomass-free)

Primary inference is the **within-between** decomposition (not the earlier blended
`log_hours_total` term). Biomass covariate excluded.

| Script | Role |
|--------|------|
| `pipeline/run_h2h3_within_between.R` | Original primary: `FP_between * phase + FP_within * phase + (1\|stat_rec)` (data-driven phases 1989/2001/2008); CAR sensitivity |
| `pipeline/run_h2h3_phase_v2_refit.R` | **Current primary:** same within-between formula with policy-anchored `phase_v2` (1992/2002/2008); comparison + BLUP Moran/Geary |
| `pipeline/run_h2h3_phase_v2_reporting.R` | phase_v2 slopes, CAR + pooled contrast, proportional effects, presentation figures |
| `pipeline/R/h2h3_within_between_helpers.R` | Decomposition helpers |
| `pipeline/run_h2h3_wb_proportional_effects.R` | Proportional / gap-change effect sizes (H2 & H3) |
| `pipeline/run_h2h3_wb_pooled_between_contrast.R` | Optional pooled-between CAR contrast (discussion) |
| `pipeline/run_h2h3_fp_between_confounding_bootstrap.R` | Supporting: rectangle-level `FP_between` spatial permutation null for H2 |
| `pipeline/build_knn_spatial_weights.R` | Spec A: k-NN (k=4) weights for lagged `FP_between` |
| `pipeline/run_h2h3_spec_a_lag_refit.R` | Spec A: primary + `FP_between_lag * phase_v2`; VIF; BLUP Moran |
| `pipeline/run_h2h3_fp_between_confounding_bootstrap_spec_a.R` | Spec A: confounding bootstrap with lag recomputed per shuffle |
| `pipeline/run_h2h3_spec_b_lag_refit.R` | Spec B / A+B: lagged neighbour `ln_B_obs`; VIF; BLUP+residual Moran; H2/H3 coef comparison |
| `pipeline/R/h2h3_knn_spatial_helpers.R` | Spec A/B helpers (k-NN, lag, VIF) |
| `pipeline/run_h2h3_presentation_figures.R` | Presentation figures from live wb outputs |
| `pipeline/run_h2h3_raw_correlation_stats.R` | Raw correlation companion stats |
| `pipeline/run_h2_n_hauls_metric_diagnostics.R` | Haul-count vs metric screen (wedge + mean-trend; sensitivity if flagged) |
| `pipeline/run_h2h3_haulcount_by_phase_map.R` | Haul-count choropleth by phase |

### H2/H3 — spatial residual diagnostics (supporting)

| Script | Role |
|--------|------|
| `pipeline/run_h2h3_primary_spatial_autocorr_check.R` | Moran/Geary on BLUPs + rectangle-mean residuals |
| `pipeline/run_h2h3_bathymetry_anisotropy_check.R` | Bathymetry anisotropy diagnostic (directional variogram + Jammalamadaka–Sarma alignment) |

Bathymetry anisotropy (2026-08-04): verdict **not_confirmed** (primary ρ_JS ≈ 0.05, p ≈ 0.55).
Outputs: `outputs/bathymetry_anisotropy_verdict.md`,
`outputs/bathymetry_by_rectangle.csv`, `outputs/bearing_alignment_test.csv`,
`outputs/directional_variogram_*.png`, `outputs/h2h3_bathymetry_anisotropy_run_log.md`.
Design: `display_discussion/Design_bathymetry_spatial_anisotropy.md`.

Shared helpers still required by the live H2/H3 scripts (also used historically by
exploratory runs): `pipeline/R/h2h3_feasibility_helpers.R`,
`pipeline/R/h2h3_results_helpers.R`.

Headline outputs: `outputs/h2h3_wb_*.csv`, `outputs/h2h3_wb_run_log.md`,
`outputs/figures/h2h3_wb_*.png`, `outputs/figures/h2h3_presentation_*.png`.
Policy-anchored primary: `outputs/primary_model_v2.rds`,
`outputs/phase_v2_vs_original_comparison.csv`,
`outputs/phase_v2_blup_diagnostic.csv`, `outputs/phase_v2_refit_run_log.md`.
Reporting stack: `outputs/phase_v2_fp_slopes_by_phase.csv`,
`outputs/phase_v2_proportional_effects_H{2,3}.csv`,
`outputs/phase_v2_pooled_between_coef.csv`,
`outputs/phase_v2_reporting_model_objects.rds`,
`outputs/figures/phase_v2_presentation_H{2,3}_gap_change_by_phase.png`,
`outputs/phase_v2_reporting_run_log.md`.
Haul-count diagnostic: `outputs/h2_n_hauls_metric_*.csv`,
`outputs/figures/h2_n_hauls_*.png`, discussion note
`display_discussion/H2_n_hauls_metric_check.md`.
FP_between spatial confounding bootstrap (supporting): 
`outputs/fp_between_confounding_bootstrap_*.csv|md|rds`,
`outputs/figures/fp_between_confounding_bootstrap_null.png`.

---

## Intermediate inputs kept at top-level `outputs/`

These files were produced by scripts now under `exploratory/`, but the live
within-between / structbreak scripts still read them (paths unchanged):

- `outputs/h2h3_designA1_year_fishing_summary.csv` (structbreak input)
- `outputs/h2h3_results_model_objects.rds` (analysis data reuse)
- `outputs/h2h3_results_fp_slopes_by_phase.csv` (blended comparison table)
- `outputs/h2h3_feasibility_round2_model_objects.rds` (CAR adjacency matrix)

---

## Not live (see `exploratory/`)

- Zone-scheme scripts (Scheme A block-merge, Scheme B pressure-tier)
- Pre-H3 feasibility / visualisation script
- Biomass-included shared-model feasibility runs
- Blended-term (undecomposed) results model and its GAM / proportional-effects tables
- Temporal-robustness check on the blended term
- Design-support zone / ICC explorations beyond the structbreak A4 outputs
