# Live pipeline — record for supervisor review

Scripts that produce **presented** results (methods model list).
Superseded work: [`exploratory/`](../exploratory/).

Updated: 2026-08-19 (H3 reported from CAR; H2/H3 share `wb_car_v2`).

---

## H1 — haul-level EEoS

| Script | Role |
|--------|------|
| `pipeline/build_datras_state_variables.R` | S, N, E_raw, m_min |
| `pipeline/build_eeos_predictions.R` | FishGlob join + EEoS B_pred + residuals |
| `pipeline/run_h1_harte_baseline.R` | Unfitted EEoS + productivity 1:1 + Harte Fig. 2 |
| `pipeline/run_h1_presentation_figures.R` | Presentation figures |
| `pipeline/fig_state_variable_distributions.R` | Methods 4-panel S/N/E/B_obs distributions |
| `pipeline/run_h1_lne_reference.R` | ln(E) OLS + model comparison |
| `pipeline/run_h1_null_model.R` | Uniform-draw null |
| `pipeline/explore_h1_catchability_scaling.R` | Quartile ratios / scalar reject |
| `pipeline/explore_h1_haul_dominance.R` | D + size_CV |
| `pipeline/explore_h1_dominance_partial_r2.R` | Partial R² beyond biomass |
| `pipeline/run_pipeline_diagnostics.R` | Audit |
| `pipeline/run_h1_dropout_diagnosis.R` | Dropout funnel |

---

## H2 inputs

| Script | Role |
|--------|------|
| `pipeline/import_couce_fishing_effort.R` | Couce hours |
| `pipeline/build_h2_rectangle_panel.R` | Rectangle panel |

---

## H2/H3 — primary (`phase_v2`) + methods robustness

| Script | Role |
|--------|------|
| `pipeline/run_h2h3_within_between.R` | Prerequisite: original-phase stack → `h2h3_wb_model_objects.rds` |
| `pipeline/run_h2h3_phase_v2_refit.R` | Companion RE `wb_primary_v2` (feeds reporting; not H3 source) |
| `pipeline/run_h2h3_phase_v2_reporting.R` | **Primary CAR** `wb_car_v2` (H2 + H3 slopes) + proportional effects + pooled contrast + figs |
| `pipeline/run_h2h3_car_reporting_diagnostics.R` | CAR supplement: leave-one-rectangle-out, Nakagawa R², signed vs unsigned audit (reads `fit_car`; does not change primary) |
| `pipeline/run_h2h3_h2_multiplicity_subsampling.R` | Bonferroni m=4 for H2 and H3 families |
| `pipeline/permutation_bootstrap_FP_between_CAR.R` | FP_between spatial permutation (CAR) |
| `pipeline/run_h2h3_rectangle_subsampling_refit.R` | Rectangle subsample RE (companion) |
| `pipeline/run_h2h3_rectangle_subsampling_car_refit.R` | Rectangle subsample CAR (H2 + H3 write-up) |
| `pipeline/build_knn_spatial_weights.R` | k-NN k=4 weights (Spec A) |
| `pipeline/run_h2h3_spec_a_car_identifiability.R` | Spec A on CAR |
| `pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R` | Spec A permutation |

Helpers: `pipeline/R/h2h3_within_between_helpers.R`, `h2h3_feasibility_helpers.R`,
`h2h3_results_helpers.R`, `h2h3_knn_spatial_helpers.R`, `h2h3_spatial_autocorr_helpers.R`.

### Headline outputs

- `outputs/primary_model_v2.rds` (`primary_model_v2`, `fit_wb_car_v2`)
- `outputs/phase_v2_fp_slopes_by_phase.csv`
- `outputs/phase_v2_proportional_effects_H{2,3}.csv`
- `outputs/phase_v2_pooled_between_coef.csv`
- `outputs/h2_car_reporting_diagnostics.md` (LOO + CAR R² + signed vs unsigned; after P2)
- `outputs/figures/phase_v2_presentation_H{2,3}_gap_change_by_phase.png`
- `outputs/h2_multiplicity_correction_*`, `outputs/h3_multiplicity_correction_*`
- `outputs/permutation_bootstrap_FP_between_CAR_*`
- `outputs/h2h3_glmm_subsampling_*`, `outputs/h2h3_car_subsampling_*`
- `outputs/primary_model_v2_spec_a_car.rds`, `outputs/spec_a_car_*`
- `outputs/permutation_bootstrap_FP_between_CAR_spec_a_*`

---

## Intermediate inputs kept at top-level `outputs/`

Live scripts still read:

- `outputs/h2h3_wb_model_objects.rds` (phase_v2 refit base)
- `outputs/h2h3_feasibility_round2_model_objects.rds` (CAR `adjMatrix`)
- `outputs/h2h3_results_model_objects.rds` / related historical inputs where referenced

---

## Not live (see `exploratory/`)

Archived 2026-08-13 (plus earlier moves):

- Early rectangle SEM (`run_h2_models.R`)
- Blended-term / zones / GAM / Bai–Perron structbreak / dose-response GAM
- RE Spec A + RE Spec A permutation; Spec B / A+B
- RE FP_between permutation
- OLS residual-proxy rectangle subsampling
- Bathymetry anisotropy / BYM smoothness / n_hauls metric / haul-count maps
- Original-phase standalone proportional / pooled / presentation / raw-correlation scripts
- Superseded results note `H2_H3_results_interpretation.md` (blended-term era)
