# Live pipeline — record for supervisor review

Scripts that produce **presented** results (methods model list).
Superseded work: [`exploratory/`](../exploratory/).

Updated: 2026-11-10 (DATRAS DataType / haul-duration fix; primary = hauls of exactly 30 min).
Previous: 2026-08-19 (H3 reported from CAR; H2/H3 share `wb_car_v2`).

---

## DATRAS DataType / haul-duration fix (branch `fix/haul-effort-datatype`, 2026-11-10)

DATRAS HL counts for DataType C are numbers per hour, but the pipeline used them as raw counts
(~2× N and E for 30-min tows) against FishGlob raw haul biomass. HL is now aligned to the HH haul
records before S, N, E are built: `HaulVal == "V"`, DataType C/R/P/S, **`HaulDur == 30` exactly
(primary; supervisor decision 10 Nov 2026)**; C counts × HaulDur/60, R/P/S × SubFactor; no duration
rescaling; `B_obs` unchanged. Sensitivity: 26–34 min window (`EEOS_HAUL_DUR_MIN=26
EEOS_HAUL_DUR_MAX=34`), outputs in `outputs/sensitivity_window_26_34/` (gitignored). Sources: ICES
(2020) NS-IBTS manual SISP 10 Rev. 11 §2.5.4 and Annex 2; DATRAS FAQ; DATRAS NS-IBTS index
procedure; FishGlob `get_datras.R` (see `CLAUDE.md`).

- Samples (exactly 30 min): H1 10,128 hauls; H2/H3 8,670 hauls, 158 rectangles (pre-fix 12,069 / 10,464).
- All live scripts in the CLAUDE.md order were re-run at exactly 30 min (19 scripts, exit 0).
- Comparison old vs 30 min vs 26–34: [`effort_fix_comparison.md`](effort_fix_comparison.md)
  (tables in `effort_fix_comparison_tables/`). Pre-fix outputs archived locally in
  `outputs/archive_pre_effort_fix/` (gitignored).
- **Results text in `display_discussion/` and CLAUDE.md is superseded pending interpretation after the
  DataType fix.**

### Stale outputs (scripts not re-run after the fix)

These scripts were not re-run; their outputs still describe the pre-fix haul set (12,069 H1 hauls,
C counts per hour) and must not be read as current:

| Script | In git? | Tracked outputs now stale |
|--------|---------|---------------------------|
| `pipeline/run_h1_presentation_figures.R` | yes | none tracked (`figures/h1_presentation_*.png` local only) |
| `pipeline/fig_state_variable_distributions.R` | no (local) | none tracked |
| `pipeline/explore_h1_catchability_scaling.R` | yes | `h1_catchability_by_quartile.csv`, `h1_catchability_scaling_summary.csv`, `figures/h1_catchability_ratio_vs_bobs.png` |
| `pipeline/explore_h1_dominance_partial_r2.R` | yes | `h1_dominance_partial_r2_{anova,coefficients,summary}.csv` |
| `pipeline/run_h2h3_car_reporting_diagnostics.R` | no (local) | none tracked (`h2_car_*` local only) |
| `pipeline/explore_h1_haul_dominance.R` | yes | other `h1_dominance_*.csv`, `figures/h1_dominance_ratio_vs_{D,sizeCV}.png`; also still builds counts from unaligned HL |

---

## H1 — haul-level EEoS

| Script | Role |
|--------|------|
| `pipeline/import_datras_hh.R` | HH haul records → `datras_hh_raw.rds` (only if missing) |
| `pipeline/validate_haul_effort_alignment.R` | Effort-alignment checks → `datras_effort_validation/` (optional) |
| `pipeline/build_datras_state_variables.R` | HL→HH effort alignment; S, N, E_raw, m_min |
| `pipeline/build_eeos_predictions.R` | FishGlob join + EEoS B_pred + residuals |
| `pipeline/run_h1_harte_baseline.R` | Unfitted EEoS + productivity 1:1 + Harte Fig. 2 |
| `pipeline/run_h1_presentation_figures.R` | Presentation figures — **stale** |
| `pipeline/fig_state_variable_distributions.R` | Methods 4-panel S/N/E/B_obs distributions — **stale** |
| `pipeline/run_h1_lne_reference.R` | ln(E) OLS + model comparison |
| `pipeline/run_h1_null_model.R` | Uniform-draw null |
| `pipeline/explore_h1_catchability_scaling.R` | Quartile ratios / scalar reject — **stale** |
| `pipeline/explore_h1_haul_dominance.R` | D + size_CV — **stale** |
| `pipeline/explore_h1_dominance_partial_r2.R` | Partial R² beyond biomass — **stale** |
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
| `pipeline/run_h2h3_car_reporting_diagnostics.R` | CAR supplement: leave-one-rectangle-out, Nakagawa R², signed vs unsigned audit (reads `fit_car`; does not change primary) — **stale** |
| `pipeline/run_h2h3_h2_multiplicity_subsampling.R` | Bonferroni m=4 for H2 and H3 families |
| `pipeline/permutation_bootstrap_FP_between_CAR.R` | FP_between spatial permutation (CAR) |
| `pipeline/run_h2h3_rectangle_subsampling_refit.R` | Rectangle subsample RE (companion) |
| `pipeline/run_h2h3_rectangle_subsampling_car_refit.R` | Rectangle subsample CAR (H2 + H3 write-up) |
| `pipeline/build_knn_spatial_weights.R` | k-NN k=4 weights (Spec A) |
| `pipeline/run_h2h3_spec_a_car_identifiability.R` | Spec A on CAR |
| `pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R` | Spec A permutation |

Helpers: `pipeline/R/h2h3_within_between_helpers.R`, `h2h3_feasibility_helpers.R`,
`h2h3_results_helpers.R`, `h2h3_knn_spatial_helpers.R`, `h2h3_spatial_autocorr_helpers.R`.

**H2/H3 data source (changed 2026-11-10, branch `fix/haul-effort-datatype`).**
`run_h2h3_within_between.R` previously read its haul-level analysis data (`residual`, Couce
`log_hours_total`, `phase`) from `outputs/h2h3_results_model_objects.rds$data`. That file is
written only by an exploratory script (`exploratory/pipeline/run_h2h3_shared_model_results.R`),
so it was **stale**: re-running the live H1 pipeline did not reach H2/H3. The script now rebuilds
the data on every run with `build_feasibility_data()` (`pipeline/R/h2h3_feasibility_helpers.R`)
from the current `haul_eeos_predictions.rds`, `h2_rectangle_panel.rds` and
`h2_couce_year_effort.rds`. Rebuilding from the pre-fix inputs reproduces the stored data exactly
(10,464 hauls; residual / FP differences 0; primary RE fixed effects identical to 7e-15). The CAR
adjacency matrix is still read from `outputs/h2h3_feasibility_round2_model_objects.rds`; the
script stops if the analysis rectangles differ from its rows. Formula, phases, priors and
fishing-pressure construction are unchanged.

The same script also did not add the project `.R_libs` to its library path (unlike the other H2/H3
scripts), so spaMM was not found and its original-phase CAR sensitivity
(`outputs/h2h3_wb_car_fixed_effects.csv`) was silently skipped. The path is now added; the CAR
sensitivity runs again (re-run 2026-11-10 at exactly 30 min; RE outputs unchanged).

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
