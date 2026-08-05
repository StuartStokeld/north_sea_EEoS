# North Sea EEoS — analysis pipeline

PhD research project testing whether the Ecological Equation of State (EEoS) predicts community biomass at individual trawl haul level in the NS-IBTS Q1 survey (1985–2015), and whether fishing pressure explains EEoS residuals spatially (H2) and temporally (H3).

**EEoS implementation:** [micbru/equation_of_state](https://github.com/micbru/equation_of_state) (`biomass.py`), called from R via `reticulate`.

**Live vs exploratory:** Scripts in this folder are the **live** pipeline (see
[`../outputs/live_pipeline_run_log.md`](../outputs/live_pipeline_run_log.md)).
Superseded zone-scheme / pre-H3 / biomass-included feasibility / blended-term work lives in
[`../exploratory/`](../exploratory/).

---

## Prerequisites

- R (4.4+) with `renv` — open `north_sea_eeos.Rproj` and allow renv to restore packages
- Python 3 with `numpy`, `scipy`, `pandas` — project uses `.venv/` at repo root  
  ```bash
  python3 -m venv .venv
  .venv/bin/pip install numpy scipy pandas
  Rscript -e 'install.packages("reticulate")'  # or renv::install("reticulate")
  ```
- Optional: `testthat` for unit tests — `install.packages("testthat")`
- Local data (not all tracked in git) — see [Required data](#required-data)

Run all scripts from the **workspace root** (`north_sea_eeos/`):

```bash
cd /path/to/north_sea_eeos
```

---

## Required data

| Path | Description |
|------|-------------|
| `FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData` | FishGlob NS-IBTS species × haul data (B_obs only) |
| `Unaggregated trawl and biological information_*/` | ICES HL CSV export (auto-imported to `datras_hl_raw.rds`) |
| `outputs/datras_hl_raw.rds` | DATRAS HL length bins (Q1 1985–2015); built from CSV if missing |
| `outputs/fishbase_lw_lookup_v2.csv` | Length–weight parameters |
| `equation_of_state/` | Cloned EEoS repo (contains `biomass.py`) |
| `data/external/couce_trawling_effort/NorthSea_trawling_effort_1985to2015_REVIEW_v2.csv` | Couce et al. (2020) fishing hours (download via Cefas Data Hub DOI 10.14466/CefasDataHub.61; see [`DATA_SOURCE.md`](../data/external/couce_trawling_effort/DATA_SOURCE.md)) |
| `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp` | ICES statistical rectangles for H2 spatial weights (tracked in repo) |

Legacy notebook: [`supplementary/fishbase_lookup.Rmd`](../supplementary/fishbase_lookup.Rmd).

---

## Pipeline — run in order

| Step | Script | Purpose | Main outputs |
|------|--------|---------|--------------|
| 0a | [`import_datras_hl_from_csv.R`](import_datras_hl_from_csv.R) | Optional: force CSV → RDS import | `outputs/datras_hl_raw.rds` |
| 0 | [`build_datras_state_variables.R`](build_datras_state_variables.R) | Imports CSV if needed; builds S, N, E_raw, m_min | `outputs/datras_haul_state.rds`, `outputs/datras_haul_E.rds` |
| 1 | [`build_eeos_predictions.R`](build_eeos_predictions.R) | Join FishGlob B_obs; call EEoS; compute residuals | `outputs/haul_state_variables.rds`, `outputs/haul_eeos_predictions.rds` |
| 2a | [`run_h1_harte_baseline.R`](run_h1_harte_baseline.R) | Harte unfitted baseline (Fig 1, Fig 2, productivity 1:1) | `outputs/h1_harte_baseline_metrics.csv`, `outputs/figures/harte_fig*.png` |
| 2b | [`run_h1_lne_reference.R`](run_h1_lne_reference.R) | Unified model comparison + ln(E) OLS (Tier 2) | `outputs/h1_model_comparison.csv`, `outputs/haul_h1_benchmarks.rds` |
| 3 | [`run_h1_null_model.R`](run_h1_null_model.R) | B-only null model (999 permutations) | `outputs/null_summary.rds`, `outputs/figures/null_r2_*.png` |
| — | [`run_pipeline_diagnostics.R`](run_pipeline_diagnostics.R) | Audit checks from `cursor_pipeline_audit.md` | `outputs/pipeline_audit_results.csv`, `outputs/h1_dropout_by_year.csv` |
| — | [`run_h1_dropout_diagnosis.R`](run_h1_dropout_diagnosis.R) | Year funnel + 1998/2013–14 spike diagnosis | `outputs/h1_dropout_funnel_by_year.csv` (local markdown report not tracked in repo) |
| — | [`explore_h1_catchability_scaling.R`](explore_h1_catchability_scaling.R) | Catchability offset exploration (**no correction applied**) | `display_discussion/H1_catchability_scaling_exploration.md` |
| — | [`explore_h1_haul_dominance.R`](explore_h1_haul_dominance.R) | Numerical dominance (D) + dominant-species size-homogeneity (size_CV) exploration (**parallel diagnostic track, no correction applied**) | `outputs/h1_dominance_*.csv`, `display_discussion/H1_dominance_size_homogeneity_exploration.md` |
| — | [`explore_h1_dominance_partial_r2.R`](explore_h1_dominance_partial_r2.R) | Partial R² of {D, size_CV} beyond log(B_obs); D×size_CV interaction test (follow-up to the dominance exploration above) | `outputs/h1_dominance_partial_r2_*.csv` |

### Hypothesis 2 (rectangle-level residuals vs fishing pressure)

Run after the H1 pipeline (`build_eeos_predictions.R` must exist).

| Step | Script | Purpose | Main outputs |
|------|--------|---------|--------------|
| H2a | [`import_couce_fishing_effort.R`](import_couce_fishing_effort.R) | Import Couce trawling hours | `outputs/h2_couce_rectangle_effort.rds`, `outputs/h2_couce_import_diagnostics.csv` |
| H2b | [`build_h2_rectangle_panel.R`](build_h2_rectangle_panel.R) | Aggregate residuals to rectangle panel | `outputs/h2_rectangle_panel.rds`, `outputs/h2_rectangle_panel.csv` |
| H2c | [`run_h2_models.R`](run_h2_models.R) | OLS, Moran's I, SEM, figures | `outputs/h2_ols_results.csv`, `outputs/h2_sem_results.csv`, `outputs/figures/h2_topline_result.png`, `outputs/figures/h2_*.png` |

```bash
Rscript pipeline/import_couce_fishing_effort.R
Rscript pipeline/build_h2_rectangle_panel.R
Rscript pipeline/run_h2_models.R
```

Requires **spdep**, **spatialreg**, and **sf** (plus **patchwork** for the topline figure). Install if not in renv:
`renv::install(c("sf", "spdep", "spatialreg", "patchwork"))`

Also requires ICES rectangle geometry: [`gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp`](../gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp).

### Exploratory / superseded H2–H3 work (moved)

Pre-H3 visualisation, Scheme A/B zone scripts, biomass-included feasibility,
blended-term results + GAM, and temporal-robustness on the blended term are under
[`../exploratory/`](../exploratory/) (see that folder’s README). They are not part of
the live pipeline.

### Structural-break check (live — defines H2/H3 phases)

| Script | Purpose | Main outputs |
|--------|---------|--------------|
| [`run_h2h3_structbreak_check.R`](run_h2h3_structbreak_check.R) | BIC structural breaks on annual mean fishing hours | `outputs/h2h3_designA4_structbreak_years.csv`, `outputs/figures/h2h3_designA4_structbreak_series.png`, `outputs/h2h3_designA4_structbreak_run_log.md` |

```bash
Rscript --vanilla pipeline/run_h2h3_structbreak_check.R
```

Requires **strucchange** (ad hoc / ambient library; not in `renv.lock`). Helpers:
`R/h2h3_structbreak_helpers.R`. Reads `outputs/h2h3_designA1_year_fishing_summary.csv`
(kept at top-level `outputs/` as a live input; producing design-support script is under
`exploratory/`).

### H2/H3 live model — within-between decomposition (biomass-free)

**This is the presented H2/H3 analysis.** Separates persistent between-rectangle fishing
pressure (**H2**) from within-rectangle year-to-year deviations (**H3**):
`residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)`.
Optional CAR sensitivity with the same decomposed terms. No biomass covariate.

| Script | Main outputs |
|--------|--------------|
| [`run_h2h3_within_between.R`](run_h2h3_within_between.R) | `outputs/h2h3_wb_primary_fixed_effects.csv`, `outputs/h2h3_wb_fp_slopes_by_phase.csv`, `outputs/h2h3_wb_wald_tests.csv`, `outputs/h2h3_wb_blended_comparison.csv`, `outputs/h2h3_wb_model_fit.csv`, `outputs/h2h3_wb_partial_pooling.csv`, `outputs/h2h3_wb_model_objects.rds`, `outputs/figures/h2h3_wb_fp_slopes_by_phase.png`, `outputs/figures/h2h3_wb_partial_pooling.png`, `outputs/h2h3_wb_run_log.md` |
| [`run_h2h3_phase_v2_refit.R`](run_h2h3_phase_v2_refit.R) | Policy-anchored primary: `outputs/primary_model_v2.rds`, `outputs/phase_v2_vs_original_comparison.csv`, `outputs/phase_v2_blup_diagnostic.csv`, `outputs/phase_v2_refit_run_log.md` |
| [`run_h2h3_phase_v2_reporting.R`](run_h2h3_phase_v2_reporting.R) | phase_v2 slopes + CAR + pooled contrast + proportional effects + presentation figures (`outputs/phase_v2_*`) |

```bash
Rscript --vanilla pipeline/run_h2h3_within_between.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R
```

Helpers: `R/h2h3_within_between_helpers.R` (plus `R/h2h3_feasibility_helpers.R`,
`R/h2h3_results_helpers.R`). Must use `--vanilla` (glmmTMB / spaMM ad hoc).
spaMM may live in project-local `.R_libs/` (gitignored); the reporting script
prepends that path when present.

`phase_v2` breakpoints (policy-anchored; original `phase` retained): 1985–1991 /
1992–2001 / 2002–2007 / 2008–2015. Formula:
`residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)`.
H2 proportional effects use CAR between slopes; H3 uses primary within slopes
(same convention as the original-phase reporting stack).

#### Pooled FP_between contrast (discussion only)

| Script | Main outputs |
|--------|--------------|
| [`run_h2h3_wb_pooled_between_contrast.R`](run_h2h3_wb_pooled_between_contrast.R) | `outputs/h2h3_wb_pooled_between_*.csv`, `outputs/h2h3_wb_pooled_between_run_log.md` |

```bash
Rscript --vanilla pipeline/run_h2h3_wb_pooled_between_contrast.R
```

#### FP_between spatial confounding bootstrap (supporting diagnostic)

Rectangle-level permutation of `FP_between` under the current primary within-between
fit: tests whether presented H2 slopes depend on the spatial arrangement of fishing
pressure relative to `(1 | stat_rec)`. Does not change the primary model.

Prefers `outputs/primary_model_v2.rds` (`phase_v2`: 1985–1991 / 1992–2001 /
2002–2007 / 2008–2015); falls back to legacy `h2h3_wb_model_objects.rds` if v2
is absent.

| Script | Main outputs |
|--------|--------------|
| [`run_h2h3_fp_between_confounding_bootstrap.R`](run_h2h3_fp_between_confounding_bootstrap.R) | `outputs/fp_between_confounding_bootstrap_results.csv`, `outputs/fp_between_confounding_bootstrap_summary.md`, `outputs/figures/fp_between_confounding_bootstrap_null.png`, `outputs/fp_between_confounding_bootstrap_run_log.md` |

```bash
Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R   # if primary_model_v2.rds missing
Rscript --vanilla pipeline/run_h2h3_fp_between_confounding_bootstrap.R
# optional: N_BOOT=1000 SEED=42 N_CORES=1
```

Requires **glmmTMB** (`--vanilla`; ad hoc / ambient library). Default `N_CORES=1`
(sequential; TMB is unsafe under forked `mclapply`). ~6 min for 1000 refits.

#### Spec A — k-NN lagged `FP_between` (mechanism probe)

Single-change refit of the primary: adds `FP_between_lag * phase_v2` using
row-standardised k-NN (k=4) weights. Does **not** replace the primary model.

```bash
Rscript --vanilla pipeline/build_knn_spatial_weights.R
Rscript --vanilla pipeline/run_h2h3_spec_a_lag_refit.R
Rscript --vanilla pipeline/run_h2h3_fp_between_confounding_bootstrap_spec_a.R
```

| Script | Main outputs |
|--------|--------------|
| [`build_knn_spatial_weights.R`](build_knn_spatial_weights.R) | `outputs/knn_listw_k4.rds`, `outputs/knn_listw_k4_audit.csv` |
| [`run_h2h3_spec_a_lag_refit.R`](run_h2h3_spec_a_lag_refit.R) | `outputs/primary_model_v2_spec_a.rds`, `outputs/fp_between_lag_rectangle.rds`, VIF + BLUP Moran |
| [`run_h2h3_fp_between_confounding_bootstrap_spec_a.R`](run_h2h3_fp_between_confounding_bootstrap_spec_a.R) | `outputs/fp_between_confounding_bootstrap_spec_a_*.{csv,md,rds}` |

#### Spec B — k-NN lagged outcome (`B_lag_neighbour`)

Contemporaneous k-NN mean of neighbours' rectangle-year mean `ln_B_obs`
(pooled main effect). Spec B alone and Spec A+B joint; Moran's I on BLUPs and
rectangle-mean residuals.

```bash
Rscript --vanilla pipeline/run_h2h3_spec_b_lag_refit.R
```

| Script | Main outputs |
|--------|--------------|
| [`run_h2h3_spec_b_lag_refit.R`](run_h2h3_spec_b_lag_refit.R) | `outputs/b_lag_neighbour_rectangle_year.rds`, `outputs/primary_model_v2_spec_b.rds`, `outputs/primary_model_v2_spec_ab.rds`, `outputs/spec_b_spatial_diagnostic.csv`, `outputs/spec_b_coefficient_comparison.csv` |

#### Proportional effect sizes (decomposed H2/H3)

| Script | Main outputs |
|--------|--------------|
| [`run_h2h3_wb_proportional_effects.R`](run_h2h3_wb_proportional_effects.R) | `outputs/h2h3_wb_proportional_effects_H2.csv`, `outputs/h2h3_wb_proportional_effects_H3.csv`, `outputs/figures/h2h3_wb_gap_change_by_phase.png`, `outputs/h2h3_wb_proportional_effects_run_log.md` |

```bash
Rscript --vanilla pipeline/run_h2h3_wb_proportional_effects.R
```

#### Presentation figures / raw correlations

```bash
Rscript --vanilla pipeline/run_h2h3_presentation_figures.R
Rscript --vanilla pipeline/run_h2h3_raw_correlation_stats.R
```

Canonical live-script list: [`../outputs/live_pipeline_run_log.md`](../outputs/live_pipeline_run_log.md).

**Unit tests:**

```bash
Rscript pipeline/tests/testthat.R
```

**Review reports:** see [`display_discussion/`](../display_discussion/).

---

## Shared R helpers

| File | Contents |
|------|----------|
| [`R/datras_constants.R`](R/datras_constants.R) | Non-fish AphiaIDs, analysis year/quarter |
| [`R/datras_hl_helpers.R`](R/datras_hl_helpers.R) | LngtCode → cm, HL cleaning, bin aggregation |
| [`R/datras_csv_import.R`](R/datras_csv_import.R) | ICES unaggregated CSV → HL RDS import |
| [`R/datras_state_helpers.R`](R/datras_state_helpers.R) | S, N, E_raw, m_min assembly and validation |
| [`R/h1_common.R`](R/h1_common.R) | Project root, log-scale R², RMSE, SS metrics |
| [`R/h1_join_helpers.R`](R/h1_join_helpers.R) | FishGlob–DATRAS join, null-E guard, filter reason codes |
| [`R/h1_harte_baseline_helpers.R`](R/h1_harte_baseline_helpers.R) | Harte unfitted baseline (Tier 1) |
| [`R/h1_lne_helpers.R`](R/h1_lne_helpers.R) | ln(E) OLS + unified model comparison |
| [`R/h1_null_helpers.R`](R/h1_null_helpers.R) | B-null distribution assessment and simulation |
| [`R/h1_dominance_helpers.R`](R/h1_dominance_helpers.R) | Haul numerical dominance (D) + dominant-species size-CV; confound checks; binned reporting |
| [`R/h2_common.R`](R/h2_common.R) | H2 constants, paths, stat_rec normalisation |
| [`R/h2_couce_helpers.R`](R/h2_couce_helpers.R) | Couce fishing effort import |
| [`R/h2_panel_helpers.R`](R/h2_panel_helpers.R) | Rectangle-level EEoS panel |
| [`R/h2_model_helpers.R`](R/h2_model_helpers.R) | OLS model fitting for H2 |
| [`R/h2_spatial_helpers.R`](R/h2_spatial_helpers.R) | Spatial weights, Moran's I, SEM, figures |
| [`R/h2h3_structbreak_helpers.R`](R/h2h3_structbreak_helpers.R) | Structural-break criterion tables / series figure helpers |
| [`R/h2h3_feasibility_helpers.R`](R/h2h3_feasibility_helpers.R) | Phase factor, haul-level analysis dataset, CAR adjacency (still required by live within-between) |
| [`R/h2h3_results_helpers.R`](R/h2h3_results_helpers.R) | H2/H3 term labelling, Nakagawa R², CAR slopes, Wald helpers (still required by live within-between) |
| [`R/h2h3_within_between_helpers.R`](R/h2h3_within_between_helpers.R) | Within-between FP decomposition helpers |

Exploratory-only helpers (`h3_pre_*`, `h3_policy_*`, `h2h3_design_support_*`,
`h2h3_temporal_robustness_*`) live under [`../exploratory/pipeline/R/`](../exploratory/pipeline/R/).

---

## Key conventions

| Variable | Units / notes |
|----------|----------------|
| `S` | `n_distinct(AphiaID)` per haul from DATRAS HL (LW-filtered, fish-only) |
| `N` | `sum(HLNoAtLngt × SubFactor)` per haul from same DATRAS rows as E |
| `B_obs` | grams (`sum(wgt) × 1000`; FishGlob `wgt` is kg) |
| `E` | normalised metabolic rate (`E_raw / m_min^0.75`; dimensionless) |
| `m_min` | grams; mass of smallest length bin in haul (fallback: min species mean mass if HL incomplete) |
| `B_pred` | grams (`B_pred_norm × m_min`) |
| `E_calibrated` | grams-equivalent productivity proxy (`E × m_min`); **headline** unfitted productivity 1:1 predictor (same E and m_min as EEoS). Replaces the former `E_raw` baseline. |
| `E_raw` | unnormalised metabolic sum (`∑ n × mass_g^0.75`); retained for diagnostics (`productivity_1to1_uncalibrated`) only — **not** the Test 1 headline baseline |
| **Primary R²** | `log_r2()` = 1 − SS_res/SS_tot (can be negative); **not** cor(log B)² |
| **Also reported** | `log_rmse()`, cor², Harte unfitted baselines (productivity 1:1 = E×m_min, E/B^(3/4) ratio) |
| **Tier 2 benchmark** | Fitted `lm(log B_obs ~ log E)` — alternative correlative model, not prior-method baseline |
| **Primary residual** | `residual = log(B_obs) − log(B_pred)` |
| `B_pred_norm`, `residual_norm` | Diagnostic only — not for primary inference |

Median `B_pred/B_obs` > 1 indicates a **systematic scale offset** (EEoS over-predicts at haul level), not a unit conversion error.

---

## Other project files

| File | Purpose |
|------|---------|
| `cursor_pipeline_audit.md` | Audit checklist for pipeline integrity |
| [`load_ns_ibts.R`](../supplementary/load_ns_ibts.R) | Quick loader for FishGlob NS-IBTS |
| [`supplementary/`](../supplementary/) | Exploratory scripts and legacy notebooks |

---

## Citation

Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). An equation of state unifies biodiversity, ecosystem functioning, and biomass in ecosystems. *Communications Biology*, 5, 957.
