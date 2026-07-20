# North Sea EEoS — Hypothesis 1 (haul-level biomass)

PhD research project testing whether the Ecological Equation of State (EEoS) predicts community biomass at individual trawl haul level in the NS-IBTS Q1 survey (1985–2015).

**EEoS implementation:** [micbru/equation_of_state](https://github.com/micbru/equation_of_state) (`biomass.py`), called from R via `reticulate`.

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
| — | [`run_h1_dropout_diagnosis.R`](run_h1_dropout_diagnosis.R) | Year funnel + 1998/2013–14 spike diagnosis | `outputs/h1_dropout_funnel_by_year.csv`, `display_discussion/H1_dropout_diagnosis.md` |
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

```bash
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
Rscript pipeline/run_pipeline_diagnostics.R
Rscript pipeline/run_h1_dropout_diagnosis.R
Rscript pipeline/explore_h1_catchability_scaling.R
Rscript pipeline/explore_h1_haul_dominance.R
Rscript pipeline/explore_h1_dominance_partial_r2.R
```

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
