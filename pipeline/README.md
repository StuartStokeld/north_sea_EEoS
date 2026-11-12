# North Sea EEoS — analysis pipeline

PhD research testing whether the Ecological Equation of State (EEoS) predicts
community biomass at haul level in NS-IBTS Q1 (1985–2015), and whether fishing
pressure explains EEoS residuals spatially (H2) and temporally (H3).

**EEoS:** [micbru/equation_of_state](https://github.com/micbru/equation_of_state)
(`biomass.py`), called from R via `reticulate`.

**Live vs exploratory:** This folder is the **live** pipeline only.
Run log / model list: [`../outputs/live_pipeline_run_log.md`](../outputs/live_pipeline_run_log.md).
Archived scripts/outputs: [`../exploratory/`](../exploratory/).

---

## Prerequisites

- R (4.4+) with `renv` — open `north_sea_eeos.Rproj` and restore packages
- Python 3 with `numpy`, `scipy`, `pandas` — `.venv/` at repo root
- Run all scripts from the **repo root**

---

## Required data

| Path | Description |
|------|-------------|
| `FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData` | FishGlob B_obs |
| `outputs/datras_hl_raw.rds` | DATRAS HL (built from CSV if missing) |
| `outputs/datras_hh_raw.rds` | DATRAS HH haul records: HaulDur, DataType, HaulVal (`import_datras_hh.R`; from the ICES CSV, fallback `icesDatras`) |
| `outputs/fishbase_lw_lookup_v2.csv` | Length–weight parameters |
| `equation_of_state/` | Cloned EEoS repo |
| `data/external/couce_trawling_effort/...REVIEW_v2.csv` | Couce fishing hours |
| `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp` | Rectangle geometry |

---

## Live pipeline — run in order

### H1

| Step | Script | Main outputs |
|------|--------|--------------|
| −1 | `import_datras_hh.R` (only if HH RDS missing) | `datras_hh_raw.rds` |
| −1′ | `validate_haul_effort_alignment.R` (optional check) | `datras_effort_validation/` |
| 0 | `build_datras_state_variables.R` | `datras_haul_state.rds`, `datras_haul_E.rds`, `datras_effort_alignment_summary.csv`, `datras_effort_raw_count_check.csv` |
| 1 | `build_eeos_predictions.R` | `haul_state_variables.rds`, `haul_eeos_predictions.rds` |
| 2a | `run_h1_harte_baseline.R` | `h1_harte_baseline_metrics.csv`, `h1_fig2_leverage_hauls.csv` |
| 2a′ | `run_h1_presentation_figures.R` | `figures/h1_presentation_*.png` — **stale** (not re-run after DataType fix) |
| 2b | `run_h1_lne_reference.R` | `h1_model_comparison.csv`, `h1_lne_coefficients.csv` |
| 3 | `run_h1_null_model.R` | `null_summary.rds`, `null_distribution_summary.csv` |
| — | `explore_h1_catchability_scaling.R` | `h1_catchability_*.csv` — **stale** |
| — | `explore_h1_haul_dominance.R` / `explore_h1_dominance_partial_r2.R` | `h1_dominance_*.csv` — **stale** (`explore_h1_haul_dominance.R` still builds counts from unaligned HL) |
| — | `run_h1_dropout_diagnosis.R` / `run_pipeline_diagnostics.R` | dropout / audit tables |

```bash
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
```

### H2/H3 inputs → primary → robustness

| Step | Script | Main outputs |
|------|--------|--------------|
| H2a | `import_couce_fishing_effort.R` | `h2_couce_rectangle_effort.rds` |
| H2b | `build_h2_rectangle_panel.R` | `h2_rectangle_panel.rds` |
| WB0 | `run_h2h3_within_between.R` | `h2h3_wb_model_objects.rds` (prerequisite) |
| P1 | `run_h2h3_phase_v2_refit.R` | `primary_model_v2.rds` (companion RE) |
| P2 | `run_h2h3_phase_v2_reporting.R` | CAR / H2+H3 slopes, proportional effects, figs |
| P2′ | `run_h2h3_car_reporting_diagnostics.R` | CAR LOO influence, R², signed vs unsigned audit (does not refit primary) — **stale** |
| R1 | `run_h2h3_h2_multiplicity_subsampling.R` | H2 + H3 Bonferroni |
| R2 | `permutation_bootstrap_FP_between_CAR.R` | CAR spatial permutation |
| R3 | `run_h2h3_rectangle_subsampling_refit.R` | RE subsample (companion) |
| R3′ | `run_h2h3_rectangle_subsampling_car_refit.R` | CAR subsample (H2 + H3) |
| S1 | `build_knn_spatial_weights.R` | `knn_listw_k4.rds` |
| S2 | `run_h2h3_spec_a_car_identifiability.R` | `primary_model_v2_spec_a_car.rds` |
| S3 | `permutation_bootstrap_FP_between_CAR_spec_a.R` | Spec A permutation |

```bash
Rscript pipeline/import_couce_fishing_effort.R
Rscript pipeline/build_h2_rectangle_panel.R
Rscript --vanilla pipeline/run_h2h3_within_between.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R
Rscript --vanilla pipeline/run_h2h3_car_reporting_diagnostics.R
Rscript --vanilla pipeline/run_h2h3_h2_multiplicity_subsampling.R
Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR.R
Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_refit.R
Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_car_refit.R
Rscript --vanilla pipeline/build_knn_spatial_weights.R
Rscript --vanilla pipeline/run_h2h3_spec_a_car_identifiability.R
Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R
```

Use `Rscript --vanilla` for glmmTMB / spaMM scripts.

**Stale after the DataType fix:** scripts marked **stale** above (and `fig_state_variable_distributions.R`)
were not re-run on branch `fix/haul-effort-datatype`; their outputs still describe the pre-fix haul set.
See `../outputs/live_pipeline_run_log.md`.

**Primary formula (`phase_v2`):**
`residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)`
(CAR; H2 and H3). Companion RE uses `(1 | stat_rec)`. Phases: 1985–1991 / 1992–2001 /
2002–2007 / 2008–2015.

---

## Shared R helpers (live)

| File | Contents |
|------|----------|
| `R/h1_*.R` | H1 metrics, join, null, dominance, Harte, ln(E) |
| `R/datras_*.R` | DATRAS import / state |
| `R/h2_common.R`, `h2_couce_helpers.R`, `h2_panel_helpers.R` | H2 inputs |
| `R/h2_spatial_helpers.R`, `h2h3_spatial_autocorr_helpers.R` | Weights / Moran helpers (Spec A) |
| `R/h2h3_within_between_helpers.R` | FP_between / FP_within |
| `R/h2h3_feasibility_helpers.R`, `h2h3_results_helpers.R` | Phase factor, CAR slopes |
| `R/h2h3_knn_spatial_helpers.R` | Spec A k-NN lag |

---

## Key conventions

| Variable | Notes |
|----------|--------|
| Haul selection | HL aligned to HH: `HaulVal == "V"`, DataType C/R/P/S, `HaulDur == 30` (primary). Sensitivity 26–34 min via `EEOS_HAUL_DUR_MIN` / `EEOS_HAUL_DUR_MAX` |
| Raw counts | `HLNoAtLngt × SubFactor × f`; `f = HaulDur/60` for DataType C (per-hour counts), 1 for R/P/S. No duration rescaling; `B_obs` unchanged |
| Hybrid / mean-length build | Disabled (cannot be effort-aligned); the state build stops if HL is incomplete |
| Primary residual | `log(B_obs) − log(B_pred)` |
| Primary R² (H1) | `log_r2()` = 1 − SS_res/SS_tot (can be negative) |
| Productivity baseline | `E × m_min` (unfitted) |
| H2 term | `FP_between` — rectangle long-run mean log(hours+1) |
| H3 term | `FP_within` — year deviation from that mean |
| Phases | Policy-anchored `phase_v2` (1992 / 2002 / 2008) |

---

## Citation

Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). An equation of state
unifies biodiversity, ecosystem functioning, and biomass in ecosystems.
*Communications Biology*, 5, 957.

DataType / haul-duration handling: ICES (2020) NS-IBTS manual, SISP 10-IBTS 10 Rev. 11
(https://doi.org/10.17895/ices.pub.7562; §2.5.4, Annex 2); ICES DATRAS FAQ For Data Users v3, 2014
(https://www.ices.dk/data/Documents/DATRAS/DATRAS_FAQs.pdf; C = adjusted to one hour, SubFactor = 1);
ICES DATRAS IBTS index calculation procedure
(https://www.ices.dk/data/Documents/DATRAS/Indices_Calculation_Steps_IBTS.pdf; non-C × 60/HaulDuration);
FishGlob `cleaning_codes/get_datras.R` (Maureaud et al. 2024,
*Scientific Data* 11, 24). Details: `../CLAUDE.md` (Key scientific conventions).

## Tests

Run the whole suite from the repo root (needs `testthat`; about 2.5 minutes):

```bash
Rscript pipeline/tests/run_all.R                     # fast tests
EEOS_SLOW_TESTS=1 Rscript pipeline/tests/run_all.R   # + slow tests (raw ICES CSV checks, solver check on all hauls)
```

- **Layout.** Tests live in `pipeline/tests/testthat/`: input contracts (`test-contracts-*.R`), hand-computed
  units (`test-units-hand.R`), the EEoS solver incl. the 42-community golden test (`test-eeos-solver.R`),
  pipeline integrity and provenance (`test-integrity.R`) and metamorphic checks (`test-metamorphic.R`).
- **Inputs.** The suite needs the committed `outputs/`, plus the local FishGlob RData and `equation_of_state/`.
  Tests skip cleanly when an input is missing.
- **Expected failures.** Tests tagged `[F-..]` encode findings in `audit/findings.md` and fail until those
  findings are fixed. `run_all.R` lists any untagged failure separately as a regression. It exits with status 1
  if anything fails.
- **Mutation testing:** `Rscript audit/scripts/mutation_test.R` (`audit/mutation_report.md`).
- **Dependency graph:** `Rscript --vanilla audit/scripts/build_dependency_graph.R` →
  `audit/dependency_graph.md` / `.csv`.

