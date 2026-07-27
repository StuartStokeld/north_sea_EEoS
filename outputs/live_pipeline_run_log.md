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
| `pipeline/run_h2h3_within_between.R` | Primary: `FP_between * phase + FP_within * phase + (1\|stat_rec)`; CAR sensitivity |
| `pipeline/R/h2h3_within_between_helpers.R` | Decomposition helpers |
| `pipeline/run_h2h3_wb_proportional_effects.R` | Proportional / gap-change effect sizes (H2 & H3) |
| `pipeline/run_h2h3_wb_pooled_between_contrast.R` | Optional pooled-between CAR contrast (discussion) |
| `pipeline/run_h2h3_presentation_figures.R` | Presentation figures from live wb outputs |
| `pipeline/run_h2h3_raw_correlation_stats.R` | Raw correlation companion stats |

Shared helpers still required by the live H2/H3 scripts (also used historically by
exploratory runs): `pipeline/R/h2h3_feasibility_helpers.R`,
`pipeline/R/h2h3_results_helpers.R`.

Headline outputs: `outputs/h2h3_wb_*.csv`, `outputs/h2h3_wb_run_log.md`,
`outputs/figures/h2h3_wb_*.png`, `outputs/figures/h2h3_presentation_*.png`.

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
