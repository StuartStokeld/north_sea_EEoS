# H2/H3 primary-model residual spatial autocorrelation check — run log

Tests whether the primary hierarchical model `(1 | stat_rec)` has absorbed the spatial clustering documented in archived H2c (`outputs/h2_spatial_diagnostics.csv`), or whether clustering remains in BLUPs and/or rectangle-collapsed model residuals.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_primary_spatial_autocorr_sessionInfo.txt
Runtime path: base R + foreign. BLUPs/residuals recovered from saved glmmTMB artifact without loading glmmTMB; queen weights from ICES DBF SOUTH/WEST (equivalent to poly2nb queen=TRUE on this regular lattice); Moran/Geary via Cliff–Ord randomization moments matching h2_global_spatial_tests() / spdep defaults (zero.policy, alternative greater).

## Step 1b — Archived residual definition (checked before producing numbers)
FINDING (from pipeline/run_h2_models.R + outputs/h2_spatial_diagnostics.csv): `ols_primary_abs_residuals` is residuals(lm(mean_abs_residual ~ mean_annual_hours_total)), i.e. signed OLS residuals from a regression whose dependent variable was already rectangle-mean absolute EEoS residual. It is NOT abs(signed residual). The 'abs' in the variable name refers to the DV (mean_abs_residual), not to an absolute-value transform of the OLS residuals themselves. Closest role-match among the new primary-model residual rows is therefore the signed rectangle-mean response residual (leftover after the fitted model), not the mean-absolute collapse — despite the archived variable name.

Primary-model residual extraction: conditional response residual = y − (Xβ + b_i), recovered from fit_wb$obj$env$last.par.best. This matches `residuals(glmmTMB_fit, type = "response")` (includes the random intercept). Pearson residuals are not used.

## Inputs
Loaded archived diagnostics verbatim from: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_spatial_diagnostics.csv
Archived rows: mean_abs_residual (N=158); log_mean_annual_hours_total (N=158); ols_primary_abs_residuals (N=158)
Archived H2 panel: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_rectangle_panel.rds (158 rectangles).
Loaded primary model artifact: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_objects.rds (fit_wb; formula: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)).
Model was NOT refit — using saved glmmTMB object and its associated analysis data.
Analysis data: 10464 hauls, 158 unique rectangles, years 1985–2015.

## Spatial weights (matched to H2c)
No persisted nb/listw RDS was found in outputs/ from run_h2_models.R. Rebuilding queen contiguity from the same ICES shapefile DBF (gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.dbf) via SOUTH/WEST grid indices (0.5° lat × 1° lon), row-standardised style W, zero.policy=TRUE — the regular-lattice equivalent of build_h2_spatial_weights()/poly2nb(queen=TRUE).
Queen contiguity weights: N=158 rectangles; n_isolated=0; mean degree=6.911.
Archived diagnostics n_isolated=0 (must match rebuild: 0).

### Weights validation against archived statistics
  mean_abs_residual: recomputed Moran I=0.59100012 vs archived 0.59100012 (|diff|=2e-15) OK
  log_mean_annual_hours_total: recomputed Moran I=0.60282359 vs archived 0.60282359 (|diff|=2e-15) OK
  ols_primary_abs_residuals: recomputed Moran I=0.55822635 vs archived 0.55822635 (|diff|=1.55e-15) OK
Weights validated: recomputed Moran I matches archived CSV to machine precision.

## Step 1 — Extract BLUPs and rectangle-collapsed residuals
BLUPs: nrow=158; setequal(panel IDs)=TRUE; aligned to panel row order for listw.
Residuals: nrow(resid_by_rect)=158; columns resid (signed mean) and abs_resid (mean |residual|).
FLAG — archived-comparable residual row: **Primary model residuals — signed mean** (`rectangle-mean residual`). Basis: script/CSV check shows ols_primary_abs_residuals = signed OLS residuals of mean_abs_residual ~ hours, not an absolute residual; role-matched leftover after the primary model is the signed response residual. Mean-absolute residual is retained as secondary (magnitude clustering).

## Step 2 — Moran / Geary on primary-model objects
  (1|stat_rec) intercepts: Moran I=0.551900 (E=-0.006369, Var=0.00182724, z=13.060, p=2.78227e-39); Geary C=0.440912 (E=1.000000, p=1.14963e-37); N=158
  rectangle-mean residual: Moran I=0.473853 (E=-0.006369, Var=0.00181052, z=11.286, p=7.69187e-30); Geary C=0.502619 (E=1.000000, p=3.16455e-30); N=158
  rectangle-mean |residual|: Moran I=0.481456 (E=-0.006369, Var=0.00183984, z=11.373, p=2.85002e-30); Geary C=0.505883 (E=1.000000, p=7.41212e-30); N=158

## Step 3 — Comparison table
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spatial_diagnostics_comparison.csv

| Stage | Variable | Moran's I | Geary's C | p (Moran) | p (Geary) | N |
|---|---|---|---|---|---|---|
| Raw DV | `mean_abs_residual` | 0.591 | 0.402 | <0.001 | <0.001 | 158 |
| Raw IV | `log_mean_annual_hours_total` | 0.603 | 0.392 | <0.001 | <0.001 | 158 |
| OLS residuals | `ols_primary_abs_residuals` | 0.558 | 0.433 | <0.001 | <0.001 | 158 |
| Primary model BLUPs | `(1|stat_rec) intercepts` | 0.552 | 0.441 | <0.001 | <0.001 | 158 |
| Primary model residuals — signed mean | `rectangle-mean residual` | 0.474 | 0.503 | <0.001 | <0.001 | 158 |
| Primary model residuals — mean absolute | `rectangle-mean |residual|` | 0.481 | 0.506 | <0.001 | <0.001 | 158 |

Footnote: archived-comparable new residual row = `Primary model residuals — signed mean` / `rectangle-mean residual`. FINDING (from pipeline/run_h2_models.R + outputs/h2_spatial_diagnostics.csv): `ols_primary_abs_residuals` is residuals(lm(mean_abs_residual ~ mean_annual_hours_total)), i.e. signed OLS residuals from a regression whose dependent variable was already rectangle-mean absolute EEoS residual. It is NOT abs(signed residual). The 'abs' in the variable name refers to the DV (mean_abs_residual), not to an absolute-value transform of the OLS residuals themselves. Closest role-match among the new primary-model residual rows is therefore the signed rectangle-mean response residual (leftover after the fitted model), not the mean-absolute collapse — despite the archived variable name.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/blups_by_rectangle.csv
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/residuals_by_rectangle.csv
Caveat: BLUPs are shrinkage estimates (attenuated vs true rectangle effects), especially for sparse rectangles — see n_years / n_hauls in blups_by_rectangle.csv. Non-significant Moran on BLUPs = no strong evidence of remaining clustering, not proof that clustering is absent.

## Step 4 — Decision-rule classification
Alpha for significance: 0.05 (Moran's I one-sided greater, matching archived).
BLUPs: Moran I=0.5519, p=2.78227e-39 → SIGNIFICANT
Primary residual input (signed mean): Moran I=0.4739, p=7.69187e-30 → SIGNIFICANT
Secondary residual (mean absolute): Moran I=0.4815, p=2.85002e-30 → SIGNIFICANT

**Decision rule 2:**
"Exchangeability directly falsified: rectangle intercepts are themselves spatially clustered. Supports prioritizing a spatial-lag covariate over further error-covariance modelling."

Signed and absolute residual Moran tests agree on significance at alpha=0.05 (both significant).
