# H2/H3 primary model — policy-anchored phase_v2 refit — run log

Refits the primary within-between model with policy-anchored phase breakpoints (1992 / 2002 / 2008) instead of the data-driven structural breaks (1989 / 2001 / 2008). Original `phase` column retained for comparison. No bootstrap or bathymetry logic.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_refit_sessionInfo.txt
glmmTMB 1.1.14

## Inputs
Loaded original primary model: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_objects.rds (formula: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec)).
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015.
Original phase levels: 1985-1988 | 1989-2000 | 2001-2007 | 2008-2015.

## phase_v2 definition (policy-anchored)
breaks = c(-Inf, 1991, 2001, 2007, Inf); right = TRUE → 1985–1991 | 1992–2001 | 2002–2007 | 2008–2015.
Anchors: pre-reform baseline; 1992 CFP reform; 2002 CFP reform (in force Jan 2003); 2008 LTMP / MSFD.

### Boundary-year check: table(year, phase_v2)
  year 1991 → phase_v2 1985-1991 (expected 1985-1991) OK
  year 1992 → phase_v2 1992-2001 (expected 1992-2001) OK
  year 2001 → phase_v2 1992-2001 (expected 1992-2001) OK
  year 2002 → phase_v2 2002-2007 (expected 2002-2007) OK
  year 2007 → phase_v2 2002-2007 (expected 2002-2007) OK
  year 2008 → phase_v2 2008-2015 (expected 2008-2015) OK

### Phase sizes (expect roughly 7 / 10 / 6 / 8 years)
  1985-1991: n_years=7, n_obs=2643 (years 1985–1991)
  1992-2001: n_years=10, n_obs=3258 (years 1992–2001)
  2002-2007: n_years=6, n_obs=2062 (years 2002–2007)
  2008-2015: n_years=8, n_obs=2501 (years 2008–2015)

## Step 1 — Refit primary model with phase_v2
Target formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec) [REML]
Note: stats::update(. ~ . - phase + phase_v2) leaves interaction terms on the old `phase` factor; refitting explicitly with the full formula.
Refit elapsed: 0.35 sec.
Fitted formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)

## Step 2 — Compare against original model (values from saved fit)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_vs_original_comparison.csv

| Quantity | Original (phase) | Revised (phase_v2) |
|---|---|---|
| FP_between_coef | -0.00595837 | -0.024729 |
| FP_between_SE | 0.0179877 | 0.0155049 |
| FP_within_coef | 0.0317424 | 0.0268781 |
| FP_within_SE | 0.0181741 | 0.0148443 |
| AIC | 13252.8 | 13256.7 |
| BIC | 13354.4 | 13358.3 |
| rectangle_RE_variance | 0.0287319 | 0.0286321 |
Delta AIC (v2 − original) = 3.938; Delta BIC = 3.938.

## Step 3 — BLUP Moran / Geary diagnostic
No persisted listw RDS exists in outputs/. Reusing the established 158-rectangle queen-contiguity construction from run_h2h3_primary_spatial_autocorr_check.R (ICES DBF SOUTH/WEST grid; style W; zero.policy=TRUE) — identical to the weights that produced the archived original BLUP Moran I ≈ 0.552 / Geary C ≈ 0.441. Not a new weights scheme.
listw: N=158; n_isolated=0; mean degree=6.911.
Archived original BLUP: Moran I=0.551900, Geary C=0.440912 (source: spatial_diagnostics_comparison.csv).
Recomputed original BLUP on same listw: Moran I=0.551900, Geary C=0.440912 (|ΔI|=0, |ΔC|=1.11e-16).
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_blup_diagnostic.csv

| Model | Moran's I | Geary's C | p (Moran) | p (Geary) | N |
|---|---|---|---|---|---|
| original_phase | 0.5519 | 0.4409 | <0.001 | <0.001 | 158 |
| phase_v2 | 0.5503 | 0.4422 | <0.001 | <0.001 | 158 |

## Outputs
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_vs_original_comparison.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_blup_diagnostic.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_refit_sessionInfo.txt
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_refit_run_log.md (this file)

Note: original `phase` and `h2h3_wb_model_objects.rds` are unchanged. `primary_model_v2.rds` is the policy-anchored primary refit.
