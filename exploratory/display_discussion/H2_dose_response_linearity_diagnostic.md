# H2 dose-response linearity diagnostic — design-review note

**Status:** archived / exploratory (not part of primary or sensitivity reporting). Diagnostic gate complete (2026-08-05). Production model unchanged.

**Script:** `exploratory/pipeline/run_h2h3_h2_dose_response_linearity.R`  
**Summary artifact:** `exploratory/outputs/h2_dose_response_linearity_summary.md`

## Decision flag

`NONLINEARITY_DETECTED` — do **not** auto-respecify the primary within-between model.

## Result (one paragraph)

H2 dose-response linearity diagnostic (mgcv, same haul data and `(1|stat_rec)` structure as `primary_model_v2`; H3 left linear): ML AIC linear = 12942.9, smooth = 12919.0 (Δ = −23.9); nested Chi-sq test p = 7.6×10⁻⁶. **1985–1991:** linear assumption held (edf = 1.00 ≈ 1); **1992–2001:** linear assumption held (edf = 1.00 ≈ 1); **2002–2007:** linearity NOT held (edf = 3.45; smooth is shallow at low–mid FP then rises more steeply above ~9, departing from the straight IQR slope); **2008–2015:** linearity NOT held (edf = 3.07; smooth shows a mid-range dip then a steeper high-FP rise). Overall: nonlinearity is confined to the two post-2002 phases; production model left unchanged pending design review of fixed-effect form before any Moran / permutation / KNN residual re-run.

## What a redesign would need to decide (not implemented)

1. Keep linear H2 in early phases and only relax 2002–2007 / 2008–2015?
2. Global smooth / spline / threshold / piecewise form for `FP_between`?
3. How that interacts with `phase_v2` and with reported IQR / per-SD contrasts.
4. Re-running downstream residual spatial checks only after the new fixed-effect form is locked.

## Unchanged

- `outputs/primary_model_v2.rds`
- H3 (`FP_within × phase_v2`)
- Phase cut points and random-intercept structure
- Moran / permutation / KNN pipelines
