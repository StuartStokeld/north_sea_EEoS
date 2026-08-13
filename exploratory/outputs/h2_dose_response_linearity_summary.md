# H2 dose-response linearity diagnostic — summary

Generated: 2026-08-05 15:42 BST

## Decision flag

**`NONLINEARITY_DETECTED`**

## One-paragraph result

H2 dose-response linearity diagnostic (mgcv, same haul data and (1|stat_rec) structure as primary_model_v2; H3 left linear): ML AIC linear = 12942.9, smooth = 12919.0 (Δ = -23.9); nested Chi-sq test p = 7.609e-06. 1985-1991: linear assumption held (edf = 1.00 ≈ 1); 1992-2001: linear assumption held (edf = 1.00 ≈ 1); 2002-2007: linearity NOT held (edf = 3.45 > 1.5; shape ≈ non-monotonic (mid-range dip / uneven rise)); 2008-2015: linearity NOT held (edf = 3.07 > 1.5; shape ≈ non-monotonic (mid-range dip / uneven rise)). Overall verdict: nonlinearity is detected in at least one phase or by the nested fit comparison — production model unchanged pending design review.

## Results table (per phase)

See `exploratory/outputs/h2_dose_response_linearity_by_phase.csv`. Gate: edf ≤ 1.5 treated as linear-adequate; nested test and ΔAIC are model-level.

| phase_v2 | edf | Ref.df | F | p (vs zero) | phase_linearity_held |
|----------|-----|--------|---|-------------|----------------------|
| 1985-1991 | 1.001 | 1.001 | 2.360 | 0.1244 | TRUE |
| 1992-2001 | 1.000 | 1.001 | 1.237 | 0.2661 | TRUE |
| 2002-2007 | 3.454 | 3.823 | 22.842 | 0 | FALSE |
| 2008-2015 | 3.069 | 3.543 | 6.331 | 0.0001242 | FALSE |

## Model comparison (ML)

- Linear GAM AIC = 12942.88
- Smooth GAM AIC = 12918.97
- ΔAIC (smooth − linear) = -23.91
- Nested anova Chi-sq: Df = 8.411, statistic = 7.689, p = 7.609e-06

## Figures

- Panel: `exploratory/outputs/figures/h2_dose_response_linearity_by_phase.png`
- Per phase: `exploratory/outputs/figures/h2_dose_response_linearity_1985_1991.png` (and siblings)

## What happens next (not implemented)

STOP — do not respecify the production model from this script. Review the phase-wise curves (threshold vs saturation vs non-monotonic) before any redesign of the fixed-effect form, phase interaction, or re-running the Moran / permutation / KNN residual pipeline.

## Explicitly unchanged

- Production model `outputs/primary_model_v2.rds`
- H3 (`FP_within × phase_v2`)
- Random-intercept structure, phase cut points
- Moran / permutation / KNN residual pipeline

