# FP_between spatial confounding bootstrap — Spec A summary

Formula: `residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + (1 | stat_rec)` (REML)

Each replicate: shuffle rectangle-level `FP_between` once, **recompute 
`FP_between_lag` from the shuffled map** (fixed k-NN structure), refit Spec A.

- Seed: `42`; n_boot: 1000; n_failed: 0; runtime: 414.0 sec

## Null distribution vs observed

| Target | Observed | Null 2.5% | Null 97.5% | Empirical p | Inside null 95%? |
|--------|----------|-----------|------------|-------------|------------------|
| FP_between (reference-phase main effect, 1985-1991) | +0.009876 | -0.030784 | 0.031979 | 0.5584 | yes |
| Phase 1985-1991 FP_between slope | +0.009876 | -0.030784 | 0.031979 | 0.5584 | yes |
| Phase 1992-2001 FP_between slope | -0.001341 | -0.031023 | 0.030852 | 0.9291 | yes |
| Phase 2002-2007 FP_between slope | -0.021929 | -0.037591 | 0.041316 | 0.2807 | yes |
| Phase 2008-2015 FP_between slope | -0.016308 | -0.030981 | 0.033009 | 0.3546 | yes |

## Headline

**Headline (Spec A):** all phase-specific `FP_between` slopes fall inside their spatial-permutation null 95% intervals under the lag-augmented model.

