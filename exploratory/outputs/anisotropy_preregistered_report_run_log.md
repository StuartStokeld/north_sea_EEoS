# Anisotropy pre-registered report — run log

**Date:** 2026-08-04 11:02 BST
**Inputs:** `outputs/bearings_by_rectangle.csv`, `outputs/bathymetry_by_rectangle.csv`
**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`

## Pre-registered rules applied (exact)

1. **TID:** exclude rectangles with pct_sounding < 50% from confirmatory tests; keep them in the bearing plot. Report 25% as labeled sensitivity only.
2. **Primary object:** signed residual-correlation bearing. BLUP = Table 2 robustness only.
3. **Along-shelf axis (Table 1):** per-rectangle local depth-gradient direction rotated +90° (along-shelf). Global compass split = Table 2 foil only.
4. **Decision:** Jammalamadaka–Sarma via `circular::cor.circular(..., test=TRUE)`. Table 1 verdict uses **asymptotic p < 0.05**. Permutation p (999 shuffles of the bearing pairing) reported alongside; coefficient magnitude is descriptive only.

### Foil operationalization (Table 2)
Global compass along-shelf (depth-independent): 90° (N–S) for lon < 4°E; 0° (E–W) for lon ≥ 4°E. Same JS test and p < 0.05 labeling as primary, but **not confirmatory**.

N rectangles: 158
TID flagged <50%: 6
TID flagged <25%: 0

## Table 1 computation
Pairing: bearing_resid × bearing_depth (= grad_bearing_along = depth-gradient + 90°). Sample: TID-OK at 50% threshold.
Primary: n=152  rho=0.0483  p_asymp=0.5523  p_perm=0.5616  confirmed=No

## Table 2 computations
BLUP robustness: n=152  rho=0.0173  p_asymp=0.8292  p_perm=0.8398  label=No
Global compass foil: n=152  rho=0.1239  p_asymp=0.1256  p_perm=0.0751  label=No
TID 25% sensitivity: n=158  rho=0.0371  p_asymp=0.6411  p_perm=0.6316  label=No

## Verdict branch
verdict_label: not_confirmed
narrative: No evidence that shelf-geometry-based direction explains the earlier unidentified decay range. Consistent with either a genuinely borderless spatial field or confounding with the fixed-effect spatial trend — does not distinguish between these two.

Wrote plot: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/anisotropy_bearing_comparison_plot.png
Plot includes all 158 rectangles; 6 TID-flagged (<50%) shown as triangles.
Wrote summary: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/anisotropy_results_summary.md
