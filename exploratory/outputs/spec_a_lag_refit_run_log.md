# Spec A — lagged FP_between refit — run log

Adds FP_between_lag * phase_v2 to the current primary model (single-change principle). No year term; (1 | stat_rec) unchanged.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_lag_refit_sessionInfo.txt

## Inputs
Loaded primary v2: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015
Loaded k-NN weights: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/knn_listw_k4.rds 
k = 4; panel N = 158 
Note: 15/158 k-th neighbour ties are lattice artifacts (equidistant opposite diagonals on the 0.5°×1° ICES grid), resolved by lowest `stat_rec` — not a distance bug.

## FP_between_lag construction
FP_between_lag range [6.2870, 10.3340]; cor(FP_between, FP_between_lag) = 0.8276
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_lag_rectangle.rds

## VIF (FP_between vs FP_between_lag, rectangle level)
VIF = 3.174 (r = 0.8276)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_vif_fp_between_lag.csv
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_vif_fp_between_lag.md

## Model refit (Spec A)
Formula: residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + (1 | stat_rec) [REML]
Refit elapsed: 0.46 sec.
Fitted formula: residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + (1 | stat_rec)

## BLUP spatial diagnostic (queen listw; baseline Moran I ≈ 0.552)
primary_v2 BLUP Moran I = 0.550265; Spec A BLUP Moran I = 0.531810 (Δ = -0.018455)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_blup_spatial_diagnostic.csv

## Outputs
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_a.rds
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_lag_rectangle.rds
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_vif_fp_between_lag.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_vif_fp_between_lag.md
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_a_blup_spatial_diagnostic.csv
