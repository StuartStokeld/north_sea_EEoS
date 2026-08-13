# Spec B / Spec A+B — lagged-outcome refit — run log

Spec B adds pooled B_lag_neighbour (neighbour mean ln_B_obs) to the primary. Spec A+B adds both FP_between_lag * phase_v2 and B_lag_neighbour. No year term; (1 | stat_rec) unchanged.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_lag_refit_sessionInfo.txt

## Inputs
Primary v2: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds
Haul predictions: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/haul_eeos_predictions.rds
k-NN weights: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/knn_listw_k4.rds (k=4)
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015
Loaded Spec A fit: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_a.rds

## B_lag_neighbour construction
Rectangle-year table: 4898 rows (158 × 31 years); with biomass: 4771
Analysis rect-years: 4768; k-NN full=4491; partial=277; none=0
Neighbour-source-only rect-years (biomass present, not in analysis data): 3
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/b_lag_neighbour_rectangle_year.rds
Joined B_lag_neighbour to hauls: range [8.6830, 14.4124]; mean n_neighbours_used=3.949

## VIF diagnostics
Spec B max VIF = 1.250 (B_lag_neighbour)
Spec A+B max VIF = 3.641 (FP_between_lag)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_vif.csv / /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_ab_vif.csv

## Model fits
Spec B formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + B_lag_neighbour +      (1 | stat_rec)
Spec B elapsed: 0.40 sec.
Spec A+B formula: residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + B_lag_neighbour + (1 | stat_rec)
Spec A+B elapsed: 0.46 sec.
B_lag_neighbour coef — Spec B: -0.135164; Spec A+B: -0.142827

## Spatial diagnostics (queen listw; BLUP + rectangle-mean residual)
| Model | Level | Moran's I | Geary's C | p (Moran) |
|---|---|---|---|---|
| primary_v2 | BLUP | 0.5503 | 0.4422 | 4.57e-39 |
| primary_v2 | rectangle_mean_residual | 0.4742 | 0.5020 | 6.84e-30 |
| spec_a | BLUP | 0.5318 | 0.4619 | 1.52e-36 |
| spec_a | rectangle_mean_residual | 0.4615 | 0.5136 | 2.76e-28 |
| spec_b | BLUP | 0.4216 | 0.5719 | 8.79e-24 |
| spec_b | rectangle_mean_residual | 0.3262 | 0.6408 | 3.74e-15 |
| spec_ab | BLUP | 0.4186 | 0.5730 | 1.75e-23 |
| spec_ab | rectangle_mean_residual | 0.3378 | 0.6251 | 4.11e-16 |
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_spatial_diagnostic.csv

## Coefficient comparison (FP_between / FP_within phase slopes)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_coefficient_comparison.csv
  spec_b H2 1985-1991: primary=-0.0247 → spec_b=-0.0624 (Δ=-0.0377)
  spec_b H2 1992-2001: primary=+0.0142 → spec_b=-0.0413 (Δ=-0.0555)
  spec_b H2 2002-2007: primary=+0.1223 → spec_b=+0.0760 (Δ=-0.0463)
  spec_b H2 2008-2015: primary=+0.0519 → spec_b=+0.0258 (Δ=-0.0261)
  spec_ab H2 1985-1991: primary=-0.0247 → spec_ab=+0.0220 (Δ=+0.0467)
  spec_ab H2 1992-2001: primary=+0.0142 → spec_ab=+0.0101 (Δ=-0.0040)
  spec_ab H2 2002-2007: primary=+0.1223 → spec_ab=-0.0174 (Δ=-0.1397)
  spec_ab H2 2008-2015: primary=+0.0519 → spec_ab=-0.0124 (Δ=-0.0643)

## Outputs
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_b.rds
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_ab.rds
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_lag_refit_summary.md
