# Spec B / Spec A+B — lagged-outcome summary

## Formulas

- Spec B: `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + B_lag_neighbour +      (1 | stat_rec)`
- Spec A+B: `residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + B_lag_neighbour + (1 | stat_rec)`

## B_lag_neighbour

- Analysis rect-years: 4768 (full k=4: 4491; partial: 277)
- Neighbour-source-only rect-years (incl. Couce gaps): 3
- Spec B `B_lag_neighbour` coef: **-0.135164**
- Spec A+B `B_lag_neighbour` coef: **-0.142827**

## Moran's I (queen listw)

| Model | BLUP I | Residual I |
|-------|--------|------------|
| primary_v2 | 0.5503 | 0.4742 |
| spec_a | 0.5318 | 0.4615 |
| spec_b | 0.4216 | 0.3262 |
| spec_ab | 0.4186 | 0.3378 |

Baselines: Finding 1 BLUP I ≈ 0.552; Spec A BLUP I ≈ 0.532.

## Outputs

- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/b_lag_neighbour_rectangle_year.rds`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_b.rds`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_ab.rds`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_spatial_diagnostic.csv`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/spec_b_coefficient_comparison.csv`

