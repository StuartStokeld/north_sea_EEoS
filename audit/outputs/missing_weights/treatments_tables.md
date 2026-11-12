## H1 metrics by treatment

| treatment | n_hauls | log_r2_1to1 | bias_mean_ln | log_r2_bias_removed | median_ratio_pred_obs | cor2_ln |
|---|---|---|---|---|---|---|
| a current (NA weight counted as 0) | 10,128 | 0.6610 | -0.6523 | 0.9282 | 1.958 | 0.9713 |
| b1 exclude hauls missing >1% biomass | 10,034 | 0.6599 | -0.6536 | 0.9281 | 1.961 | 0.9712 |
| b5 exclude hauls missing >5% biomass | 10,098 | 0.6604 | -0.6529 | 0.9282 | 1.959 | 0.9713 |
| b10 exclude hauls missing >10% biomass | 10,113 | 0.6608 | -0.6526 | 0.9282 | 1.959 | 0.9713 |
| c exclude all affected hauls | 8,979 | 0.6693 | -0.6518 | 0.9297 | 1.960 | 0.9716 |
| d fill missing taxa from DATRAS lengths x LW | 10,128 | 0.6611 | -0.6517 | 0.9280 | 1.957 | 0.9712 |

## CAR slopes by treatment (slope (SE); * raw p<0.05; ** Bonferroni p<0.05)

| treatment | response | component | 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015 | n_hauls_h2h3 | n_rect |
|---|---|---|---|---|---|---|---|---|
| a current (NA weight counted as 0) | residual | FP_between | +0.004 (0.014) | +0.026 (0.013) * | +0.013 (0.013) | -0.025 (0.013) * | 8,670 | 158 |
| a current (NA weight counted as 0) | residual | FP_within | -0.011 (0.011) | +0.034 (0.012) ** | -0.032 (0.012) ** | +0.032 (0.009) ** | 8,670 | 158 |
| a current (NA weight counted as 0) | abs_residual | FP_between | -0.006 (0.012) | -0.030 (0.011) ** | -0.016 (0.012) | +0.021 (0.011) | 8,670 | 158 |
| a current (NA weight counted as 0) | abs_residual | FP_within | +0.012 (0.009) | -0.035 (0.011) ** | +0.027 (0.011) * | -0.031 (0.008) ** | 8,670 | 158 |
| b1 exclude hauls missing >1% biomass | residual | FP_between | +0.004 (0.014) | +0.025 (0.013) | +0.012 (0.013) | -0.026 (0.013) * | 8,589 | 158 |
| b1 exclude hauls missing >1% biomass | residual | FP_within | -0.014 (0.011) | +0.033 (0.012) ** | -0.031 (0.012) ** | +0.032 (0.009) ** | 8,589 | 158 |
| b1 exclude hauls missing >1% biomass | abs_residual | FP_between | -0.006 (0.012) | -0.029 (0.011) ** | -0.016 (0.012) | +0.021 (0.011) | 8,589 | 158 |
| b1 exclude hauls missing >1% biomass | abs_residual | FP_within | +0.015 (0.010) | -0.034 (0.011) ** | +0.027 (0.011) * | -0.031 (0.008) ** | 8,589 | 158 |
| b5 exclude hauls missing >5% biomass | residual | FP_between | +0.004 (0.014) | +0.025 (0.013) | +0.012 (0.013) | -0.026 (0.013) * | 8,640 | 158 |
| b5 exclude hauls missing >5% biomass | residual | FP_within | -0.012 (0.011) | +0.035 (0.012) ** | -0.032 (0.012) ** | +0.032 (0.009) ** | 8,640 | 158 |
| b5 exclude hauls missing >5% biomass | abs_residual | FP_between | -0.006 (0.012) | -0.029 (0.011) ** | -0.016 (0.012) | +0.021 (0.011) | 8,640 | 158 |
| b5 exclude hauls missing >5% biomass | abs_residual | FP_within | +0.013 (0.010) | -0.036 (0.011) ** | +0.027 (0.011) * | -0.031 (0.008) ** | 8,640 | 158 |
| b10 exclude hauls missing >10% biomass | residual | FP_between | +0.004 (0.014) | +0.025 (0.013) * | +0.012 (0.013) | -0.026 (0.013) * | 8,655 | 158 |
| b10 exclude hauls missing >10% biomass | residual | FP_within | -0.011 (0.011) | +0.035 (0.012) ** | -0.032 (0.012) ** | +0.032 (0.009) ** | 8,655 | 158 |
| b10 exclude hauls missing >10% biomass | abs_residual | FP_between | -0.006 (0.012) | -0.029 (0.011) ** | -0.016 (0.012) | +0.021 (0.011) | 8,655 | 158 |
| b10 exclude hauls missing >10% biomass | abs_residual | FP_within | +0.012 (0.010) | -0.036 (0.011) ** | +0.027 (0.011) * | -0.031 (0.008) ** | 8,655 | 158 |
| c exclude all affected hauls | residual | FP_between | +0.001 (0.014) | +0.025 (0.013) | +0.011 (0.013) | -0.025 (0.013) | 8,070 | 158 |
| c exclude all affected hauls | residual | FP_within | -0.020 (0.011) | +0.039 (0.013) ** | -0.030 (0.013) * | +0.032 (0.009) ** | 8,070 | 158 |
| c exclude all affected hauls | abs_residual | FP_between | -0.001 (0.012) | -0.029 (0.011) ** | -0.015 (0.012) | +0.020 (0.011) | 8,070 | 158 |
| c exclude all affected hauls | abs_residual | FP_within | +0.023 (0.010) * | -0.040 (0.011) ** | +0.025 (0.011) * | -0.031 (0.008) ** | 8,070 | 158 |
| d fill missing taxa from DATRAS lengths x LW | residual | FP_between | +0.004 (0.014) | +0.026 (0.013) * | +0.013 (0.013) | -0.026 (0.013) * | 8,670 | 158 |
| d fill missing taxa from DATRAS lengths x LW | residual | FP_within | -0.011 (0.011) | +0.034 (0.012) ** | -0.032 (0.012) ** | +0.032 (0.009) ** | 8,670 | 158 |
| d fill missing taxa from DATRAS lengths x LW | abs_residual | FP_between | -0.006 (0.012) | -0.030 (0.011) ** | -0.016 (0.012) | +0.021 (0.011) | 8,670 | 158 |
| d fill missing taxa from DATRAS lengths x LW | abs_residual | FP_within | +0.012 (0.009) | -0.035 (0.011) ** | +0.027 (0.011) * | -0.031 (0.008) ** | 8,670 | 158 |
