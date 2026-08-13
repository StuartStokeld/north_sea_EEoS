## H2 haul-count vs metric diagnostics

Generated: 2026-07-30 11:02

Rectangle panel: 158 rectangles, n_hauls range 16–111 (median 66).
Mean-trend escalate threshold: |Spearman rho| >= 0.25 (central tendency only).

Live WB universe: 10464 hauls in 158 rectangles.

### Partial pooling (existing artefact)
Spearman(log n_hauls, shrinkage_ratio) = 0.351 (p = 6.25e-06, n = 158).
Shrinkage rising with n is the model answer to unequal precision;
it is separate from whether the DV mean trends with n.

### Decision rule
Wedge / convergence (SE or SD vs n declining): YES — expected pattern present
Linear/monotone mean-metric trend (|rho| >= 0.25): YES — escalate (sensitivity)

Flagged mean-trend rows:
  - wb_phase / 2001-2007 / mean_abs_residual: Spearman rho = -0.412 (p = 7.24e-08)
  - wb_phase / 2001-2007 / mean_residual: Spearman rho = 0.415 (p = 5.89e-08)

### Sensitivity (escalation)
Rectangle-level OLS of residual metric ~ FP_between: unweighted, n-weighted, inverse-SE^2 weighted, and lowest-n-decile dropped. Descriptive robustness only — not a replacement for the haul-level WB primary.

Phase 2001-2007: n_hauls decile cut = 12.0 (drop n <= cut); retained 140 / 158 rectangles.
  After drop — mean_abs vs log n: Spearman rho = -0.516 (flag = TRUE)
  After drop — mean_residual vs log n: Spearman rho = 0.519 (flag = TRUE)

Sensitivity slopes (residual metric ~ FP_between):
  panel_unweighted | all | mean_abs_residual: slope = -0.0319 (SE = 0.0123, p = 0.0102, n = 158)
  panel_weight_n | all | mean_abs_residual: slope = -0.0308 (SE = 0.0123, p = 0.0137, n = 158)
  panel_weight_inv_se2 | all | mean_abs_residual: slope = -0.0370 (SE = 0.0125, p = 0.00346, n = 158)
  panel_unweighted | all | mean_residual: slope = 0.0332 (SE = 0.0129, p = 0.0109, n = 158)
  panel_weight_n | all | mean_residual: slope = 0.0320 (SE = 0.0130, p = 0.0147, n = 158)
  panel_weight_inv_se2 | all | mean_residual: slope = 0.0380 (SE = 0.0130, p = 0.00403, n = 158)
  phase_unweighted | 2001-2007 | mean_abs_residual: slope = -0.1196 (SE = 0.0172, p = 8.66e-11, n = 158)
  phase_weight_n | 2001-2007 | mean_abs_residual: slope = -0.1208 (SE = 0.0170, p = 3.69e-11, n = 158)
  phase_weight_inv_se2 | 2001-2007 | mean_abs_residual: slope = -0.0830 (SE = 0.0165, p = 1.36e-06, n = 158)
  phase_drop_lowest_n_decile | 2001-2007 | mean_abs_residual: slope = -0.1272 (SE = 0.0176, p = 3.03e-11, n = 140)
  phase_unweighted | 2001-2007 | mean_residual: slope = 0.1205 (SE = 0.0179, p = 3.29e-10, n = 158)
  phase_weight_n | 2001-2007 | mean_residual: slope = 0.1221 (SE = 0.0177, p = 1.26e-10, n = 158)
  phase_weight_inv_se2 | 2001-2007 | mean_residual: slope = 0.0840 (SE = 0.0170, p = 1.9e-06, n = 158)
  phase_drop_lowest_n_decile | 2001-2007 | mean_residual: slope = 0.1300 (SE = 0.0181, p = 3.61e-11, n = 140)

Robustness: weighted / drop-low-n slopes keep the same sign and similar magnitude as unweighted phase OLS — H2 FP_between association for flagged phase(s) is not an artefact of sparse rectangles. No change to primary WB model.

### Key Spearman results (central tendency)
  rectangle_panel | all | mean_abs_residual: rho = 0.024 (p = 0.763), OLS slope = 0.0800, flag = FALSE
  rectangle_panel | all | mean_residual: rho = -0.018 (p = 0.819), OLS slope = -0.0796, flag = FALSE
  wb_universe | all | mean_abs_residual: rho = 0.027 (p = 0.737), OLS slope = 0.0795, flag = FALSE
  wb_universe | all | mean_residual: rho = -0.021 (p = 0.793), OLS slope = -0.0790, flag = FALSE
  wb_phase | 1985-1988 | mean_abs_residual: rho = 0.210 (p = 0.00853), OLS slope = 0.0650, flag = FALSE
  wb_phase | 1985-1988 | mean_residual: rho = -0.204 (p = 0.0107), OLS slope = -0.0643, flag = FALSE
  wb_phase | 1989-2000 | mean_abs_residual: rho = 0.021 (p = 0.793), OLS slope = 0.0537, flag = FALSE
  wb_phase | 1989-2000 | mean_residual: rho = -0.019 (p = 0.814), OLS slope = -0.0498, flag = FALSE
  wb_phase | 2001-2007 | mean_abs_residual: rho = -0.412 (p = 7.24e-08), OLS slope = -0.1499, flag = TRUE
  wb_phase | 2001-2007 | mean_residual: rho = 0.415 (p = 5.89e-08), OLS slope = 0.1530, flag = TRUE
  wb_phase | 2008-2015 | mean_abs_residual: rho = -0.085 (p = 0.291), OLS slope = 0.0106, flag = FALSE
  wb_phase | 2008-2015 | mean_residual: rho = 0.086 (p = 0.284), OLS slope = -0.0136, flag = FALSE

### Key Spearman results (dispersion / wedge)
  rectangle_panel | sd_abs_residual: rho = 0.243 (p = 0.00206)
  rectangle_panel | se_abs_residual: rho = -0.301 (p = 0.000121)
  wb_universe | sd_abs_residual: rho = 0.241 (p = 0.00227)
  wb_universe | sd_residual: rho = 0.213 (p = 0.00713)
  wb_universe | se_abs_residual: rho = -0.305 (p = 9.71e-05)

### Outputs
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_n_hauls_metric_diagnostics.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_n_hauls_metric_bins.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_n_hauls_metric_sensitivity.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_n_hauls_vs_metric.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_n_hauls_wedge_sd.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_n_hauls_vs_metric_by_phase.png

### Verdict summary
Unequal haul counts are handled by partial pooling for inference. Overall rectangle metrics show no mean trend with n_hauls and the expected SE wedge. A phase-specific mean trend (see flagged rows) triggered weighted / drop-low-n sensitivity of residual ~ FP_between; see sensitivity table and robustness note above. Primary haul-level WB model unchanged.
