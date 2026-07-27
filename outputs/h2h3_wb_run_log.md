# H2/H3 within-between fishing-pressure decomposition — run log

Separates between-rectangle (H2) and within-rectangle (H3) fishing-pressure variation. Biomass excluded. Same panel / phase / (1|stat_rec) as the biomass-free primary results run. No ecological interpretation beyond labelling which coefficients answer which hypothesis.

## Session
sessionInfo written to: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_sessionInfo.txt
glmmTMB 1.1.14; spaMM available = TRUE (4.6.1)

## Data and FP within-between transformation
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015 (reused from biomass-free results RDS).
FP_between: time-invariant rectangle mean of log(hours+1); var = 1.1908; range = [6.223, 10.881].
FP_within: haul-level deviation from rectangle mean; var = 0.7043; range = [-5.538, 2.795].

### Sanity check: FP_between vs FP_within correlation
Correlation(FP_between, FP_within) = -7.42842e-16 (expected ≈ 0 by construction).
Max |within-rectangle mean of FP_within| = 1.385e-15; mean abs = 4.092e-16 (expected ≈ 0).
Sanity check passed: between and within components are effectively uncorrelated.

## Primary model fit
Formula: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec) [REML]
H2 answered by FP_between (+ phase interactions / phase-specific between slopes). H3 answered by FP_within (+ phase interactions / phase-specific within slopes).
Primary within-between fit time: 0.88 sec.

## Sensitivity: CAR with decomposed terms
Fitting CAR: residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)
CAR within-between fit time: 1.16 sec.
Saved CAR FE table: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_car_fixed_effects.csv

## Deferred sensitivity
GAM / continuous-year version of the within-between decomposition is DEFERRED (timeline). Not silently dropped — flag here for a follow-up if needed.
Saved model objects: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_model_objects.rds

## 1. Primary fixed-effects table (H2/H3 labelled)
Fixed effects:
  - (Intercept)                         est=-1.2768 SE=0.1646 z=-7.76 p=8.593e-15  [control]
  - FP_between                          est=-0.0060 SE=0.0180 z=-0.33 p=7.405e-01  [H2_spatial_between]
  - phase1989-2000                      est=+0.0110 SE=0.1390 z=+0.08 p=9.366e-01  [control]
  - phase2001-2007                      est=-1.1232 SE=0.1471 z=-7.64 p=2.210e-14  [control]
  - phase2008-2015                      est=-0.3202 SE=0.1463 z=-2.19 p=2.858e-02  [control]
  - FP_within                           est=+0.0317 SE=0.0182 z=+1.75 p=8.071e-02  [H3_temporal_within]
  - FP_between:phase1989-2000           est=-0.0074 SE=0.0149 z=-0.49 p=6.212e-01  [H2_spatial_between]
  - FP_between:phase2001-2007           est=+0.1315 SE=0.0159 z=+8.27 p=1.366e-16  [H2_spatial_between]
  - FP_between:phase2008-2015           est=+0.0580 SE=0.0157 z=+3.69 p=2.217e-04  [H2_spatial_between]
  - phase1989-2000:FP_within            est=+0.0097 SE=0.0232 z=+0.42 p=6.745e-01  [H3_temporal_within]
  - phase2001-2007:FP_within            est=-0.0433 SE=0.0249 z=-1.74 p=8.211e-02  [H3_temporal_within]
  - phase2008-2015:FP_within            est=+0.0254 SE=0.0228 z=+1.11 p=2.653e-01  [H3_temporal_within]
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_primary_fixed_effects.csv

Phase-specific slopes (primary):
  - [H2_spatial_between] FP_between / 1985-1988: slope=-0.0060 SE=0.0180 CI=[-0.0412, +0.0293] p=7.405e-01
  - [H2_spatial_between] FP_between / 1989-2000: slope=-0.0133 SE=0.0139 CI=[-0.0406, +0.0139] p=3.382e-01
  - [H2_spatial_between] FP_between / 2001-2007: slope=+0.1256 SE=0.0150 CI=[+0.0963, +0.1549] p=4.655e-17
  - [H2_spatial_between] FP_between / 2008-2015: slope=+0.0520 SE=0.0147 CI=[+0.0233, +0.0807] p=3.869e-04
  - [H3_temporal_within] FP_within / 1985-1988: slope=+0.0317 SE=0.0182 CI=[-0.0039, +0.0674] p=8.071e-02
  - [H3_temporal_within] FP_within / 1989-2000: slope=+0.0415 SE=0.0154 CI=[+0.0113, +0.0716] p=6.997e-03
  - [H3_temporal_within] FP_within / 2001-2007: slope=-0.0116 SE=0.0158 CI=[-0.0425, +0.0193] p=4.619e-01
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0571 SE=0.0122 CI=[+0.0332, +0.0811] p=2.947e-06
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_fp_slopes_by_phase.csv

## 2. Joint Wald tests
  - H2_all_between_phase_slopes_joint_zero: chi2(4)=178.288, p=1.739e-37
  - H3_all_within_phase_slopes_joint_zero: chi2(4)=35.556, p=3.571e-07
  - H3_within_phase_interactions_joint_zero: chi2(3)=12.656, p=5.443e-03
  - H2_between_phase_interactions_joint_zero: chi2(3)=169.334, p=1.772e-36
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_wald_tests.csv

## 3. Comparison to blended log_hours_total model
Per-phase blended vs decomposed slopes:
  - 1985-1988: blended=-0.0008 | between(H2)=-0.0060 (p=7.405e-01) | within(H3)=+0.0317 (p=8.071e-02)
  - 1989-2000: blended=-0.0056 | between(H2)=-0.0133 (p=3.382e-01) | within(H3)=+0.0415 (p=6.997e-03)
  - 2001-2007: blended=+0.0885 | between(H2)=+0.1256 (p=4.655e-17) | within(H3)=-0.0116 (p=4.619e-01)
  - 2008-2015: blended=+0.0580 | between(H2)=+0.0520 (p=3.869e-04) | within(H3)=+0.0571 (p=2.947e-06)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_blended_comparison.csv

## 5. Diagnostics (primary within-between model)
Convergence: code=0, Hessian PD=TRUE, max|grad|=9.717e-08, any NA SE=FALSE
Nakagawa R2: marginal=0.0513, conditional=0.1711; rect SD=0.1695; resid SD=0.4457
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_model_fit.csv
Partial pooling: Spearman cor(log n_hauls, shrinkage_ratio)=0.351 across 158 rectangles.
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_partial_pooling.csv
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_partial_pooling.png
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_fp_slopes_by_phase.png

## Coefficient → hypothesis mapping (no further interpretation)
H2 (spatial / persistent pressure): FP_between main effect, FP_between × phase interactions, and the four phase-specific FP_between slopes.
H3 (temporal / own-rectangle deviation): FP_within main effect, FP_within × phase interactions, and the four phase-specific FP_within slopes. The joint test of FP_within × phase interactions is the direct parallel to the previous blended interaction test for change across phases.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_primary_fixed_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_fp_slopes_by_phase.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_wald_tests.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_blended_comparison.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_partial_pooling.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_model_fit.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_car_fixed_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_partial_pooling.png
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_fp_slopes_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_sessionInfo.txt
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_run_log.md (this file)

---

## Live pipeline record

This within-between run is part of the **live** H2/H3 pipeline. Canonical live-script list: [`live_pipeline_run_log.md`](live_pipeline_run_log.md). Exploratory / superseded material: [`../exploratory/`](../exploratory/).
