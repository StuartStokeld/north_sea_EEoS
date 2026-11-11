# H2/H3 within-between fishing-pressure decomposition — run log

Separates between-rectangle (H2) and within-rectangle (H3) fishing-pressure variation. Biomass excluded. Same panel / phase / (1|stat_rec) as the biomass-free primary results run. No ecological interpretation beyond labelling which coefficients answer which hypothesis.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_sessionInfo.txt
glmmTMB 1.1.14; spaMM available = TRUE (4.6.65)

## Data and FP within-between transformation
Analysis data: 8670 hauls, 158 rectangles, years 1985–2015 (rebuilt via build_feasibility_data() from current haul_eeos_predictions.rds; rectangles match CAR adjacency).
FP_between: time-invariant rectangle mean of log(hours+1); var = 1.2652; range = [5.974, 10.824].
FP_within: haul-level deviation from rectangle mean; var = 0.6793; range = [-5.205, 2.798].

### Sanity check: FP_between vs FP_within correlation
Correlation(FP_between, FP_within) = 3.34165e-16 (expected ≈ 0 by construction).
Max |within-rectangle mean of FP_within| = 1.197e-15; mean abs = 4.883e-16 (expected ≈ 0).
Sanity check passed: between and within components are effectively uncorrelated.

## Primary model fit
Formula: residual ~ FP_between * phase + FP_within * phase + (1 | stat_rec) [REML]
H2 answered by FP_between (+ phase interactions / phase-specific between slopes). H3 answered by FP_within (+ phase interactions / phase-specific within slopes).
Primary within-between fit time: 0.30 sec.

## Sensitivity: CAR with decomposed terms
Fitting CAR: residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)
CAR within-between fit time: 0.29 sec.
Saved CAR FE table: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_car_fixed_effects.csv

## Deferred sensitivity
GAM / continuous-year version of the within-between decomposition is DEFERRED (timeline). Not silently dropped — flag here for a follow-up if needed.
Saved model objects: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_objects.rds

## 1. Primary fixed-effects table (H2/H3 labelled)
Fixed effects:
  - (Intercept)                         est=-0.8094 SE=0.1309 z=-6.18 p=6.332e-10  [control]
  - FP_between                          est=+0.0253 SE=0.0143 z=+1.76 p=7.796e-02  [H2_spatial_between]
  - phase1989-2000                      est=-0.1952 SE=0.1037 z=-1.88 p=5.966e-02  [control]
  - phase2001-2007                      est=-0.2168 SE=0.1072 z=-2.02 p=4.307e-02  [control]
  - phase2008-2015                      est=+0.2451 SE=0.1071 z=+2.29 p=2.211e-02  [control]
  - FP_within                           est=-0.0165 SE=0.0132 z=-1.25 p=2.115e-01  [H3_temporal_within]
  - FP_between:phase1989-2000           est=+0.0141 SE=0.0111 z=+1.27 p=2.033e-01  [H2_spatial_between]
  - FP_between:phase2001-2007           est=+0.0138 SE=0.0116 z=+1.20 p=2.316e-01  [H2_spatial_between]
  - FP_between:phase2008-2015           est=-0.0295 SE=0.0115 z=-2.56 p=1.034e-02  [H2_spatial_between]
  - phase1989-2000:FP_within            est=+0.0321 SE=0.0167 z=+1.92 p=5.474e-02  [H3_temporal_within]
  - phase2001-2007:FP_within            est=-0.0226 SE=0.0179 z=-1.26 p=2.065e-01  [H3_temporal_within]
  - phase2008-2015:FP_within            est=+0.0497 SE=0.0166 z=+2.99 p=2.758e-03  [H3_temporal_within]
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_primary_fixed_effects.csv

Phase-specific slopes (primary):
  - [H2_spatial_between] FP_between / 1985-1988: slope=+0.0253 SE=0.0143 CI=[-0.0028, +0.0534] p=7.796e-02
  - [H2_spatial_between] FP_between / 1989-2000: slope=+0.0394 SE=0.0115 CI=[+0.0169, +0.0619] p=6.014e-04
  - [H2_spatial_between] FP_between / 2001-2007: slope=+0.0391 SE=0.0118 CI=[+0.0159, +0.0623] p=9.522e-04
  - [H2_spatial_between] FP_between / 2008-2015: slope=-0.0042 SE=0.0117 CI=[-0.0272, +0.0188] p=7.218e-01
  - [H3_temporal_within] FP_within / 1985-1988: slope=-0.0165 SE=0.0132 CI=[-0.0423, +0.0094] p=2.115e-01
  - [H3_temporal_within] FP_within / 1989-2000: slope=+0.0157 SE=0.0111 CI=[-0.0060, +0.0373] p=1.563e-01
  - [H3_temporal_within] FP_within / 2001-2007: slope=-0.0391 SE=0.0112 CI=[-0.0610, -0.0171] p=4.788e-04
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0332 SE=0.0089 CI=[+0.0157, +0.0508] p=2.034e-04
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_fp_slopes_by_phase.csv

## 2. Joint Wald tests
  - H2_all_between_phase_slopes_joint_zero: chi2(4)=44.508, p=5.031e-09
  - H3_all_within_phase_slopes_joint_zero: chi2(4)=30.744, p=3.452e-06
  - H3_within_phase_interactions_joint_zero: chi2(3)=30.150, p=1.283e-06
  - H2_between_phase_interactions_joint_zero: chi2(3)=39.003, p=1.733e-08
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_wald_tests.csv

## 3. Comparison to blended log_hours_total model
Per-phase blended vs decomposed slopes:
  - 1985-1988: blended=-0.0008 | between(H2)=+0.0253 (p=7.796e-02) | within(H3)=-0.0165 (p=2.115e-01)
  - 1989-2000: blended=-0.0056 | between(H2)=+0.0394 (p=6.014e-04) | within(H3)=+0.0157 (p=1.563e-01)
  - 2001-2007: blended=+0.0885 | between(H2)=+0.0391 (p=9.522e-04) | within(H3)=-0.0391 (p=4.788e-04)
  - 2008-2015: blended=+0.0580 | between(H2)=-0.0042 (p=7.218e-01) | within(H3)=+0.0332 (p=2.034e-04)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_blended_comparison.csv

## 5. Diagnostics (primary within-between model)
Convergence: code=0, Hessian PD=TRUE, max|grad|=3.739e-05, any NA SE=FALSE
Nakagawa R2: marginal=0.0170, conditional=0.2220; rect SD=0.1500; resid SD=0.2923
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_fit.csv
Partial pooling: Spearman cor(log n_hauls, shrinkage_ratio)=0.369 across 158 rectangles.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_partial_pooling.csv
Saved figure: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2h3_wb_partial_pooling.png
Saved figure: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2h3_wb_fp_slopes_by_phase.png

## Coefficient → hypothesis mapping (no further interpretation)
H2 (spatial / persistent pressure): FP_between main effect, FP_between × phase interactions, and the four phase-specific FP_between slopes.
H3 (temporal / own-rectangle deviation): FP_within main effect, FP_within × phase interactions, and the four phase-specific FP_within slopes. The joint test of FP_within × phase interactions is the direct parallel to the previous blended interaction test for change across phases.

## Outputs
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_primary_fixed_effects.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_fp_slopes_by_phase.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_wald_tests.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_blended_comparison.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_partial_pooling.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_fit.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_model_objects.rds
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_car_fixed_effects.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2h3_wb_partial_pooling.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2h3_wb_fp_slopes_by_phase.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_sessionInfo.txt
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_wb_run_log.md (this file)
