# H2 dose-response linearity diagnostic — run log

Additive diagnostic only. Does not modify `primary_model_v2` or any downstream Moran / permutation / KNN pipeline. H3 left parametric.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_sessionInfo.txt
mgcv 1.9.4

## Inputs
Loaded /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds — 10464 hauls, 158 rectangles, years 1985–2015.
phase_v2 levels: 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015
FP_between (rectangle means): min=6.223, p25=8.100, median=8.966, p75=9.855, max=10.881 (n_unique=158).
Production formula (unchanged): residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)
Diagnostic linear (mgcv): residual ~ FP_between * phase_v2 + FP_within * phase_v2 + s(stat_rec, bs="re")
Diagnostic smooth (mgcv): residual ~ phase_v2 + s(FP_between, by=phase_v2, k=5) + FP_within * phase_v2 + s(stat_rec, bs="re")

## Model fitting
Both models fit in mgcv so AIC and the nested comparison share the same likelihood / EDF penalty construction. Production glmmTMB AIC is recorded for reference only (not used for the linearity gate).
Fitting linear GAM (ML) ...
  elapsed: 2.8 sec
Fitting smooth GAM (ML) ...
  elapsed: 7.9 sec
Fitting REML copies for marginal-effect plots ...
  elapsed: 9.0 sec

## Fit comparison (ML)
Linear GAM ML AIC = 12942.88; BIC = 14055.27
Smooth GAM ML AIC = 12918.97; BIC = 14087.29
Delta AIC (smooth − linear) = -23.91  (negative favours smooth)
Nested anova(linear, smooth, test="Chisq"):
Analysis of Deviance Table

Model 1: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + s(stat_rec, 
    bs = "re")
Model 2: residual ~ phase_v2 + s(FP_between, by = phase_v2, k = K_SMOOTH) + 
    FP_within * phase_v2 + s(stat_rec, bs = "re")
  Resid. Df Resid. Dev     Df Deviance  Pr(>Chi)    
1     10298     2049.6                              
2     10289     2041.9 8.4108    7.689 7.609e-06 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
Extracted: Df=8.411, deviance/Chi=7.689, p=7.609e-06
Production glmmTMB REML AIC (reference only) = 13256.71
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_model_comparison.csv

## Per-phase smooth diagnostics
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_by_phase.csv
Per-phase edf:
  1985-1991: edf=1.001 (Ref.df=1.001); F=2.360; approx p(vs zero)=0.1244; phase_linearity_held=TRUE
  1992-2001: edf=1.000 (Ref.df=1.001); F=1.237; approx p(vs zero)=0.2661; phase_linearity_held=TRUE
  2002-2007: edf=3.454 (Ref.df=3.823); F=22.842; approx p(vs zero)=0; phase_linearity_held=FALSE
  2008-2015: edf=3.069 (Ref.df=3.543); F=6.331; approx p(vs zero)=0.0001242; phase_linearity_held=FALSE
Note: summary.gam p-values test the smooth against a zero function, not against a linear null. The linearity gate uses edf + nested AIC/anova.

## Marginal-effect curves
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_curves.csv
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_dose_response_linearity_by_phase.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_dose_response_linearity_1985_1991.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_dose_response_linearity_1992_2001.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_dose_response_linearity_2002_2007.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/h2_dose_response_linearity_2008_2015.png

## Decision gate
Decision flag: NONLINEARITY_DETECTED
any_edf_high=TRUE; nested_sig=TRUE; aic_favours_smooth=TRUE

### Summary paragraph
H2 dose-response linearity diagnostic (mgcv, same haul data and (1|stat_rec) structure as primary_model_v2; H3 left linear): ML AIC linear = 12942.9, smooth = 12919.0 (Δ = -23.9); nested Chi-sq test p = 7.609e-06. 1985-1991: linear assumption held (edf = 1.00 ≈ 1); 1992-2001: linear assumption held (edf = 1.00 ≈ 1); 2002-2007: linearity NOT held (edf = 3.45 > 1.5; shape ≈ non-monotonic (mid-range dip / uneven rise)); 2008-2015: linearity NOT held (edf = 3.07 > 1.5; shape ≈ non-monotonic (mid-range dip / uneven rise)). Overall verdict: nonlinearity is detected in at least one phase or by the nested fit comparison — production model unchanged pending design review.

### Recommended next step (do not implement here)
STOP — do not respecify the production model from this script. Review the phase-wise curves (threshold vs saturation vs non-monotonic) before any redesign of the fixed-effect form, phase interaction, or re-running the Moran / permutation / KNN residual pipeline.
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_summary.md
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2_dose_response_linearity_model_objects.rds
