# H2/H3 shared hierarchical model — temporal robustness check run log
## Categorical phase vs continuous linear year vs smooth year

FEASIBILITY/ROBUSTNESS CHECK ONLY. No fishing-pressure or fishing-pressure x time coefficient below is treated as the final H2/H3 answer. This task does not finalise the time structure as the committed approach — that follows supervisor discussion, not this script.

## Tool and formula choices
Base specification: the plain non-spatial (1 | stat_rec) random intercept recommended as the primary model after the spatial feasibility checks (Rounds 1–2) — NOT the CAR or continuous-exp spatial models. Categorical-phase simple_nonspatial fit is REUSED from Round 1 (outputs/h2h3_feasibility_model_objects.rds), not refit under REML.
Continuous linear year: glmmTMB, residual ~ log_hours_total * year_c + mean_ln_B_obs + (1 | stat_rec), REML = TRUE. year_c = year - mean(year) over the analysis hauls (mean ≈ 1999.5) for numerical stability; the centring constant is recorded and used when recovering the fishing-pressure slope as a function of calendar year.
Smooth year (GAM): mgcv::gam (version 1.9.1), residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) + mean_ln_B_obs + s(stat_rec, bs = "re"), method = "REML". This is a VARYING-COEFFICIENT form: the by-smooth contributes log_hours_total * f(year) to the linear predictor, so f(year) IS the fishing-pressure effect as a function of year (recovered via predict(type = "terms") at log_hours_total = 1). Chosen over the brief's ti(log_hours_total, year) example for that direct recoverability; ti() would give a more general 2-D interaction surface that is harder to read as "the FP effect by year" without additional slicing. mgcv chosen over gamm4 because s(stat_rec, bs = "re") already implements the random intercept inside gam(), fits cleanly here, and exposes gam.check() diagnostics. mgcv ships with base R — NO new ad-hoc dependency (unlike glmmTMB / spaMM / strucchange).
MUST be run with `Rscript --vanilla` (renv not activated) — same environment note as Round 1/2.

## Data
Reused Round 1 simple_nonspatial formula: residual ~ log_hours_total + phase + mean_ln_B_obs + (1 | stat_rec) + log_hours_total:phase
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015 (identical to Round 1 / Round 2 universe).
year_c centre (mean calendar year over analysis hauls): 1999.4839

## Model fitting
Fitting continuous_linear_year (glmmTMB, REML) ...
continuous_linear_year REML fit time: 0.82 sec.
Fitting smooth_year_gam (mgcv::gam, REML; expect ~1–2 min) ...
smooth_year_gam REML fit time: 87.77 sec.
Refitting all three models under ML for a fair AIC/BIC comparison ...
ML refits (all three) total time: 103.66 sec.
Saved model objects: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_model_objects.rds

## 1. Fit-quality comparison (AIC / BIC)
CAVEAT on information criteria: REML AIC/BIC are NOT strictly comparable across models with different fixed-effect structures (the REML criterion profiles out fixed effects, so the "likelihood" being compared is not the same quantity). The PRIMARY comparison below therefore uses the ML refits. A further caveat specific to the GAM: mgcv's AIC uses an effective-degrees-of-freedom penalty based on the smooth's EDF, which is a DIFFERENT penalty construction from glmmTMB's parametric AIC — so even the ML AIC comparison between the GAM and the two glmmTMB models is approximate, not an exact nested-model likelihood-ratio comparison. Reported with that caveat, not as a formal model-selection test.
ML AIC/BIC (primary comparison):
  - smooth_year_gam: AIC = 12750.3, BIC = 13809.4, logLik = -6229.2, n_params_reported = 145.968635573656
  - categorical_phase: AIC = 13133.7, BIC = 13213.5, logLik = -6555.9, n_params_reported = 11
  - continuous_linear_year: AIC = 13306.1, BIC = 13356.8, logLik = -6646.0, n_params_reported = 7
REML AIC/BIC (reported for completeness; NOT the primary comparison — see caveat above):
  - smooth_year_gam: AIC = 12749.8, BIC = 13810.8, logLik = -6228.7, n_params_reported = 146.226979877202
  - categorical_phase: AIC = 13199.9, BIC = 13279.7, logLik = -6588.9, n_params_reported = 11
  - continuous_linear_year: AIC = 13356.0, BIC = 13406.8, logLik = -6671.0, n_params_reported = 7
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_ic_comparison.csv

## 2. Does the fishing-pressure × time interaction still show up?
Continuous linear year — interaction log_hours_total:year_c: estimate = 0.002370, SE = 0.000424, z = 5.591, p = 2.254e-08. (Positive estimate means the fishing-pressure slope increases with year; NOT interpreted here as an H2/H3 finding.)
Smooth year GAM — by-smooth s(year):log_hours_total: edf = 7.267, Ref.df = 7.693, F = 44.050, approx p = 0.000e+00. (Highly significant approx p indicates a non-flat fishing-pressure effect across years; NOT interpreted here as an H2/H3 finding.)
Categorical phase — fishing-pressure × phase interaction terms (reused Round 1 REML fit):
  - log_hours_total:phase1989-2000: estimate = -0.0102, SE = 0.0146, z = -0.696, p = 4.867e-01
  - log_hours_total:phase2001-2007: estimate = 0.0872, SE = 0.0146, z = 5.957, p = 2.576e-09
  - log_hours_total:phase2008-2015: estimate = 0.0579, SE = 0.0145, z = 4.003, p = 6.248e-05
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_effects.csv

## 3. Visual comparison: fishing-pressure effect vs year (overlay)
Per-phase fishing-pressure slopes (categorical phase model):
  - 1985-1988 (1985–1988): slope = -0.0149, SE = 0.0139, 95% CI = [-0.0422, 0.0125]
  - 1989-2000 (1989–2000): slope = -0.0251, SE = 0.0089, 95% CI = [-0.0425, -0.0076]
  - 2001-2007 (2001–2007): slope = 0.0724, SE = 0.0078, 95% CI = [0.0571, 0.0876]
  - 2008-2015 (2008–2015): slope = 0.0430, SE = 0.0071, 95% CI = [0.0291, 0.0568]
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_fp_slope_by_year.csv
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_temporal_robustness_fp_effect_overlay.png
How to read the figure: agreement between the red step function and the green smooth curve (within CIs, and especially whether the step jumps land near where the smooth bends) would support the categorical phase structure as a reasonable simplification of a smoother underlying pattern; material disagreement (e.g. the smooth staying flat across a phase boundary the step treats as a level shift, or vice versa) would argue the discrete phases capture something a continuous term misses — or overstate a smooth change as a break. Verdict below.

## 4. Convergence diagnostics
  - categorical_phase (glmmTMB): converged = TRUE; hessian_positive_definite = TRUE; max|gradient| = 9.756e-04; optimizer message = 'relative convergence (4)'
  - continuous_linear_year (glmmTMB): converged = TRUE; hessian_positive_definite = TRUE; max|gradient| = 1.312e-02; optimizer message = 'relative convergence (4)'
  - smooth_year_gam (mgcv::gam): converged = TRUE; hessian_positive_definite = TRUE; max|gradient| = 2.685e-03; optimizer message = 'full convergence'
    gam.check() (basis-dimension / k check):  | Method: REML   Optimizer: outer newton | full convergence after 7 iterations. | Gradient range [-0.002684628,0.0006997747] | (score 6489.185 & scale 0.195247). | Hessian positive definite, eigenvalue range [1.711191,5230.305]. | Model rank =  175 / 175  |  | Basis dimension (k) checking results. Low p-value (k-index<1) may | indicate that k is too low, especially if edf is close to k'. |  |                             k'    edf k-index p-value | s(year)                   7.00   6.34    1.03    0.99 | s(year):log_hours_total   8.00   7.27    1.03    0.98 | s(stat_rec)             158.00 128.72      NA      NA
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_convergence.csv

## Feasibility / robustness verdict
Year-by-year correlation of the fitted fishing-pressure effect: phase vs linear = 0.742; phase vs GAM = 0.778; linear vs GAM = 0.555. Fraction of years with the same sign of the FP effect: phase vs linear = 83.9%; phase vs GAM = 77.4%.
Lowest-ML-AIC model: smooth_year_gam (see caveat above — not a formal selection test).
PLAIN STATEMENT (per the brief — categorical phase is a reasonable simplification / materially disagrees with the continuous alternatives):
  CATEGORICAL PHASE IS A REASONABLE SIMPLIFICATION. The fishing-pressure × time interaction shows up under all three time structures (phase interactions, linear interaction p = 2.25e-08, GAM by-smooth approx p = 0.00e+00), and the phase model's step-function FP effect tracks the continuous alternatives at the year-by-year level (phase–GAM correlation = 0.78; same-sign years = 77%). The GAM is more flexible and preferred by ML AIC, but it does not tell a qualitatively different story from the discrete phases — the step function is a coarse but recognisable summary of the smoother curve. RECOMMENDATION for supervisor discussion: keep the categorical phase structure as the primary (interpretable, break-aligned) specification; treat the linear and GAM fits as supporting robustness checks, not replacements.
This verdict is about TIME-STRUCTURE AGREEMENT only — it does not interpret, and is not based on treating, the sign or magnitude of any fishing-pressure coefficient as an H2/H3 finding.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_ic_comparison.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_fp_slope_by_year.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_convergence.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_temporal_robustness_fp_effect_overlay.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_temporal_robustness_run_log.md (this file)
