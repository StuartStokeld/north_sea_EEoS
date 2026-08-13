# CAR global-smoothness false-negative check (BYM) — run log

Single fit only. Relaxes CAR's forced-global-smoothness assumption by adding an unstructured rectangle intercept alongside adjacency(1|stat_rec). Not a primary-model candidate. No permutation layer.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_sessionInfo.txt
spaMM 4.6.65

## Inputs
Adjacency: loaded from h2h3_feasibility_round2_model_objects.rds$adjMatrix; confirmed identical to phase_v2_reporting_model_objects.rds$adjMatrix
Data: 10464 hauls, 158 rectangles; formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 |      stat_rec) + (1 | stat_rec)
Estimation framework: frequentist REML via spaMM::fitme (same family as primary CAR — not Bayesian CARBayes/INLA).

## Runtime check
Primary-CAR-equivalent refit: 0.35 sec
Projected BYM (adjacency + IID) cost: previously timed ~0.7 sec on this panel — comparable to CAR; proceeding with single fit.

## Fit BYM-style mixing model
BYM fit time: 0.82 sec
lambda_spatial (adjacency term): 0.0124203
lambda_iid (unstructured term): 1.98608e-07
Mixing proportion φ = λ_spatial / (λ_spatial + λ_iid) = 0.999984
CAR rho = 0.131292 (admissible [-0.263426, 0.131705]; fraction of upper bound = 0.9990)
Interpretation of φ: near 1 → estimated structure is almost fully spatial (global-smoothness end of the spectrum); little support for local IID rectangle intercepts beyond the CAR term.

## Phase-specific FP_between slopes (BYM)
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_fp_slopes.csv
  1985-1991: -0.062797 (SE 0.015917; 95% CI [-0.093994, -0.031599]; p = 7.973e-05)
  1992-2001: -0.024352 (SE 0.014719; 95% CI [-0.053202, +0.004497]; p = 0.09804)
  2002-2007: +0.084112 (SE 0.015919; 95% CI [+0.052910, +0.115314]; p = 1.267e-07)
  2008-2015: +0.013256 (SE 0.015114; 95% CI [-0.016368, +0.042879]; p = 0.3805)
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_vs_car.csv

## CAR vs BYM comparison
  1985-1991: CAR -0.0628 (confounded=TRUE) | BYM -0.0628 | Δ=+1.41e-07
  1992-2001: CAR -0.0244 (confounded=TRUE) | BYM -0.0244 | Δ=+1.40e-07
  2002-2007: CAR +0.0841 (confounded=TRUE) | BYM +0.0841 | Δ=+1.40e-07
  2008-2015: CAR +0.0133 (confounded=FALSE) | BYM +0.0133 | Δ=+1.43e-07

## Interpretation
Relaxing CAR's forced global smoothness via a BYM (adjacency + IID) term does not change the H2 conclusion. The estimated mixing proportion φ = 1.0000 sits at the fully-spatial end of the spectrum (IID variance numerically zero), and phase-specific FP_between coefficients match primary CAR to < 1.0e-04. The data do not support extra local unstructured intercepts beyond the CAR term, so the three-confounded / one-clean CAR pattern is not an artifact of forcing global smoothness. Caveat: this check is REML/spaMM BYM, not Bayesian Leroux/BYM2; confounding flags were not re-tested by permutation (per scope).
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_summary.md
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_objects.rds
