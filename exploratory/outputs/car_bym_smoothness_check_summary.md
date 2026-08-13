# CAR global-smoothness false-negative check — summary

## Model choice

**Option A analogue: BYM (adjacency + unstructured IID) via `spaMM::fitme`.**

- Preferred Leroux/BYM2 via `CARBayes` / R-INLA not installed in this stack.
- spaMM already used for primary CAR; adding `(1 | stat_rec)` beside
  `adjacency(1 | stat_rec)` gives a BYM-style mix without a new dependency.
- Mixing diagnostic: φ = λ_spatial / (λ_spatial + λ_iid)
  (φ → 1 fully spatial/global; φ → 0 fully local/IID).
- DAGAR not used (Option A path available).

## Estimation framework

**Frequentist REML** (`spaMM::fitme`), same framework as primary CAR.
Not Bayesian MCMC (CARBayes) or INLA. Coefficients below are REML point
estimates with Wald-style CIs — directly comparable to primary CAR on
estimation footing (unlike a CARBayes/INLA posterior).

- Formula: `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 |      stat_rec) + (1 | stat_rec)`
- Adjacency: loaded from h2h3_feasibility_round2_model_objects.rds$adjMatrix; confirmed identical to phase_v2_reporting_model_objects.rds$adjMatrix
- Fit time: 0.82 sec (CAR reference refit 0.35 sec)

## Mixing / dependence parameter

- λ_spatial (adjacency): **0.0124203**
- λ_iid (unstructured): **1.98608e-07**
- **φ = 0.999984** (fraction of RE variance that is spatial)
- CAR ρ = 0.131292 (admissible [-0.2634, 0.1317]; 99.9% of upper bound)

φ near 1 → data favour **global spatial** structure; little support for local IID intercepts beyond CAR.

## Phase-specific FP_between (BYM)

| Phase | H2 (`FP_between`) coefficient — BYM | 95% CI |
|---|---|---|
| 1985-1991 | -0.063 | [-0.094, -0.032] |
| 1992-2001 | -0.024 | [-0.053, +0.004] |
| 2002-2007 | +0.084 | [+0.053, +0.115] |
| 2008-2015 | +0.013 | [-0.016, +0.043] |

## Direct comparison vs primary CAR

Primary CAR confounding flags come from the existing spatial-permutation
bootstrap (`permutation_bootstrap_FP_between_CAR`). **BYM confounding was
not re-tested** (brief: no permutation layer). When FE match CAR, the same
confounded pattern is implied.

| Phase | CAR coef | CAR confounded? | BYM coef | Δ (BYM−CAR) | BYM status |
|---|---|---|---|---|---|
| 1985-1991 | -0.063 | Yes | -0.063 | +1.41e-07 | not re-tested; FE ≈ CAR ⇒ same pattern implied |
| 1992-2001 | -0.024 | Yes | -0.024 | +1.40e-07 | not re-tested; FE ≈ CAR ⇒ same pattern implied |
| 2002-2007 | +0.084 | Yes | +0.084 | +1.40e-07 | not re-tested; FE ≈ CAR ⇒ same pattern implied |
| 2008-2015 | +0.013 | No | +0.013 | +1.43e-07 | not re-tested; FE ≈ CAR ⇒ same pattern implied |

## Interpretation

Relaxing CAR's forced global smoothness via a BYM (adjacency + IID) term does not change the H2 conclusion. The estimated mixing proportion φ = 1.0000 sits at the fully-spatial end of the spectrum (IID variance numerically zero), and phase-specific FP_between coefficients match primary CAR to < 1.0e-04. The data do not support extra local unstructured intercepts beyond the CAR term, so the three-confounded / one-clean CAR pattern is not an artifact of forcing global smoothness. Caveat: this check is REML/spaMM BYM, not Bayesian Leroux/BYM2; confounding flags were not re-tested by permutation (per scope).

## Scope

- Single fit only; no permutation bootstrap.
- Not a candidate replacement for primary CAR.
- No further modelling without check-back.

## Outputs

- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_fp_slopes.csv`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_vs_car.csv`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_summary.md`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_run_log.md`
- `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/car_bym_smoothness_check_objects.rds`

