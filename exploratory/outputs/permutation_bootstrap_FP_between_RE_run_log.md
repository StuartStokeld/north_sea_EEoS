# FP_between spatial confounding bootstrap — run log (ARCHIVED RE)

> Relocated from `outputs/fp_between_confounding_bootstrap_*`. Paths below
> record the original run destination; current archive lives under
> `exploratory/outputs/permutation_bootstrap_FP_between_RE_*`.
> Current CAR version: `pipeline/permutation_bootstrap_FP_between_CAR.R`.

# FP_between spatial confounding bootstrap — run log

Rectangle-level permutation of FP_between under the current primary within-between glmmTMB model. Destroys the spatial arrangement of the H2 covariate while leaving residual, FP_within, the phase factor, and (1 | stat_rec) untouched.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_sessionInfo.txt
glmmTMB 1.1.14
N_BOOT = 1000; SEED = 42; N_CORES = 1

## Inputs
Loaded: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds
Phase scheme: phase_v2 (column `phase_v2`)
Phase levels: 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015
Formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)
FP_between unique per rectangle: 158 values; range [6.223, 10.881]
Permutation sanity check passed: only FP_between changes; FP_within / residual / stat_rec / phase_v2 unchanged.

## Step 1 — Observed H2 coefficients (baseline)
Observed fixef FP_between (reference phase 1985-1991): -0.024729
  SE = 0.015505; 95% CI [-0.055119, 0.005661]; z = -1.595; p = 0.1107
  Phase 1985-1991 slope = -0.024729 (SE 0.015505; 95% CI [-0.055119, 0.005661]; p = 0.1107)
  Phase 1992-2001 slope = +0.014173 (SE 0.014137; 95% CI [-0.013536, 0.041881]; p = 0.3161)
  Phase 2002-2007 slope = +0.122331 (SE 0.015465; 95% CI [0.092020, 0.152642]; p = 2.567e-15)
  Phase 2008-2015 slope = +0.051887 (SE 0.014636; 95% CI [0.023201, 0.080574]; p = 0.0003924)

## Steps 2–3 — Rectangle-level FP_between permutation + refit
Each replicate: sample(FP_between) across the 158 rectangles once, reassign to all hauls in that rectangle, refit primary model with REML = TRUE.
Running sequentially (replicate)
Finished: n_boot = 1000; n_failed = 0; n_ok = 1000; runtime = 335.3 sec (0.34 sec/rep)

### Matched permutation procedure (explicit)
For each of the 1000 replicates the procedure is:

1. **One global shuffle** of the 158 rectangle-level `FP_between` values (whole-series long-run means; not recomputed within phase).
2. **One refit** of the full primary formula
   `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)`.
3. **One `fixef()` extract** from that single fitted object, recording:
   - reference-phase main effect `FP_between` (1985–1991), and
   - the three interaction increments `FP_between:phase_v2…` (1992–2001, 2002–2007, 2008–2015),
   which are then combined to phase-specific slopes as `main + interaction`.

There are **not** four independent shuffle-and-refit procedures. The four null distributions (each of length 1000) are columns from the same 1000 × K coefficient matrix — one row per replicate, all coefficients matched to the same randomization. The four empirical p-values therefore come from a coherent, matched permutation procedure.

## Step 4 — Observed vs permutation null
FP_between (reference-phase main effect, 1985-1991): obs=-0.024729; null mean=+0.000165 sd=0.015810; 95% null [-0.029975, 0.031692]; p_emp=0.1199; inside_95=TRUE
  Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
  Interpretation: H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.
Phase 1985-1991 FP_between slope: obs=-0.024729; null mean=+0.000165 sd=0.015810; 95% null [-0.029975, 0.031692]; p_emp=0.1199; inside_95=TRUE
  Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
  Interpretation: H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.
Phase 1992-2001 FP_between slope: obs=+0.014173; null mean=+0.000284 sd=0.015232; 95% null [-0.030904, 0.030465]; p_emp=0.3487; inside_95=TRUE
  Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
  Interpretation: H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.
Phase 2002-2007 FP_between slope: obs=+0.122331; null mean=+0.000840 sd=0.019782; 95% null [-0.038040, 0.040760]; p_emp=0.000999 (< 0.001); inside_95=FALSE
  Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
  Interpretation: The rectangle intercept and FP_between are competing for the same spatial signal when FP_between's real spatial arrangement is present. The reported H2 effect is at least partly attributable to this overlap, not purely to the fishing-pressure effect itself.
Phase 2008-2015 FP_between slope: obs=+0.051887; null mean=+0.000146 sd=0.016808; 95% null [-0.031329, 0.033119]; p_emp=0.002997; inside_95=FALSE
  Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
  Interpretation: The rectangle intercept and FP_between are competing for the same spatial signal when FP_between's real spatial arrangement is present. The reported H2 effect is at least partly attributable to this overlap, not purely to the fishing-pressure effect itself.

## Deliverables
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_results.csv
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/fp_between_confounding_bootstrap_null.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_summary.md
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_objects.rds

## Headline
**Headline (phase-specific H2 slopes, phase_v2):** at least one presented phase-specific `FP_between` slope falls outside its spatial-permutation null 95% interval, indicating that the reported H2 effect depends on the real spatial arrangement of fishing pressure and is at least partly entangled with the rectangle intercept's spatial signal.
