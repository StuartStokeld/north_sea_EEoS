# FP_between spatial confounding bootstrap — Spec A — run log

Rectangle-level permutation of FP_between under the Spec A model. FP_between_lag recomputed from shuffled FP_between after each permutation; k-NN neighbour structure fixed. Targets: same FP_between × phase_v2 coefficients as the original bootstrap for comparison.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_spec_a_sessionInfo.txt
glmmTMB 1.1.14
N_BOOT = 1000; SEED = 42; N_CORES = 1

## Inputs
Loaded Spec A: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2_spec_a.rds
k-NN weights: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/knn_listw_k4.rds
Phase levels: 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015
Analysis data: 10464 hauls, 158 rectangles, years 1985–2015
Formula: residual ~ FP_between * phase_v2 + FP_between_lag * phase_v2 +      FP_within * phase_v2 + (1 | stat_rec)
Permutation sanity check passed: FP_between and FP_between_lag change; FP_within / residual / stat_rec / phase_v2 unchanged.

## Step 1 — Observed H2 coefficients (Spec A baseline)

## Steps 2–3 — Permutation + refit (lag recomputed each replicate)
Running sequentially
  progress 1/1000 (0.5 sec elapsed, ~0.53 sec/rep)
  progress 100/1000 (41.5 sec elapsed, ~0.42 sec/rep)
  progress 200/1000 (83.3 sec elapsed, ~0.42 sec/rep)
  progress 300/1000 (124.3 sec elapsed, ~0.41 sec/rep)
  progress 400/1000 (165.9 sec elapsed, ~0.41 sec/rep)
  progress 500/1000 (208.4 sec elapsed, ~0.42 sec/rep)
  progress 600/1000 (249.8 sec elapsed, ~0.42 sec/rep)
  progress 700/1000 (291.3 sec elapsed, ~0.42 sec/rep)
  progress 800/1000 (332.0 sec elapsed, ~0.41 sec/rep)
  progress 900/1000 (372.7 sec elapsed, ~0.41 sec/rep)
  progress 1000/1000 (414.0 sec elapsed, ~0.41 sec/rep)
Finished: n_boot = 1000; n_failed = 0; n_ok = 1000; runtime = 414.0 sec

## Step 4 — Observed vs permutation null (Spec A)
FP_between (reference-phase main effect, 1985-1991): obs=+0.009876; null 95% [-0.030784, 0.031979]; p_emp=0.5584; inside_95=TRUE
Phase 1985-1991 FP_between slope: obs=+0.009876; null 95% [-0.030784, 0.031979]; p_emp=0.5584; inside_95=TRUE
Phase 1992-2001 FP_between slope: obs=-0.001341; null 95% [-0.031023, 0.030852]; p_emp=0.9291; inside_95=TRUE
Phase 2002-2007 FP_between slope: obs=-0.021929; null 95% [-0.037591, 0.041316]; p_emp=0.2807; inside_95=TRUE
Phase 2008-2015 FP_between slope: obs=-0.016308; null 95% [-0.030981, 0.033009]; p_emp=0.3546; inside_95=TRUE

## Comparison vs original primary bootstrap
Phase 1985-1991: original inside_95=yes; Spec A inside_95=yes
Phase 1992-2001: original inside_95=yes; Spec A inside_95=yes
Phase 2002-2007: original inside_95=no; Spec A inside_95=yes
Phase 2008-2015: original inside_95=no; Spec A inside_95=yes

## Deliverables
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_spec_a_results.csv
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_spec_a_summary.md
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/fp_between_confounding_bootstrap_spec_a_null.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/fp_between_confounding_bootstrap_spec_a_objects.rds

**Headline (Spec A):** all phase-specific `FP_between` slopes fall inside their spatial-permutation null 95% intervals under the lag-augmented model.
