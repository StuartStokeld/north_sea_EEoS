# DATRAS haul-effort fix: old vs exactly 30 min (primary) vs 26–34 min (sensitivity)

Generated 2026-11-10 by `pipeline/compare_effort_fix.R` on branch `fix/haul-effort-datatype`.

- **old**: `outputs/archive_pre_effort_fix`. Pipeline before the fix (DataType C per-hour counts used as raw counts; all hauls, any duration).
- **30 min (primary)**: `outputs`. Effort fix: valid hauls (HaulVal V), DataType C/R/P/S, **HaulDur == 30 exactly**; C counts × HaulDur/60.
- **26–34 (sensitivity)**: `outputs/sensitivity_window_26_34`. Same fix with the window 26 ≤ HaulDur ≤ 34 (run first as the primary, then moved there unchanged).
- Model formulas, phases, priors and fishing-pressure construction are unchanged. `B_obs` (FishGlob raw haul weight) is unchanged.

## H1 — haul-level EEoS

### Haul sets

H1 hauls: old 12,069; 30 min 10,128; 26–34 11,027. The 30-min set nests in both others (0 of its hauls missing from old, 0 from 26–34).

| DataType | old | 30 min | 26–34 |
|---|---|---|---|
| C | 10220 | 8493 | 9225 |
| P |   102 |   92 |   98 |
| R |  1738 | 1543 | 1704 |
| (no unique HH) |     9 |    0 |    0 |

By recorded HaulDur (minutes):

| HaulDur | old | 30 min | 26–34 |
|---|---|---|---|
| (no unique HH) |     9 |     0 |     0 |
| 26–34, not 30 |   899 |     0 |   899 |
| 30 | 10128 | 10128 | 10128 |
| other |   336 |     0 |     0 |
| ≥ 55 |   697 |     0 |     0 |

Definitions (log scale; obs = ln B_obs, pred = ln B_pred):

- **1:1 log_r2** = 1 − Σ(obs − pred)² / Σ(obs − mean obs)². Penalises pattern *and* offset; can be negative.
- **Pearson r**: pattern only (invariant to offset and scale).
- **OLS slope** of obs on pred: 1 means predicted spread matches observed spread; < 1 means predictions too spread out.
- **median ratio** = median(B_pred / B_obs); **mean log bias** = mean(pred − obs) (> 0 = over-prediction; exp(bias) = geometric-mean ratio).
- **log_r2 bias-removed** = 1:1 log_r2 after subtracting the mean log bias from pred: what 1:1 log_r2 would be with the offset gone.

### Benchmarks, each run on its own haul set (as each run would be reported)

| model | run | n | 1:1 log_r2 | Pearson r | OLS slope | median ratio | mean log bias | log_r2 bias-removed |
|---|---|---|---|---|---|---|---|---|
| EEoS (S, N, E), unfitted | old | 12,069 | -0.223 | 0.962 | 0.781 | 4.081 | 1.322 | 0.853 |
| EEoS (S, N, E), unfitted | 30 min | 10,128 | 0.661 | 0.986 | 0.826 | 1.958 | 0.652 | 0.928 |
| EEoS (S, N, E), unfitted | 26–34 | 11,027 | 0.658 | 0.985 | 0.827 | 1.957 | 0.652 | 0.928 |
| Productivity 1:1, log(E × m_min), headline baseline | old | 12,069 | 0.736 | 0.948 | 0.870 | 0.667 | -0.480 | 0.878 |
| Productivity 1:1, log(E × m_min), headline baseline | 30 min | 10,128 | 0.256 | 0.968 | 0.911 | 0.369 | -1.033 | 0.927 |
| Productivity 1:1, log(E × m_min), headline baseline | 26–34 | 11,027 | 0.245 | 0.967 | 0.911 | 0.368 | -1.036 | 0.927 |
| Productivity 1:1 uncalibrated, log E_raw | old | 12,069 | 0.619 | 0.955 | 0.927 | 0.536 | -0.683 | 0.906 |
| Productivity 1:1 uncalibrated, log E_raw | 30 min | 10,128 | -0.002 | 0.975 | 0.970 | 0.293 | -1.231 | 0.950 |
| Productivity 1:1 uncalibrated, log E_raw | 26–34 | 11,027 | -0.013 | 0.974 | 0.971 | 0.292 | -1.231 | 0.949 |
| ln(E) OLS, fitted | old | 12,069 | 0.595 | 0.771 | 1.000 | 0.962 | 0.000 | 0.595 |
| ln(E) OLS, fitted | 30 min | 10,128 | 0.603 | 0.777 | 1.000 | 0.960 | 0.000 | 0.603 |
| ln(E) OLS, fitted | 26–34 | 11,027 | 0.595 | 0.772 | 1.000 | 0.961 | -0.000 | 0.595 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | old | 12,069 | -0.434 | 0.911 | 1.048 | 0.348 | -0.992 | 0.828 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | 30 min | 10,128 | 0.593 | 0.957 | 1.061 | 0.604 | -0.489 | 0.914 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | 26–34 | 11,027 | 0.595 | 0.958 | 1.062 | 0.604 | -0.489 | 0.915 |

### Benchmarks on the same 10,128 hauls (the 30-min set) in every run

Separates the effect of the count conversion from the change of haul set. On these hauls the 30-min and 26–34
inputs are identical (max relative difference in S, N, E, B_pred: 0.0e+00), so the 26–34 rows differ from the 30-min rows
only for ln(E) OLS, which is fitted to each run's own full set.

| model | run | n | 1:1 log_r2 | Pearson r | OLS slope | median ratio | mean log bias | log_r2 bias-removed |
|---|---|---|---|---|---|---|---|---|
| EEoS (S, N, E), unfitted | old | 10,128 | -0.282 | 0.970 | 0.776 | 4.199 | 1.350 | 0.863 |
| EEoS (S, N, E), unfitted | 30 min | 10,128 | 0.661 | 0.986 | 0.826 | 1.958 | 0.652 | 0.928 |
| EEoS (S, N, E), unfitted | 26–34 | 10,128 | 0.661 | 0.986 | 0.826 | 1.958 | 0.652 | 0.928 |
| Productivity 1:1, log(E × m_min), headline baseline | old | 10,128 | 0.757 | 0.953 | 0.863 | 0.688 | -0.452 | 0.886 |
| Productivity 1:1, log(E × m_min), headline baseline | 30 min | 10,128 | 0.256 | 0.968 | 0.911 | 0.369 | -1.033 | 0.927 |
| Productivity 1:1, log(E × m_min), headline baseline | 26–34 | 10,128 | 0.256 | 0.968 | 0.911 | 0.369 | -1.033 | 0.927 |
| Productivity 1:1 uncalibrated, log E_raw | old | 10,128 | 0.653 | 0.962 | 0.921 | 0.546 | -0.650 | 0.919 |
| Productivity 1:1 uncalibrated, log E_raw | 30 min | 10,128 | -0.002 | 0.975 | 0.970 | 0.293 | -1.231 | 0.950 |
| Productivity 1:1 uncalibrated, log E_raw | 26–34 | 10,128 | -0.002 | 0.975 | 0.970 | 0.293 | -1.231 | 0.950 |
| ln(E) OLS, fitted | old | 10,128 | 0.611 | 0.783 | 0.996 | 1.012 | 0.045 | 0.613 |
| ln(E) OLS, fitted | 30 min | 10,128 | 0.603 | 0.777 | 1.000 | 0.960 | 0.000 | 0.603 |
| ln(E) OLS, fitted | 26–34 | 10,128 | 0.603 | 0.777 | 1.005 | 0.953 | -0.008 | 0.603 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | old | 10,128 | -0.497 | 0.917 | 1.039 | 0.341 | -1.013 | 0.840 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | 30 min | 10,128 | 0.593 | 0.957 | 1.061 | 0.604 | -0.489 | 0.914 |
| Productivity ratio E/B^(3/4) (Harte Fig 2; different target) | 26–34 | 10,128 | 0.593 | 0.957 | 1.061 | 0.604 | -0.489 | 0.914 |

Notes: the productivity ratio (Harte Fig 2) compares E/B^(3/4) predicted vs observed, a different target from ln B.
ln(E) OLS is fitted (to each run's full H1 set), so its mean bias is ~0 on its own set by construction.
Full-set log_r2 values for EEoS, productivity 1:1 (calibrated and uncalibrated) and ln(E) OLS were checked to reproduce
each run's `h1_model_comparison.csv`.

### EEoS residual summary

| hauls | run | n_hauls | median B_pred/B_obs | log_r2 | mean residual (= −bias) | residual SD |
|---|---|---|---|---|---|---|
| own set | old | 12,069 | 4.081 | -0.223 | -1.322 | 0.489 |
| own set | 30 min | 10,128 | 1.958 | 0.661 | -0.652 | 0.338 |
| own set | 26–34 | 11,027 | 1.957 | 0.658 | -0.652 | 0.337 |
| same (30-min) hauls | old | 10,128 | 4.199 | -0.282 | -1.350 | 0.467 |
| same (30-min) hauls | 30 min | 10,128 | 1.958 | 0.661 | -0.652 | 0.338 |
| same (30-min) hauls | 26–34 | 10,128 | 1.958 | 0.661 | -0.652 | 0.338 |

### Same hauls, old → 30 min: what changed in the inputs (by DataType)

| DataType | n_hauls | median N new/old | median E new/old | S unchanged | median B_pred new/old | mean resid old | mean resid 30 min |
|---|---|---|---|---|---|---|---|
| C | 8493 | 0.500 | 0.500 | 100.0% | 0.435 | -1.488 | -0.655 |
| P |   92 | 1.000 | 1.000 | 100.0% | 1.000 | -0.575 | -0.575 |
| R | 1543 | 1.000 | 1.000 | 100.0% | 1.000 | -0.640 | -0.640 |

### Uniform-draw null model (`null_distribution_summary.csv`, each run's own H1 set)

| run | n | log_median | log_iqr | log_range | iqr_range_ratio | skewness | excess_kurtosis | hist_peak_ratio | log_q025 | log_q975 |
|---|---|---|---|---|---|---|---|---|---|---|
| old | 12069 | 11.897 | 1.609 | 13.456 | 0.120 | -0.193 | 1.016 | 4.957 | 9.314 | 14.378 |
| 30 min | 10128 | 11.850 | 1.588 | 13.456 | 0.118 | -0.232 | 1.122 | 5.006 | 9.293 | 14.269 |
| 26–34 | 11027 | 11.844 | 1.586 | 13.456 | 0.118 | -0.214 | 1.082 | 5.063 | 9.306 | 14.244 |

## Residual means by country × period (DataType C → R switches)

Before the fix, DataType C hauls had ~2× N and E, so B_pred was too high and residuals too low while a
country reported C. A C → R switch then shows as an upward step in that country's mean residual.
`step` = mean residual in the 4 years from the switch minus the 4 years before (GB-SCT: 2014–2015 after);
`net` subtracts the same-years change in all other countries (removes North-Sea-wide year effects).
Each run on its own H1 haul set.

Switch years confirmed from HH (share of valid hauls with DataType C, 4 years before → after):

| Country | switch_year | share C before | share C after |
|---|---|---|---|
| DK | 2000 | 1.00 | 0.00 |
| NL | 2000 | 0.75 | 0.00 |
| GB-SCT | 2014 | 1.00 | 0.00 |

| Country | switch year | step old | step 30 min | step 26–34 | net old | net 30 min | net 26–34 | n before/after old | n before/after 30 min | n before/after 26–34 |
|---|---|---|---|---|---|---|---|---|---|---|
| DK | 2000 | 0.822 | -0.044 | -0.010 | 0.880 | 0.049 | 0.085 | 176 / 154 | 158 / 121 | 175 / 154 |
| NL | 2000 | 0.547 | -0.073 | -0.071 | 0.620 | 0.018 | 0.015 | 202 / 257 | 202 / 251 | 202 / 252 |
| GB-SCT | 2014 | 0.902 | 0.062 | 0.061 | 1.021 | 0.113 | 0.111 | 214 / 102 | 211 / 92 | 213 / 98 |

### Year by year around each switch (H1 hauls)

| Country | year | share C | n old | n 30 min | n 26–34 | resid old | resid 30 min | resid 26–34 |
|---|---|---|---|---|---|---|---|---|
| DK | 1996 | 1.00 | 46 | 46 | 46 | -1.519 | -0.686 | -0.686 |
| DK | 1997 | 1.00 | 47 | 47 | 47 | -1.577 | -0.743 | -0.743 |
| DK | 1998 | 1.00 | 50 | 32 | 49 | -1.606 | -0.754 | -0.780 |
| DK | 1999 | 1.00 | 33 | 33 | 33 | -1.557 | -0.723 | -0.723 |
| DK | 2000 | 0.00 | 34 | 20 | 34 | -0.620 | -0.668 | -0.620 |
| DK | 2001 | 0.00 | 41 | 35 | 41 | -0.734 | -0.733 | -0.734 |
| DK | 2002 | 0.00 | 39 | 29 | 39 | -0.808 | -0.837 | -0.808 |
| DK | 2003 | 0.00 | 40 | 37 | 40 | -0.799 | -0.805 | -0.799 |
| GB-SCT | 2010 | 1.00 | 52 | 49 | 51 | -1.487 | -0.653 | -0.652 |
| GB-SCT | 2011 | 1.00 | 50 | 50 | 50 | -1.516 | -0.685 | -0.685 |
| GB-SCT | 2012 | 1.00 | 56 | 56 | 56 | -1.502 | -0.668 | -0.668 |
| GB-SCT | 2013 | 1.00 | 56 | 56 | 56 | -1.383 | -0.550 | -0.550 |
| GB-SCT | 2014 | 0.00 | 49 | 42 | 46 | -0.573 | -0.578 | -0.581 |
| GB-SCT | 2015 | 0.00 | 53 | 50 | 52 | -0.566 | -0.572 | -0.573 |
| NL | 1996 | 1.00 | 45 | 45 | 45 | -1.351 | -0.528 | -0.528 |
| NL | 1997 | 0.00 | 50 | 50 | 50 | -0.305 | -0.305 | -0.305 |
| NL | 1998 | 1.00 | 54 | 54 | 54 | -1.526 | -0.692 | -0.692 |
| NL | 1999 | 1.00 | 53 | 53 | 53 | -1.340 | -0.510 | -0.510 |
| NL | 2000 | 0.00 | 60 | 55 | 55 | -0.652 | -0.634 | -0.634 |
| NL | 2001 | 0.00 | 71 | 71 | 71 | -0.692 | -0.692 | -0.692 |
| NL | 2002 | 0.00 | 58 | 58 | 58 | -0.321 | -0.321 | -0.321 |
| NL | 2003 | 0.00 | 68 | 67 | 68 | -0.652 | -0.662 | -0.652 |

### All countries × `phase_v2` period

| Country | period | share C (old hauls) | n old | n 30 min | n 26–34 | resid old | resid 30 min | resid 26–34 |
|---|---|---|---|---|---|---|---|---|
| DE | 1985-1991 | 1.00 | 641 | 636 | 637 | -1.510 | -0.673 | -0.673 |
| DE | 1992-2001 | 1.00 | 751 | 746 | 747 | -1.592 | -0.755 | -0.756 |
| DE | 2002-2007 | 1.00 | 418 | 354 | 410 | -1.541 | -0.719 | -0.712 |
| DE | 2008-2015 | 1.00 | 523 | 405 | 486 | -1.514 | -0.639 | -0.644 |
| DK | 1985-1991 | 1.00 | 221 | 198 | 213 | -1.533 | -0.706 | -0.702 |
| DK | 1992-2001 | 0.82 | 427 | 317 | 420 | -1.459 | -0.763 | -0.772 |
| DK | 2002-2007 | 0.00 | 250 | 227 | 249 | -0.706 | -0.705 | -0.705 |
| DK | 2008-2015 | 0.00 | 326 | 252 | 314 | -0.709 | -0.717 | -0.711 |
| FR | 1985-1991 | 1.00 | 258 | 250 | 253 | -1.361 | -0.532 | -0.530 |
| FR | 1992-2001 | 1.00 | 587 | 491 | 566 | -1.438 | -0.602 | -0.606 |
| FR | 2002-2007 | 1.00 | 453 | 333 | 439 | -1.480 | -0.652 | -0.648 |
| FR | 2008-2015 | 1.00 | 630 | 481 | 615 | -1.453 | -0.625 | -0.624 |
| GB | 1985-1991 | 1.00 | 336 | 325 | 335 | -1.500 | -0.665 | -0.669 |
| GB | 1992-2001 | 1.00 |  53 |  52 |  52 | -1.565 | -0.733 | -0.733 |
| GB | 2002-2007 | 1.00 |  57 |  54 |  56 | -1.524 | -0.681 | -0.684 |
| GB-SCT | 1985-1991 | 1.00 | 371 |   4 |   4 | -0.772 | -0.771 | -0.771 |
| GB-SCT | 1992-2001 | 1.00 | 503 | 182 | 184 | -1.089 | -0.750 | -0.751 |
| GB-SCT | 2002-2007 | 1.00 | 312 | 300 | 306 | -1.520 | -0.684 | -0.684 |
| GB-SCT | 2008-2015 | 0.76 | 425 | 409 | 418 | -1.257 | -0.627 | -0.626 |
| NL | 1985-1991 | 1.00 | 512 | 506 | 506 | -1.414 | -0.589 | -0.589 |
| NL | 1992-2001 | 0.65 | 524 | 514 | 514 | -1.143 | -0.595 | -0.595 |
| NL | 2002-2007 | 0.00 | 383 | 374 | 379 | -0.593 | -0.598 | -0.597 |
| NL | 2008-2015 | 0.00 | 480 | 425 | 470 | -0.630 | -0.630 | -0.628 |
| NO | 1985-1991 | 1.00 | 302 | 270 | 283 | -1.233 | -0.392 | -0.392 |
| NO | 1992-2001 | 1.00 | 413 | 383 | 402 | -1.473 | -0.640 | -0.640 |
| NO | 2002-2007 | 1.00 | 199 | 109 | 145 | -1.559 | -0.629 | -0.604 |
| NO | 2008-2015 | 0.86 | 306 | 212 | 268 | -1.260 | -0.503 | -0.499 |
| SE | 1985-1991 | 1.00 | 296 | 268 | 269 | -1.475 | -0.690 | -0.688 |
| SE | 1992-2001 | 1.00 | 463 | 437 | 455 | -1.625 | -0.778 | -0.782 |
| SE | 2002-2007 | 1.00 | 285 | 275 | 279 | -1.537 | -0.684 | -0.686 |
| SE | 2008-2015 | 1.00 | 364 | 339 | 353 | -1.466 | -0.609 | -0.605 |

## H2 / H3 — primary within-between model (`phase_v2`)

`residual ~ FP_between * phase_v2 + FP_within * phase_v2`; CAR `adjacency(1|stat_rec)` (`wb_car_v2`, primary for H2 and H3)
and companion RE `(1|stat_rec)` (`wb_primary_v2`). Each run's own haul panel.

| run | n_hauls | n_rectangles | n_1985-1991 | n_1992-2001 | n_2002-2007 | n_2008-2015 |
|---|---|---|---|---|---|---|
| old | 10464 | 158 | 2643 | 3258 | 2062 | 2501 |
| 30 min | 8670 | 158 | 2192 | 2685 | 1745 | 2048 |
| 26–34 | 9482 | 158 | 2233 | 2885 | 1975 | 2389 |

Same rectangle set in all three runs: yes (158 rectangles).

| run | CAR R²m | CAR R²c | RE R²m | RE R²c |
|---|---|---|---|---|
| old | 0.043 | 0.157 | 0.051 | 0.171 |
| 30 min | 0.016 | 0.209 | 0.019 | 0.225 |
| 26–34 | 0.016 | 0.210 | 0.018 | 0.223 |

Nakagawa R² as in `nakagawa_r2_spamm()` / `nakagawa_r2_glmmtmb()` (`pipeline/R/h2h3_results_helpers.R`).

### Phase-specific slopes — CAR `wb_car_v2` (presented H2 and H3 effects)

Each cell: slope (p; Bonferroni p). Bonferroni p = min(1, 4 × p) within each m = 4 family (H2, H3 separately);
† = Bonferroni p < 0.05. SEs and 95% CIs are in `effort_fix_comparison_tables/h2h3_slopes_by_phase.csv`.

| hypothesis | phase | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|---|
| H2 FP_between | 1985-1991 | -0.0628 (8.0e-05; 3.2e-04) † | 0.0042 (0.758; 1.000) | 0.0057 (0.669; 1.000) |
| H2 FP_between | 1992-2001 | -0.0244 (0.098; 0.392) | 0.0257 (0.044; 0.175) | 0.0244 (0.052; 0.208) |
| H2 FP_between | 2002-2007 | 0.0841 (1.3e-07; 5.1e-07) † | 0.0129 (0.331; 1.000) | 0.0121 (0.347; 1.000) |
| H2 FP_between | 2008-2015 | 0.0133 (0.380; 1.000) | -0.0255 (0.048; 0.191) | -0.0250 (0.046; 0.185) |
| H3 FP_within | 1985-1991 | 0.0273 (0.065; 0.261) | -0.0114 (0.283; 1.000) | -0.0094 (0.361; 1.000) |
| H3 FP_within | 1992-2001 | 0.0462 (0.007; 0.028) † | 0.0343 (0.005; 0.021) † | 0.0270 (0.024; 0.097) |
| H3 FP_within | 2002-2007 | -0.0043 (0.806; 1.000) | -0.0316 (0.011; 0.043) † | -0.0327 (0.006; 0.023) † |
| H3 FP_within | 2008-2015 | 0.0567 (3.2e-06; 1.3e-05) † | 0.0321 (3.1e-04; 0.001) † | 0.0292 (5.3e-04; 0.002) † |

Sign and significance, old / 30 min / 26–34:

| hypothesis | phase | sign | p < 0.05 | Bonferroni p < 0.05 |
|---|---|---|---|---|
| H2 | 1985-1991 | − / + / + | yes / no / no | yes / no / no |
| H2 | 1992-2001 | − / + / + | no / yes / no | no / no / no |
| H2 | 2002-2007 | + / + / + | yes / no / no | yes / no / no |
| H2 | 2008-2015 | + / − / − | no / yes / yes | no / no / no |
| H3 | 1985-1991 | + / − / − | no / no / no | no / no / no |
| H3 | 1992-2001 | + / + / + | yes / yes / yes | yes / yes / no |
| H3 | 2002-2007 | − / − / − | no / yes / yes | no / yes / yes |
| H3 | 2008-2015 | + / + / + | yes / yes / yes | yes / yes / yes |

### Fixed effects — CAR `wb_car_v2`

Each cell: estimate (SE; p).

| term | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|
| FP_between | -0.0628 (0.0159; p 8.0e-05) | 0.0042 (0.0136; p 0.758) | 0.0057 (0.0134; p 0.669) |
| phase_v21992-2001 | -0.4084 (0.1134; p 3.2e-04) | -0.2969 (0.0835; p 3.8e-04) | -0.2719 (0.0818; p 8.9e-04) |
| phase_v22002-2007 | -1.2397 (0.1263; p 9.9e-23) | -0.1501 (0.0900; p 0.095) | -0.1279 (0.0864; p 0.139) |
| phase_v22008-2015 | -0.4803 (0.1201; p 6.4e-05) | 0.2494 (0.0868; p 0.004) | 0.2543 (0.0834; p 0.002) |
| FP_within | 0.0273 (0.0148; p 0.065) | -0.0114 (0.0106; p 0.283) | -0.0094 (0.0103; p 0.361) |
| FP_between:phase_v21992-2001 | 0.0384 (0.0121; p 0.002) | 0.0215 (0.0089; p 0.015) | 0.0186 (0.0087; p 0.032) |
| FP_between:phase_v22002-2007 | 0.1469 (0.0136; p 4.4e-27) | 0.0087 (0.0097; p 0.370) | 0.0064 (0.0093; p 0.492) |
| FP_between:phase_v22008-2015 | 0.0761 (0.0128; p 2.8e-09) | -0.0297 (0.0092; p 0.001) | -0.0308 (0.0089; p 5.3e-04) |
| phase_v21992-2001:FP_within | 0.0189 (0.0225; p 0.400) | 0.0457 (0.0161; p 0.004) | 0.0364 (0.0156; p 0.020) |
| phase_v22002-2007:FP_within | -0.0315 (0.0237; p 0.183) | -0.0202 (0.0169; p 0.231) | -0.0233 (0.0163; p 0.151) |
| phase_v22008-2015:FP_within | 0.0294 (0.0201; p 0.142) | 0.0435 (0.0145; p 0.003) | 0.0386 (0.0139; p 0.006) |

### Phase-specific slopes — companion RE `wb_primary_v2`

| hypothesis | phase | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|---|
| H2 FP_between | 1985-1991 | -0.0247 (0.111; 0.443) | 0.0245 (0.052; 0.208) | 0.0247 (0.047; 0.189) |
| H2 FP_between | 1992-2001 | 0.0142 (0.316; 1.000) | 0.0467 (5.8e-05; 2.3e-04) † | 0.0439 (1.2e-04; 4.8e-04) † |
| H2 FP_between | 2002-2007 | 0.1223 (2.6e-15; 1.0e-14) † | 0.0339 (0.005; 0.021) † | 0.0317 (0.007; 0.029) † |
| H2 FP_between | 2008-2015 | 0.0519 (3.9e-04; 0.002) † | -0.0042 (0.723; 1.000) | -0.0052 (0.652; 1.000) |
| H3 FP_within | 1985-1991 | 0.0269 (0.070; 0.281) | -0.0126 (0.235; 0.938) | -0.0106 (0.306; 1.000) |
| H3 FP_within | 1992-2001 | 0.0437 (0.011; 0.044) † | 0.0316 (0.010; 0.042) † | 0.0244 (0.042; 0.167) |
| H3 FP_within | 2002-2007 | -0.0018 (0.919; 1.000) | -0.0297 (0.017; 0.067) | -0.0310 (0.009; 0.036) † |
| H3 FP_within | 2008-2015 | 0.0577 (2.4e-06; 9.6e-06) † | 0.0333 (1.9e-04; 7.5e-04) † | 0.0304 (3.3e-04; 0.001) † |

### Fixed effects — companion RE `wb_primary_v2`

| term | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|
| FP_between | -0.0247 (0.0155; p 0.111) | 0.0245 (0.0126; p 0.052) | 0.0247 (0.0125; p 0.047) |
| phase_v21992-2001 | -0.4108 (0.1136; p 3.0e-04) | -0.3011 (0.0836; p 3.1e-04) | -0.2756 (0.0819; p 7.6e-04) |
| phase_v22002-2007 | -1.2395 (0.1265; p 1.1e-22) | -0.1558 (0.0900; p 0.084) | -0.1335 (0.0865; p 0.123) |
| phase_v22008-2015 | -0.4836 (0.1204; p 5.9e-05) | 0.2420 (0.0869; p 0.005) | 0.2470 (0.0835; p 0.003) |
| FP_within | 0.0269 (0.0148; p 0.070) | -0.0126 (0.0106; p 0.235) | -0.0106 (0.0103; p 0.306) |
| FP_between:phase_v21992-2001 | 0.0389 (0.0122; p 0.001) | 0.0221 (0.0089; p 0.013) | 0.0191 (0.0087; p 0.028) |
| FP_between:phase_v22002-2007 | 0.1471 (0.0136; p 4.5e-27) | 0.0093 (0.0097; p 0.336) | 0.0070 (0.0093; p 0.450) |
| FP_between:phase_v22008-2015 | 0.0766 (0.0128; p 2.4e-09) | -0.0287 (0.0093; p 0.002) | -0.0299 (0.0089; p 7.8e-04) |
| phase_v21992-2001:FP_within | 0.0168 (0.0226; p 0.455) | 0.0442 (0.0161; p 0.006) | 0.0350 (0.0156; p 0.025) |
| phase_v22002-2007:FP_within | -0.0287 (0.0238; p 0.229) | -0.0171 (0.0169; p 0.312) | -0.0205 (0.0163; p 0.209) |
| phase_v22008-2015:FP_within | 0.0308 (0.0202; p 0.127) | 0.0459 (0.0145; p 0.002) | 0.0409 (0.0140; p 0.003) |

### Robustness — FP_between spatial permutation (CAR, 1,000 shuffles)

Each cell: two-sided p_perm, recomputed from the null draws.

| phase | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|
| 1985-1991 | 1.0e-03 | 0.673 | 0.618 |
| 1992-2001 | 0.033 | 0.008 | 0.011 |
| 2002-2007 | 1.0e-03 | 0.115 | 0.141 |
| 2008-2015 | 0.284 | 0.011 | 0.010 |

### Robustness — CAR rectangle subsampling (80% retention, 1,000 draws)

Each cell: % of draws matching the full-data sign and significance (% matching sign only).

| hypothesis | phase | old | 30 min (primary) | 26–34 (sensitivity) |
|---|---|---|---|---|
| H2 | 1985-1991 | 93.2 (100.0) | 87.7 (91.3) | 88.2 (93.4) |
| H2 | 1992-2001 | 88.3 (91.2) | 76.7 (100.0) | 25.7 (100.0) |
| H2 | 2002-2007 | 100.0 (100.0) | 77.3 (98.6) | 77.8 (98.7) |
| H2 | 2008-2015 | 75.2 (98.3) | 10.8 (96.5) | 10.2 (96.6) |
| H3 | 1985-1991 | 73.6 (99.1) | 88.2 (95.5) | 87.4 (92.3) |
| H3 | 1992-2001 | 79.3 (100.0) | 83.5 (100.0) | 53.6 (99.8) |
| H3 | 2002-2007 | 60.9 (60.9) | 76.0 (100.0) | 85.8 (100.0) |
| H3 | 2008-2015 | 100.0 (100.0) | 99.7 (100.0) | 99.2 (100.0) |

### Robustness — Spec A (CAR + k-NN spatial lag) identifiability gate

| run | gate | rho_primary | max_se_ratio_vs_fp_between | n_warnings |
|---|---|---|---|---|
| old | CLEAR | 0.131 | 1.487 | 0.000 |
| 30 min | CLEAR | 0.131 | 1.550 | 0.000 |
| 26–34 | CLEAR | 0.131 | 1.552 | 0.000 |

## Run notes

Exactly-30-min live pipeline run (CLAUDE.md order): 19 scripts, total 16.7 min; exit codes: 0.

| script | exit | seconds |
|---|---|---|
| pipeline/build_datras_state_variables.R | 0 |   7 |
| pipeline/build_eeos_predictions.R | 0 |  38 |
| pipeline/run_h1_harte_baseline.R | 0 |   1 |
| pipeline/run_h1_lne_reference.R | 0 |   1 |
| pipeline/run_h1_null_model.R | 0 |   1 |
| pipeline/run_pipeline_diagnostics.R | 0 |   4 |
| pipeline/run_h1_dropout_diagnosis.R | 0 |   4 |
| pipeline/import_couce_fishing_effort.R | 0 |   0 |
| pipeline/build_h2_rectangle_panel.R | 0 |   1 |
| pipeline/run_h2h3_within_between.R | 0 |   2 |
| pipeline/run_h2h3_phase_v2_refit.R | 0 |   1 |
| pipeline/run_h2h3_phase_v2_reporting.R | 0 |   4 |
| pipeline/run_h2h3_h2_multiplicity_subsampling.R | 0 |   9 |
| pipeline/permutation_bootstrap_FP_between_CAR.R | 0 | 118 |
| pipeline/run_h2h3_rectangle_subsampling_refit.R | 0 | 504 |
| pipeline/run_h2h3_rectangle_subsampling_car_refit.R | 0 | 268 |
| pipeline/build_knn_spatial_weights.R | 0 |   0 |
| pipeline/run_h2h3_spec_a_car_identifiability.R | 0 |   2 |
| pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R | 0 |  38 |

- `run_h2h3_within_between.R` did not add the project `.R_libs` to its library path (unlike the other H2/H3 scripts), so
  spaMM was not found and its original-phase CAR sensitivity was skipped in the pipeline runs above. The path was added and
  the script re-run on its own at exactly 30 min, regenerating `outputs/h2h3_wb_car_fixed_effects.csv`; its log is otherwise
  identical to the pipeline run (same 8,670 hauls; same RE fit read by `run_h2h3_phase_v2_refit.R`), so nothing downstream
  changes. The primary `wb_car_v2` fit (in `run_h2h3_phase_v2_reporting.R`) was never affected. The 26–34 sensitivity
  folder does not contain this original-phase CAR table.

**Package versions differ from `renv.lock`.** `renv.lock` records R 4.4.3; old run: R version 4.6.1; R version 4.4.3; 30 min: R version 4.6.1; 26–34: R version 4.6.1.
`renv::restore()` fails on R 4.6 (pinned versions such as rlang 1.1.6 do not compile), so the installed
current versions were used and `renv.lock` was not modified. Versions below are read from the `*_sessionInfo.txt`
files the H2/H3 scripts write (H1 scripts do not write one); several values mean different scripts saw different libraries.

| package | renv_lock | old | 30 min | 26–34 |
|---|---|---|---|---|
| spaMM | not in lock | 4.6.1, 4.6.65 | 4.6.65 | 4.6.65 |
| glmmTMB | not in lock | 1.1.14 | 1.1.14 | 1.1.14 |
| dplyr | 1.1.4 | 1.1.4, 1.2.1 | 1.2.1 | 1.2.1 |
| rlang | 1.1.6 | 1.1.6, 1.3.0 | 1.3.0 | 1.3.0 |
| Matrix | 1.7-2 | 1.7-2, 1.7-5 | 1.7-5 | 1.7-5 |
| TMB | not in lock | 1.9.20, 1.9.23 | 1.9.23 | 1.9.23 |

