# H1 haul-level dominance and size-homogeneity exploration

Generated: 2026-07-18 — **parallel diagnostic track; does not modify the primary H1 pipeline (log_r2/cor2, null model, dropout diagnostics, catchability scaling).**

## Question

Does EEoS prediction failure vary systematically with how "typical" a haul's
community structure is (numerically dominated by one species; that species
caught at a near-uniform size), independent of / conditional on the already-
reported biomass-magnitude bias (median `B_pred/B_obs` rising from ~3x in the
lowest `B_obs` quartile to >5x in the highest)?

N hauls in this diagnostic: **12069** (of 12069 EEoS predictions).

## Metric definitions

- **D** (Berger-Parker dominance) = `n_max_species / N_haul`, using the same
  SubFactor-raised DATRAS HL abundance already used for `N` in the EEoS call
  (LW-matched bins only — a species with no FishBase LW parameters is not
  N-eligible under the production pipeline's own convention, so it cannot be
  numerically dominant here either).
- **size_CV** = raised-count-weighted coefficient of variation of bin-level
  LW mass (`W = a*L^b` per 1cm bin, not on a single mean length) for the
  dominant species' own length bins in that haul. `n_bins_dominant_species`
  is reported alongside so mechanically-low CV from sparse bins can be told
  apart from genuine single-cohort catches (see data-quality check below).

## Confound check 1 — correlation matrix

| var1 | var2 | n | pearson_r | spearman_r | flag_gt_0.4 |
|---|---|---|---|---|---|
| D | size_CV | 12069 | -0.149 | -0.150 | FALSE |
| D | N | 12069 | 0.262 | 0.363 | FALSE |
| D | B_obs | 12069 | 0.156 | 0.137 | FALSE |
| D | n_bins_dominant_species | 12069 | -0.153 | -0.162 | FALSE |
| size_CV | N | 12069 | -0.105 | -0.242 | FALSE |
| size_CV | B_obs | 12069 | -0.134 | -0.207 | FALSE |
| size_CV | n_bins_dominant_species | 12069 | 0.430 | 0.429 | TRUE |
| N | B_obs | 12069 | 0.524 | 0.746 | TRUE |
| N | n_bins_dominant_species | 12069 | -0.029 | -0.100 | FALSE |
| B_obs | n_bins_dominant_species | 12069 | 0.096 | 0.106 | FALSE |

**2 of 10 pairs** exceed |r| > 0.4 (either Pearson or Spearman) and are
 therefore read via the conditional (within-B_obs-quartile) check below rather than marginally.

## Confound check 2 — conditional on B_obs quartile

Correlation of `D` and `size_CV` with `ln_ratio` (`log(B_pred/B_obs)`),
recomputed separately within each existing `B_obs` quartile:

| b_obs_quartile | predictor | n | pearson_r | spearman_r |
|---|---|---|---|---|
| 1 | D | 3018 | 0.040 | 0.067 |
| 1 | size_CV | 3018 | -0.052 | -0.070 |
| 2 | D | 3017 | 0.130 | 0.145 |
| 2 | size_CV | 3017 | -0.057 | -0.064 |
| 3 | D | 3017 | 0.098 | 0.104 |
| 3 | size_CV | 3017 | -0.048 | -0.067 |
| 4 | D | 3017 | 0.168 | 0.171 |
| 4 | size_CV | 3017 | -0.108 | -0.148 |

If these within-quartile correlations are materially weaker than the marginal
correlations above, dominance/size-homogeneity are largely re-describing the
biomass-magnitude bias rather than an independent axis of failure.

## Confound check 3 — data-quality cross-check

Bottom decile of `size_CV` (n = 1207 of 12069 hauls). Fraction falling in
 flagged years (1985 discovery-file era, 1998, 2013-2014 dropout spikes):

| flagged | n_low_cv_hauls | pct_of_low_cv_total |
|---|---|---|
| TRUE | 161 | 13.3 |
| FALSE | 1046 | 86.7 |

Same check against ICES rectangles with any recorded EEoS-filter dropout
(`h1_dropout_by_stat_rec.csv`, bonus check beyond the year-based ask):

| flagged | n_low_cv_hauls | pct_of_low_cv_total |
|---|---|---|
| TRUE | 310 | 25.7 |
| FALSE | 897 | 74.3 |

Per `H1_dropout_diagnosis.md`, the 1998/2013-2014 spikes and 1985-only HL
coverage are historical artefacts already fixed upstream (NA-LW propagation);
the current `datras_hl_raw.rds` build (`hl_bins_full` mode) has 0 mean-length-
fallback hauls. A concentration of low-CV hauls in these years would still be
worth flagging as a residual data-quality signal rather than pure ecology.

## Confound check 4 — taxonomic breakdown (descriptive only)

Dominant-species identity in the extreme decile of each metric — descriptive
only, **not** used to construct a species-based flag or filter:

| metric | extreme | dominant_species | n_hauls | pct_of_decile |
|---|---|---|---|---|
| D | top | Trisopterus esmarkii | 329 | 27.3 |
| D | top | Clupea harengus | 326 | 27.0 |
| D | top | Sprattus sprattus | 318 | 26.4 |
| D | top | Limanda limanda | 104 | 8.6 |
| D | top | Merlangius merlangus | 71 | 5.9 |
| D | top | Melanogrammus aeglefinus | 21 | 1.7 |
| D | top | Eutrigla gurnardus | 15 | 1.2 |
| D | top | Ammodytes marinus | 6 | 0.5 |
| D | top | Scomber scombrus | 6 | 0.5 |
| D | top | Ammodytes | 5 | 0.4 |
| D | top | Ammodytes tobianus | 1 | 0.1 |
| D | top | Gadus morhua | 1 | 0.1 |
| D | top | Hippoglossoides platessoides | 1 | 0.1 |
| D | top | Pollachius virens | 1 | 0.1 |
| D | top | Trisopterus minutus | 1 | 0.1 |
| size_CV | top | Melanogrammus aeglefinus | 391 | 32.4 |
| size_CV | top | Limanda limanda | 284 | 23.5 |
| size_CV | top | Merlangius merlangus | 177 | 14.7 |
| size_CV | top | Sprattus sprattus | 100 | 8.3 |
| size_CV | top | Clupea harengus | 88 | 7.3 |
| size_CV | top | Trisopterus esmarkii | 65 | 5.4 |
| size_CV | top | Hippoglossoides platessoides | 56 | 4.6 |
| size_CV | top | Trisopterus minutus | 11 | 0.9 |
| size_CV | top | Gadus morhua | 5 | 0.4 |
| size_CV | top | Trachurus trachurus | 4 | 0.3 |
| size_CV | top | Eutrigla gurnardus | 3 | 0.2 |
| size_CV | top | Pleuronectes platessa | 3 | 0.2 |
| size_CV | top | Scomber scombrus | 3 | 0.2 |
| size_CV | top | Gadiculus argenteus | 2 | 0.2 |
| size_CV | top | Pollachius virens | 2 | 0.2 |
| size_CV | top | Scyliorhinus canicula | 2 | 0.2 |
| size_CV | top | Amblyraja radiata | 1 | 0.1 |
| size_CV | top | Ammodytes | 1 | 0.1 |
| size_CV | top | Ammodytes marinus | 1 | 0.1 |
| size_CV | top | Buglossidium luteum | 1 | 0.1 |
| size_CV | top | Merluccius merluccius | 1 | 0.1 |
| size_CV | top | Myoxocephalus scorpius | 1 | 0.1 |
| size_CV | top | Pomatoschistus | 1 | 0.1 |
| size_CV | top | Raja clavata | 1 | 0.1 |
| size_CV | top | Spondyliosoma cantharus | 1 | 0.1 |
| size_CV | top | Trisopterus luscus | 1 | 0.1 |
| size_CV | bottom | Clupea harengus | 626 | 51.9 |
| size_CV | bottom | Trisopterus esmarkii | 270 | 22.4 |
| size_CV | bottom | Sprattus sprattus | 182 | 15.1 |
| size_CV | bottom | Scomber scombrus | 29 | 2.4 |
| size_CV | bottom | Merlangius merlangus | 28 | 2.3 |
| size_CV | bottom | Melanogrammus aeglefinus | 21 | 1.7 |
| size_CV | bottom | Pollachius virens | 8 | 0.7 |
| size_CV | bottom | Limanda limanda | 6 | 0.5 |
| size_CV | bottom | Sardina pilchardus | 5 | 0.4 |
| size_CV | bottom | Engraulis encrasicolus | 4 | 0.3 |
| size_CV | bottom | Ammodytes | 3 | 0.2 |
| size_CV | bottom | Micromesistius poutassou | 3 | 0.2 |
| size_CV | bottom | Scyliorhinus canicula | 3 | 0.2 |
| size_CV | bottom | Trachurus trachurus | 3 | 0.2 |
| size_CV | bottom | Ammodytes marinus | 2 | 0.2 |
| size_CV | bottom | Ammodytes tobianus | 2 | 0.2 |
| size_CV | bottom | Eutrigla gurnardus | 2 | 0.2 |
| size_CV | bottom | Trisopterus minutus | 2 | 0.2 |
| size_CV | bottom | Agonus cataphractus | 1 | 0.1 |
| size_CV | bottom | Arnoglossus laterna | 1 | 0.1 |
| size_CV | bottom | Callionymus lyra | 1 | 0.1 |
| size_CV | bottom | Hyperoplus immaculatus | 1 | 0.1 |
| size_CV | bottom | Platichthys flesus | 1 | 0.1 |
| size_CV | bottom | Pomatoschistus | 1 | 0.1 |
| size_CV | bottom | Spinachia spinachia | 1 | 0.1 |
| size_CV | bottom | Taurulus bubalis | 1 | 0.1 |

## Binned reporting — D (Berger-Parker dominance)

Quartile table, formatted analogously to the existing `B_obs`-quartile
magnitude-bias result (`h1_catchability_by_quartile.csv`), with IQR added:

| bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|
| 1 | 3018 | 0.160 | 0.459 | 3.704 | 2.557 | 4.704 |
| 2 | 3017 | 0.459 | 0.601 | 4.120 | 2.784 | 5.215 |
| 3 | 3017 | 0.602 | 0.780 | 4.200 | 2.914 | 5.394 |
| 4 | 3017 | 0.780 | 1.000 | 4.409 | 2.994 | 5.830 |

Decile table (finer resolution; sample size per bin ≈ n/10):

| bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|
| 1 | 1207 | 0.160 | 0.363 | 3.454 | 2.351 | 4.447 |
| 2 | 1207 | 0.363 | 0.430 | 3.905 | 2.690 | 4.857 |
| 3 | 1207 | 0.430 | 0.484 | 3.948 | 2.718 | 5.028 |
| 4 | 1207 | 0.484 | 0.539 | 4.064 | 2.742 | 5.188 |
| 5 | 1207 | 0.539 | 0.601 | 4.238 | 2.877 | 5.328 |
| 6 | 1207 | 0.602 | 0.669 | 4.226 | 2.894 | 5.410 |
| 7 | 1207 | 0.669 | 0.741 | 4.233 | 2.968 | 5.398 |
| 8 | 1207 | 0.741 | 0.822 | 4.157 | 2.836 | 5.410 |
| 9 | 1207 | 0.822 | 0.901 | 4.309 | 2.992 | 5.720 |
| 10 | 1206 | 0.901 | 1.000 | 4.570 | 3.106 | 6.223 |

Broken out within each existing B_obs quartile (D quartile x B_obs quartile):

| b_obs_quartile | bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 755 | 0.160 | 0.440 | 2.906 | 1.850 | 3.642 |
| 1 | 2 | 755 | 0.441 | 0.582 | 3.154 | 2.128 | 3.921 |
| 1 | 3 | 754 | 0.582 | 0.765 | 3.152 | 2.056 | 4.013 |
| 1 | 4 | 754 | 0.765 | 0.997 | 3.108 | 1.935 | 4.158 |
| 2 | 1 | 755 | 0.186 | 0.433 | 3.604 | 2.678 | 4.331 |
| 2 | 2 | 754 | 0.433 | 0.565 | 3.935 | 2.699 | 4.703 |
| 2 | 3 | 754 | 0.565 | 0.735 | 3.984 | 3.044 | 4.824 |
| 2 | 4 | 754 | 0.735 | 0.998 | 4.037 | 2.945 | 5.017 |
| 3 | 1 | 755 | 0.193 | 0.459 | 4.367 | 3.103 | 5.033 |
| 3 | 2 | 754 | 0.459 | 0.595 | 4.558 | 3.493 | 5.355 |
| 3 | 3 | 754 | 0.596 | 0.771 | 4.562 | 3.434 | 5.483 |
| 3 | 4 | 754 | 0.772 | 1.000 | 4.607 | 3.333 | 5.745 |
| 4 | 1 | 755 | 0.222 | 0.508 | 5.205 | 3.356 | 6.146 |
| 4 | 2 | 754 | 0.508 | 0.667 | 5.430 | 4.068 | 6.445 |
| 4 | 3 | 754 | 0.668 | 0.836 | 5.490 | 3.755 | 6.609 |
| 4 | 4 | 754 | 0.836 | 1.000 | 5.942 | 4.240 | 7.795 |

## Binned reporting — size_CV (dominant-species size homogeneity)

Quartile table:

| bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|
| 1 | 3018 | 0.000 | 0.325 | 4.458 | 3.061 | 5.798 |
| 2 | 3017 | 0.325 | 0.454 | 4.176 | 2.853 | 5.316 |
| 3 | 3017 | 0.454 | 0.623 | 4.062 | 2.689 | 5.108 |
| 4 | 3017 | 0.623 | 4.432 | 3.726 | 2.626 | 4.755 |

Decile table:

| bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|
| 1 | 1207 | 0.000 | 0.240 | 4.480 | 3.123 | 6.003 |
| 2 | 1207 | 0.240 | 0.299 | 4.396 | 3.009 | 5.660 |
| 3 | 1207 | 0.299 | 0.350 | 4.421 | 3.057 | 5.627 |
| 4 | 1207 | 0.350 | 0.402 | 4.220 | 2.848 | 5.362 |
| 5 | 1207 | 0.402 | 0.454 | 4.106 | 2.739 | 5.161 |
| 6 | 1207 | 0.454 | 0.515 | 4.126 | 2.653 | 5.140 |
| 7 | 1207 | 0.515 | 0.585 | 4.038 | 2.743 | 5.110 |
| 8 | 1207 | 0.585 | 0.671 | 3.905 | 2.623 | 5.008 |
| 9 | 1207 | 0.671 | 0.808 | 3.788 | 2.630 | 4.741 |
| 10 | 1206 | 0.808 | 4.432 | 3.611 | 2.612 | 4.656 |

Broken out within each existing B_obs quartile (size_CV quartile x B_obs quartile):

| b_obs_quartile | bin | n | metric_min | metric_max | median_ratio | iqr_ratio_low | iqr_ratio_high |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 755 | 0.000 | 0.370 | 3.228 | 2.179 | 4.033 |
| 1 | 2 | 755 | 0.370 | 0.508 | 3.072 | 1.896 | 4.031 |
| 1 | 3 | 754 | 0.509 | 0.687 | 3.018 | 1.879 | 3.920 |
| 1 | 4 | 754 | 0.687 | 4.432 | 2.932 | 2.007 | 3.736 |
| 2 | 1 | 755 | 0.000 | 0.344 | 3.955 | 2.937 | 4.781 |
| 2 | 2 | 754 | 0.344 | 0.484 | 4.013 | 2.785 | 4.882 |
| 2 | 3 | 754 | 0.484 | 0.648 | 3.851 | 2.904 | 4.673 |
| 2 | 4 | 754 | 0.648 | 3.672 | 3.723 | 2.784 | 4.461 |
| 3 | 1 | 755 | 0.000 | 0.326 | 4.664 | 3.442 | 5.557 |
| 3 | 2 | 754 | 0.326 | 0.446 | 4.546 | 3.348 | 5.406 |
| 3 | 3 | 754 | 0.446 | 0.613 | 4.593 | 3.458 | 5.404 |
| 3 | 4 | 754 | 0.613 | 2.688 | 4.306 | 3.186 | 5.138 |
| 4 | 1 | 755 | 0.096 | 0.283 | 5.976 | 4.532 | 7.213 |
| 4 | 2 | 754 | 0.283 | 0.386 | 5.474 | 3.790 | 6.627 |
| 4 | 3 | 754 | 0.386 | 0.539 | 5.344 | 3.673 | 6.433 |
| 4 | 4 | 754 | 0.542 | 1.763 | 5.188 | 3.909 | 6.229 |

## Notes

- `D` and `size_CV` are reported and binned **separately** throughout — no
  composite dominance score is constructed at this stage.
- This diagnostic reuses `outputs/datras_hl_raw.rds`, `outputs/fishbase_lw_lookup_v2.csv`,
  `outputs/haul_state_variables.rds`, and `outputs/haul_eeos_predictions.rds` exactly
  as built by the production pipeline; it does not recompute or alter N, E, B_pred,
  or the primary residual.

![Ratio vs D](../outputs/figures/h1_dominance_ratio_vs_D.png)

![Ratio vs size_CV](../outputs/figures/h1_dominance_ratio_vs_sizeCV.png)

*Outputs: `outputs/h1_dominance_haul_table.csv`, `outputs/h1_dominance_correlation_matrix.csv`,
`outputs/h1_dominance_conditional_correlation.csv`, `outputs/h1_dominance_dataquality_*.csv`,
`outputs/h1_dominance_taxonomic_breakdown.csv`, `outputs/h1_dominance_by_D_bins.csv`,
`outputs/h1_dominance_by_sizeCV_bins.csv`, `outputs/h1_dominance_by_D_x_bobs_quartile.csv`,
`outputs/h1_dominance_by_sizeCV_x_bobs_quartile.csv`.*
