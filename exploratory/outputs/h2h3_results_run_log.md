# H2/H3 shared hierarchical model — RESULTS run log

RESULTS RUN (not a feasibility check). Fishing-pressure and fishing-pressure x phase coefficients ARE interpreted below as the substantive H2 and H3 findings, within the limits of what the coefficients directly support.

CORRECTION / RE-RUN: biomass (`mean_ln_B_obs`) has been REMOVED from the primary model and both sensitivity models. Corrected fixed effects = `log_hours_total * phase` only. This re-run overwrites the previous with-biomass results outputs; a before/after comparison of headline numbers is logged below.

## Session and package versions
Full sessionInfo() written to: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_sessionInfo.txt
Key packages: R 4.4.3; glmmTMB 1.1.14; spaMM 4.6.1; mgcv 1.9.1; DHARMa 0.4.7 (available = TRUE).
Environment note (unchanged from feasibility): glmmTMB and spaMM were installed ad hoc into the ambient/user library, NOT added to renv.lock. This script MUST be run with `Rscript --vanilla pipeline/run_h2h3_shared_model_results.R`.

## Data and analysis universe
Rebuilt analysis data via build_feasibility_data() (identical pipeline to feasibility checks): 10464 hauls across 158 rectangles, years 1985–2015. Universe = established H2 panel (>= 10 hauls AND Couce coverage) = 158 rectangles. Dropped 6 hauls lacking rectangle-year Couce.
Phase definition (break year = first year of new phase; 1997 excluded):
  - 1985-1988 (years 1985–1988): n_hauls = 1583, n_rectangles = 156
  - 1989-2000 (years 1989–2000): n_hauls = 3933, n_rectangles = 156
  - 2001-2007 (years 2001–2007): n_hauls = 2447, n_rectangles = 158
  - 2008-2015 (years 2008–2015): n_hauls = 2501, n_rectangles = 157

## Model fitting (corrected specification — no biomass)
Round 2 data stored vs rebuilt checksum: n=10464/10464, n_rect=158/158, sum_resid match=TRUE
Prior feasibility / temporal-robustness fitted objects included mean_ln_B_obs and are therefore not reused. All three models are REFIT under the corrected formula. Round 2 adjMatrix is still reused (same queen-adjacency neighbour definition).
Corrected formulas used for this re-run:
  - Primary: residual ~ log_hours_total * phase + (1 | stat_rec)  [glmmTMB REML]
  - CAR:     residual ~ log_hours_total * phase + adjacency(1 | stat_rec)  [spaMM REML]
  - GAM:     residual ~ s(year, k=8) + s(year, by=log_hours_total, k=8) + s(stat_rec, bs="re")  [mgcv REML]
Fitting primary model ...
Primary fit time: 0.70 sec.
CAR adjacency matrix: REUSED from Round 2 RDS.
Fitting CAR sensitivity model ...
CAR fit time: 1.09 sec.
Fitting GAM sensitivity model (expect ~1–2 min) ...
GAM fit time: 111.23 sec.
Fitted formulas:
  - primary: residual ~ log_hours_total * phase + (1 | stat_rec)
  - car: residual ~ log_hours_total * phase + adjacency(1 | stat_rec)
  - gam: residual ~ s(year, k = 8) + s(year, by = log_hours_total, k = 8) + s(stat_rec, bs = "re")
Confirmed: mean_ln_B_obs absent from all three fitted formulas.

### Reuse summary
  - Primary (1|stat_rec): REFIT (biomass removed)
  - CAR adjacency model:  REFIT (biomass removed); adjMatrix REUSED from Round 2
  - GAM smooth year:      REFIT (biomass removed)

## Specification note
Biomass removed from primary and both sensitivities relative to the immediately preceding results run. Otherwise unchanged: plain RE primary; CAR and varying-coefficient GAM sensitivities; phase breaks 1989/2001/2008; canonical residual; 158-rectangle H2 panel.
Saved model objects: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_model_objects.rds

## 1. Primary model — full fixed-effects table
Primary formula: residual ~ log_hours_total * phase + (1 | stat_rec) [REML; no biomass]
Fixed effects (estimate, SE, z, p) with H2/H3 labels:
  - (Intercept)                               est = -1.3135  SE = 0.1324  z = -9.920  p = 3.399e-23  [control (intercept)]
  - log_hours_total                           est = -0.0008  SE = 0.0141  z = -0.056  p = 9.554e-01  [H2 (FP effect, reference phase 1985-1988)]
  - phase1989-2000                            est = +0.0020  SE = 0.1379  z = +0.015  p = 9.883e-01  [control (phase main effect)]
  - phase2001-2007                            est = -0.7394  SE = 0.1365  z = -5.417  p = 6.065e-08  [control (phase main effect)]
  - phase2008-2015                            est = -0.3364  SE = 0.1338  z = -2.515  p = 1.191e-02  [control (phase main effect)]
  - log_hours_total:phase1989-2000            est = -0.0048  SE = 0.0147  z = -0.329  p = 7.424e-01  [H3 (FP x phase interaction vs reference)]
  - log_hours_total:phase2001-2007            est = +0.0893  SE = 0.0147  z = +6.076  p = 1.232e-09  [H3 (FP x phase interaction vs reference)]
  - log_hours_total:phase2008-2015            est = +0.0588  SE = 0.0145  z = +4.051  p = 5.109e-05  [H3 (FP x phase interaction vs reference)]
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_primary_fixed_effects.csv

## Phase-specific fishing-pressure slopes (derived; used for H2 reporting)
Per-phase FP slopes (primary):
  - 1985-1988: slope = -0.0008, SE = 0.0141, 95% CI = [-0.0285, +0.0269]
  - 1989-2000: slope = -0.0056, SE = 0.0091, 95% CI = [-0.0234, +0.0122]
  - 2001-2007: slope = +0.0885, SE = 0.0079, 95% CI = [+0.0731, +0.1039]
  - 2008-2015: slope = +0.0580, SE = 0.0071, 95% CI = [+0.0441, +0.0719]
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_fp_slopes_by_phase.csv
Joint Wald H3 (all FP x phase interactions = 0): chi2(3) = 113.045, p = 2.427e-24
Joint Wald H2-style (FP slope = 0 in every phase): chi2(4) = 200.576, p = 2.825e-42
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_wald_tests.csv

## 2. Side-by-side sensitivity comparison
Comparison table has two blocks: (1) parametric fixed effects shared by primary and CAR (GAM columns filled only where the GAM has a matching parametric term — typically Intercept only under the no-biomass GAM; phase terms are NA for GAM by design); (2) phase-specific FP slopes for all three models (GAM = mean of year-by-year by-smooth within each phase).
Parametric block — primary vs CAR estimate shift:
  - (Intercept): primary -1.3135 -> CAR -1.1707 (+10.9%)
  - log_hours_total: primary -0.0008 -> CAR -0.0076 (-865.2%)
  - phase1989-2000: primary +0.0020 -> CAR +0.0154 (+660.9%)
  - phase2001-2007: primary -0.7394 -> CAR -0.7439 (-0.6%)
  - phase2008-2015: primary -0.3364 -> CAR -0.3517 (-4.5%)
  - log_hours_total:phase1989-2000: primary -0.0048 -> CAR -0.0062 (-29.0%)
  - log_hours_total:phase2001-2007: primary +0.0893 -> CAR +0.0893 (-0.0%)
  - log_hours_total:phase2008-2015: primary +0.0588 -> CAR +0.0596 (+1.3%)
Phase-slope block:
  - fp_slope:1985-1988: primary -0.0008 | CAR -0.0076 | GAM(mean) +0.0097
  - fp_slope:1989-2000: primary -0.0056 | CAR -0.0138 | GAM(mean) -0.0043
  - fp_slope:2001-2007: primary +0.0885 | CAR +0.0817 | GAM(mean) +0.1084
  - fp_slope:2008-2015: primary +0.0580 | CAR +0.0520 | GAM(mean) +0.0582
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_sensitivity_comparison.csv

## 3. Effect visualisation: FP effect by phase (primary) + GAM overlay
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_fp_effect_by_phase.png
Saved year-grid slopes: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_fp_slopes_by_year.csv

## 4. Random-effects summary (partial pooling)
Rectangle random-intercept SD and partial pooling for the primary model: Spearman cor log(n_hauls) vs shrinkage_ratio = 0.440 across 158 rectangles (positive = more hauls, less shrinkage — expected pattern).
Rectangle intercept SD = 0.1679 (95% CI [0.1484, 0.1899]); residual SD = 0.4470.
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_partial_pooling.csv
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_partial_pooling.png

## 5. Model diagnostics
Primary model fit: Nakagawa R2 marginal = 0.0469, conditional = 0.1647 (manual Nakagawa–Schielzeth implementation — performance/MuMIn not installed). logLik = -6630.5, AIC = 13280.9, BIC = 13353.5.
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_model_fit.csv
Saved residual diagnostics figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_residual_diagnostics.png
DHARMa simulated-residual screen (primary model; n = 250 simulations):
  - Uniformity (KS): D = 0.0737, p = 7.698e-50
  - Dispersion: ratio/stat = 0.9976, p = 9.200e-01
  - NOTE: with n = 10,464, the KS uniformity test is almost always significant even for small departures; read alongside the QQ plot rather than as a hard reject of the model.
Influence screen (approximate Cook's D using FE hat matrix): flag if (Cook's D >= sample 99th percentile = 0.00088 AND |Pearson| > 3) OR (hat > 3 * mean hat = 0.0023). Flagged 352 / 10464 hauls. Top 5 by Cook's D among flags:
  - row 390, 42F2, year 1985: cooks_approx = 0.0114, hat = 0.0033, pearson = -5.289
  - row 73, 46E6, year 1985: cooks_approx = 0.0050, hat = 0.0310, pearson = -1.105
  - row 1056, 44E7, year 1987: cooks_approx = 0.0044, hat = 0.0025, pearson = +3.758
  - row 7934, 42E8, year 2007: cooks_approx = 0.0040, hat = 0.0079, pearson = +1.990
  - row 648, 49E8, year 1986: cooks_approx = 0.0031, hat = 0.0130, pearson = +1.357
Saved flagged observations: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_influence_flags.csv
NOTE: Cook's D here is approximate (ignores random-effect structure in the hat matrix). Use as a screen for write-up caveats, not as a formal deletion diagnostic.
Convergence (primary): code = 0, Hessian PD = TRUE, max|grad| = 1.633e-07, any NA SE = FALSE.

## 6. Plain-language summary of H2 and H3 findings

### Framing
Outcome = haul-level residual = log(B_obs) - log(B_pred). A positive fishing-pressure slope means higher Couce fishing hours associate with larger (more positive) residuals (observed biomass above EEOS prediction); a negative slope means the opposite. No biomass covariate. Effect sizes are on the log-residual scale per unit log(hours+1).

### H2 — Does fishing pressure predict residual?
In the reference phase (1985–1988), the fishing-pressure coefficient is -0.0008 (SE 0.0141, z = -0.056, p = 9.554e-01). Per-phase slopes (primary model):
  - 1985-1988: -0.0008 (95% CI [-0.0285, +0.0269]; z = -0.06, p = 9.554e-01)
  - 1989-2000: -0.0056 (95% CI [-0.0234, +0.0122]; z = -0.62, p = 5.372e-01)
  - 2001-2007: +0.0885 (95% CI [+0.0731, +0.1039]; z = +11.26, p = 2.027e-29)
  - 2008-2015: +0.0580 (95% CI [+0.0441, +0.0719]; z = +8.17, p = 3.031e-16)
Joint Wald test that the FP slope is zero in every phase: chi2(4) = 200.576, p = 2.825e-42. This is the global H2-style test under the primary model; phase-specific slopes and CIs above are the effect-size detail.
Factual H2 statement (coefficients only):
  Fishing-pressure slopes differ from zero at p < 0.05 in: 2001-2007 (+0.0885); 2008-2015 (+0.0580). Joint test that all phase slopes are zero: p = 2.825e-42. Magnitudes remain small on the residual scale (see slopes above); significance and effect size should be read together.

### H3 — Does the fishing-pressure–residual relationship change across phases?
  - FP x phase 1989-2000 (change in FP slope vs 1985-1988): -0.0048 (SE 0.0147, z = -0.329, p = 7.424e-01)
  - FP x phase 2001-2007 (change in FP slope vs 1985-1988): +0.0893 (SE 0.0147, z = +6.076, p = 1.232e-09)
  - FP x phase 2008-2015 (change in FP slope vs 1985-1988): +0.0588 (SE 0.0145, z = +4.051, p = 5.109e-05)
Joint Wald test that all FP x phase interactions are zero: chi2(3) = 113.045, p = 2.427e-24.
Factual H3 statement (coefficients only):
  At least one FP x phase interaction differs from zero at p < 0.05: log_hours_total:phase2001-2007 (+0.0893, p = 1.232e-09); log_hours_total:phase2008-2015 (+0.0588, p = 5.109e-05). Joint interaction test p = 2.427e-24. This supports a phase-dependent change in the fishing-pressure–residual relationship (H3) under the primary model; read alongside the per-phase slopes and the sensitivity overlay (GAM), which should agree if the design justification holds.

### Sensitivity agreement (primary vs CAR vs GAM)
Phase-slope correlation across the four phases: primary vs CAR = 1.000; primary vs GAM-within-phase mean = 0.988. Primary vs CAR parametric FE shifts are logged in section 2; the overlay figure shows primary step vs GAM smooth year-by-year.
No conclusions beyond the coefficients: this summary does not claim mechanism, management causality, or that non-significance equals evidence of no effect.

## Before/after: with-biomass previous run vs this no-biomass re-run
Source of previous numbers: hardcoded snapshot from the completed with-biomass results run (on-disk CSVs had already been overwritten by a no-biomass re-run).
Primary phase slopes (previous with biomass -> this re-run without):
  - 1985-1988: -0.0149 -> -0.0008 (delta +0.0141)
  - 1989-2000: -0.0251 -> -0.0056 (delta +0.0194)
  - 2001-2007: +0.0724 -> +0.0885 (delta +0.0162)
  - 2008-2015: +0.0430 -> +0.0580 (delta +0.0150)
Primary H3 interaction terms (previous -> this re-run):
  - log_hours_total:phase1989-2000: -0.0102 -> -0.0048 (delta +0.0054)
  - log_hours_total:phase2001-2007: +0.0872 -> +0.0893 (delta +0.0021)
  - log_hours_total:phase2008-2015: +0.0579 -> +0.0588 (delta +0.0010)
Nakagawa R2 (previous -> this): marginal 0.1042 -> 0.0469; conditional 0.1630 -> 0.1647.
Sensitivity phase-slope correlations this re-run: primary-CAR = 1.000; primary-GAM = 0.988 (previous with-biomass run: primary-CAR ≈ 1.000, primary-GAM ≈ 0.989).

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_primary_fixed_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_sensitivity_comparison.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_fp_slopes_by_phase.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_fp_slopes_by_year.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_partial_pooling.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_model_fit.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_influence_flags.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_wald_tests.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_fp_effect_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_partial_pooling.png
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_residual_diagnostics.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_sessionInfo.txt
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_run_log.md (this file)
