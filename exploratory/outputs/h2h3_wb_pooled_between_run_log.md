# H2 contrast: pooled FP_between (no phase interaction) — run log

ADDITIONAL CONTRAST ONLY — not a replacement for the primary phase-specific H2 model (`FP_between * phase`). Quantifies a single time-stable between-rectangle fishing-pressure slope for discussion. No figures produced.

## Session
sessionInfo: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_sessionInfo.txt
spaMM 4.6.1

## Data and specification
Data reused from h2h3_wb_model_objects.rds: 10464 hauls, 158 rectangles.
adjMatrix reused from Round 2 RDS (same queen adjacency as prior CAR fits).

Contrast formula (this run):
  residual ~ FP_between + FP_within * phase + adjacency(1 | stat_rec)  [CAR REML]
Primary within-between CAR (unchanged; for reference only):
  residual ~ FP_between * phase + FP_within * phase + adjacency(1 | stat_rec)

Confirmed changes: ONLY FP_between × phase dropped (FP_between enters as a single pooled main effect).
Confirmed unchanged: FP_within * phase; phase factor; CAR adjacency; no biomass; same panel.

## Model fit
Fit time: 0.85 sec.
Saved model object: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_model_objects.rds

## Convergence
warnings=0 (none); robust_to_starting_value=TRUE; fitted_rho=0.131288; any non-finite FE SE=FALSE
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_convergence.csv

## 1. Pooled FP_between coefficient
FP_between (pooled): estimate=+0.000586, SE=0.013363, 95% CI=[-0.025605, +0.026778], z=+0.044, p=0.965002
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_coef.csv
FP_within-related terms still in this contrast fit: FP_within; FP_within:phase1989-2000; FP_within:phase2001-2007; FP_within:phase2008-2015
Confirmed: no FP_between × phase interaction terms in the contrast fit.

## 2. Pooled vs phase-specific FP_between (CAR)
Side-by-side:
  - all_years_pooled    est=+0.000586  SE=0.013363  CI=[-0.025605, +0.026778]  p=0.965002
  - 1985-1988           est=-0.043421  SE=0.018251  CI=[-0.079193, -0.007649]  p=0.0173539
  - 1989-2000           est=-0.051798  SE=0.014488  CI=[-0.080193, -0.023402]  p=0.000349806
  - 2001-2007           est=+0.087210  SE=0.015424  CI=[+0.056978, +0.117442]  p=1.56732e-08
  - 2008-2015           est=+0.013313  SE=0.015115  CI=[-0.016312, +0.042938]  p=0.37843
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_comparison.csv

## 3. Plain numeric observation (not an interpretation)
The pooled FP_between estimate (+0.0006) lies between the most negative phase-specific CAR slope (-0.0518 in 1989-2000) and the most positive (+0.0872 in 2001-2007); |pooled| (0.0006) is smaller than the largest |phase-specific| slope (0.0872), so pooling partially cancels the phase-level reversal in magnitude.

## Outputs (no figures)
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_coef.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_comparison.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_convergence.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_sessionInfo.txt
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_run_log.md (this file)
