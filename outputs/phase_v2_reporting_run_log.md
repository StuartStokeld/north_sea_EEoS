# phase_v2 reporting — slopes / CAR / proportional effects / figures — run log

Builds the reporting stack for the policy-anchored primary model (`phase_v2`: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015). Does not overwrite original-phase `h2h3_wb_*` artifacts.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_sessionInfo.txt
glmmTMB 1.1.14; spaMM available = TRUE (4.6.65)

## Inputs
Loaded /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds: 10464 hauls, 158 rectangles; formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)
phase_v2 levels: 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015

## Primary (RE) fixed effects and phase-specific slopes
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_primary_fixed_effects.csv

## CAR sensitivity (phase_v2)
Fitting CAR: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
CAR fit time: 0.36 sec.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_car_fixed_effects.csv

### Pooled FP_between CAR contrast (H2 figure overlay)
Fitting: residual ~ FP_between + FP_within * phase_v2 + adjacency(1 | stat_rec)
Pooled CAR fit time: 0.27 sec.
Pooled FP_between = +0.000033 (SE 0.013364; 95% CI [-0.026161, +0.026227]; p = 0.998)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_pooled_between_coef.csv

## Phase-specific slopes
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_fp_slopes_by_phase.csv
### wb_car_v2
  - [H2_spatial_between] FP_between / 1985-1991: slope=-0.0628 SE=0.0159 CI=[-0.0940, -0.0316] p=7.97e-05
  - [H2_spatial_between] FP_between / 1992-2001: slope=-0.0244 SE=0.0147 CI=[-0.0532, +0.0045] p=0.098
  - [H2_spatial_between] FP_between / 2002-2007: slope=+0.0841 SE=0.0159 CI=[+0.0529, +0.1153] p=1.27e-07
  - [H2_spatial_between] FP_between / 2008-2015: slope=+0.0133 SE=0.0151 CI=[-0.0164, +0.0429] p=0.38
  - [H3_temporal_within] FP_within / 1985-1991: slope=+0.0273 SE=0.0148 CI=[-0.0017, +0.0562] p=0.0652
  - [H3_temporal_within] FP_within / 1992-2001: slope=+0.0462 SE=0.0172 CI=[+0.0126, +0.0798] p=0.00707
  - [H3_temporal_within] FP_within / 2002-2007: slope=-0.0043 SE=0.0174 CI=[-0.0383, +0.0298] p=0.806
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0567 SE=0.0122 CI=[+0.0328, +0.0806] p=3.19e-06
### wb_primary_v2
  - [H2_spatial_between] FP_between / 1985-1991: slope=-0.0247 SE=0.0155 CI=[-0.0551, +0.0057] p=0.111
  - [H2_spatial_between] FP_between / 1992-2001: slope=+0.0142 SE=0.0141 CI=[-0.0135, +0.0419] p=0.316
  - [H2_spatial_between] FP_between / 2002-2007: slope=+0.1223 SE=0.0155 CI=[+0.0920, +0.1526] p=2.57e-15
  - [H2_spatial_between] FP_between / 2008-2015: slope=+0.0519 SE=0.0146 CI=[+0.0232, +0.0806] p=0.000392
  - [H3_temporal_within] FP_within / 1985-1991: slope=+0.0269 SE=0.0148 CI=[-0.0022, +0.0560] p=0.0702
  - [H3_temporal_within] FP_within / 1992-2001: slope=+0.0437 SE=0.0172 CI=[+0.0100, +0.0774] p=0.011
  - [H3_temporal_within] FP_within / 2002-2007: slope=-0.0018 SE=0.0174 CI=[-0.0360, +0.0324] p=0.919
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0577 SE=0.0122 CI=[+0.0337, +0.0816] p=2.4e-06

## Proportional effect sizes
H2 source: `wb_car_v2` FP_between phase slopes.
H3 source: primary (`wb_primary_v2`) FP_within phase slopes.
Max |primary − CAR| FP_within phase slope = 0.00249.
Baselines (phase_v2 median residual → ratio):
  - 1985-1991: median residual=-1.4030; ratio=0.2459 (24.6%); IQR(FP_within)=0.967
  - 1992-2001: median residual=-1.4799; ratio=0.2277 (22.8%); IQR(FP_within)=0.476
  - 2002-2007: median residual=-1.3630; ratio=0.2559 (25.6%); IQR(FP_within)=0.793
  - 2008-2015: median residual=-1.2747; ratio=0.2795 (28.0%); IQR(FP_within)=0.890
H2 IQR(FP_between) across 158 rectangles = 1.7546.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H2.csv
  H2 1985-1991: slope=-0.0628; IQR gap widened by 3.40% [-4.96, -1.76]
  H2 1992-2001: slope=-0.0244; IQR gap widened by 1.23% [-2.63, 0.23]
  H2 2002-2007: slope=+0.0841; IQR gap closed by 5.47% [3.35, 7.71]
  H2 2008-2015: slope=+0.0133; IQR gap closed by 0.91% [-1.10, 3.03]
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H3.csv
  H3 1985-1991: slope=+0.0269; IQR gap closed by 0.86% [-0.07, 1.81]
  H3 1992-2001: slope=+0.0437; IQR gap closed by 0.62% [0.14, 1.11]
  H3 2002-2007: slope=-0.0018; IQR gap widened by 0.05% [-0.97, 0.90]
  H3 2008-2015: slope=+0.0577; IQR gap closed by 2.04% [1.18, 2.92]

## Figures
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_wb_gap_change_by_phase.png
H2 pooled overlay (IQR scale): +0.002% [-1.513%, +1.589%]
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H2_gap_change_by_phase.png
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H3_gap_change_by_phase.png
  H2 fig 1985-1991: gap -3.40% [-4.96, -1.76]; significant
  H2 fig 1992-2001: gap -1.23% [-2.63, +0.23]; not significant
  H2 fig 2002-2007: gap +5.47% [3.35, +7.71]; significant
  H2 fig 2008-2015: gap +0.91% [-1.10, +3.03]; not significant
  H3 fig 1985-1991: gap +0.86% [-0.07, +1.81]; not significant
  H3 fig 1992-2001: gap +0.62% [0.14, +1.11]; significant
  H3 fig 2002-2007: gap -0.05% [-0.97, +0.90]; not significant
  H3 fig 2008-2015: gap +2.04% [1.18, +2.92]; significant
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_model_objects.rds
Updated /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds with CAR / pooled fits.

## Outputs
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_fp_slopes_by_phase.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_primary_fixed_effects.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_car_fixed_effects.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H2.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H3.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_pooled_between_coef.csv
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_model_objects.rds
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_wb_gap_change_by_phase.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H2_gap_change_by_phase.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H3_gap_change_by_phase.png
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_sessionInfo.txt
- /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_run_log.md (this file)

display_discussion drafts not rewritten by this script — use these phase_v2_* artifacts when updating One page / H2 / H3 drafts.
