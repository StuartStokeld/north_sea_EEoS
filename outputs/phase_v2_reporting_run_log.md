# phase_v2 reporting — slopes / CAR / proportional effects / figures — run log

Builds the reporting stack for the policy-anchored primary model (`phase_v2`: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015). Does not overwrite original-phase `h2h3_wb_*` artifacts.

## Session
sessionInfo: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_reporting_sessionInfo.txt
glmmTMB 1.1.14; spaMM available = TRUE (4.6.65)

## Inputs
Loaded /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds: 8670 hauls, 158 rectangles; formula: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)
phase_v2 levels: 1985-1991 | 1992-2001 | 2002-2007 | 2008-2015

## Primary (RE) fixed effects and phase-specific slopes
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_primary_fixed_effects.csv

## CAR sensitivity (phase_v2)
Fitting CAR: residual ~ FP_between * phase_v2 + FP_within * phase_v2 + adjacency(1 | stat_rec)
CAR fit time: 0.30 sec.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_car_fixed_effects.csv

### Pooled FP_between CAR contrast (H2 figure overlay)
Fitting: residual ~ FP_between + FP_within * phase_v2 + adjacency(1 | stat_rec)
Pooled CAR fit time: 0.21 sec.
Pooled FP_between = +0.003750 (SE 0.011985; 95% CI [-0.019740, +0.027240]; p = 0.7544)
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_pooled_between_coef.csv

## Phase-specific slopes
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_fp_slopes_by_phase.csv
### wb_car_v2
  - [H2_spatial_between] FP_between / 1985-1991: slope=+0.0042 SE=0.0136 CI=[-0.0225, +0.0309] p=0.758
  - [H2_spatial_between] FP_between / 1992-2001: slope=+0.0257 SE=0.0128 CI=[+0.0007, +0.0507] p=0.0438
  - [H2_spatial_between] FP_between / 2002-2007: slope=+0.0129 SE=0.0132 CI=[-0.0131, +0.0388] p=0.331
  - [H2_spatial_between] FP_between / 2008-2015: slope=-0.0255 SE=0.0129 CI=[-0.0507, -0.0003] p=0.0477
  - [H3_temporal_within] FP_within / 1985-1991: slope=-0.0114 SE=0.0106 CI=[-0.0321, +0.0094] p=0.283
  - [H3_temporal_within] FP_within / 1992-2001: slope=+0.0343 SE=0.0123 CI=[+0.0102, +0.0584] p=0.00529
  - [H3_temporal_within] FP_within / 2002-2007: slope=-0.0316 SE=0.0124 CI=[-0.0558, -0.0073] p=0.0107
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0321 SE=0.0089 CI=[+0.0146, +0.0496] p=0.000312
### wb_primary_v2
  - [H2_spatial_between] FP_between / 1985-1991: slope=+0.0245 SE=0.0126 CI=[-0.0002, +0.0493] p=0.0519
  - [H2_spatial_between] FP_between / 1992-2001: slope=+0.0467 SE=0.0116 CI=[+0.0239, +0.0694] p=5.78e-05
  - [H2_spatial_between] FP_between / 2002-2007: slope=+0.0339 SE=0.0122 CI=[+0.0100, +0.0577] p=0.00536
  - [H2_spatial_between] FP_between / 2008-2015: slope=-0.0042 SE=0.0118 CI=[-0.0272, +0.0189] p=0.723
  - [H3_temporal_within] FP_within / 1985-1991: slope=-0.0126 SE=0.0106 CI=[-0.0334, +0.0082] p=0.235
  - [H3_temporal_within] FP_within / 1992-2001: slope=+0.0316 SE=0.0123 CI=[+0.0074, +0.0557] p=0.0104
  - [H3_temporal_within] FP_within / 2002-2007: slope=-0.0297 SE=0.0124 CI=[-0.0540, -0.0054] p=0.0166
  - [H3_temporal_within] FP_within / 2008-2015: slope=+0.0333 SE=0.0089 CI=[+0.0158, +0.0508] p=0.000187

## Proportional effect sizes
H2 source: `wb_car_v2` FP_between phase slopes.
H3 source: `wb_car_v2` FP_within phase slopes.
Max |CAR − RE| FP_within phase slope = 0.00274 (RE retained as companion).
Baselines (phase_v2 median residual → ratio):
  - 1985-1991: median residual=-0.6359; ratio=0.5295 (52.9%); IQR(FP_within)=0.950
  - 1992-2001: median residual=-0.7154; ratio=0.4890 (48.9%); IQR(FP_within)=0.481
  - 2002-2007: median residual=-0.6752; ratio=0.5091 (50.9%); IQR(FP_within)=0.777
  - 2008-2015: median residual=-0.6441; ratio=0.5251 (52.5%); IQR(FP_within)=0.896
H2 IQR(FP_between) across 158 rectangles = 1.8156.
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H2.csv
  H2 1985-1991: slope=+0.0042; IQR gap closed by 0.86% [-4.51, 6.49]
  H2 1992-2001: slope=+0.0257; IQR gap closed by 4.58% [0.13, 9.23]
  H2 2002-2007: slope=+0.0129; IQR gap closed by 2.45% [-2.43, 7.57]
  H2 2008-2015: slope=-0.0255; IQR gap widened by 5.00% [-9.72, -0.05]
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/phase_v2_proportional_effects_H3.csv
  H3 1985-1991: slope=-0.0114; IQR gap widened by 1.21% [-3.38, 1.01]
  H3 1992-2001: slope=+0.0343; IQR gap closed by 1.59% [0.47, 2.73]
  H3 2002-2007: slope=-0.0316; IQR gap widened by 2.51% [-4.40, -0.59]
  H3 2008-2015: slope=+0.0321; IQR gap closed by 3.23% [1.46, 5.02]

## Figures
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_wb_gap_change_by_phase.png
H2 pooled overlay (IQR scale): +0.720% [-3.711%, +5.344%]
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H2_gap_change_by_phase.png
Saved: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/figures/phase_v2_presentation_H3_gap_change_by_phase.png
  H2 fig 1985-1991: gap +0.86% [-4.51, +6.49]; not significant
  H2 fig 1992-2001: gap +4.58% [0.13, +9.23]; significant
  H2 fig 2002-2007: gap +2.45% [-2.43, +7.57]; not significant
  H2 fig 2008-2015: gap -5.00% [-9.72, -0.05]; significant
  H3 fig 1985-1991: gap -1.21% [-3.38, +1.01]; not significant
  H3 fig 1992-2001: gap +1.59% [0.47, +2.73]; significant
  H3 fig 2002-2007: gap -2.51% [-4.40, -0.59]; significant
  H3 fig 2008-2015: gap +3.23% [1.46, +5.02]; significant
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

Supplementary CAR diagnostics (leave-one-rectangle-out, CAR R², signed vs unsigned residual audit) are a separate step that reads this fit_car and does not refit the signed primary model: Rscript --vanilla pipeline/run_h2h3_car_reporting_diagnostics.R
