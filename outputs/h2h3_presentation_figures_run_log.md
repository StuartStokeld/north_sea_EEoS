# Final presentation figures (H2 / H3) — run log

Standalone supervisor figures. No models refit. All phase bar values read directly from proportional-effects CSVs (IQR convention).

## Source files
- H2 phase bars: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H2.csv
- H3 phase bars: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H3.csv
- H2 pooled overlay: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_pooled_between_coef.csv

## H2 pooled overlay conversion (same IQR scale as phase bars)
Pooled slope = +0.000586 (95% CI [-0.025605, +0.026778]) from pooled coef CSV.
Delta-x = IQR(FP_between) = 1.7546 from H2 proportional CSV; baseline ratio = mean of 4 phase baselines = 0.2525 from H2 proportional CSV.
Converted gap change (IQR): +0.035% [-1.484%, +1.625%] (transform from run_h2h3_wb_proportional_effects.R; not a refit).

## Figure 1 (H2)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_presentation_H2_gap_change_by_phase.png
  - 1985-1988: gap change -2.39% [-4.23%, -0.43%]; significant
  - 1989-2000: gap change -2.65% [-4.01%, -1.23%]; significant
  - 2001-2007: gap change +5.54% [3.52%, +7.66%]; significant
  - 2008-2015: gap change +0.92% [-1.09%, +3.04%]; not significant

## Figure 2 (H3)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_presentation_H3_gap_change_by_phase.png
No pooled overlay (FP_within pooled contrast not computed).
  - 1985-1988: gap change +1.13% [-0.14%, +2.44%]; IQR(FP_within)=1.072; not significant
  - 1989-2000: gap change +0.68% [0.18%, +1.18%]; IQR(FP_within)=0.530; significant
  - 2001-2007: gap change -0.30% [-1.10%, +0.51%]; IQR(FP_within)=0.787; not significant
  - 2008-2015: gap change +2.02% [1.16%, +2.90%]; IQR(FP_within)=0.890; significant

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_presentation_H2_gap_change_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_presentation_H3_gap_change_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_presentation_figures_run_log.md (this file)
