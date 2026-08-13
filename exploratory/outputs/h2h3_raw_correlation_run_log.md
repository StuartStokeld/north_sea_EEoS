# H2/H3 supplementary raw correlation & range statistics — run log

## Framing (raw vs adjusted)
These correlations are raw and unadjusted: they do not control for the other decomposed fishing-pressure term, CAR spatial structure, or partial pooling the way the fitted within-between models do. They are expected to be directionally consistent with, but not numerically identical to, the adjusted phase slopes already reported.

## Data sources (reused; no new data prep)
- Haul-level panel: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_model_objects.rds (10464 hauls, 158 rectangles)
- Model slopes for sign comparison: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_fp_slopes_by_phase.csv
  - H2 reference slopes: wb_car FP_between (same as primary H2 reporting)
  - H3 reference slopes: wb_primary FP_within (same as primary H3 reporting)

## H2 — rectangle-level correlation (mean residual vs FP_between), per phase
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_raw_correlation_H2.csv
  - 1985-1988: n=156 rectangles; r=-0.0666 (r²=0.0044, p=0.4086); model slope=-0.0434; min=-1.856 (rect 36F0, FP_between=8.554); max=-0.545 (rect 43F3, FP_between=7.972)
  - 1989-2000: n=156 rectangles; r=-0.0852 (r²=0.0073, p=0.2901); model slope=-0.0518; min=-1.780 (rect 47E8, FP_between=6.957); max=-0.713 (rect 40F4, FP_between=9.173)
  - 2001-2007: n=158 rectangles; r=+0.4736 (r²=0.2243, p=3.293e-10); model slope=+0.0872; min=-1.792 (rect 48F0, FP_between=8.483); max=-0.375 (rect 34F1, FP_between=7.317)
  - 2008-2015: n=157 rectangles; r=+0.2928 (r²=0.0857, p=0.0001981); model slope=+0.0133; min=-1.742 (rect 49F0, FP_between=9.270); max=-0.435 (rect 34F1, FP_between=7.317)

## H3 — rectangle-year-level correlation (mean residual vs FP_within), per phase
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_raw_correlation_H3.csv
  - 1985-1988: n=605 rectangle-years; r=+0.0488 (r²=0.0024, p=0.2311); model slope=+0.0317; min=-2.525 (rect 44F0, year 1988, FP_within=+0.509); max=-0.095 (rect 35F2, year 1985, FP_within=-0.632)
  - 1989-2000: n=1825 rectangle-years; r=+0.0507 (r²=0.0026, p=0.03018); model slope=+0.0415; min=-2.304 (rect 39F0, year 1995, FP_within=+0.715); max=+1.062 (rect 40F4, year 1997, FP_within=+0.086)
  - 2001-2007: n=1091 rectangle-years; r=+0.1357 (r²=0.0184, p=6.854e-06); model slope=-0.0116; min=-2.331 (rect 45F1, year 2002, FP_within=+0.228); max=+0.216 (rect 34F1, year 2006, FP_within=-0.411)  ** SIGN DISAGREEMENT **
  - 2008-2015: n=1247 rectangle-years; r=+0.1639 (r²=0.0269, p=5.787e-09); model slope=+0.0571; min=-2.890 (rect 46F2, year 2015, FP_within=+1.048); max=+0.182 (rect 38F8, year 2010, FP_within=+1.323)

## Sign-disagreement flags (raw r vs adjusted model slope)
  - H3 2001-2007: raw r=+0.1357 vs model slope=-0.0116

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_raw_correlation_H2.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_raw_correlation_H3.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_raw_correlation_run_log.md (this file)
