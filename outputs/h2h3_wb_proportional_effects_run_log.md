# H2/H3 within-between proportional effect sizes — run log

Reporting only — no models refit. Re-expresses decomposed FP_between (H2) and FP_within (H3) phase slopes on proportional / gap-change scales.

SUPERSEDES the earlier blended-term proportional table (`outputs/h2h3_results_proportional_effects.csv` / `h2h3_results_gap_closed_by_phase.png`). That table used the undecomposed `log_hours_total` slopes and is no longer the primary proportional reporting artefact. Use the H2/H3 tables produced here instead.

## Coefficient sources
H2: CAR model (`wb_car`) FP_between phase slopes (per design).
H3: primary model (`wb_primary`) FP_within phase slopes. CAR within slopes agree closely (max |primary − CAR| phase-slope difference logged below); primary used as the reported H3 source.
Max |primary − CAR| FP_within phase slope = 0.00209 (phases agree for reporting).

## Baseline and Δx benchmarks
Baseline for percentage-point ratio shift and gap-change: each phase's MEDIAN residual; starting ratio = exp(median residual). Same convention as the (now superseded) blended-term proportional table.
  - 1985-1988: median residual=-1.4032; baseline ratio=0.2458 (24.6%); phase IQR(FP_within)=1.072
  - 1989-2000: median residual=-1.4533; baseline ratio=0.2338 (23.4%); phase IQR(FP_within)=0.530
  - 2001-2007: median residual=-1.3827; baseline ratio=0.2509 (25.1%); phase IQR(FP_within)=0.787
  - 2008-2015: median residual=-1.2747; baseline ratio=0.2795 (28.0%); phase IQR(FP_within)=0.890
H2 IQR benchmark: IQR(FP_between) across 158 rectangles = 1.7546 (applied to all phases; FP_between is time-invariant).
H3 IQR benchmark: phase-specific IQR(FP_within) among hauls in that phase (realistic within-rectangle year-to-year spread in that period).
H3 doubling convention: NOT APPLICABLE. FP_within is a mean-zero deviation (can be negative); 'doubling fishing hours' is not a meaningful contrast on this scale. Only the IQR-of-deviation benchmark is reported for H3.

## H2 table (FP_between, CAR)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H2.csv
Gap labelling: positive pct_gap_change = gap closed/narrowed; negative = gap widened.
  - 1985-1988: slope=-0.0434; doubling: ratio change=-2.96%, gap widened by 0.97% [-1.74, -0.17]; IQR(Δx=1.755): ratio change=-7.34%, gap widened by 2.39%
  - 1989-2000: slope=-0.0518; doubling: ratio change=-3.53%, gap widened by 1.08% [-1.65, -0.49]; IQR(Δx=1.755): ratio change=-8.69%, gap widened by 2.65%
  - 2001-2007: slope=+0.0872; doubling: ratio change=+6.23%, gap closed by 2.09% [1.35, 2.84]; IQR(Δx=1.755): ratio change=+16.53%, gap closed by 5.54%
  - 2008-2015: slope=+0.0133; doubling: ratio change=+0.93%, gap closed by 0.36% [-0.44, 1.17]; IQR(Δx=1.755): ratio change=+2.36%, gap closed by 0.92%

## H3 table (FP_within, primary) — IQR-of-deviation only
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H3.csv
NOTE: No doubling-convention columns for H3. FP_within is mean-zero by construction; a doubling benchmark does not apply. See doubling_note column in the CSV.
  - 1985-1988: slope=+0.0317; IQR(Δx=1.072): ratio change=+3.46% [-0.42, +7.49]; ratio 24.6% → 25.4% (+0.85 pp); gap closed by 1.13% [-0.14, 2.44]
  - 1989-2000: slope=+0.0415; IQR(Δx=0.530): ratio change=+2.23% [+0.60, +3.87]; ratio 23.4% → 23.9% (+0.52 pp); gap closed by 0.68% [0.18, 1.18]
  - 2001-2007: slope=-0.0116; IQR(Δx=0.787): ratio change=-0.91% [-3.29, +1.53]; ratio 25.1% → 24.9% (-0.23 pp); gap widened by 0.30% [-1.10, 0.51]
  - 2008-2015: slope=+0.0571; IQR(Δx=0.890): ratio change=+5.21% [+3.00, +7.48]; ratio 28.0% → 29.4% (+1.46 pp); gap closed by 2.02% [1.16, 2.90]

## Figure
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_gap_change_by_phase.png
Two series: H2 gap-change under doubling; H3 gap-change under phase IQR(FP_within). Sign convention labelled on the y-axis (closed vs widened).

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H2.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_H3.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_wb_gap_change_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_wb_proportional_effects_run_log.md (this file)
