# H2/H3 proportional effect-size reporting — run log

Reporting/formatting only. No models refit. Re-expresses primary-model phase slopes from the biomass-free results run on proportional and gap-closure scales.
Confirmed source model is biomass-free primary specification.

## Phase-specific baseline choice
Baseline for percentage-point ratio shift and percent-gap-closed: each phase's MEDIAN residual (ratio = exp(median residual)). Mean residual is recorded but not used.
Why median: residual is left-skewed within phases; median is the typical haul and matches the interpretation note's worked-example convention.
  - 1985-1988: n=1583; median residual=-1.4032 (ratio=0.2458); mean residual=-1.3453 (ratio=0.2605); log(hours+1) IQR=1.156 (hours Q25–Q75: 6610–21011)
  - 1989-2000: n=3933; median residual=-1.4533 (ratio=0.2338); mean residual=-1.3711 (ratio=0.2538); log(hours+1) IQR=1.561 (hours Q25–Q75: 6292–29981)
  - 2001-2007: n=2447; median residual=-1.3827 (ratio=0.2509); mean residual=-1.2703 (ratio=0.2808); log(hours+1) IQR=2.153 (hours Q25–Q75: 2662–22933)
  - 2008-2015: n=2501; median residual=-1.2747 (ratio=0.2795); mean residual=-1.1791 (ratio=0.3076); log(hours+1) IQR=2.209 (hours Q25–Q75: 1222–11131)

## Table
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_proportional_effects.csv

Column groups:
  - pct_ratio_change_* : percent change in B_obs/B_pred (= (exp(slope*Δx)-1)*100)
  - pp_ratio_shift_*   : percentage-point change in the ratio (e.g. 25.0% → 26.4% = +1.4 pp)
  - pct_gap_closed_*   : percent of remaining gap (1 - ratio) closed; negative = gap widens
  - *_doubling / *_iqr : Δx = log(2) vs phase-specific IQR of log(hours+1)
  - *_lo / *_hi        : same transforms of the slope 95% CI bounds

### Summary (doubling of fishing hours)
  - 1985-1988: slope=-0.0008; % ratio change=-0.05% [-1.95, +1.88]; ratio 24.6% → 24.6% (-0.01 pp); gap closed=-0.02% [-0.64, +0.61]
  - 1989-2000: slope=-0.0056; % ratio change=-0.39% [-1.61, +0.85]; ratio 23.4% → 23.3% (-0.09 pp); gap closed=-0.12% [-0.49, +0.26]
  - 2001-2007: slope=+0.0885; % ratio change=+6.33% [+5.20, +7.47]; ratio 25.1% → 26.7% (+1.59 pp); gap closed=+2.12% [+1.74, +2.50]
  - 2008-2015: slope=+0.0580; % ratio change=+4.10% [+3.10, +5.11]; ratio 28.0% → 29.1% (+1.15 pp); gap closed=+1.59% [+1.20, +1.98]

### Summary (phase log-hours IQR)
  - 1985-1988: Δx(IQR)=1.156; % ratio change=-0.09% [-3.24, +3.16]; ratio 24.6% → 24.6% (-0.02 pp); gap closed=-0.03% [-1.06, +1.03]
  - 1989-2000: Δx(IQR)=1.561; % ratio change=-0.87% [-3.59, +1.92]; ratio 23.4% → 23.2% (-0.20 pp); gap closed=-0.27% [-1.10, +0.59]
  - 2001-2007: Δx(IQR)=2.153; % ratio change=+20.99% [+17.05, +25.08]; ratio 25.1% → 30.4% (+5.27 pp); gap closed=+7.03% [+5.71, +8.40]
  - 2008-2015: Δx(IQR)=2.209; % ratio change=+13.67% [+10.23, +17.22]; ratio 28.0% → 31.8% (+3.82 pp); gap closed=+5.30% [+3.97, +6.68]

## Figure
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_gap_closed_by_phase.png
Percent-of-gap-closed (doubling) by phase with CI bars; year-midpoint x-axis and break-year markers match the absolute-scale FP-effect figure layout.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_proportional_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_results_gap_closed_by_phase.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_results_proportional_effects_run_log.md (this file)
