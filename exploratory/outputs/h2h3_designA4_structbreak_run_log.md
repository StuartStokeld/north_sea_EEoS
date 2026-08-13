# Structural-break check — fishing-pressure time series — run log

Checks whether Section A's eyeballed 4-phase pattern is statistically supported. Reports the comparison as numbers only — does not force a match, does not adjust the visual boundaries, does not recommend which set of boundaries to use.

## Method / package
Package: strucchange (version 1.5.4, Bai & Perron dynamic-programming least-squares breakpoint estimation, `breakpoints()`). NOTE: strucchange was NOT part of this project's renv-managed dependency set at the time this script was written — it was installed ad hoc (`install.packages('strucchange', repos = 'https://cloud.r-project.org')`) into the ambient R library, not added to renv.lock. Flagged as an environment note, not resolved here (out of scope for this task to restructure the project's dependency management).
Model: `mean_hours ~ year` — breaks allow BOTH the intercept and the linear slope (hours/year) to change at each breakpoint, matching how Section A's within-phase trends (A3) are themselves framed (separate linear fits per phase), rather than a pure mean-shift ('~1') model. This is a method/formula choice, not a value dictated by the briefing, so it is surfaced explicitly here.

## Data
Read /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA1_year_fishing_summary.csv directly (31 year rows, 1985-2015), NOT recomputed from raw Couce data — same full 215-rectangle Couce universe as Section A. Series used: `mean_hours` (mean annual Couce fishing hours across all rectangles with coverage that year).

## Break-count model-selection criterion (item 2)
`h` = 0.15 (package default; minimum segment length = h*n = 4.65 ~= 5 observations) determines the maximum number of breaks the dynamic program can evaluate for n=31 observations: max feasible = 6 breaks (>= the 5 the briefing asked for; all 7 candidate counts from 0 to 6 are reported below, not just 0-5).
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA4_structbreak_criterion.csv
  - m=0 breaks: RSS=1.499e+08, BIC=575.4027
  - m=1 breaks: RSS=1.717e+07, BIC=518.5399
  - m=2 breaks: RSS=8.784e+06, BIC=508.0669
  - m=3 breaks: RSS=3.421e+06, BIC=489.1390
  - m=4 breaks: RSS=2.359e+06, BIC=487.9123  <- BIC-optimal
  - m=5 breaks: RSS=1.980e+06, BIC=492.7946
  - m=6 breaks: RSS=1.826e+06, BIC=500.5751
BIC-optimal break count = 4. NOTE: BIC at m=4 (487.9123) is only marginally lower than at m=3 (489.1390) — reported as a numeric fact about how decisive the BIC minimum is, not interpreted further.

## Selected break year(s) at the BIC-optimal count, with 95% CIs (item 3)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA4_structbreak_years.csv
  - break 1: year = 1989, 95% CI = [1988, 1990]
  - break 2: year = 1997, 95% CI = [1996, 1998]
  - break 3: year = 2001, 95% CI = [2000, 2002]
  - break 4: year = 2008, 95% CI = [2007, 2009]

## Comparison: statistically located break years vs visual boundaries (item 4)
Visual boundaries (from the briefing, eyeballed from h3_pre_D1/D3 and Section A): 1990, 2000, 2010. Statistically located break years (BIC-optimal, m=4): 1989, 1997, 2001, 2008.
  - [visual_boundary_to_nearest_statistical_break] 1990 -> nearest 1989 (|diff| = 1 years)
  - [visual_boundary_to_nearest_statistical_break] 2000 -> nearest 2001 (|diff| = 1 years)
  - [visual_boundary_to_nearest_statistical_break] 2010 -> nearest 2008 (|diff| = 2 years)
  - [statistical_break_to_nearest_visual_boundary] 1989 -> nearest 1990 (|diff| = 1 years)
  - [statistical_break_to_nearest_visual_boundary] 1997 -> nearest 2000 (|diff| = 3 years)
  - [statistical_break_to_nearest_visual_boundary] 2001 -> nearest 2000 (|diff| = 1 years)
  - [statistical_break_to_nearest_visual_boundary] 2008 -> nearest 2010 (|diff| = 2 years)

## Figure
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designA4_structbreak_series.png

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA4_structbreak_criterion.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA4_structbreak_years.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_designA4_structbreak_series.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_designA4_structbreak_run_log.md (this file)
