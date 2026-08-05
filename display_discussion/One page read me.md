Our hypotheses move from establishing whether the EEoS baseline is reachable in this system, to detecting a spatial disturbance signal, to detecting a temporal one. **The challenge is isolated these signals from noise.**

- H1 - EEoS predicts observed biomass B from (S, N, E) better than chance at the haul level
- H2 - Rectangles with higher fishing pressure show larger EEoS residuals
- H3 - Within rectangles, periods of higher fishing pressure are associated with larger EEoS residuals

H1 was tested at individual haul level (n = 12,069), comparing EEoS's parameter-free biomass prediction from observed S, N, and E directly against observed catch biomass with no fitting.

Results:

- EEoS does not predict absolute haul-level biomass in the North Sea trawl survey, a simple productivity map outperforms the full equation of state by a wide margin.
- H1 quantitative result: log_r² = −0.223 (worse than mean-prediction) against a 1:1 productivity baseline of 0.736, but cor² = 0.926: High correlation with a systematic ~3.7–4.1× overprediction.
- EEoS overprediction is not a fixed offset, but **grows with observed biomass magnitude**, from roughly 3× in the lowest quartile to over 5× in the highest, indicating the catchability-type bias scales with haul size rather than acting as a simple constant correction factor.
- However, the **relative productivity-ratio structure Harte et al. used to validate EEoS is reproduced, in the form Harte associated with disturbed systems.**
- “Big shoals” (hauls at the extreme of the dominance/size-homogeneity distribution) show larger overprediction than hauls with a more typical spread across species and sizes. This gradient is smaller than the biomass-magnitude effect, but it persists after conditioning on biomass magnitude, indicating a second, largely independent source of prediction failure.
- Detail: methods: [H1_methods_draft.md], headline results: [H1_results_summary.md]

**H2/H3 were tested using a single mixed-effects model at haul level** (n = 10,464 hauls, 158 rectangles). 
Rectangle-year fishing pressure was decomposed into a between-rectangle term (a rectangle's own long-run mean fishing pressure) and a within-rectangle term (that year's deviation from the rectangle's own mean). 
Rectangle entered as a random intercept with a CAR spatial-correlation structure. 
Temporal periods were defined using policy-anchored breakpoints (1992, 2002, 2008), giving phases 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015 (1992 and 2002 CFP reforms; 2008 LTMP / MSFD).

**H2 results:**
- Fishing pressure predicts residual, but the direction is not stable over time.
- Comparing a rectangle at the 25th percentile of persistent fishing pressure to one at the 75th percentile (within the same phase), the overprediction gap (the shortfall between EEoS's predicted biomass and the observed biomass) widens by 3.4% in 1985–1991 (significant; higher-fishing rectangles predict *worse*, consistent with the original disturbance hypothesis). In 1992–2001 the point estimate is still a 1.2% widening but is not significant. The association then reverses to close by 5.5% in 2002–2007 (higher-fishing rectangles predict *better*), before showing no detectable effect in 2008–2015 (0.9% closed, not significant).
- *sensitivity analysis* if fishing pressure's spatial effect on residual were assumed stable across the full study period (no phase interaction), the estimated effect is +0.00003 (95% CI [−0.026, +0.026], p ≈ 1.00) statistically indistinguishable from zero. The phase-specific analysis was necessary to detect a real, substantial effect that exists within specific periods.
- These percentages reflect the fitted relationship across all 158 rectangles in each phase, adjusted for spatial autocorrelation. Raw, unadjusted correlations (still from the original-phase run; not yet recomputed for `phase_v2`) were weak in both early windows (r = −0.07 and −0.09, r² under 1%, not significant) but stronger later (r = +0.47, r² = 22% in the old 2001–2007 window; r = +0.29, r² = 9% in 2008–2015), directionally consistent with the adjusted post-2002 pattern.


**H3 results:**
- The within-rectangle relationship is intermittent.
- Comparing a rectangle's own lower-fishing year to its own higher-fishing year (25th vs. 75th percentile of that rectangle's year-to-year fishing-pressure fluctuation), the overprediction gap closes by 0.9% in 1985–1991 (not significant) and by 0.6% in 1992–2001 (significant) — a rectangle's own fishing-pressure rise tracks a smaller gap that same year — then widens by 0.05% in 2002–2007 (not significant), before closing again by 2.0% in 2008–2015 (significant).

This partially confirms and partially contradicts the working hypothesis stated above, neither H2 nor H3 shows the single, uniform disturbance signal originally anticipated.
**Fishing pressure and time period, together, explain only about 5% of why haul-level prediction error varies (marginal R² = 0.051).** The remaining ~95% is driven by something else entirely.

H2/ H3 limitation:
- The low R² (0.051) is unlikely to improve through further refinement of the spatial or temporal structure alone. Three spatial specifications (plain random intercept, continuous distance-decay, CAR) were already tested and R² barely moved (0.047–0.051). Switching from data-driven phases to policy-anchored `phase_v2` likewise left marginal R² ≈ 0.051 — same ceiling, because the model is still built from the same two predictors (fishing pressure, phase). Raising R² meaningfully would need new predictors already flagged in H1: the dominance/size-CV diagnostics, environmental covariates, or a catchability correction, rather than a different model of fishing pressure and time.
- Couce data has some areas missing (Skagerrak/Kattegat and eastern English Channel) (highlighted/ map in H2 methods)
- Couce fishing-hours data may not adequately capture disturbance - SAR data may offer a better signal (highlighted in H2 methods)
