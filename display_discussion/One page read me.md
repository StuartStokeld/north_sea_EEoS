Our hypotheses move from establishing whether the EEoS baseline is reachable in this system, to detecting a spatial disturbance signal, to detecting a temporal one. **The challenge is isolated these signals from noise.**

- H1 - Can EEoS be used to predict biomass using the state variables from the haul data.
- H2 - Do failures in EEoS prediction correlate with disturbance spatially.
- H3 - Do failures in EEoS prediction correlate with disturbance temporally.

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
Temporal periods were defined using a structural-break analysis (1989, 2001, 2008), here a break also corresponds with the 2002 Common Fisheries Policy (introduced stricter fishing effort controls).

**H2 results:**
- Fishing pressure predicts residual, but the direction is not stable over time.
- Comparing a rectangle at the 25th percentile of persistent fishing pressure to one at the 75th percentile (within the same phase), the overprediction gap (the shortfall between EEoS's predicted biomass and the observed biomass) widens by 2.4% in 1985–1988 and 2.7% in 1989–2000 (higher-fishing rectangles predict *worse*, consistent with the original disturbance hypothesis), then reverses to close by 5.5% in 2001–2007 (higher-fishing rectangles predict *better*), before showing no detectable effect in 2008–2015 (0.9% closed, not significant).
- *sensitivity analysis* if fishing pressure's spatial effect on residual were assumed stable across the full study period (no phase interaction), the estimated effect is +0.0006 (95% CI [−0.026, +0.027], p = 0.97) statistically indistinguishable from zero. The phase-specific analysis was necessary to detect a real, substantial effect that exists within specific periods.
- These percentages reflect the fitted relationship across all 158 rectangles in each phase, adjusted for spatial autocorrelation; the raw, unadjusted correlation between a rectangle's persistent fishing pressure and its mean residual is weak in both early phases (r = −0.07 and −0.09, r² under 1%, not significant) but strengthens substantially post-2001 (r = +0.47, r² = 22% in 2001–2007; r = +0.29, r² = 9% in 2008–2015), consistent in direction with the adjusted slopes throughout. 


**H3 results:**
- The within-rectangle relationship is intermittent.
- Comparing a rectangle's own lower-fishing year to its own higher-fishing year (25th vs. 75th percentile of that rectangle's year-to-year fishing-pressure fluctuation), the overprediction gap closes by 1.1% in 1985–1988 (not significant) and by 0.7% in 1989–2000 (significant) — a rectangle's own fishing-pressure rise tracks a smaller gap that same year — then widens by 0.3% in 2001–2007 (not significant), before closing again by 2.0% in 2008–2015 (significant).

This partially confirms and partially contradicts the working hypothesis stated above, neither H2 nor H3 shows the single, uniform disturbance signal originally anticipated.
**Fishing pressure and time period, together, explain only about 5% of why haul-level prediction error varies (marginal R² = 0.051).** The remaining ~95% is driven by something else entirely.

H2/ H3 limitation:
- The low R² (0.051) is unlikely to improve through further refinement of the spatial or temporal structure alone. Three spatial specifications (plain random intercept, continuous distance-decay, CAR) were already tested and R² barely moved (0.047–0.051), and a different time-window choice faces the same ceiling, since it's still built from the same two predictors (fishing pressure, phase). Raising R² meaningfully would need new predictors already flagged in H1: the dominance/size-CV diagnostics, environmental covariates, or a catchability correction, rather than a different model of fishing pressure and time.
- Couce data has some areas missing (Skagerrak/Kattegat and eastern English Channel) (highlighted/ map in H2 methods)
- Couce fishing-hours data may not adequately capture disturbance - SAR data may offer a better signal (highlighted in H2 methods)

**Separating noise and bias from theory failure.**

METE fails when the community structure is not as the theory assumes: IE “intact/ equilibrium”. We can expect METE to fail in our project for two key reasons:

A) **Due to disturbance** (good failure, we are looking for this).

- The North Sea is a fundamentally disturbed system. METE would fail (for example) due to the selective removal of larger or longer-lived taxa under sustained fishing pressure, resulting in the aggregate state variables (S, N, E) no longer describe a community in the configuration METE assumes. High fishing pressure across the study area is our working expectation/ explanation for H1 failure.
- H2/H3 aim to isolate this, testing for correlation between magnitude of disturbance and magnitude theory prediction failure.

B) **Due to catchability or haul level unrepresentativeness** (bad failure: we may assume this is just theory failure, OR associate it spuriously with disturbance).

- Initial research design assumed that haul data accurately represents the ecological assemblage within the North Sea (because we have enough hauls and they are well distributed spatially and temporally). This is challenged by a key confounder *“catchability”* at the **haul level** (in fisheries = haul systematically underestimates true standing biomass by a roughly consistent factor).
- Alternatively, A **single haul** may not be an entirely accurate representation of wider community structure. Fish are highly mobile, often travel in shoals, and the state variable distribution in any given haul is driven by the ability for a haul to capture a random sample of the North Sea assemblage.
- The key example I started with and am using currently is a large haul of many similar sized individuals of a single species (big shoal). This is an unrepresentative snapshot, and we would expect METE to fail here, because it falls outside METE’s assumption that individuals are distributed across species and body sizes without systematic aggregation. There are other examples of these (’unexpected METE haul distributions’) but not yet explored.
- Detail: Disturbance vs sampling logic: [Methods_note_for_discussion.md]. Catchability / scale offset: [H1_catchability_scaling_exploration.md]. Dominance & size homogeneity (“big shoals”) → [H1_dominance_results_summary.md]. Full exploration → [H1_dominance_size_homogeneity_exploration.md].

We have two “extremes” of research design. 1) assume that we have enough hauls to not control for catchability at all 2) attempt to fully control for catchability.

1. The problem here is that we cannot say whether failure is due to A or B.
2. The problem here is that if we control too strongly, we might obscure the signal that we are looking for or give us spurious results.

I am attempting to find an informed balance between the two, accept that we cannot fully control for catchability but that there are low hanging fruit (starting with big shoals), and take reasonable steps to ensure that controlling does not end up circularly impacting our results (IE giving us a spurious result that we are “looking for”).

Broken down above / in **methods note for discussion:**

- Atypical-haul filtering strategy - second source of *rectangle-level covariate: could be tested for their own confound relationship with fishing pressure before being trusted as clean predictors.*

DynaMETE

- Predicts how disturbance shifts the state variables (S, N, E) over time, rather than treating them as fixed. This could explain why the within-rectangle relationship is intermittent rather than continuous: a rectangle's own fishing-pressure change may only track its residual once that rectangle has moved far enough along its own disturbance trajectory, not simply whenever fishing pressure itself moves.
- This would change what `B_pred` *is -* Instead of adding a term to the H2/H3 regression, it replaces the static, disturbance-blind EEoS prediction with one that lets S, N, E evolve dynamically under fishing-driven demographic change.
