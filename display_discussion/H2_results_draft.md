# H2 results – draft

*Updated July 2026: primary H2 narrative below uses the within-between decomposition (persistent between-rectangle fishing pressure, CAR-adjusted, phase-specific slopes). An earlier cross-sectional rectangle-scale analysis (OLS / SEM) is retained at the end for reference.*

---

## H2 — does a rectangle's persistent fishing-pressure level predict its prediction error?

| Phase | Effect | Significant? | Gap change (typical low- vs high-fishing rectangle) |
| --- | --- | --- | --- |
| 1985–1988 | Negative — persistently higher-fishing rectangles have a *larger* overprediction gap | Yes (p ≈ 0.017) | widened 2.4% [0.4, 4.2] |
| 1989–2000 | Negative, same direction | Yes (p ≈ 3.5×10⁻⁴) | widened 2.7% [1.2, 4.0] |
| 2001–2007 | Positive — reverses; higher-fishing rectangles now have a *smaller* gap | Yes (p ≈ 1.6×10⁻⁸) | closed 5.5% [3.5, 7.7] |
| 2008–2015 | No detectable effect | No (p ≈ 0.38) | closed 0.9% [−1.1, 3.0] (CI crosses zero) |

"Typical low- vs high-fishing rectangle" compares a rectangle at the 25th percentile of persistent fishing pressure to one at the 75th percentile, within the same phase (IQR of `FP_between` across rectangles).

![H2: fishing pressure's spatial effect reverses direction over time — proportional gap change by phase (CAR model, IQR benchmark). Dashed red line = pooled FP_between estimate (all years, no phase interaction).](../outputs/figures/h2h3_presentation_H2_gap_change_by_phase.png)

These percentages reflect the fitted relationship across all 158 rectangles in each phase, adjusted for spatial autocorrelation; the raw, unadjusted correlation between a rectangle's persistent fishing pressure and its mean residual is weak in both early phases (r = −0.07 and −0.09, r² under 1%, not significant) but strengthens substantially post-2001 (r = +0.47, r² = 22% in 2001–2007; r = +0.29, r² = 9% in 2008–2015), consistent in direction with the adjusted slopes throughout.

Notably, the raw correlation remains statistically significant in 2008–2015 (p < 0.001) even though the spatially-adjusted model slope for that phase is not (p ≈ 0.38) — this is expected given the between-rectangle term's known sensitivity to spatial confounding, and is a direct illustration of why the CAR-adjusted, not raw, estimate is used as the primary H2 result.

Individual rectangles vary considerably around these fitted relationships (e.g. mean residual across rectangles spans −1.79 to −0.38 in 2001–2007), so the reported percentage describes a central tendency, not a bound on any one rectangle's actual shift.

*Sources: `outputs/h2h3_wb_proportional_effects_H2.csv`, `outputs/h2h3_raw_correlation_H2.csv`, `outputs/figures/h2h3_presentation_H2_gap_change_by_phase.png`.*

---

## Earlier analysis: cross-sectional rectangle panel (OLS / SEM)

Fishing pressure was not associated with EEoS residual magnitude at the rectangle scale in this earlier specification, and this null result was robust to spatial model specification.

### Key figure (earlier analysis)

![H2 topline result: OLS negative association between fishing hours and mean |EEoS residual| is not robust after spatial error correction (SEM). N = 158 ICES rectangles, 1985–2015.](../outputs/figures/h2_topline_result.png)

The primary OLS regression of mean absolute EEoS residual on mean annual fishing hours (Couce et al. 2020) across 158 ICES rectangles returned a significant negative coefficient (β = −2.74 × 10⁻⁶, SE = 1.06 × 10⁻⁶, t = −2.59, p = 0.010, R² = 0.041), opposite in direction to the pre-registered hypothesis that heavier fishing pressure would be associated with larger residuals. However, OLS residuals showed strong positive spatial autocorrelation (Moran's I = 0.56, Geary's C = 0.43, both p ≪ 0.001). After accounting for this via a spatial error model, the fishing-pressure effect was no longer significant (β = −6.1 × 10⁻⁷, p = 0.541, λ = 0.858), indicating the naive OLS association was not robust to spatial structure in the data.

This null result held across sample definitions. Re-running the spatial models on panels built with lower and higher haul-count thresholds (≥5 hauls, N = 161; ≥20 hauls, N = 156) showed no consistent, significant relationship at any threshold (SEM β range: −6.1 × 10⁻⁷ to +1.6 × 10⁻⁷, all p > 0.5), and the sign of the (non-significant) coefficient reversed at the lowest threshold. The apparent stability of the OLS result across thresholds (p = 0.009–0.017) therefore reflected an unmodelled spatial artefact rather than a genuine, threshold-independent signal. The OLS result was also not attributable to high-leverage rectangles: excluding the seven rectangles exceeding the conventional Cook's distance threshold strengthened, rather than weakened, the negative association (β = −3.68 × 10⁻⁶, p = 0.00096, N = 151), ruling out an influential-point explanation for the effect that spatial correction subsequently removed.

Adding mean log observed biomass as a covariate reproduced a pattern seen in the original catchability-confound hypothesis under OLS, that fishing pressure became positively associated with residual magnitude once biomass was controlled for (β = +2.82 × 10⁻⁶ equivalent flip, R² = 0.55) but this too failed to survive spatial correction. Under SEM, fishing hours remained non-significant (β = −1.1 × 10⁻⁷, p = 0.861), while mean log biomass was strongly and independently significant (β = 0.187, p ≪ 0.001, λ = 0.887). This pattern was confirmed by formal model comparison: Lagrange Multiplier tests favoured a spatial lag specification for the primary model (robust LM-lag p = 0.039, robust LM-error p = 0.844), but the reverse for the biomass-covariate model (robust LM-error p ≪ 0.001, robust LM-lag p = 0.773), consistent with an AIC comparison strongly favouring the spatial error specification for the biomass model (ΔAIC = 54.8 in favour of SEM) while the two specifications were statistically indistinguishable for the primary model (ΔAIC = 0.08). Under the spatial lag specification for the biomass-covariate model, fishing pressure again appeared positive and significant (β = +1.91 × 10⁻⁶, p = 0.0016), reproducing the apparent effect seen under OLS. Given this specification fit the data substantially worse by both AIC and formal LM diagnostics, this result is interpreted as an artefact of spatial misspecification rather than evidence for a genuine catchability-driven relationship between fishing pressure and EEoS residual magnitude.

Overall, this cross-sectional specification found no support for a time-averaged association between fishing pressure and EEoS residual magnitude at the ICES rectangle scale. The within-between decomposition above separates persistent between-rectangle variation from within-rectangle temporal change and is the primary H2 result going forward.

---

*Limitations and potential expansions of H2: see [H2 methods draft](H2_methods_draft.md).*
