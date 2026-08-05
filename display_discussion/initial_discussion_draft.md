The thesis moves in three steps of increasing specificity: is the EEoS baseline reachable at all (H1)? if not achievable in absolute terms, does its residual still carry a *spatial* disturbance signal (H2)? and does that signal also appear *temporally*, within a site (H3)? The core challenge throughout is isolating a disturbance signal from noise, artifact, and geography.

---

## H1 fails in absolute terms, partially holds in structure

EEoS does **not** predict absolute haul-level biomass, a simple productivity map outperforms it outright (log_r² = −0.223 vs. a 0.736 baseline). But correlation is high (cor² = 0.926): EEoS systematically **overpredicts** by ~3.7–4.1×, and that overprediction is not a fixed offset, it grows with haul biomass (3× lowest quartile → >5× highest) and independently with shoal homogeneity ("big shoal" hauls overpredict more, even after conditioning on biomass).

The theory fails as a stand-alone absolute predictor, but the *relative productivity-ratio structure* Harte used to validate EEoS elsewhere is reproduced, in the specific form Harte associates with disturbed systems. H1 fails on its literal terms while reproducing the qualitative signature the framework predicts for a disturbed community. Two candidate explanations remain: a scale-dependent catchability artifact, and/or a genuine disturbance signature, H1 alone cannot separate them.

## H2 “real” but very fragile, and partly a mechanism artifact

Fishing pressure predicts EEoS residual, but the *direction* is not stable: a significant disturbance-consistent effect in 1985–1991, weakening to non-significance by 2001, then **reversing** post-2002 before disappearing again by 2008–2015. A phase-invariant model finds nothing (effect ≈ 0), the signal only exists when the policy-era structure is modelled explicitly.

Extensive robustness testing complicates this further: the rectangle random intercept does not remove underlying spatial autocorrelation from the data, and a spatial-permutation bootstrap shows the post-2002 reversal is **not distinguishable from geographic clustering alone**. A specific mechanism, fishing-effort displacement to neighbouring rectangles, explains that confounding once modelled directly, but even this only nudges residual spatial clustering by ~3%.

Only the earliest phase (1985–1991) offers a fishing-pressure effect independent of geography. Everything after is either non-significant or explainable by spatial effort displacement rather than disturbance acting within each rectangle. H2 should be read/ interpreted with caution - a signal that exists in a narrow historical window and is otherwise substantially entangled with unmodelled spatial structure.

## H3 intermittent, no consistent direction

The within-rectangle effect flips sign and significance across every phase (non-significant closing → significant closing → non-significant widening → significant closing), with no monotonic or policy-consistent pattern. Unlike H2, these estimates have not been tested against a spatial-confounding null, so even the phases reported as "significant" carry an open question H2's checks would raise if applied here.

**Read:** H3 does not show the within-cell temporal signature. This is the hypothesis with the weakest support of the three. A temporal lag mechanism would have been an improvement — to mention in the discussion and not built given the time.

---

Explanatory power

Fishing pressure and policy phase together explain only ~5% of haul-level prediction-error variance (marginal R² = 0.051), a ceiling that held across every spatial specification tested. **~95% of why EEoS gets it wrong is unexplained by disturbance intensity as currently measured.**

## Why the hypotheses fail — synthesis

1. **A confound sits underneath all three hypotheses**: EEoS overprediction scales with biomass and shoal homogeneity (H1), in a way that looks like catchability/sampling bias as much as ecological disturbance. Until this is isolated, H2 and H3's residuals are testing a partly-contaminated signal.
2. **Geography is doing more work than fishing pressure**: H2's headline result depends heavily on spatial structure the model doesn't fully absorb; the strongest "disturbance" period is explainable by effort displacement between neighbours, not disturbance intensity per se.
3. **Low variance explained (5%) suggests fishing pressure/phase are the wrong, or at least insufficient, predictors**, not that the underlying EEoS-residual signal is absent. The H1 diagnostics (dominance/size-CV) (potentially) point to where real signal may be hiding.

## Avenues for further research

- **Catchability correction**: test whether a size- or biomass-scaled correction collapses the H1 overprediction gradient, isolating artifact from disturbance.
- **Better disturbance proxy**: SAR-derived effort data in place of Couce fishing-hours, which may not capture disturbance intensity well.
- **Extend spatial-confounding checks to H3**: bootstrap/neighbour-lag tests were only run on H2's between-rectangle term.
- **DynaMETE disturbance-type discrimination**: test whether residual *shape*, not just magnitude, distinguishes disturbance types, rather than relying on fishing pressure as a single intensity axis.
- **Dominance/size-CV as a direct covariate**: H1's second, independent overprediction gradient (shoal homogeneity) was never incorporated into the H2/H3 model, a natural next step given it persists after conditioning on biomass.
