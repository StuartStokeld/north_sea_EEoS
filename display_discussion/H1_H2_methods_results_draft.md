# Methods and results — draft (Rule 7 structure)

**Central claim (working):** EEoS does not predict absolute haul-level biomass in the North Sea trawl survey, but haul-specific residuals carry structured information about community configuration; when aggregated to rectangles, that residual magnitude is **not** robustly associated with fishing pressure after spatial correction — biomass scale, not disturbance, is the dominant cross-rectangle pattern in |residual|.

**Survey:** NS-IBTS Q1, 1985–2015 · **Primary units:** haul (H1), ICES statistical rectangle (H2)

---

# Methods

## Overview

We tested whether the Ecological Equation of State (EEoS; Harte et al. 2022) predicts community biomass at individual trawl haul resolution (H1), then whether mean absolute EEoS residual magnitude at rectangle scale covaries with reconstructed fishing effort (H2). State variables (species richness *S*, raised abundance *N*, normalised metabolic rate *E*) were derived from ICES DATRAS length-frequency data; observed biomass *B*<sub>obs</sub> came from FishGlob — deliberately separate sources so the predictor–outcome comparison is not circular. EEoS predictions *B*<sub>pred</sub> used the authors' reference Python solver (`biomass.py`) via R/reticulate, rescaled to grams with the same minimum-individual mass *m*<sub>min</sub> used to normalise *E*. Length–weight relationships (FishBase, genus/family fallback) were applied per 1 cm bin, not to haul-mean length. Hauls failing partition-function requirements (*S* < 2, *E* ≤ *N*, non-finite inputs, *B*<sub>obs</sub> ≤ 0) were excluded with logged reasons (*N* = 12,069 hauls retained).

The headline H1 metric is log-scale coefficient of determination (log_r² = 1 − SS<sub>res</sub>/SS<sub>tot</sub> on log *B*<sub>obs</sub> vs log *B*<sub>pred</sub>), which penalises systematic scale error; squared Pearson correlation (cor²) is reported as a diagnostic only. H1 tests replicate Harte et al.'s unfitted validations before any fitted alternative. H2 uses mean |log *B*<sub>obs</sub> − log *B*<sub>pred</sub>| per rectangle (pooled 1985–2015 hauls) as the dependent variable and mean annual reconstructed otter + beam trawling hours (Couce et al. 2020) as the independent variable (*N* = 158 rectangles with ≥10 hauls and Couce coverage). Spatial structure was modelled with queen-contiguity weights on ICES rectangles (row-standardised; spatial error model as pre-registered primary inference).

---

## H1: haul-level prediction and residual characterisation

### Unfitted baseline tests

**Test 1 (Harte Fig 1):** log *B*<sub>obs</sub> vs log *B*<sub>pred</sub> on a strict unfitted 1:1 line, compared to the same unfitted 1:1 test for a productivity baseline log(*E* × *m*<sub>min</sub>) using identical normalised *E* and *m*<sub>min</sub>.

**Test 2 (Harte Fig 2):** predicted vs observed productivity ratio *E*/*B*^0.75 via unweighted linear regression on raw ratios (Pearson R² per Harte source code), with a leverage trim excluding the top two hauls by predicted×observed ratio magnitude (trim applies to this regression only).

**Extended comparison (separate):** fitted OLS log *B*<sub>obs</sub> ~ log *E* — reported after unfitted tests; not treated as prior-method validation.

### Null model and pipeline integrity

*B*<sub>pred</sub> was held fixed; *B*<sub>obs</sub> was randomised over 999 permutations under two schemes (uniform draw on central 95% of log *B*<sub>obs</sub>; shuffle of observed *B*<sub>obs</sub> across hauls), with scheme choice driven by distribution shape diagnostics. A year-by-rectangle processing funnel tracked sample attrition. Systematic overprediction (median *B*<sub>pred</sub>/*B*<sub>obs</sub> ≈ 4×) was examined across *B*<sub>obs</sub> quartiles; no scalar catchability correction was applied.

### Compositional typicality (dominance and size homogeneity)

To separate disturbance-driven from sampling-driven deviation from the METE well-mixed assumption, each haul received Berger–Parker numerical dominance *D* (share of *N* in the single most abundant species) and dominant-species size homogeneity (size_CV: raised-count-weighted CV of bin-level mass within that species), computed from the same HL bins as *N* and *E*. Before interpreting any *D*/size_CV–residual association, four pre-specified confound checks were run: correlation matrix (flag |*r*| > 0.4); association within each *B*<sub>obs</sub> quartile; cross-tab of low-size_CV decile vs known dropout-flagged years/rectangles; descriptive taxonomic breakdown in metric extremes. Partial R² of {*D*, size_CV} beyond log *B*<sub>obs</sub> was estimated with nested OLS (including *D* × size_CV interaction test). No dominance-based filter was applied to the primary H1 metrics; this track is diagnostic and sets design constraints for H2/H3.

---

## H2: rectangle-level residuals vs fishing pressure

### Primary models

Baseline OLS: mean_abs_residual ~ mean_annual_hours_total. Moran's *I* and Geary's *C* on OLS residuals assessed spatial autocorrelation. Primary inference used a spatial error model (`errorsarlm`, same weights, `zero.policy = TRUE`).

### Robustness and alternative specifications

Models were re-run on panels with ≥5 and ≥20 hauls per rectangle (*N* = 161, 156). Cook's distance and leverage identified high-influence rectangles (threshold 4/*N*); models were refit excluding flagged units. A biomass-covariate specification added mean log haul biomass per rectangle (mean of log *B*<sub>obs</sub>, not log of the mean). Lagrange Multiplier tests (classic and robust) compared spatial error vs spatial lag specifications on OLS residuals (primary and biomass-covariate models). Spatial lag models (`lagsarlm`) were fit for both specifications and compared to SEM by coefficient, AIC, and log-likelihood.

### Limitations

Couce effort does not cover Skagerrak/Kattegat and eastern English Channel approaches (~23 rectangles excluded for missing effort; additional exclusions for <10 hauls). Effort is in hours, not swept-area ratio.

---

# Results

*Each subsection states a conclusion in its title; paragraphs follow question → evidence → answer.*

---

## H1 results

### The analysis pipeline yields an auditable haul set suitable for unfitted EEoS testing

**Figure 0 (conceptual funnel, optional):** *Sample funnel from raw DATRAS HL to 12,069 EEoS-passing hauls with independent FishGlob biomass.*

After joining DATRAS length data to FishGlob catch records and applying structural filters, **12,069 NS-IBTS Q1 hauls (1985–2015)** entered analysis. Predictors (*S*, *N*, *E*) and outcome (*B*<sub>obs</sub>) were built from separate aggregations; exclusions were logged by reason. This establishes that the H1 test is not circular and that sample loss is traceable rather than opaque.

---

### EEoS fails the unfitted absolute-biomass test that prior validations used

**Figure 1:** *EEoS does not predict absolute haul biomass on an unfitted 1:1 map (log_r² = −0.22), while unfitted productivity 1:1 outperforms it (log_r² = +0.74).*

On the primary unfitted comparison — log *B*<sub>obs</sub> vs log *B*<sub>pred</sub>, no fitted intercept or slope — EEoS achieved **log_r² = −0.22** (worse than predicting the mean). Under the identical unfitted footing, productivity 1:1 (log[*E* × *m*<sub>min</sub>] vs log *B*<sub>obs</sub>) achieved **log_r² = +0.74**. EEoS left **4.6×** more unexplained variance than the productivity baseline (Harte criterion not met). The EEoS cloud sits systematically above the 1:1 line (median *B*<sub>pred</sub>/*B*<sub>obs</sub> = **4.1×**); cor² = **0.93** shows tight rank association but does not indicate accurate absolute prediction. **Conclusion:** at North Sea haul scale, EEoS does not reproduce the unfitted biomass–productivity equivalence reported in prior EEoS validations.

| Model | Fitted | log_r² | cor² | Median pred/obs |
|-------|--------|-------:|-----:|----------------:|
| EEoS (*S*, *N*, *E*) | No | −0.22 | 0.93 | 4.1× |
| Productivity 1:1 (*E* × *m*<sub>min</sub>) | No | +0.74 | 0.90 | 0.67× |

*Legend note for Figure 1:* Three panels — left: Harte Fig 1 reproduction (EEoS); centre: productivity 1:1 baseline; right: log_r² comparison bar chart. All points are individual hauls; no post-hoc calibration.

---

### High cor² masks a systematic scale offset that grows with observed biomass

**Figure 2 (optional / supplementary):** *Overprediction ratio increases with observed biomass quartile (3.1× to 5.5×), ruling out a single catchability multiplier.*

Median *B*<sub>pred</sub>/*B*<sub>obs</sub> rose from **3.1×** (lowest *B*<sub>obs</sub> quartile) to **5.5×** (highest). A uniform scalar correction would therefore misrepresent the error structure and inflate apparent fit without mechanistic justification. **Conclusion:** H2/H3 should use residual **magnitude**, not bias-corrected absolute biomass, as the condition-relevant quantity.

---

### EEoS residuals are structured beyond chance, satisfying a minimal prerequisite for downstream hypothesis tests

**Figure 3 (optional):** *Observed log_r² lies far above both null permutation distributions (p < 0.001).*

Under 999 permutations (uniform 95% log *B*<sub>obs</sub> draw and *B*<sub>obs</sub> shuffle), observed log_r² = −0.22 exceeded both null medians (−2.02 and −2.60; **p < 0.001** each). **Conclusion:** haul-specific (*S*, *N*, *E*) → *B*<sub>pred</sub> mapping carries information not explained by randomising *B*<sub>obs</sub> alone — residuals are structured enough to support H2/H3, even though absolute prediction fails.

---

### Numerical dominance and dominant-species size homogeneity explain modest residual variation beyond biomass magnitude

**Figure 4:** *EEoS overprediction rises with Berger–Parker dominance and falls with dominant-species size heterogeneity, independently of biomass quartile.*

Hauls high in *D* and low in size_CV (numerically dominated, near-uniform dominant-species catch — the shoal/cohort configuration METE treats as atypical) showed larger median *B*<sub>pred</sub>/*B*<sub>obs</sub> than less dominated, more size-heterogeneous hauls (*D* quartiles: 3.70× → 4.41×; size_CV quartiles: 4.46× → 3.73×). Confound checks supported interpreting this as a second axis, not a re-labelling of the biomass-magnitude bias: |*r*| between *D*/size_CV and *B*<sub>obs</sub>/*N* ≤ 0.26; associations persisted within every *B*<sub>obs</sub> quartile; low-size_CV deciles were not enriched in dropout-flagged years/rectangles. Partial R² of {*D*, size_CV} beyond log *B*<sub>obs</sub> ≈ **1%**; no *D* × size_CV interaction (ΔR² < 0.0001, *p* = 0.79). **Conclusion:** compositional atypicality contributes to EEoS failure in the direction METE's well-mixed assumption predicts, but explains far less variance than biomass scale; dominance and size homogeneity are tracked as separate diagnostic axes, not combined into a single correction applied to H1 primary metrics.

![Ratio vs D](../outputs/figures/h1_dominance_ratio_vs_D.png)

*Figure 4a legend:* log(*B*<sub>pred</sub>/*B*<sub>obs</sub>) vs Berger–Parker dominance *D*; loess smooth; *N* = 12,069 hauls.

![Ratio vs size_CV](../outputs/figures/h1_dominance_ratio_vs_sizeCV.png)

*Figure 4b legend:* log(*B*<sub>pred</sub>/*B*<sub>obs</sub>) vs dominant-species size_CV; loess smooth.

---

### Fitted ln(*E*) OLS outperforms EEoS but answers a different question

Fitted log *B*<sub>obs</sub> ~ log *E* yielded log_r² = **+0.60** vs EEoS **−0.22** (SS ratio 3.0×). **Conclusion:** allowing the model to absorb systematic offset improves fit, but this is not equivalent to unfitted mechanistic validation; it is reported only as extended context after Tests 1–2.

---

## H2 results

### Rectangle-level aggregation links H1 residuals to a spatially explicit disturbance test

**Figure 5 (conceptual / map panel, optional):** *Mean |EEoS residual| and mean annual fishing hours across 158 ICES rectangles, 1985–2015.*

H1 established structured haul-level residuals with scale-dependent bias. H2 aggregates mean |log *B*<sub>obs</sub> − log *B*<sub>pred</sub>| per rectangle and regresses it on mean annual reconstructed trawling hours (Couce et al. 2020), testing whether heavier fishing associates with **larger** residual magnitude (pre-registered positive β). Default panel: **158 rectangles** (≥10 hauls, Couce coverage).

---

### Naive OLS suggests higher fishing pressure associates with *smaller* residuals — opposite to prediction

**Figure 6:** *The significant negative OLS association between fishing hours and mean |EEoS residual| does not survive spatial error correction.*

Primary OLS: β = **−2.74 × 10⁻⁶** h⁻¹ (SE = 1.06 × 10⁻⁶, *t* = −2.59, **p = 0.010**, R² = 0.041, *N* = 158) — opposite the pre-registered direction. **Conclusion:** a naive cross-sectional reading does not support the disturbance hypothesis; if anything, it points the wrong way.

![H2 topline result](../outputs/figures/h2_topline_result.png)

*Figure 6 legend:* Panel A — scatter of mean annual hours vs mean |residual| with OLS trend (Moran's *I* on OLS residuals = 0.56, *p* ≪ 0.001). Panel B — fishing-pressure coefficient per 10,000 annual hours: OLS significant and negative; SEM not significant (coefficients scaled for readability; * p < 0.05).

---

### OLS residuals are strongly spatially autocorrelated, invalidating uncorrected inference

Moran's *I* = **0.56**, Geary's *C* = **0.43** (both *p* ≪ 0.001 on OLS residuals). **Conclusion:** standard OLS standard errors and *p*-values for the fishing-pressure term cannot be trusted without a spatial model.

---

### After spatial error correction, fishing pressure is not associated with residual magnitude

SEM: β = **−6.1 × 10⁻⁷** (*p* = **0.541**, λ = 0.858). **Conclusion:** the OLS association is not robust to spatial structure; pre-registered success criterion (positive significant β in SEM) is **not met**.

---

### The null result is not an artefact of panel definition, leverage, or influential rectangles

Re-running SEM at ≥5 and ≥20 haul thresholds (*N* = 161, 156) gave SEM β from −6.1 × 10⁻⁷ to +1.6 × 10⁻⁷ (all *p* > 0.5); OLS significance across thresholds ( *p* = 0.009–0.017) reflected unmodelled spatial structure, not a stable signal. Excluding seven high-Cook's rectangles **strengthened** the negative OLS association (β = −3.68 × 10⁻⁶, *p* = 0.00096, *N* = 151). **Conclusion:** neither sample threshold nor leverage explains away the spatially corrected null.

---

### Biomass-adjusted models show fishing pressure only under spatial misspecification; biomass scale dominates

Adding mean log *B*<sub>obs</sub> flipped the OLS hours coefficient positive (*p* ≈ 0.019; R² = 0.55) — mirroring H1 catchability confounding — but SEM kept hours non-significant (β = −1.1 × 10⁻⁷, *p* = 0.861) while mean log biomass was strongly significant (β = 0.187, *p* ≪ 0.001). LM tests favoured spatial lag for the hours-only model (robust LM-lag *p* = 0.039) but spatial error for the biomass model (robust LM-error *p* ≪ 0.001; ΔAIC = 54.8 favouring SEM). SAR on the biomass model reproduced a positive significant hours effect (β = +1.91 × 10⁻⁶, *p* = 0.0016) with substantially worse fit. **Conclusion:** apparent positive fishing-pressure effects under OLS or misspecified SAR are interpreted as spatial/biomass confounding artefacts, not evidence for disturbance-driven residual magnitude.

---

### Fishing pressure does not predict EEoS residual magnitude at rectangle scale; biomass scale does

Across specifications, **no robust positive association** between fishing hours and mean |EEoS residual| was found. **Mean absolute residual magnitude scaled positively with mean observed biomass** independent of fishing pressure — consistent with H1's scale-dependent bias as the dominant cross-rectangle driver of residual structure. **Conclusion:** H2 does not support the hypothesis that heavier fishing is associated with larger EEoS deviation magnitude; the more consistent community-level pattern is biomass-related prediction error, not fishing disturbance per se.

---

# Synthesis table (logical progression)

| Step | Question | Answer | Key figure |
|------|----------|--------|------------|
| 1 | Does EEoS predict absolute haul biomass unfitted? | **No** (log_r² = −0.22) | Fig 1 |
| 2 | Is failure a fixed scale offset? | **No** — grows with *B*<sub>obs</sub> | Fig 2 |
| 3 | Are residuals structured? | **Yes** (null *p* < 0.001) | Fig 3 |
| 4 | Does compositional atypicality add explanation? | **Modestly**, beyond biomass | Fig 4 |
| 5 | Does fishing pressure predict \|residual\| (OLS)? | **Significant negative** (wrong sign) | Fig 6A |
| 6 | Is OLS trustworthy spatially? | **No** (Moran's *I* = 0.56) | Fig 6A |
| 7 | Does SEM support disturbance signal? | **No** (*p* = 0.54) | Fig 6B |
| 8 | Robust to threshold/leverage/specification? | **Yes — null throughout** | — |
| 9 | What *does* predict \|residual\|? | **Mean log biomass** | — |

---

# Figure list for submission (conclusion-forward titles)

1. **EEoS underpredicts observed biomass on an unfitted 1:1 map less well than unfitted productivity 1:1** — Harte Fig 1 reproduction + productivity baseline + log_r² bars (`outputs/figures/harte_fig*.png`, productivity 1:1).
2. **Overprediction ratio increases with observed biomass** — quartile bar chart or scatter (from catchability exploration outputs).
3. **Haul-level EEoS residuals exceed null permutation expectations** — null density + observed log_r² (`outputs/figures/null_r2_*.png`).
4. **Dominance and size homogeneity associate with EEoS overprediction beyond biomass magnitude** — `h1_dominance_ratio_vs_D.png`, `h1_dominance_ratio_vs_sizeCV.png`.
5. **Fishing pressure is not associated with mean |EEoS residual| after spatial correction** — `h2_topline_result.png`.

---

*Sources: [`H1_methods_draft.md`](H1_methods_draft.md), [`H2_methods_draft.md`](H2_methods_draft.md), [`H1_results_summary.md`](H1_results_summary.md), [`H1_dominance_results_summary.md`](H1_dominance_results_summary.md), [`H2_results_draft.md`](H2_results_draft.md), [`Methods_note_for_discussion.md`](Methods_note_for_discussion.md).*
