# H2/H3 results — interpretation note (biomass-free primary model)

> **SUPERSEDED (Aug 2026).** This note describes the earlier **blended** model
> `residual ~ log(hours+1) × phase + (1 | stat_rec)` with data-driven phases
> (1989 / 2001 / 2008). It is **not** the presented primary.
>
> **Current primary:** within-between decomposition with policy-anchored `phase_v2`
> (1992 / 2002 / 2008) — see
> [`H2_results_draft.md`](H2_results_draft.md),
> [`H3_results_draft.md`](H3_results_draft.md),
> [`One page read me.md`](One%20page%20read%20me.md), and
> `outputs/phase_v2_*` artifacts from `pipeline/run_h2h3_phase_v2_reporting.R`.
>
> Kept below as a historical plain-language reading of the blended-term run only.
> Sign-convention material in §0 remains useful; phase-specific slopes and R² in
> later sections do not match the current primary.

**Audience:** Stuart + supervisors  
**Source run:** `pipeline/run_h2h3_shared_model_results.R` (no-biomass re-run) — superseded  
**Primary model (this note only):** `residual ~ log(hours+1) × phase + (1 | stat_rec)`  
**Data:** 10,464 hauls, 158 ICES rectangles, 1985–2015  

This note translates the coefficient tables into plain language. It does **not** offer ecological mechanism.

---

## 0. Sign convention (checked against this run’s data)

### Definition used

The outcome is the pipeline’s canonical residual:

\[
\text{residual} = \log(B_{\text{obs}}) - \log(B_{\text{pred}}) = \log\!\left(\frac{B_{\text{obs}}}{B_{\text{pred}}}\right)
\]

- Residual **= 0** → observed biomass matches EEOS prediction.  
- Residual **< 0** → observed biomass is **below** prediction (theory **overpredicts**).  
- Residual **> 0** → observed biomass is **above** prediction (theory **underpredicts**).

### What the data actually look like in this run

| Check | Value in this run |
|-------|-------------------|
| Share of hauls with residual < 0 | **99.0%** |
| Median residual | **−1.39** |
| Median \(B_{\text{obs}}/B_{\text{pred}}\) | \(\exp(-1.39) \approx\) **0.25** |

So the typical haul has observed biomass about **one-quarter** of what EEOS predicts — a large, systematic overprediction gap. Interpretation of fishing-pressure slopes should be read against that background, not against a residual centred at zero.

### What a *positive* fishing-pressure coefficient means

The model’s fishing-pressure slope is:

\[
\frac{\partial\,\widehat{\text{residual}}}{\partial\,\log(\text{hours}+1)}
\]

A **positive** slope means: hauls (or rectangle-years) with **higher** fishing hours are predicted to have **higher** residual — i.e. less negative, **closer to zero**.

In plain terms, given that residual is almost always negative:

> **A positive fishing-pressure coefficient means the overprediction gap shrinks as fishing pressure rises** (observed biomass is a larger fraction of predicted biomass where fishing hours are higher).

A **negative** slope would mean the opposite: higher fishing hours → more negative residual → **larger** overprediction gap.

### Worked example (using the 2001–2007 primary slope)

Primary slope in 2001–2007: **+0.0885**.

Take two otherwise comparable hauls that differ by **one unit** of \(\log(\text{hours}+1)\). The model predicts their residuals differ by:

\[
0.0885 \times 1 = 0.0885
\]

That is a shift of the log ratio \(B_{\text{obs}}/B_{\text{pred}}\) by +0.0885, i.e. the biomass ratio is multiplied by \(\exp(0.0885) \approx 1.093\) (**about +9%**).

More tangibly, a **doubling** of fishing hours is a change of \(\log 2 \approx 0.693\) on the predictor scale:

\[
\Delta\text{residual} = 0.0885 \times 0.693 \approx 0.061
\]

At this run’s **median** residual (−1.39):

| | Residual | \(B_{\text{obs}}/B_{\text{pred}}\) | Distance from zero (\(|\text{residual}|\)) |
|--|----------|--------------------------------------|-------------------------------------------|
| Typical haul | −1.392 | 0.249 | 1.392 |
| Same haul after +1 doubling of hours (model prediction) | −1.331 | 0.264 | 1.331 |

The residual moves **toward zero** by 0.061 log units; the observed/predicted ratio rises from ~25% to ~26% of predicted. That is a **shrinkage of the overprediction gap**, not a move into underprediction.

*(Same logic for 2008–2015 slope +0.058: doubling hours → Δ residual ≈ +0.040 → ratio ×1.041, about +4%.)*

---

## 1. H2 — Does fishing pressure predict residual?

**Short answer:** Yes in some phases, no in others. There is **no** single fishing-pressure effect that holds across 1985–2015. The finding is phase-specific.

Primary-model fishing-pressure slopes (change in residual per unit \(\log(\text{hours}+1)\)):

| Phase | Slope | 95% CI | Differs from zero? | Plain-language direction |
|-------|------:|--------|--------------------|--------------------------|
| 1985–1988 | −0.001 | [−0.028, +0.027] | No (p ≈ 0.96) | No detectable link between fishing hours and the overprediction gap |
| 1989–2000 | −0.006 | [−0.023, +0.012] | No (p ≈ 0.54) | Same — gap does not reliably change with fishing hours |
| 2001–2007 | **+0.089** | [+0.073, +0.104] | **Yes** (p ≈ 2×10⁻²⁹) | Higher fishing hours → **smaller** overprediction gap |
| 2008–2015 | **+0.058** | [+0.044, +0.072] | **Yes** (p ≈ 3×10⁻¹⁶) | Same direction, somewhat **weaker** than 2001–2007 |

Joint test that the slope is zero in *every* phase: χ²(4) = 201, p ≈ 3×10⁻⁴² — driven by the two post-2001 phases, not by a uniform effect.

**Do not collapse this into one headline slope.** Pre-2001: no detectable FP–residual association. Post-2001: positive association (overprediction gap shrinks where fishing pressure is higher).

---

## 2. H3 — Does that relationship change over time?

**Short answer:** Yes. Relative to 1985–1988, the fishing-pressure slope is **unchanged** through 1989–2000, then **rises** in 2001–2007 and remains higher in 2008–2015.

H3 is answered by the FP × phase **interaction** terms (change in slope vs the 1985–1988 reference):

| Contrast vs 1985–1988 | Interaction | p | Plain reading |
|-----------------------|------------:|---|---------------|
| 1989–2000 | −0.005 | 0.74 | No detectable change in the FP–residual link |
| 2001–2007 | **+0.089** | 1×10⁻⁹ | FP–residual link becomes clearly more positive |
| 2008–2015 | **+0.059** | 5×10⁻⁵ | Still more positive than 1985–1988 |

Joint test that all interactions are zero: χ²(3) = 113, p ≈ 2×10⁻²⁴.

So H3’s phase pattern matches H2’s: the relationship is flat/null early, then positive after 2001. (No mechanism claimed here.)

---

## 3. Effect size in tangible units

A full four-phase table (percent ratio change, percentage-point ratio shift, and percent of remaining gap closed — each under both a doubling of hours and that phase’s own fishing-hours IQR, with CIs) is in `outputs/h2h3_results_proportional_effects.csv`, with a gap-closure companion figure at `outputs/figures/h2h3_results_gap_closed_by_phase.png`. Worked examples for the two significant phases follow.

Predictor is \(\log(\text{hours}+1)\). Useful conversions:

| Change in fishing hours | Change on predictor | Residual shift (2001–2007, slope 0.089) | Residual shift (2008–2015, slope 0.058) |
|-------------------------|--------------------:|----------------------------------------:|----------------------------------------:|
| **Doubling** | +0.693 | **+0.061** | **+0.040** |
| ≈ IQR of \(\log(\text{hours}+1)\) in this data (~1.89) | +1.89 | +0.167 | +0.110 |

**In words (post-2001 only):**

- **2001–2007:** a doubling of Couce fishing hours is associated with residual higher by ~0.06 — about a **6%** higher \(B_{\text{obs}}/B_{\text{pred}}\) ratio, i.e. the overprediction gap shrinks modestly.  
- **2008–2015:** same doubling → residual higher by ~0.04 — about a **4%** higher ratio.

These are small moves relative to the typical gap (median residual −1.39). They are **detectable**, not **large** relative to how far EEOS already sits from the data.

---

## 4. R² — statistical robustness ≠ dominant explanation

| Quantity | Value | Plain meaning |
|----------|------:|---------------|
| **Marginal R²** | **0.047** | Fixed effects alone (fishing pressure × phase) explain about **5%** of residual variance |
| **Conditional R²** | **0.165** | Fixed effects **plus** rectangle random intercepts explain about **17%** |

Most of the lift from 0.047 → 0.165 is rectangle-level pooling, not fishing pressure.

**State both of these together:**

- With ~10k hauls, the post-2001 slopes and interactions have **very small p-values** — the associations are estimated precisely and are unlikely to be sampling noise.  
- At the same time, **marginal R² ≈ 0.05** means fishing pressure × phase is **not** the main story in haul-level residual variation. Most residual scatter remains unexplained by these fixed effects.

Significance says “reliably non-zero in this sample.” Marginal R² says “still a small share of total variation.” Both are true; neither replaces the other.

---

## 5. Sensitivity — result does not hinge on model choice

Phase-specific fishing-pressure slopes agree closely across the three settled specifications:

| Comparison | Correlation of the four phase slopes |
|------------|--------------------------------------:|
| Primary vs CAR (spatial sensitivity) | **1.000** |
| Primary vs GAM (continuous-year sensitivity) | **0.988** |

In plain terms: whether you use a plain rectangle intercept, a CAR spatial structure, or a smooth year-varying coefficient, you get the **same phase pattern** — near-zero early, positive after 2001. The H2/H3 reading above is not an artefact of one random-effect or time structure.

---

## 6. Primary figure

![Fishing-pressure effect by phase](../outputs/figures/h2h3_results_fp_effect_by_phase.png)

**Figure.** Estimated fishing-pressure effect on residual by year. Red steps = primary model (constant slope within each phase, 95% CI bands). Green = GAM sensitivity (smooth over year). Dotted lines = structural breaks 1989, 2001, 2008.

**How to read the vertical axis in plain terms:** values above zero mean higher fishing hours → **higher** residual → **smaller** overprediction gap (see §0). Values near zero with CI crossing zero mean no detectable association. The post-2001 red steps sit clearly above zero; the pre-2001 steps do not. The green smooth tracks that same broad pattern.

---

## Appendix — biomass removed (re-run note)

An earlier results draft included `mean_ln_B_obs` as a covariate. That was incorrect for the intended design. The numbers in this note are from the **biomass-free** re-run.

Relative to the with-biomass draft (full detail in `outputs/h2h3_results_run_log.md`, “Before/after” section):

- Phase slopes shifted by ~0.01–0.02 (early slopes moved toward zero; post-2001 slopes slightly larger).  
- H3 interactions barely changed.  
- Marginal R² fell from 0.104 → **0.047** (expected after dropping biomass); conditional R² stayed ~0.16.  
- Primary–CAR / primary–GAM correlations remained ~1.00 / ~0.99.

The qualitative H2/H3 pattern (null early, positive post-2001) is the same with or without biomass; the biomass-free run is the one to cite.
