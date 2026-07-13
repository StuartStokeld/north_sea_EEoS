# H1 metrics review — statistical best practice and guidance for H2/H3

**Purpose:** Rank the H1 statistics by ecological and statistical utility, distinguish metrics appropriate for *predictive validation* (H1) from those appropriate for *residual inference* (H2/H3), and recommend a minimal reporting set going forward.

**Context:** H1 tested whether EEoS predicts haul-level biomass from (S, N, E). H2 asks whether EEoS *residuals* covary with fishing pressure at StatRec × decade scale. H3 asks whether residual dynamics over time at StatRec scale reflect disturbance. These are different statistical questions requiring different metrics.

---

## 1. What H1 established (and what it did not)

| Finding | Implication for H2/H3 |
|---------|----------------------|
| log_r2 = −0.22 at haul level | EEoS is a poor absolute predictor at fine scale — **do not re-test this at H2/H3** |
| cor² = 0.93 | Rank agreement is high but scale is wrong — **misleading if reused as success metric** |
| Null model p < 0.001 | Residuals are structured, not random — **necessary precondition for H2/H3** |
| Median B_pred/B_obs = 4×, increasing with biomass quartile | Systematic calibration error — **must be controlled or modelled in H2/H3** |
| Productivity 1:1 log_r2 = +0.74 > EEoS −0.22 | Simple E×m_min map beats EEoS at haul scale — **context for interpreting residual signal** |

H1 failure at haul level does not invalidate H2/H3. Aggregating hauls to StatRec × decade reduces stochastic noise; the question shifts from "does EEoS predict biomass?" to "do deviations from EEoS track disturbance?"

---

## 2. Metric ranking for H2/H3

### Tier A — Primary (use in H2/H3 analysis and reporting)

| Metric | H1 value (EEoS) | Role in H2/H3 | Ecology rationale |
|--------|----------------:|---------------|-------------------|
| **Mean absolute log residual** `mean\|log(B_obs) − log(B_pred)\|` | mean 1.328 | **Primary dependent variable (H2/H3)** | EEoS is a *static equilibrium* theory: under disturbance a community can deviate from the equilibrium prediction in either direction, so the **magnitude** of deviation — not its sign — is the quantity of interest. Signed residuals would net out directionally meaningful deviation when aggregated. This is the established design principle and the DV implemented in the H2 pipeline (`mean_abs_residual ~ fishing pressure`). |
| **Signed log residual** `log(B_obs) − log(B_pred)` | median −1.41 | **Diagnostic / calibration only** | Log-scale residuals are additive and appropriate for right-skewed biomass (Burnham & Anderson 2002; standard in fisheries CPUE modelling). Retained to characterise the systematic overprediction offset (feeds the catchability framing), **not** as the H2/H3 response — it would cancel opposing deviations on aggregation. |
| **Median B_pred/B_obs (or pred/obs ratio)** | 4.1× overall; 3.1–5.5× by quartile | **Calibration diagnostic; covariate or stratifier** | Multiplicative bias is interpretable for catchability and survey selectivity. If fishing pressure correlates with haul biomass quartile, the quartile-dependent offset (§6 of H1 summary) can generate spurious spatial residual patterns — test and adjust. |
| **Null-model significance (p < 0.001)** | Beat both schemes | **Precondition, not repeated test** | Establishes that haul-specific (S,N,E) → B pairing carries information. H2/H3 inherit this: residuals are not exchangeable noise. Report once; do not re-permute at aggregate scale without new justification. |
| **Sample size per stratum (n hauls per StatRec × decade)** | 12,069 total | **Design and inference** | Ecological inference at aggregate scale requires adequate replication. Report n explicitly; down-weight or exclude sparse cells (standard practice in meta-analysis and mixed models). |

### Tier B — Supporting (report in methods/supplement; use for diagnostics)

| Metric | Role | Notes |
|--------|------|-------|
| **log RMSE** | Compare model skill if testing alternative residual definitions | Useful for comparing EEoS vs productivity 1:1 at a given scale. Less interpretable than signed median residual for disturbance inference. |
| **Median \|log residual\|** | Dispersion of errors | Complements signed median. High |residual| with stable signed median suggests heteroscedasticity — relevant for mixed-model specification. |
| **SS_res ratio vs productivity 1:1** | Theoretical benchmark | Harte criterion is for *predictive* comparison, not residual–pressure association. Report in H1/H2 methods as context; do not use as H2 test statistic. |
| **Productivity ratio (E/B^3/4)** | Theoretical scaling check | Harte Fig 2 replication failed at haul scale. Could be recomputed at StatRec aggregate if theory predicts scaling at community level — but not the primary H2/H3 metric. |
| **Slice log_r2 within strata** | **Avoid as primary** | Unstable when B_obs range is narrow (Q2/Q3 artifacts in H1). If reported at StatRec level, use only where SS_tot is sufficient. |

### Tier C — Retire for H2/H3 (haul-level only; do not carry forward as test statistics)

| Metric | Why retire |
|--------|------------|
| **Global haul-level cor² = 0.93** | Measures rank covariation driven largely by E–B association. High cor² coexists with failed prediction. Using cor² at aggregate scale would repeat the same confounding. |
| **Global haul-level log_r2** | Answers "does EEoS predict?" — already answered (no). Recomputing at StatRec level changes the question to local predictive skill, not disturbance detection. |
| **ln(E) OLS R² = 0.60** | Fitted benchmark; unfair comparison with unfitted EEoS. Irrelevant once analysis shifts to residuals. |
| **Raw RMSE (grams)** | Dominated by large hauls; poor ecological interpretability for community condition. |
| **Harte criterion vs ln(E) OLS** | Tier 2 comparison; OLS minimises SS_res by construction. Not informative for residual–pressure inference. |

---

## 3. Statistical best practice — ecology and fisheries context

### 3.1 Predictive validation vs residual inference

H1 is **predictive validation**: fixed theory, no parameters fitted to B_obs, evaluate against observations (Roberts et al. 2017; Dormann et al. 2018 — distinction between validation and inference).

H2/H3 is **residual inference**: treat EEoS as a null model for expected biomass given (S, N, E); test whether deviations covary with fishing pressure or time. This follows the **calibration–validation split** common in species distribution and ecosystem models: a model may fail calibration (absolute scale) yet retain useful structure in residuals if errors are structured and understood.

**Recommendation:** Frame H2/H3 explicitly as residual analysis against an EEoS null, not as a second prediction contest.

### 3.2 Why log scale

Biomass spans orders of magnitude. Log transformation:

- Stabilises variance (common in fisheries CPUE and biomass surveys)
- Makes multiplicative errors additive (4× overprediction = constant log offset)
- Aligns with EEoS theory, which operates in log-metabolic space

**Recommendation:** Primary H2/H3 response = **mean absolute residual** `mean|log(B_obs) − log(B_pred)|` aggregated by StatRec (and decade for H2), because EEoS is a static equilibrium null and disturbance can push a community either side of it (see Tier A). Report the signed residual and `B_obs/B_pred` alongside for calibration/communication, but do not use signed residual as the aggregated response — opposing deviations would cancel.

### 3.3 R² is not enough (and cor² is worse)

Ecological model evaluation best practice (Araujo et al. 2005; Pineiro et al. 2008; Thiele et al. 2022):

- **R² alone** does not assess calibration (systematic bias vs random scatter)
- **cor²** ignores intercept and scale — explicitly misleading when EEoS overpredicts uniformly
- Report **multiple metrics**: bias (median residual), precision (IQR or RMSE of residuals), and **graphical calibration** (residual vs fitted, residual vs pressure)

H1 already demonstrates why: cor² = 0.93 with log_r2 = −0.22.

**Recommendation for H2/H3:** Report effect sizes (slope of residual ~ pressure, or median residual difference between pressure groups), confidence intervals, and **spatial/temporal residual maps** — not R² of EEoS at aggregate scale.

### 3.4 Aggregation and scale

Haul-level noise includes:

- Stochastic trawl catchability
- Local patchiness
- Measurement error in S, N, E

StatRec × decade aggregation (H2) or StatRec time series (H3) follows standard survey analysis: **variance decreases with √n** for independent hauls, allowing signal to emerge.

**Recommendation:**

- H2: median (or GAM/mixed-model) residual per StatRec × decade; relate to mean annual fishing hours (Couce et al. data)
- H3: StatRec-level residual time series; test for trends or breakpoints
- Always report **n hauls per cell** and sensitivity to minimum-n threshold

### 3.5 Confounding and the quartile offset

H1 shows median B_pred/B_obs increases from 3.1× (smallest hauls) to 5.5× (largest). If:

- High fishing pressure → depleted communities → smaller catches → lower biomass quartile
- Lower quartile → smaller EEoS overprediction ratio

…then a **negative correlation between fishing pressure and EEoS residual could arise from scale-dependent bias alone**, not from ecological condition.

**Recommendation for H2:**

1. Include **log(B_obs)** or biomass quartile as covariate
2. Test residual ~ pressure with and without adjustment
3. Report sensitivity analysis stratified by biomass quartile
4. Consider **relative residual** (residual − StatRec mean residual) to remove spatial calibration offsets

### 3.6 Spatial and temporal autocorrelation

ICES StatRec rectangles are spatially contiguous; consecutive years are autocorrelated. Standard OLS on StatRec-level residuals **inflates significance** (Legendre 1993; Dormann et al. 2007).

**Recommendation:**

- Use mixed models with StatRec random effect, or
- Spatial models (CAR/ICAR on StatRec adjacency), or
- At minimum: block bootstrap by StatRec or year when testing pressure–residual association

### 3.7 Permutation and null models

H1 null (randomise B_obs, hold B_pred fixed) is appropriate for "is there haul-specific signal?" It is **not** appropriate for H2 ("does residual covary with pressure?") without redesign.

**Recommendation for H2 null:** Permute fishing pressure labels within StatRec or within spatial blocks, preserving the marginal distribution of pressure. Test whether observed residual–pressure association exceeds chance.

---

## 4. Recommended minimal reporting set

### For H1 (complete — now in `H1_results_summary`)

- Tier 1 (primary H1 benchmark): log_r2, cor², log RMSE, median pred/obs for EEoS and **unfitted productivity 1:1**
- Extended: ln(E) OLS (clearly labelled fitted; not prior-method baseline)
- Null model p-value
- Quartile table of median pred/obs
- Signed median residual

### For H2 (proposed)

| Statistic | Definition |
|-----------|------------|
| **Primary response** | Mean `\|log(B_obs) − log(B_pred)\|` per StatRec × decade (magnitude of deviation from the EEoS equilibrium) |
| **Primary predictor** | Mean annual reconstructed fishing hours (otter + beam) per ICES rectangle, 1985–2015 |
| **Effect size** | Slope or group difference in mean absolute residual vs pressure (with CI) |
| **Diagnostics** | n hauls per cell; median pred/obs per cell; map of residuals and pressure |
| **Sensitivity** | Model with biomass quartile covariate; minimum-n threshold (e.g. ≥ 5 hauls) |

### For H3 (proposed)

| Statistic | Definition |
|-----------|------------|
| **Primary response** | StatRec-level mean absolute residual time series (annual or decadal) |
| **Primary test** | Trend, changepoint, or correlation with cumulative pressure index |
| **Diagnostics** | Autocorrelation (ACF/PACF); n hauls per year per StatRec |

---

## 5. Summary judgement

**Most useful H1 statistics looking forward:**

1. **Mean absolute log residual** — the H2/H3 response variable after aggregation (magnitude of deviation from the EEoS equilibrium; signed residual retained only as a calibration diagnostic)
2. **Median B_pred/B_obs by stratum** — calibration diagnostic and confound control
3. **Null model result** — justifies that residuals are informative, not noise
4. **Quartile-dependent offset** — critical confound for spatial pressure analysis

**Most misleading if reused:**

1. **cor² = 0.93** — suggests success where prediction failed
2. **Global log_r2** — wrong question for H2/H3
3. **ln(E) OLS comparison** — unfair fitted benchmark; wrong frame for residual inference

**Best ecological framing for supervisors:**

> H1 shows EEoS does not predict absolute haul biomass in a trawl survey (consistent with catchability and scale mismatch). The residuals are structured (p < 0.001) and biased in a scale-dependent way. H2/H3 ask whether that structured deviation tracks fishing pressure — a standard null-model inference problem — not whether EEoS beats ln(E) on R².

---

## References (selected)

- Araujo, M. B., et al. (2005). Evaluation of species distribution models. *Ecography*, 28, 693–705.
- Burnham, K. P., & Anderson, D. R. (2002). *Model selection and multimodel inference*. Springer.
- Dormann, C. F., et al. (2007). Methods to account for spatial autocorrelation. *Ecography*, 30, 609–619.
- Dormann, C. F., et al. (2018). Model selection in ecology and evolution. *Journal of Animal Ecology*, 87, 343–349.
- Harte, J., et al. (2022). An equation of state unifies biodiversity, ecosystem functioning, and biomass. *Communications Biology*, 5, 957.
- Legendre, P. (1993). Spatial autocorrelation: trouble or new paradigm? *Ecology*, 74, 1659–1673.
- Pineiro, G., et al. (2008). How to evaluate models: Observed vs predicted or predicted vs observed? *Ecological Modelling*, 216, 316–322.
- Roberts, D. R., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40, 913–929.
- Thiele, J. C., et al. (2022). Best practices for calibrating models. *Ecological Modelling*, 471, 110038.

---

*Companion to [`H1_results_summary.md`](H1_results_summary.md). Pipeline outputs: `outputs/h1_model_comparison.csv`, `outputs/haul_h1_benchmarks.rds`.*
