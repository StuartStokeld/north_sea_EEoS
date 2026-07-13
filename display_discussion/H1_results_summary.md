# H1 Results Summary

**NS-IBTS Q1 (1985–2015) · 12,069 hauls · EEoS haul-level biomass prediction**

---

## Results

### H1: Does EEoS predict haul-level community biomass?

Hypothesis 1 asks whether the Ecological Equation of State (EEoS) predicts community biomass at individual trawl haul level from three state variables — species richness (S), abundance (N), and normalised metabolic rate (E) — in the NS-IBTS Q1 survey (1985–2015). After joining DATRAS length data to FishGlob catch biomass and applying standard filters, **12,069 hauls** entered the analysis. For each haul, EEoS returns a parameter-free prediction B_pred from (S, N, E); this is compared to observed catch biomass B_obs without any post-hoc calibration.

The **primary comparison** replicates prior EEoS validation work (Harte et al. 2022): predicted biomass is evaluated against productivity **directly**, via an **unfitted 1:1 map** — not via a statistical model fit to the same data. Prior findings that “EEoS-predicted biomass is equivalent to productivity” refer to equivalence under a one-to-one map, **with productivity, not to a fit model**. I compare:

1. **EEoS (Harte Fig 1):** log(B_obs) vs log(B_pred), 1:1 line, no fitting.
2. **Productivity 1:1 (prior baseline):** log(B_obs) vs log(E × m_min), same identity map — same normalised E as `biomass()`, same m_min as B_pred grams conversion. (Former E_raw map retained as diagnostic only.)

Fitted ln(E) OLS is reported separately as an **extended comparison** — a different comparison point; R² values are not directly comparable.

---

### Main result: EEoS vs unfitted productivity 1:1

EEoS **does not** predict absolute haul-level biomass. log_r2 = **−0.22** (worse than the mean). Under the same unfitted footing, productivity 1:1 (E × m_min) achieves log_r2 = **+0.74**. EEoS leaves **4.6×** more unexplained variance than the productivity baseline — Harte criterion **NOT MET**.

The EEoS scatter sits above the 1:1 line (median B_pred/B_obs = **4.1×**); cor² = **0.93** (diagnostic only). At haul scale, EEoS **underperforms** a simple unfitted productivity map.

**Figure 1** (see HTML): *Left* — Harte Fig 1, EEoS log(B_pred) vs log(B_obs). *Centre* — productivity 1:1, log(E × m_min) vs log(B_obs). *Right* — log_r2 bar chart (both unfitted).

| Model | Fitted | log_r2 | cor² | log RMSE | Median pred/obs | SS_res ratio vs prod1 |
|-------|--------|-------:|-----:|---------:|----------------:|----------------------:|
| EEoS (S, N, E) | No | −0.22 | 0.93 | 1.41 | 4.1× | 4.6× |
| Productivity 1:1 (E × m_min) | No | +0.74 | 0.90 | 0.65 | 0.67× | — |

---

### Does NS-IBTS behave like prior EEoS datasets?

Prior validations report EEoS biomass equivalent to productivity under a direct 1:1 comparison. On NS-IBTS hauls, that equivalence **does not hold** for EEoS (log_r2 = −0.22; cloud above 1:1 throughout). Unfitted productivity 1:1 alone outperforms EEoS (+0.74 vs −0.22) despite EEoS using S, N, and E through the full partition function. The marine trawl data **shows something different** from prior EEoS datasets at haul scale — established on the same unfitted metric footing as Harte et al. before any fitted alternative test.

---

### Supporting statistics

**Two R² types:** **log_r2** (coefficient of determination) penalises scale error; **cor²** (correlation squared) measures rank agreement only. EEoS is a mechanistic post hoc check; productivity 1:1 is unfitted the same way. Neither is directly comparable to fitted OLS R² without noting OLS absorbs intercept/shift.

#### log_r2 = −0.22 (EEoS) vs +0.74 (productivity 1:1)

**What it measures / how calculated:** 1 − SS_res/SS_tot on log scale. EEoS: no parameters adjusted. Productivity 1:1: E × m_min as direct predictor (unfitted; m_min not fit to B_obs).

**Why I am using it:** Only unfitted metric penalising absolute scale error. EEoS fails to beat the mean; productivity 1:1 succeeds. Primary evidence this dataset differs from prior EEoS validations.

#### cor² = 0.93 (diagnostic only)

**What it measures:** cor(log B_obs, log B_pred)² — rank agreement, not coefficient of determination.

**Why I am using it:** High cor² with negative log_r2 explains tight scatter above 1:1 line. Not the H1 test statistic.

#### log RMSE = 1.41 (EEoS) vs 0.79 (productivity 1:1)

**What it measures:** Typical log-scale prediction error.

**Why I am using it:** Readable magnitude of failure; EEoS ≈ log(4) consistent with ~4× overprediction.

#### Median B_pred/B_obs = 4.1×

**What it measures:** Systematic multiplicative overprediction.

**Why I am using it:** Direct read on upward shift; catchability interpretation; not uniform across quartiles.

#### Median signed residual = −1.41

**What it measures:** Median log(B_obs) − log(B_pred).

**Why I am using it:** H2/H3 aggregation quantity.

#### SS_res ratio = 4.6× vs productivity 1:1

**What it measures:** Harte criterion on unfitted baseline (E × m_min). NOT MET.

**Why I am using it:** Like-for-like unfitted comparison Jake requested.

#### Null model: p < 0.001

| Null scheme | Null median log_r2 | Observed | p |
|-------------|---------------------:|---------:|---|
| Uniform 95% log(B_obs) | −2.02 | −0.22 | **< 0.001** |
| B shuffle | −2.60 | −0.22 | **< 0.001** |

**Why I am using it:** Residuals structured, not random — minimal condition for H2/H3. Does not mean accurate prediction.

#### Quartile overprediction

| Quartile | B range (log g) | Median B_pred/B_obs |
|----------|-----------------|---------------------|
| Q1 | 3.9–11.1 | 3.1× |
| Q2 | 11.1–11.9 | 3.9× |
| Q3 | 11.9–12.7 | 4.5× |
| Q4 | 12.7–17.4 | 5.5× |

**Why I am using it:** Not a fixed catchability factor; H2 confound if pressure correlates with biomass.

---

### Extended comparison: fitted ln(E) OLS (not prior-method baseline)

**Comparability warning:** EEoS = mechanistic post hoc check. ln(E) OLS = fit to minimise residuals. R² values test different things; OLS almost always wins on coefficient of determination because it unfolds the shift. Follow-on context only — **after** unfitted baseline.

| Model | Fitted | log_r2 | SS ratio vs ln(E) |
|-------|--------|-------:|------------------:|
| EEoS | No | −0.22 | 3.0× |
| ln(E) OLS | Yes | +0.60 | — |

Equation: log(B_obs) = 4.32 + 0.714 × log(E). Harte criterion vs ln(E): NOT MET.

---

### Brief conclusion

Unfitted replication: EEoS fails (log_r2 = −0.22); productivity 1:1 (E × m_min) wins (+0.74); Harte criterion NOT MET (4.6×). NS-IBTS differs from prior datasets on this test. Fitted ln(E) comparison follows but is not like-for-like on R². Structured residuals → H2/H3. B_pred = B_pred_norm × m_min (grams).

---

## Supplementary results

- **Harte Fig 2:** E/B^(3/4) ratio (normalised E, raw axes) — report **both** `fig2_r2_pearson_all` (0.964) and `fig2_r2_pearson_trimmed` (0.853, N=2); trimmed is the primary on-plot number. Extended cod R² in table footnote only.
- **Full model table:** see HTML / `outputs/h1_model_comparison.csv` (includes `productivity_1to1_uncalibrated` diagnostic)
- **Fig 2 leverage hauls:** `outputs/h1_fig2_leverage_hauls.csv`
- **Note:** productivity 1:1 uses E × m_min (unfitted); ln(E) OLS fits on normalised E — not like-for-like
- Raw RMSE = 7,668,817 g (supplementary)

*Full report: `H1_results_summary.html`*
