# H1 null model — statistics and conclusions

Generated from `run_h1_null_model.R` on haul-level predictions (`outputs/haul_eeos_predictions.rds`, **n = 12,069** hauls). Refresh: `Rscript run_h1_null_model.R`.

## 1. What question does this test?

**Not implemented (and rejected):** shuffle S, N, E independently across hauls. That breaks EEoS constraints (`E > N`, integer state variables, soft hierarchy `N >> S`, `E >> N`).

**Implemented:** keep **observed S, N, E** per haul → fixed **`B_pred`**. Randomise **`B_obs`** only, then recompute log-scale R²:

\[
R^2_{\text{null}} = R^2\big(\log B_{\text{obs,null}},\; \log B_{\text{pred}}\big)
\]

Compare to observed \(R^2_{\text{obs}} = R^2(\log B_{\text{obs}},\; \log B_{\text{pred}})\).

Primary residual for H1: **`residual = log(B_obs) − log(B_pred)`** (both in grams).

---

## 2. Null sampling scheme

Automated assessment of **`log(B_obs)`** shape (`R/h1_null_helpers.R`):

| Metric | Value | Peaked threshold |
|--------|-------|------------------|
| IQR / range | **0.12** | < 0.25 → peaked |
| Histogram peak / uniform | **4.96** | > 2.5 → peaked |
| Excess kurtosis | **1.02** | > 0.75 → peaked |

All three criteria flagged **peaked**. Primary null = **uniform on central 95% of log(B_obs)**; robustness = **B shuffle**.

Full metrics: `outputs/null_distribution_summary.csv`  
Decision record: `outputs/null_sampling_decision.rds`

---

## 3. Observed H1 performance

| Statistic | Value |
|-----------|-------|
| \(R^2_{\text{obs}}\) (log scale, log_r2) | **−0.223** |
| cor(log B)² (diagnostic) | **0.926** |
| Median \|residual\| | **1.41** |
| Median \(B_{\text{pred}} / B_{\text{obs}}\) | **4.08×** (EEoS **over**predicts) |

Negative log_r2: EEoS predictions are further from observed biomass than predicting the mean log(B_obs) would be — poor absolute haul-level prediction.

---

## 4. Null simulation results (999 permutations, seed 42)

### Primary: uniform on central 95% of log(B_obs)

| Quantity | Value |
|----------|-------|
| Null \(R^2\) median | **−2.02** |
| Null \(R^2\) 2.5%–97.5% | **−2.08** to **−1.96** |
| Observed \(R^2\) | **−0.223** |
| \(p\) one-sided: P(null \(R^2\) ≥ observed) | **< 0.001** |
| Median \|residual\| under null | **1.78** |

Observed performance is **substantially better** than random biomass from the uniform null envelope.

### Robustness: B shuffle

| Quantity | Value |
|----------|-------|
| Null \(R^2\) median | **−2.60** |
| Null \(R^2\) 2.5%–97.5% | **−2.64** to **−2.56** |
| \(p\) one-sided | **< 0.001** |
| Median \|residual\| under null | **1.67** |

Observed pairing of (S, N, E) → B_pred with haul-specific B_obs beats permuted biomass labels.

---

## 5. Conclusions (supervisor-ready)

1. **S/N/E permutation null was not run** — correctly omitted; violates EEoS structural constraints.

2. **Distribution-driven choice:** log(B_obs) is peaked → primary null uses uniform draws on central 95% of log(B), not empirical resampling.

3. **Absolute predictive failure:** log_r2 ≈ −0.22; median overprediction ~4×. EEoS does not recover haul-level biomass from S, N, E alone in NS-IBTS Q1.

4. **Relative to null:** observed fit beats both null schemes (p < 0.001) — weak non-random co-structure, not useful prediction.

5. **vs ln(E) benchmark:** ln(E) log_r2 ≈ 0.60; Harte criterion **not met** (EEoS SS ~3× ln(E) SS). See `outputs/h1_model_comparison.csv`.

6. **Suggested phrasing:**

   > “EEoS does not predict absolute haul biomass (negative log_r2, ~4× median overprediction). The observed B_obs–B_pred pairing exceeds random biomass labels under both null schemes (p < 0.001), but ln(E) alone outperforms EEoS on the Harte criterion — S and N do not improve on metabolic rate alone at haul scale.”

---

## 6. Related files

| File | Contents |
|------|----------|
| [`review_supervisor_h1_briefing.Rmd`](../display_discussion/review_supervisor_h1_briefing.Rmd) | Supervisor briefing (methods + results; share HTML) |
| [`review_h1_talk_track.md`](../display_discussion/review_h1_talk_track.md) | Personal talk track (not for sharing) |
| `outputs/null_summary.rds` | Combined null summary |
| `outputs/figures/null_r2_*.png` | Null R² histograms |
