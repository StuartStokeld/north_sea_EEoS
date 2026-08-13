Aim - self critique my models approach to spatial autocorrelation

1. Diagnose the spatial structure of my data 
    - rectangles are highly similar
2. Review how my model attempts to address this 
    - a plain rectangle random intercept, treating rectangles as independent of one another (no distance or neighbour weighting: which were tested but not incorporated).
3. Are the H2 results biased? use bootstrapping sensitivity analysis to destroy the spatial structure
    - Yes but not uniformly - first two phases are not confounded, second two are
4. Test using a more ecologically relevant spatial structure
    - Use GEBCO bathymetry data - test whether spatial correlation run stronger along the shelf than across it? Could this be used as a more environmentally informed distance metric?
5. Test an approach used in the same study zone / data

#### 1. Diagnostics: raw spatial structure of my inputs (rectangles are highly similar)

Spatial-similarity, measuring **neighbourhood autocorrelation on the ICES rectangles.**

- “Similarity” here means: neighbouring rectangles share similar values

| **Variable** | **Moran’s I** | **Geary’s C** |
| --- | --- | --- |
| log fishing hours | 0.60 | 0.39 |
| Mean |residual| | 0.59 | 0.40 |
| OLS residuals | **0.56** | 0.43 |

These three are testing different things at the same spatial scale: the raw dependent variable, the raw predictor (fishing hours), and the residuals of the simplest possible model (OLS, no rectangle effects).

All p ≪ 0.001. Neighbouring rectangles are **strongly clustered** in both fishing pressure and residuals. These are diagnostics of the raw spatial structure in the inputs - the key question to answer and defend is how the model approach absorbed that structure.

#### 2. How does my model attempt to address this? (initial approach assumed rectangle differences are exchangeable not distance or neighbour weighted)

Two alternatives were tried before settling on the model.

**A) Distance decay (correlation falling off with distance)** 

where - `exp(pos + 0 | dummy)` is an **exponential** spatial correlation over lon/lat centroids:

correlation between two rectangles ≈ exp(−distance / range)

- **range** = correlation-decay distance (in **degrees**, because centroids are lon/lat)

Result: fitted value from that check: **~17.7°**, with a huge CI, and larger than the max inter-centroid distance (~13°)

There finite decay scale "correlation dies off after X distance" visible in the data/ study area - so this was **not adopted**

**B) CAR (correlation only between direct neighbours) - used as sensitivity**

The *adjacency structure* (queen contiguity) is close to explaining all the between-rectangle variance on its own: 

- **rectangle intercepts (BLUPs) at 0.55**, and **model residuals at 0.47–0.48**

knowing which rectangles are next-door neighbours explains almost as much of the between-rectangle differences as just letting every rectangle have its own separate intercept `(1|rectangle)`

**Primary model**

The primary model uses a plain rectangle random intercept, (1|rectangle), treating rectangles as independent of one another (no distance or neighbour weighting).

- **Primary spatial assumption:** rectangles differ (random intercepts), but those differences are **exchangeable**, not distance-weighted.

### Does the primary model actually address spatial autocorrelation? (No)

I ran the same Moran's I / Geary's C test on the rectangle intercepts and the model's residuals.

| Stage | Moran's I | Geary's C |
| --- | --- | --- |
| Raw inputs (above) | 0.56–0.60 | 0.39–0.43 |
| **Rectangle intercepts from our model** | **0.55** | **0.44** |
| **Model residuals** | **0.47–0.48** | **0.50–0.51** |

All still p ≪ 0.001.

The rectangle intercepts are barely smaller than the raw signal, meaning the model isn't really treating rectangles as independent/exchangeable in practice, just re-expressing the same spatial pattern.

#### 3. Are the H2 results biased? use bootstrapping sensitivity analysis to destroy the spatial structure

is the reported `FP_between` (H2, spatial fishing-pressure) effect inflated or distorted because the rectangle intercept and `FP_between` are both tracking the same spatial pattern (spatial confounding)?

**Test:** spatial permutation bootstrap. `FP_between`'s rectangle-level values were
randomly reassigned across the 158 rectangles (breaking the covariate's real spatial
arrangement while leaving everything else in the data untouched), the full model refit
1000 times, and the observed coefficient(s) compared against the resulting null
distribution.

Originally run against the exchangeable RE model `(1 | stat_rec)`. Re-run (2026-08-06)
against the **CAR** primary used for reported H2 spatial contrasts
(`adjacency(1 | stat_rec)`, same queen adjMatrix as `run_h2h3_phase_v2_reporting.R`).
RE results retained below for comparison; current check is the CAR table.

**Result (CAR — current):**

| Phase | H2 (`FP_between`) coefficient — CAR | Spatially confounded? |
| --- | --- | --- |
| 1985–1991 (reference) | −0.063 | **Yes** |
| 1992–2001 | −0.024 | **Yes** |
| 2002–2007 | +0.084 | **Yes** |
| 2008–2015 | +0.013 | No |

**Result (RE — archived sensitivity):**

| Phase | H2 (`FP_between`) coefficient — RE | Spatially confounded? |
| --- | --- | --- |
| 1985–1991 (reference) | −0.025 | No |
| 1992–2001 | +0.014 | No |
| 2002–2007 | +0.122 | **Yes** |
| 2008–2015 | +0.052 | **Yes** |

**spatially confounded** = YES if the observed coefficient falls outside the null distribution's 95% interval

Null-width comparison (CAR null is modestly narrower than RE; SD ratio CAR/RE ≈ 0.70–0.93):

| Phase | RE null SD | CAR null SD | RE null IQR | CAR null IQR |
| --- | --- | --- | --- | --- |
| 1985–1991 | 0.0158 | 0.0147 | 0.0214 | 0.0198 |
| 1992–2001 | 0.0152 | 0.0111 | 0.0199 | 0.0153 |
| 2002–2007 | 0.0198 | 0.0139 | 0.0254 | 0.0190 |
| 2008–2015 | 0.0168 | 0.0126 | 0.0228 | 0.0164 |

Scripts/outputs: `pipeline/permutation_bootstrap_FP_between_CAR.R`,
`outputs/permutation_bootstrap_FP_between_CAR_summary.md`; archived RE at
`exploratory/pipeline/permutation_bootstrap_FP_between_RE.R`.

This is not uniform across periods (and the phase pattern differs between RE and CAR).

#### **4. Next step** - Use something more ecologically relevant.

**Test:** GEBCO bathymetry (depth, depth gradient)
extracted per rectangle

does spatial correlation in our model residuals run stronger along the shelf than across it? Could this be used as a more environmentally informed distance metric? (rather than a more complex covariance structure)

The goal is to decide between two explanations for the unidentified decay range:

1. Correlation genuinely has no finite spatial range at this resolution (confounded with
the fixed-effect spatial trend, or a borderless field).
2. Correlation has a real, finite range, but geographic (Euclidean lon/lat) distance is
the wrong metric to detect it, shelf geometry is a better one.

**Test - targets model residual spatial structure**

For each of 158 ICES rectangles, we derived **local along-shelf direction** from GEBCO bathymetry (depth gradient rotated 90°), and **residual-correlation direction** from where nearby pairs showed the strongest similarity in residuals. The confirmatory check was a **circular–circular correlation** (Jammalamadaka–Sarma): do those two directions align? Confirmed if **p < 0.05**.

Result:

**Anisotropy not confirmed** (ρ = 0.05, p = 0.55).

No evidence that shelf-geometry-based direction explains the earlier unidentified decay range. Consistent with either a genuinely borderless spatial field or confounding with the fixed-effect spatial trend, does not distinguish between these two.

Residual-correlation direction does **not** align with local shelf geometry. That argues against “geographic distance was the wrong metric” as the explanation for the earlier unidentified decay range. The result is instead consistent with a **borderless spatial field** or **confounding from the fixed-effect spatial trend,** but it does not distinguish between those two.

#### 5. Test an approach used in the same study zone / data

https://www.researchgate.net/publication/283210413_Spatio-temporal_Bayesian_network_models_with_latent_variables_for_revealing_trophic_dynamics_and_functional_networks_in_fisheries_ecology#pf2 

"we enforce three parent nodes that represent the average biomass... from the spatial neighbourhood (the three or four nearest neighbours) of the current area” : citing Aderhold et al. (2012) for the technique itself

for each area, compute the mean of a variable across its k nearest neighbours, and add that mean as an ordinary covariate on the right-hand side.

**What the model does:**

The primary H2 specification is CAR (`adjacency(1 | stat_rec)`). Spec A adds one term: `FP_between_lag` — the average long-run fishing pressure of each rectangle's four nearest neighbours (k-NN, k=4) — interacted with `phase_v2`, alongside the existing `FP_between` / `FP_within` phase structure.

**What it tests:** whether the spatial confounding already found in `FP_between` under CAR (phases 1985–1991, 1992–2001, and 2002–2007 outside the spatial-permutation null; 2008–2015 clean) is explained by fishing effort spilling over between adjacent rectangles.

**Identifiability (Task 1):** a single CAR+lag fit cleared the gate — ρ essentially unchanged vs primary CAR (0.131292 → 0.131287), `|cor(FP_between_lag, CAR BLUPs)| ≈ 0.32` (below the 0.7 flag), lag SEs not inflated vs RE Spec A or vs `FP_between` in the same model, clean convergence.

**Results (Task 2, CAR+lag permutation; lag recomputed each shuffle):**

| Phase | CAR-alone confounded? | CAR+lag confounded? | CAR+lag `FP_between` |
|-------|------------------------|---------------------|----------------------|
| 1985–1991 | Yes | **No** | +0.010 |
| 1992–2001 | Yes | **No** | −0.002 |
| 2002–2007 | Yes | **No** | −0.022 |
| 2008–2015 | No | No | −0.018 |

All three CAR-flagged phases fall inside the null 95% once neighbour fishing pressure is in the mean structure. Spec A is reported as a **confirmed mechanism / sensitivity check on the primary CAR H2 model**: fishing-effort displacement between neighbouring rectangles explains the **1985–2007** confounding in H2 under CAR (not the earlier RE-era “2002–2015” claim).

Historical note: Spec A was first tested against the superseded RE `(1 | stat_rec)` baseline, where it cleared the RE-flagged 2002–2007 / 2008–2015 pattern. That RE-scoped result is retained as the development path; the **current claim** is the CAR-scoped result above.

Moran's I on RE BLUPs moved only negligibly when lag was added there (0.552 → 0.532). Whatever drives general residual spatial clustering is *not* primarily about fishing pressure's spatial arrangement. Spec B (lagged neighbour biomass) was explored for that question and is **archived** under `exploratory/` (not in the live methods model list).

- the primary CAR model's residual spatial structure is not erased by `FP_between_lag`; Spec A is a mechanism account for H2 coefficient confounding, not a cure for exchangeability failure
