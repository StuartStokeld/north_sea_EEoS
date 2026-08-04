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

**Result:**

| Phase | H2 (`FP_between`) coefficient | Spatially confounded? |
| --- | --- | --- |
| 1985–1991 (reference) | −0.025 | No |
| 1992–2001 | +0.014 | No |
| 2002–2007 | +0.122 | **Yes** |
| 2008–2015 | +0.052 | **Yes** |

**spatially confounded** = YES if the observed coefficient falls outside the null distribution's 95% interval

This is not uniform across periods.

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

a mixed-effects regression predicting EEoS prediction error (`residual`, the gap between predicted and observed biomass) from a rectangle's own long-run fishing pressure (`FP_between`), that same pressure's year-to-year fluctuation (`FP_within`), and policy-era phase — with one addition: a new term, `FP_between_lag`, representing the average long-run fishing pressure of each rectangle's four nearest neighbours.

**What it tests:** whether the spatial confounding already found in `FP_between` (its coefficient behaving inconsistently across phases, especially post-2002, in a way tied to rectangles' spatial arrangement rather than fishing pressure itself) is explained by fishing effort spilling over between adjacent rectangles

Results:

Both previously-flagged phases (2002–2007, 2008–2015) now fall inside the null 95% under the new step 5 model

Moran's I moved negligibly **0.552** to ****0.532 is a small, likely negligible movement (Δ = −0.018, ~3% relative reduction), residual spatial structure remains strong. 

(1 | rectangle) remains the primary model. Spec A is reported as a confirmed mechanism finding and robustness check on H2

fishing-effort displacement between neighbouring rectangles explains the 2002–2015 confounding in H2.

whatever's driving the general spatial clustering in residuals is *not* primarily about fishing pressure's spatial arrangement

- the primary model's core assumption (rectangle differences are exchangeable) remains empirically false after every fix attempted so far (distance decay unidentified, CAR non-identifiable, `FP_between_lag` moves BLUPs by only ~3% relative)
