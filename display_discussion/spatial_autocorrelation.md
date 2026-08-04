Aim - self critique my models approach to spatial autocorrelation

1. Diagnose the spatial structure of my data 
    - rectangles are highly similar
2. Review how my model attempts to address this 
    - a plain rectangle random intercept, treating rectangles as independent of one another (no distance or neighbour weighting: which were tested but not incorporated).
3. Are the H2 results biased? use bootstrapping sensitivity analysis to destroy the spatial structure
    - Yes but not uniformly - first two phases are not confounded, second two are
4. Test using a more ecologically relevant spatial structure
    - Use GEBCO bathymetry data - test whether spatial correlation run stronger along the shelf than across it? Could this be used as a more environmentally informed distance metric?

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

The *adjacency structure* (queen contiguity) is close to explaining all the between-rectangle variance on its own: knowing which rectangles are next-door neighbours explains almost as much of the between-rectangle differences as just letting every rectangle have its own separate intercept `(1|rectangle)`

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

This is not uniform across periods. This is interesting - there is a potential 

#### **4. Next step** - Use something more ecologically relevant.

**Planned test:** GEBCO bathymetry (depth, depth gradient)
extracted per rectangle

does spatial correlation in our model residuals run stronger along the shelf than across it? Could this be used as a more environmentally informed distance metric? (rather than a more complex covariance structure)

The goal is to decide between two explanations for the unidentified decay range:

1. Correlation genuinely has no finite spatial range at this resolution (confounded with
the fixed-effect spatial trend, or a borderless field).
2. Correlation has a real, finite range, but geographic (Euclidean lon/lat) distance is
the wrong metric to detect it, shelf geometry is a better one.
