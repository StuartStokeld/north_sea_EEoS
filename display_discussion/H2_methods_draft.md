# H2 methods – draft

## H2: Fishing pressure and EEoS residual magnitude

### 1. Study design and variables

**Implemented:** Cross-sectional analysis at the ICES statistical rectangle level (N = 158, default panel), using NS-IBTS Q1 hauls from 1985 to 2015. The dependent variable is the mean absolute log-residual between observed and EEoS-predicted biomass per rectangle (mean_abs_residual), pooled across all hauls in that rectangle over the study period. The independent variable is mean annual fishing hours per rectangle (Couce et al. 2020), averaged across 1985 to 2015.

**Why:** Key comparison metric for H2, absolute residuals are used, disturbance can push observed biomass above or below the prediction, so only the magnitude of deviation is theoretically informative for a disturbance hypothesis.

### 2. Ordinary least squares (OLS) with spatial diagnostics

**Implemented:** A baseline OLS regression of mean_abs_residual on mean annual fishing hours, followed by Moran's I and Geary's C tests on the OLS residuals using queen contiguity spatial weights (row-standardised).

**Why:** Ecological and fishing-pressure data are both spatially structured, so residual spatial autocorrelation was checked before treating the OLS result as reliable. Significant autocorrelation would indicate that standard errors and significance from OLS cannot be trusted.

### 3. Spatial error model (SEM)

**Implemented:** A spatial error model (errorsarlm()) fit on the same specification, using the same spatial weights.

**Why:** SEM was the pre-registered approach for correcting standard errors and coefficients for spatial autocorrelation, and was treated as the primary inferential model for testing the H2 hypothesis.

### 4. Threshold sensitivity across haul-count panels

**Implemented:** The core model set (OLS, Moran's I/Geary's C, SEM) was re-run on two additional panels built with lower and higher minimum haul-count thresholds per rectangle (≥5 hauls, N = 161; ≥20 hauls, N = 156), with spatial weights rebuilt independently for each panel.

**Why:** The original robustness check only re-ran OLS across thresholds. Spatial models were added at each threshold to test whether the "stable" OLS result held once spatial autocorrelation was accounted for, rather than assuming it would.

### 5. Leverage and influence diagnostics

**Implemented:** Cook's distance and leverage (hat values) computed for the primary OLS model. Rectangles exceeding the conventional threshold (4/N) were identified and the model was refit excluding them.

**Why:** With a low OLS R² (approximately 0.04), it was necessary to check whether the significant coefficient was driven by a small number of high-leverage rectangles rather than reflecting a broad pattern.

### 6. Biomass covariate specification

**Implemented:** A second model specification added mean log observed biomass per rectangle (mean_ln_B_obs, the mean of log-transformed haul-level biomass, not the log of the mean) as a covariate alongside fishing hours, fit under both OLS and SEM.

**Why:** H1 established that EEoS systematically over-predicts biomass with an offset that may vary by biomass level. This specification tests whether that scale-dependent bias confounds the fishing-pressure relationship at rectangle scale.

### 7. Lagrange Multiplier (LM) specification tests

**Implemented:** Classic and robust LM tests for spatial error versus spatial lag dependence (lm.LMtests()), run separately on the residuals of the primary OLS model and the biomass-covariate OLS model, using the same spatial weights as the SEM fits.

**Why:** SEM was chosen a priori without a formal test against alternative spatial specifications. This step tests empirically whether a spatial error or spatial lag structure better matches the data, and whether that answer differs once the biomass covariate is included.

### 8. Spatial lag model (SAR) as an alternative specification

**Implemented:** A spatial lag model (lagsarlm()) fit for both the primary and biomass-covariate specifications, on the same panel and weights as the corresponding SEM fits, compared by coefficient, spatial parameter, log-likelihood, and AIC.

**Why:** The LM tests indicated some support for a lag specification in the primary model. This step tests whether model choice (error vs lag) changes the substantive conclusion for either specification, rather than relying on SEM by assumption.

---

## Critiques / limitations

Couce et al.'s fishing-pressure reconstruction does not cover two edge regions of the survey area (the Skagerrak/Kattegat approach and the eastern English Channel approach), excluding around 23 rectangles from the H2 model on this basis alone (plus a further few for insufficient haul count (N=10).

---

## Potential extensions

### Potential H2 expansion: adding SAR data

The current H2 model uses Couce et al. (2020) fishing hours as the pressure covariate, which offers strong temporal depth (1985–2015) but measures effort rather than physical seabed disturbance directly. The ICES WGSFD annual swept area ratio (SAR) product provides a gear-disaggregated, more mechanistically direct disturbance measure at fine spatial resolution, aggregable to ICES rectangles, but covers only 2009–2015. Adding SAR would allow three linked model comparisons rather than one: the existing Couce-based model over the full period, a SAR-based model restricted to the 2009–2015 overlap, and a Couce-based model restricted to the same overlap window as a like-for-like check. This structure directly tests whether the current H2 result (fishing hours showing no robust association with residual magnitude across sample threshold, leverage, and spatial specification checks) reflects a genuine absence of a fishing-pressure effect, or whether fishing hours understate a disturbance signal that a more direct SAR-based measure would detect.

### Dominance-covariate extension

Berger-Parker dominance (D) and dominant-species size homogeneity (size_CV) were tested as haul-level correlates of EEoS residual magnitude, computed from the same raised abundance and length data underlying N and E rather than from an external diversity index. A nested regression on the H1 haul-level panel showed both terms contribute significant explanatory power beyond biomass magnitude alone (baseline R² = 0.237, R² with D and size_CV added = 0.245; D: β = 0.197, p < 10⁻²²; size_CV: β = −0.055, p < 10⁻³), with no interaction between them (ΔR² < 0.0001, p = 0.79). Before incorporating D into the H2 rectangle-level model as a control covariate, its association with the Couce et al. (2020) fishing-pressure covariate was tested directly at the rectangle level, independent of the residual variable. This showed a weak but significant positive correlation (Pearson r = 0.183, Spearman ρ = 0.246, p < 0.05, n = 161), robust to removal of high-leverage rectangles, while size_CV showed no such association (p = 0.298). D is therefore included as a covariate in the H2 model on this basis; size_CV is retained in the H1 analysis but not required as an H2 control.
