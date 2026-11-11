# H1 methods – draft

> **Updated Nov 2026 for the DATRAS DataType / haul-duration fix** (§1b; counts in §1 and §4).
> Numerical results quoted in the "Why" paragraphs (§5, §9, §10, §13: log_r2, cor², ~4× overprediction,
> quartile ratios, partial R²) are **pre-fix values, superseded pending interpretation after the DataType fix**;
> see `outputs/effort_fix_comparison.md`.

## H1: Does the Ecological Equation of State predict haul-level community biomass?

### 1. Study design and variables

**Implemented:** Haul-level analysis of NS-IBTS Q1 trawls, 1985–2015 (N = 10,128 valid 30-min hauls after filtering; 11,027 in the 26–34 min sensitivity run). For each haul, three state variables — species richness (S), abundance (N), and normalised community metabolic rate (E) — are derived from DATRAS length-frequency (HL) data and passed to the parameter-free Ecological Equation of State (EEoS) solver to generate predicted biomass (B_pred). This is compared against observed catch biomass (B_obs) from FishGlob, an independently processed data source, for the same hauls.

**Why:** Predictor (S, N, E) and outcome (B_obs) are deliberately built from two separately-aggregated sources rather than one, so the test is not circular. Haul level (rather than rectangle- or year-aggregated) is the finest resolution at which EEoS makes a prediction, and matches the resolution at which Harte et al. (2022) validated the theory.

### 1b. Haul selection and effort alignment (DataType, haul duration)

**Implemented:** DATRAS length-frequency (HL) records are joined to the survey's haul (HH) records on survey, year, quarter, country, ship and haul number (100% of HL hauls matched). Hauls are kept if the haul is valid (HaulVal = V), the DataType is C, R, P or S, and the recorded tow duration is **exactly 30 min** (primary analysis); a 26–34 min window (30 ± 4 min) is run as a sensitivity analysis. For DataType C, DATRAS reports numbers at length per hour of fishing (with SubFactor = 1), so counts are converted back to the haul: n = HLNoAtLngt × SubFactor × HaulDur/60; for R, P and S, n = HLNoAtLngt × SubFactor. No haul is rescaled to a common duration. Primary sample: 10,128 H1 hauls (DataType C 8,493; R 1,543; P 92).

**Why:** S, N, E and B_obs must describe the same tow. FishGlob's observed biomass is the raw haul weight (its `get_datras.R` converts DataType C back to the haul with HaulDur/60), so leaving C counts per hour would roughly double N and E for a 30-min tow relative to B_obs; after conversion, DATRAS-implied biomass sits at a consistent ~0.87–0.91 of B_obs for every DataType. Duration is restricted rather than rescaled because S cannot be rescaled: a 60-min tow (Scottish tows to 1998) roughly doubles N, E and biomass but raises S much less, so it is a different sampling unit, confounded with country and period. The NS-IBTS standard tow is 30 min, and exactly 30 min is the unit agreed for the primary analysis (supervisor decision, 10 Nov 2026). The current FishGlob `get_datras.R` excludes DataType P, but our local FishGlob build includes the 92 P hauls (likely predating their resubmission as P), and their biomass ratio (0.87) is consistent with raw haul weight, so they are kept. Sources: ICES (2020) *Manual for the North Sea International Bottom Trawl Surveys*, Series of ICES Survey Protocols SISP 10-IBTS 10, Revision 11, https://doi.org/10.17895/ices.pub.7562 (§2.5.4 tow duration; Annex 2 survey history); ICES DATRAS FAQ For Data Users, version 3 (2014), https://www.ices.dk/data/Documents/DATRAS/DATRAS_FAQs.pdf (for DataType C, HLNoAtLngt is adjusted to one hour of catching and SubFactor = 1); ICES DATRAS IBTS index calculation procedure, https://www.ices.dk/data/Documents/DATRAS/Indices_Calculation_Steps_IBTS.pdf (non-C counts × 60/HaulDuration); FishGlob `cleaning_codes/get_datras.R` (Maureaud et al., 2024).

### 2. Length–weight parameterisation and mass-at-length

**Implemented:** Species-level length–weight coefficients (W = a·Lᵇ) are obtained from FishBase (rfishbase::length_weight()), with a genus-mean → family-mean fallback for species lacking a direct record, and dropped otherwise. These coefficients are applied per 1 cm length bin (not to a single haul-mean length) to compute individual mass, which underlies every downstream quantity (N, E, m_min, dominance).

**Why:** Bin-level (not mean-length) mass application prevents one mis-measured or rare-length bin from distorting the whole haul's mass total. Consistently dropping species without a length–weight relationship — rather than imputing an untested value — ensures N, E and B_pred are all built from the same, internally consistent species pool.

### 3. Haul-level state variables (S, N, E)

**Implemented:** For each haul, S = number of distinct species, N = summed raw abundance (per-haul counts after the effort alignment in §1b), E_raw = Σ(abundance × mass^0.75), and E = E_raw normalised by the metabolic rate of the smallest individual in the haul (m_min^0.75), so the minimum individual's metabolic rate equals 1. Hauls are retained only if N ≥ S and E_raw, m_min are finite and positive.

**Why:** The normalisation sets the boundary condition required by the METE partition function that EEoS is built on; without it, E is not on the scale the solver expects. The N ≥ S and positivity checks remove internally inconsistent hauls before they can produce a numerically meaningless B_pred.

### 4. EEoS prediction and input filtering

**Implemented:** B_pred is generated by calling the reference Python implementation of the EEoS solver (biomass.py, from the original authors' codebase) via reticulate, one haul at a time, then rescaled to grams using the same m_min used to normalise E. Hauls failing structural requirements of the partition function (S < 2, E ≤ N, missing/non-finite inputs, B_obs ≤ 0) are excluded and logged with a specific exclusion reason before the solver is called.

**Why:** Using the original authors' solver, rather than a re-implementation, removes numerical-coding error as a possible confound when comparing against their published results. Logging exclusions by reason (rather than a single dropped-row count) makes the sample funnel from raw hauls to final N = 10,128 fully auditable (the funnel now includes an effort-alignment stage: invalid, wrong DataType, or duration ≠ 30 min).

### 5. Primary evaluation metric

**Implemented:** The headline statistic is the log-scale coefficient of determination, log_r2 = 1 − SS_res/SS_tot on log(B_obs) vs log(B_pred) (can be negative). Squared Pearson correlation of the same log values (cor²) and log-RMSE are reported alongside as diagnostics only.

**Why:** log_r2 penalises systematic scale error, whereas cor² only measures rank/linear association and stays high even under a large constant multiplicative bias. This distinction matters here: EEoS shows cor² ≈ 0.93 but log_r2 ≈ −0.22, because predictions sit tightly but consistently above the 1:1 line — reporting cor² alone would hide that failure mode.

### 6. Test 1 — unfitted 1:1 reproduction (Harte Fig 1)

**Implemented:** EEoS's B_pred is compared to B_obs on a strict, unfitted 1:1 line (no slope/intercept estimated), exactly reproducing Harte et al.'s (2022) original validation figure. The identical unfitted 1:1 test is run against a naïve productivity baseline (E × m_min, using the same normalised E and m_min as EEoS) for comparison.

**Why:** Reproducing the original authors' own unfitted, parameter-free comparison means any failure cannot be attributed to setting a harder bar than the source paper. Testing the identical comparison against a simple E × m_min proxy controls for the possibility that any apparent EEoS performance is not actually coming from the added structure of S and N.

### 7. Test 2 — E/B^0.75 ratio reproduction (Harte Fig 2)

**Implemented:** Predicted ratio (E/B_pred^0.75) is compared to observed ratio (E/B_obs^0.75) via unweighted linear regression on raw (unlogged) ratios, with Pearson R² reported to match the convention in Harte et al.'s source code exactly. A leverage diagnostic ranks hauls by predicted×observed ratio magnitude and reports Pearson R² with the top 2 highest-leverage hauls excluded, alongside the untrimmed value; these hauls are excluded only from this ratio regression, not from the main dataset.

**Why:** Matching the exact statistic (raw-ratio Pearson R², not a log-scale version) is required for a like-for-like comparison to Harte's reported value. The leverage trim exists because NS-IBTS ratios span 2–3 orders of magnitude (versus ~1–9 in the original 42 plots), so raw-ratio Pearson R² here is disproportionately sensitive to a few extreme hauls; reporting both values prevents an outlier-driven number being presented as equivalently stable to the original.

### 8. Extended comparison — fitted ln(E) OLS

**Implemented:** A fitted regression, log(B_obs) ~ log(E), is estimated on the same hauls and reported strictly after, and separately from, the two unfitted tests above.

**Why:** This isolates what happens when a model is allowed to absorb the systematic offset that the unfitted Test 1 penalises. Keeping it clearly labelled as a different, easier question prevents a good fitted R² from being read as validation of the unfitted mechanistic claim.

### 9. Null model (permutation test)

**Implemented:** B_pred is held fixed and B_obs is randomised across 999 permutations (seed fixed) under one of two schemes: a uniform draw over the central 95% of log(B_obs), or a shuffle of observed B_obs across hauls. The primary scheme is chosen automatically from the shape of the log(B_obs) distribution (peakedness diagnosed via IQR/range, histogram peak ratio, and excess kurtosis), with the other retained as a robustness check. Observed log_r2 is compared to both null distributions (one-sided p-value).

**Why:** This tests whether the haul-specific (S, N, E) → B_pred mapping carries any real information, independent of whether it is accurate. Choosing the null scheme based on distribution shape, rather than defaulting to one, avoids an empirical resample that would trivially reproduce an already-peaked biomass distribution and understate how informative EEoS residual structure is (observed log_r2 sits far above both nulls, p < 0.001).

### 10. Pipeline-integrity diagnostics: dropout and catchability scaling

**Implemented:** A year-by-year and rectangle-by-rectangle funnel tracks haul counts through each processing stage (raw length data → state variables → FishGlob join → EEoS-passing), flagging known dropout-spike years. Separately, the systematic B_pred/B_obs overprediction (median ≈ 4×) is tested for stability as a function of observed biomass magnitude by fitting log(B_obs) ~ log(B_pred) and comparing the ratio across biomass quartiles; no correction is applied.

**Why:** The funnel distinguishes genuine survey-coverage gaps from pipeline errors at each stage rather than leaving an unexplained sample-size discrepancy. The catchability check exists to explicitly test, and reject, a post-hoc scalar correction: because the overprediction ratio grows with biomass magnitude (~3× to ~5.5× across quartiles) rather than being constant, a single multiplier would inflate log_r2 without validating a mechanism, and is instead carried forward as a design constraint for H2 (residual magnitude, not bias-corrected absolute biomass, as the dependent variable).

### 11. Dominance and size-homogeneity diagnostic

**Implemented:** A parallel diagnostic track computes, per haul, Berger–Parker numerical dominance (D = abundance of the single most numerous species ÷ total N) and a population-weighted coefficient of variation of mass within that dominant species' catch (size_CV), using the same length–weight-matched bins as the production N. D is joined against the production N (not independently recomputed) and validated against an independent recount for consistency.

**Why:** This tests a distinct component of the broader hypothesis — whether numerical dominance and within-species size-homogeneity are associated with EEoS prediction error, a pattern consistent with disturbance-driven biotic homogenisation. Tying D to the production N, with an independent validation check, guards against the dominance and biomass figures silently referring to different underlying bin sets.

### 12. Confound checks on the dominance metrics

**Implemented:** Before treating any D/size_CV–residual association as meaningful, four checks are run: a Pearson/Spearman correlation matrix among D, size_CV, N, B_obs and bin count (flagging |r| > 0.4); the same association recomputed within each B_obs quartile rather than marginally; a cross-tabulation of the lowest size_CV decile against known dropout-flagged years and rectangles; and a descriptive-only breakdown of dominant-species identity in the extreme deciles of each metric.

**Why:** Each check targets a specific way the association could be spurious rather than ecological: whether D/size_CV are just re-expressing biomass magnitude or abundance; whether the association survives after conditioning on biomass magnitude; whether it is driven by known bad-data years/rectangles; and whether it is an artefact of one or two heavily-sampled species. This is a pre-specified screening protocol run before, not after, the headline result.

### 13. Partial R² regression follow-up

**Implemented:** Three nested OLS models are fit on the per-haul dominance table: log-residual ratio ~ log(B_obs); ~ log(B_obs) + D + size_CV; and ~ log(B_obs) + D × size_CV. Partial R² of {D, size_CV} beyond log(B_obs) alone is calculated as (R²₁ − R²₀)/(1 − R²₀), and the interaction term is tested via nested-model ANOVA.

**Why:** This converts descriptive binned gradients into a single effect-size statistic that explicitly nets out biomass magnitude, since log(B_obs) is included in the baseline model — the resulting ≈1% partial R² is variance explained beyond what biomass already accounts for, not the raw (confounded) marginal association. The interaction test separately checks whether dominance and size-homogeneity only matter jointly rather than additively.
