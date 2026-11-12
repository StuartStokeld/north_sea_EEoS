# Methods reference for the production audit

Source: methods text supplied by Stuart Stokeld, 12 Nov 2026 (verbatim). **Note:** this text predates the DataType /
haul-duration fix (PR #1): its sample sizes (12,069 H1 hauls; 10,464 H2/H3 hauls) and data-processing description do
not yet include the HH alignment, the exactly-30-min filter or the C → raw conversion. The updated methods drafts on
the fix branch (`display_discussion/methods_draft.md`, `H1_methods_draft.md`) describe those steps. Audit against
both; treat count differences caused by the fix as expected, and report every other mismatch.

> **Decisions after GATE A (Stuart Stokeld, 12 Nov 2026), recorded by the audit.** The decisions box described in
> the GATE A review was not in this file when the audit resumed, so these are transcribed from the review message.
>
> 1. **Residual:** the methods text is correct and the code is in error. The H2/H3 response is the absolute log
>    difference |ln B_obs − ln B_pred|. Tests are written to this specification (they fail on the current code).
> 2. **FP_between:** the code's haul-weighted definition is correct; the methods wording was the error. Revised
>    wording: *FP_between is the haul-weighted mean of ln(hours + 1) over the rectangle-years in which the rectangle
>    has H2/H3 hauls; FP_within is the rectangle-year's deviation from it, so FP_within averages to zero across each
>    rectangle's hauls.* The year-weighted mean over all Couce years is a possible later sensitivity analysis (not run).
> 3. **Downloads:** DATRAS HL CSV 9 April 2026; DATRAS HH from the 9 November 2026 re-download.
> 4. **Interpretations 1–7 at GATE A accepted:** "valid AphiaID" = non-missing; genus/family fallback = separate
>    arithmetic means of a and b; "smallest individual" = midpoint mass of the lightest bin; "log" = natural log;
>    the ≥ 10-H1-haul rule is applied before the Couce rectangle-year filter; the CAR conditional mean is ρ × the
>    sum of neighbouring intercepts (binary adjacency); non-fish taxa = the hard-coded exclusion list.
> 5. **Five figure scripts** with hard-coded pre-fix counts are archived (out of scope).

---

## Data processing

### Trawl Survey Data

We obtained community data from a long term standardised bottom haul survey: the International Bottom Trawl Survey (NS-IBTS), which has collected data annually in the first quarter (January–February) since 1963. Trawl data came from two independently processed products of the same survey (NS-IBTS): the FishGlob cleaned NS-IBTS file, NS-IBTS_clean.RData (Maureaud et al., 2024; github.com/fishglob/FishGlob_data), and the ICES Database of Trawl Surveys (ICES, 2023), specifically the product "Unaggregated trawl and biological information, survey NS-IBTS, Quarter 1, 1985–2015" (downloaded 9 April 2026; datras.ices.dk/Data_products/Download/Download_Data_public.aspx). Both datasets are collected from standardised hauls following consistent protocols for duration, swept area, and geolocation. All analyses used Quarter 1 records: by survey design, this controls for seasonal variation in catchability and species movement, giving a stable analysis panel and a consistent sampling window across years (Maureaud et al., 2024). Each data point represents a haul nested within an ICES statistical rectangle.

DATRAS length-frequency records contain, for each haul, the number of individuals of each species in each length class, together with a length-unit code, taxonomic identifier, and a raising factor for subsampled measurements. Length rows were dropped if they (a) lacked a valid length code or length class, (b) lacked a valid WoRMS AphiaID (numeric code for identifying taxa), (c) had missing or non-positive length counts, or (d) belonged to non-fish taxa (cephalopods, crustaceans and bivalves). Counts were converted to haul-raised numbers using the DATRAS subsampling factor so that length measurements taken on a subsample represent the full catch.

From DATRAS HL we obtain species abundance (N), richness (S) and size structure. From FishGlob we obtain total haul biomass (B_obs), by summing individual mass per haul, and the ICES rectangle ID. The two products were matched on survey, year, quarter, country, ship and haul number. All 12,117 FishGlob hauls matched a DATRAS haul; 272 DATRAS hauls with no FishGlob B_obs were excluded. A further 48 matched hauls were then dropped (44 with missing or non-positive B_obs, and 4 with fewer than two species), leaving 12,069 hauls for the EEoS analysis.

### Fishing pressure data processing

Direct disturbance was approximated using fishing pressure, represented by reconstructed demersal trawling effort from Couce et al. (2020), downloaded from the Cefas Data Hub (Reconstruction of North Sea trawling effort 1985–2015, version 2; DOI 10.14466/CefasDataHub.61). Couce et al. report hours of otter and beam trawling per ICES statistical rectangle per calendar year for 1985–2015. For each rectangle-year we summed otter and beam hours (treating missing gear as zero) to give total trawling hours, yielding a rectangle-year time series of fishing pressure. Both the survey and fishing pressure datasets are therefore aggregated to ICES statistical rectangles.

Decomposing this was necessary for spatial and temporal analysis, specifically to avoid the within-between conflation problem (Bell & Jones, 2015). Two decomposed parameters were calculated. FP_between, the spatial component, is each rectangle's mean fishing pressure of log(hours + 1) over 1985–2015, calculated once per rectangle. FP_within, the temporal component, is that rectangle-year's deviation from the same mean; all hauls in a rectangle-year share one value. In practical terms, FP_between describes how heavily a place is fished in general: a rectangle with high FP_between has been chronically trawled throughout the 30-year record, relative to other rectangles. FP_within describes whether a site in a given year was fished more heavily than average.

We also defined a set of fishing pressure "phases" to allow comparisons across policy eras. Phases were assigned by haul year independently of the FP decomposition, at the 1992 (European Commission, 1992) and 2002 (European Commission, 2002) Common Fisheries reforms and the 2008 Long-Term Management Plan and Marine Strategy Framework Directive (European Union, 2008), giving four periods (1985–1991, 1992–2001, 2002–2007, 2008–2015). All data processing and analysis were carried out in R version 4.6.1 (R Core Team, 2026), except the EEoS function, run via the reference Python implementation called from R using reticulate.

## Derivation of state variables

Species richness (S) and abundance (N) were calculated from the filtered DATRAS length-frequency data, and observed biomass (B_obs), calculated by summing all individual mass per haul, from FishGlob, after the two products were merged. S and N include only taxa that receive a length–weight coefficient below, so they share the same species pool as E.

EEoS also requires community metabolic rate (E), which the survey does not record. DATRAS gives the number of individuals of each species in each length class. Individual mass in each class was estimated from FishBase length–weight coefficients (Froese & Pauly, 2026; accessed May 2026 via rfishbase) using W = aL^b (Froese, 2006), where W is mass and L is the DATRAS length-class midpoint in cm. For each species, FishBase supplies two coefficients: a, which sets mass for a fish of a given length, and b, which sets how steeply mass increases with length. Where a species had no usable FishBase record, genus-mean then family-mean a and b were used (n = 31 taxa; 0.13% of individuals). Remaining unresolved taxa were dropped (n = 32 names; 0.33% of individuals), leaving 177 taxa with coefficients.

Raw metabolic rate was defined as Σ(abundance × mass^0.75) (Harte, 2011; Harte et al., 2022); however, this is incomparable across communities. Therefore, E was normalised by the metabolic rate of the smallest individual in the haul (m_min^0.75), the boundary condition required by the METE partition function underlying EEoS (Harte, 2011).

## Establishing Ecological Equation of State

The ecological equation of state (EEoS; Harte et al., 2022) predicts total community biomass (B) from three state variables: species richness (S), abundance (N) and community metabolic rate (E). It combines the Maximum Entropy Theory of Ecology (METE; Harte, 2011), which derives an ecosystem structure function, R(n, ε | S, N, E), giving the most likely distribution of species abundance (n) and individual metabolic rate (ε), with the metabolic scaling ε ∝ m^(3/4). With units set so the smallest individual has mass and metabolic rate of 1, each individual's mass is ε^(4/3), and biomass is obtained by averaging nε^(4/3) over R and multiplying by S:

B = S Σₙ n ∫ ε^(4/3) R(n, ε | S, N, E) dε

The EEoS has no free parameters, so any mismatch between predicted and observed biomass reflects community structure that S, N and E alone do not capture.

### Predicting biomass from state variables

Predicted biomass (B_pred) was generated with the reference Python EEoS implementation (https://github.com/micbru/equation_of_state). For each haul, the function received S, N and normalised E and returned a dimensionless biomass, which was converted to grams by multiplying by the m_min value used to normalise E and compared against observed biomass (B_obs) on a strict, unfitted 1:1 line (no slope or intercept estimated), reproducing Harte et al. (2022). This tests whether EEoS predicts biomass from state variables alone, with zero free parameters. EEoS residual is defined as the magnitude of failure of the predicted against observed biomass, log(B_obs) compared to log(B_pred), irrespective of sign.

## Modelling spatial and temporal fishing pressure

We tested whether EEoS residuals carry a systematic fishing-pressure signal, spatially (H2) and temporally (H3). This uses a subset of the EEoS haul and rectangle set, requiring the haul's ICES rectangle-year to have Couce trawling-hours coverage and rectangles to contain at least 10 H1 hauls over 1985–2015, yielding 158 analysis rectangles and 10,464 hauls.

### Primary CAR model

The primary model was a generalised linear mixed model with a Gaussian (normal) error distribution and identity link, fitted by restricted maximum likelihood. Spatial dependence among ICES rectangles was included as a conditional autoregressive (CAR) random intercept: each rectangle's intercept, given its neighbours, depends on the mean intercept of neighbouring rectangles (queen contiguity: shared edge or corner). In equation form:

EEoS residual ~ FP_between × phase + FP_within × phase + CAR(rectangle)

An interaction between FP_between and phase, which is anchored to policy reforms at 1992, 2002 and 2008 (1985–1991, 1992–2001, 2002–2007, 2008–2015), as well as FP_within and phase, was included in the model to account for different effects of fishing pressure across different policy eras, and we report the phase-specific slope of haul EEoS residual on FP_between and on FP_within from this model.

## Robustness and sensitivity tests

### Bonferroni correction

The critical outputs for hypothesis testing are the four phase-specific slopes of EEoS residual on FP_between (H2) and the four phase-specific slopes of EEoS residual on FP_within (H3), from the primary model. As such, there is a risk of multiple testing and Type I errors. Therefore, we applied a Bonferroni correction within each hypothesis family (m = 4; α = 0.05/4 = 0.0125), separately, as they test distinct disturbance pathways.

### Rectangle sub-sampling

Stability was assessed by refitting the primary model on 1,000 random subsamples of rectangles, without replacement, at 80%, 60% and 40% of the 158 rectangles. For each phase-specific spatial (H2) and temporal (H3) slope we recorded the proportion of converged subsamples that matched the full-sample sign, and the proportion that matched both sign and significance (raw p < 0.05).

### Spatial permutation

To test whether the H2 signal depends on the spatial layout of fishing, we ran 1,000 permutations that reassigned the observed FP_between values among rectangles, leaving EEoS residual, FP_within, phase, and the CAR neighbour structure unchanged, and refitted the primary CAR model each time.

### Spatial lag

To isolate whether the spatial signal reflects spillover from neighbouring rectangles rather than local pressure alone, a distinct question from whether the signal is present, we refitted the primary CAR model adding a k-nearest-neighbour (k = 4) lag of neighbouring rectangles' FP_between, interacted with phase, and report this alongside the primary local-effect coefficients as evidence on the spatial scale over which fishing pressure operates.
