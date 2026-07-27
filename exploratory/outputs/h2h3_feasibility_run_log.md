# H2/H3 shared hierarchical model — feasibility check run log

FEASIBILITY CHECK ONLY. No fishing-pressure or fishing-pressure x phase coefficient below is interpreted as answering H2 or H3. This task does not finalise the model as the committed approach — that follows supervisor discussion, not this script.

## Package choice and new dependency
Package: glmmTMB (version 1.1.14). Chosen over brms (would need Stan toolchain + MCMC compile time, much slower for a feasibility check) and spaMM (less commonly used/documented in this ecosystem) because glmmTMB supports the required exp() spatial covariance structure over arbitrary point coordinates natively, fits via fast ML/REML (TMB automatic differentiation) rather than MCMC, and integrates with the same glm-style formula interface already familiar from this project's lme4/nlme-adjacent H2 work.
NEW DEPENDENCY: glmmTMB was NOT part of this project's renv-managed set (absent from renv.lock and from the renv project library). Installed ad hoc via install.packages('glmmTMB', repos = 'https://cloud.r-project.org') into the ambient/user R library — same pattern as strucchange for the h2h3_designA4 structural-break task (see that run log). NOT added to renv.lock. A version mismatch warning is emitted at load time (glmmTMB was built against TMB 1.9.19; TMB 1.9.20 was installed as its dependency) — verified with a trivial test model that this does not prevent fitting or produce incorrect results; noted here as an environment fact, not resolved.
CONSEQUENCE FOR HOW THIS SCRIPT MUST BE RUN: because glmmTMB (and, per an ambient-library check run for this task, also spdep/spatialreg/lme4 used by the existing H2 spatial scripts) are not visible under a renv-activated R session, this script must be run WITHOUT renv activation: `Rscript --vanilla pipeline/run_h2h3_shared_model_feasibility.R` from the project root. Flagged as an environment inconsistency already present in the project (run_h2_models.R has the same unstated requirement), not introduced or resolved here.

## Data and analysis universe
Analysis universe: the established H2 rectangle panel (outputs/h2_rectangle_panel.rds), i.e. rectangles with >= 10 hauls AND Couce fishing-pressure coverage — 158 rectangles. No new rectangle-inclusion rule introduced for this task.
Outcome is haul-level (not rectangle-year aggregated), per the brief. 10470 hauls fall in this 158-rectangle universe with a finite `residual`; of these, 6 lack a Couce record for their SPECIFIC rectangle-year (rectangle has overall Couce coverage, but not every year) and are dropped by the inner join to rectangle-year fishing pressure, leaving 10464 hauls for the model.
Fixed effects: log(hours_total + 1) [Couce rectangle-YEAR fishing pressure, same log1p transform as the established H2 log_mean_annual_hours_total convention] x phase [4-level, see below] + mean_ln_B_obs [rectangle-level mean log observed biomass, the established H2 confound, from h2_rectangle_panel.rds — constant within rectangle, does not vary by year].
Phase definition and haul counts (break year = first year of new phase; see build_phase_factor()):
  - 1985-1988 (years 1985-1988): n_hauls = 1583, n_rectangles = 156
  - 1989-2000 (years 1989-2000): n_hauls = 3933, n_rectangles = 156
  - 2001-2007 (years 2001-2007): n_hauls = 2447, n_rectangles = 158
  - 2008-2015 (years 2008-2015): n_hauls = 2501, n_rectangles = 157
Spatial structure: 158 unique rectangle centroids (= number of analysis rectangles; one spatial random-effect level per rectangle, shared by all hauls in that rectangle, matching the H2 SEM's rectangle-level spatial error term in spirit).

## Model fitting
Fitting full spatial model: residual ~ log_hours_total * phase + mean_ln_B_obs + exp(pos + 0 | dummy) ...
Full spatial model fit time: 1.6 sec elapsed (1.2 sec user, 0.2 sec system).
Fitting simpler non-spatial comparison model: residual ~ log_hours_total * phase + mean_ln_B_obs + (1 | stat_rec) ...
Simpler non-spatial model fit time: 0.7 sec elapsed (0.5 sec user, 0.1 sec system).
Practicality: the full spatial model takes 2.3x as long to fit as the simpler comparison model.
Saved fitted model objects (+ analysis data.frame): /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_model_objects.rds

## 1. Convergence diagnostics
  - full_spatial: optimizer_convergence_code = 0 (0 = converged per nlminb), message = 'relative convergence (4)', Hessian positive-definite = TRUE, max|gradient| = 2.994e-06, any NA std. errors = FALSE
    glmmTMB::diagnose() output: Unusually large Z-statistics (|x|>5): |  | disp~(Intercept)  |        -115.8878  |  | Large Z-statistics (estimate/std err) suggest a *possible* failure of | the Wald approximation - often also associated with parameters that are | at or near the edge of their range (e.g. random-effects standard | deviations approaching 0).  (Alternately, they may simply represent | very well-estimated parameters; intercepts of non-centered models may | fall in this category.) While the Wald p-values and standard errors | listed in summary() may be unreliable, profile confidence intervals | (see ?confint.glmmTMB) and likelihood ratio test p-values derived by | comparing models (e.g. ?drop1) are probably still OK.  (Note that the | LRT is conservative when the null value is on the boundary, e.g. a | variance or zero-inflation value of 0 (Self and Liang 1987; Stram and | Lee 1994; Goldman and Whelan 2000); in simple cases these p-values are | approximately twice as large as they should be.) |  | 
  - simple_nonspatial: optimizer_convergence_code = 0 (0 = converged per nlminb), message = 'relative convergence (4)', Hessian positive-definite = TRUE, max|gradient| = 9.756e-04, any NA std. errors = FALSE
    glmmTMB::diagnose() output: Unusually large Z-statistics (|x|>5): |  |   disp~(Intercept) theta_1|stat_rec.1  |         -115.48934          -30.06111  |  | Large Z-statistics (estimate/std err) suggest a *possible* failure of | the Wald approximation - often also associated with parameters that are | at or near the edge of their range (e.g. random-effects standard | deviations approaching 0).  (Alternately, they may simply represent | very well-estimated parameters; intercepts of non-centered models may | fall in this category.) While the Wald p-values and standard errors | listed in summary() may be unreliable, profile confidence intervals | (see ?confint.glmmTMB) and likelihood ratio test p-values derived by | comparing models (e.g. ?drop1) are probably still OK.  (Note that the | LRT is conservative when the null value is on the boundary, e.g. a | variance or zero-inflation value of 0 (Self and Liang 1987; Stram and | Lee 1994; Goldman and Whelan 2000); in simple cases these p-values are | approximately twice as large as they should be.) |  | 
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_convergence_diagnostics.csv

## 2. Variance components
  - full_spatial / rectangle_spatial_sd: estimate = 0.3117, 95% CI = [0.0106, 9.1761], residual SD = 0.4467
  - full_spatial / spatial_range: estimate = 17.7051, 95% CI = [0.0187, 16742.8317], residual SD = 0.4467
  - simple_nonspatial / rectangle_intercept_sd: estimate = 0.1185, 95% CI = [0.1032, 0.1362], residual SD = 0.4471
Boundary check (spatial model): range estimate = 17.705 (degrees, since centroids are in lon/lat); range CI upper bound = 16742.832. Flagged as degenerate if the range is pinned near 0, near/above the maximum inter-centroid distance in the data, or if its CI is many orders of magnitude wide (all checked explicitly below).
Maximum inter-centroid distance in the data (degrees): 13.000
FLAG: spatial range estimate EXCEEDS the maximum inter-centroid distance in the data — the fitted correlation structure is close to a single shared field across all rectangles, not a locally-decaying one. Reported as a degeneracy flag, not resolved.
Simple non-spatial model's rectangle random-intercept SD (0.1185) is not near zero — not flagged as singular.
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_variance_components.csv

## 3. Partial pooling check
Shrinkage measure: shrinkage_ratio = model random intercept / unpooled_intercept, where unpooled_intercept is the per-rectangle mean of (observed residual - the model's OWN population-level [fixed-effects-only, re.form=NA] prediction) — i.e. what each rectangle's intercept would be under NO pooling at all, on the same scale as the random intercept (both are deviations from the SAME fixed-effects prediction). Expected pattern: shrinkage_ratio closer to 0 (heavy shrinkage toward the population mean) for low-haul-count rectangles, closer to 1 (little shrinkage) for high-haul-count rectangles.
Spearman correlation of log(n_hauls) with shrinkage_ratio across all 158 rectangles: full_spatial rho = 0.251; simple_nonspatial rho = 0.338. Both are POSITIVE, consistent with the expected shrinkage pattern (more hauls -> shrinkage_ratio closer to 1), but only weak-to-moderate and noisy at the individual-rectangle level (see sample below — a few rectangles have shrinkage_ratio outside [0, 1], including sign flips). For the full_spatial model this is expected to be noisier than a plain random intercept: the exp() structure pools each rectangle toward a SPATIALLY-WEIGHTED combination of its neighbours, not just the flat population mean, so haul count alone does not fully determine the amount or direction of shrinkage. Reported as observed, not assumed, and not smoothed over.
Sample of rectangles spanning the haul-count range (full spatial model):
  - 47E6: n_hauls = 16, unpooled_intercept = -0.0418, random_intercept = 0.0318, shrinkage_ratio = -0.761
  - 51F0: n_hauls = 53, unpooled_intercept = 0.2195, random_intercept = 0.1336, shrinkage_ratio = 0.609
  - 43F0: n_hauls = 60, unpooled_intercept = -0.0440, random_intercept = -0.0062, shrinkage_ratio = 0.140
  - 33F1: n_hauls = 62, unpooled_intercept = 0.2774, random_intercept = 0.1289, shrinkage_ratio = 0.465
  - 48F1: n_hauls = 65, unpooled_intercept = -0.2016, random_intercept = -0.1650, shrinkage_ratio = 0.818
  - 42F2: n_hauls = 67, unpooled_intercept = 0.0107, random_intercept = 0.0109, shrinkage_ratio = 1.019
  - 43E8: n_hauls = 69, unpooled_intercept = 0.1716, random_intercept = 0.1033, shrinkage_ratio = 0.602
  - 40F7: n_hauls = 74, unpooled_intercept = -0.1408, random_intercept = -0.2560, shrinkage_ratio = 1.819
  - 33F3: n_hauls = 82, unpooled_intercept = 0.0411, random_intercept = -0.0760, shrinkage_ratio = -1.851
  - 38F7: n_hauls = 111, unpooled_intercept = -0.1957, random_intercept = -0.3164, shrinkage_ratio = 1.617
Same rectangles, simple non-spatial model:
  - 47E6: n_hauls = 16, unpooled_intercept = -0.0418, random_intercept = -0.0175, shrinkage_ratio = 0.419
  - 51F0: n_hauls = 53, unpooled_intercept = 0.2195, random_intercept = 0.1787, shrinkage_ratio = 0.814
  - 43F0: n_hauls = 60, unpooled_intercept = -0.0440, random_intercept = -0.0153, shrinkage_ratio = 0.347
  - 33F1: n_hauls = 62, unpooled_intercept = 0.2774, random_intercept = 0.2367, shrinkage_ratio = 0.853
  - 48F1: n_hauls = 65, unpooled_intercept = -0.2016, random_intercept = -0.1432, shrinkage_ratio = 0.710
  - 42F2: n_hauls = 67, unpooled_intercept = 0.0107, random_intercept = 0.0279, shrinkage_ratio = 2.610
  - 43E8: n_hauls = 69, unpooled_intercept = 0.1716, random_intercept = 0.1505, shrinkage_ratio = 0.877
  - 40F7: n_hauls = 74, unpooled_intercept = -0.1408, random_intercept = -0.1399, shrinkage_ratio = 0.994
  - 33F3: n_hauls = 82, unpooled_intercept = 0.0411, random_intercept = 0.0008, shrinkage_ratio = 0.019
  - 38F7: n_hauls = 111, unpooled_intercept = -0.1957, random_intercept = -0.1921, shrinkage_ratio = 0.982
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_partial_pooling.csv
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_feasibility_partial_pooling.png

## 4. Fixed effects: full spatial vs simpler non-spatial comparison
Descriptive only — NOT a formal model-selection test, and NOT an interpretation of the fishing-pressure or interaction terms as answering H2/H3. Reporting how much fixed-effect estimates and SEs shift with vs without the spatial term:
  - (Intercept): estimate 0.6135 -> 1.1656 (+90.0%); SE 0.2164 -> 0.3395 (+56.9%)
  - log_hours_total: estimate -0.0149 -> -0.0107 (+27.8%); SE 0.0139 -> 0.0137 (-1.7%)
  - phase1989-2000: estimate 0.0583 -> 0.0309 (-47.1%); SE 0.1377 -> 0.1371 (-0.5%)
  - phase2001-2007: estimate -0.7240 -> -0.7322 (-1.1%); SE 0.1360 -> 0.1349 (-0.8%)
  - phase2008-2015: estimate -0.3408 -> -0.3433 (-0.8%); SE 0.1330 -> 0.1318 (-1.0%)
  - mean_ln_B_obs: estimate -0.1518 -> -0.1974 (-30.0%); SE 0.0136 -> 0.0152 (+11.7%)
  - log_hours_total:phase1989-2000: estimate -0.0102 -> -0.0078 (+23.2%); SE 0.0146 -> 0.0146 (-0.5%)
  - log_hours_total:phase2001-2007: estimate 0.0872 -> 0.0877 (+0.6%); SE 0.0146 -> 0.0145 (-0.8%)
  - log_hours_total:phase2008-2015: estimate 0.0579 -> 0.0581 (+0.4%); SE 0.0145 -> 0.0143 (-0.9%)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_fixed_effects.csv

## Feasibility verdict
Optimizer-level convergence (both models converge per nlminb, positive-definite Hessian): TRUE.
Spatial-parameter identifiability flag (full spatial model's range estimate EXCEEDS the maximum inter-centroid distance, and/or its CI spans > 4 orders of magnitude): TRUE.
PLAIN STATEMENT (per the brief — feasible / feasible only in simplified form / not feasible as specified): 
  FEASIBLE ONLY IN SIMPLIFIED FORM. The full model as specified (rectangle intercept + exp() spatial correlation over centroids) reaches a numerically converged optimum (optimizer code 0, positive-definite Hessian, finite standard errors) and does NOT throw a hard error or singular-fit warning — so it is not a flat non-convergence failure. However its spatial covariance parameters are weakly identified at this rectangle density/extent: the range estimate (17.7 degrees) exceeds the maximum inter-centroid distance in the data (13.0 degrees), and both the range and sill 95% CIs span several orders of magnitude (range CI [0.0187, 16742.8]; sill CI [0.0106, 9.18]). This is consistent with a likelihood surface that is nearly flat in the range direction — the data (158 point locations spanning 13 degrees, many clustered) do not contain enough spatial contrast to pin down a continuous exponential-decay range parameter, in contrast to the original H2 SEM which used a coarser, fixed queen-adjacency spatial structure (no range parameter to estimate). The plain (1|stat_rec) comparison model, by contrast, is fully well-behaved (rectangle-intercept SD 95% CI is narrow: [0.103, 0.136]). RECOMMENDATION for supervisor discussion: either (a) fit the spatial term with a FIXED or externally-constrained range/adjacency structure (e.g. reusing the H2 SEM's queen-contiguity weights rather than a freely-estimated continuous range), or (b) proceed with the plain (1|stat_rec) random intercept, which converges cleanly with well-identified variance components, as the practical fallback.
This verdict is about model FIT BEHAVIOUR only — it does not interpret, and is not based on, the sign, magnitude, or significance of the fishing-pressure or fishing-pressure x phase coefficients.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_fixed_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_variance_components.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_partial_pooling.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_convergence_diagnostics.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_feasibility_partial_pooling.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_run_log.md (this file)
