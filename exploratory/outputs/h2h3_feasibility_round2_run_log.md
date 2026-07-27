# H2/H3 shared hierarchical model — feasibility check ROUND 2 run log
## Adjacency-based (CAR) spatial structure

FEASIBILITY CHECK ONLY (Round 2). No fishing-pressure or fishing-pressure x phase coefficient below is interpreted as answering H2 or H3. This task does not finalise the model as the committed approach — that follows supervisor discussion, not this script.
Parallel in structure to the Round 1 run log (outputs/h2h3_feasibility_run_log.md) for direct side-by-side comparison; only the spatial random-effect term changes.

## Adjacency/weights object: reuse confirmation
No serialized nb/listw object exists anywhere in this repository (checked: outputs/*.rds, outputs/*listw*, outputs/*weight* — none found). Every existing H2 spatial script (run_h2_sar_lag_models.R, run_h2_biomass_fishing_sar_diagnostics.R, run_h2_models.R) calls R/h2_spatial_helpers.R::build_h2_spatial_weights(panel, rectangles_sf, project_root) fresh each time it is needed. That function is a DETERMINISTIC pure function of (panel$stat_rec, the ICES rectangle shapefile at gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp): spdep::poly2nb(., queen = TRUE) then spdep::nb2listw(nb, style = "W", zero.policy = TRUE). This script calls that EXACT SAME function, unmodified, on the EXACT SAME panel object (outputs/h2_rectangle_panel.rds, the established 158-rectangle H2 universe) that the original H2 SEM/SAR analysis used — so the resulting neighbour list is PROVABLY identical (same code path, same inputs), not an independently-built approximation. No rebuild-and-match step was needed because there was nothing separate to match against other than re-deriving via the identical function.
Resulting neighbour list: n_rectangles = 158, n_isolated (zero-neighbour rectangles) = 0.
Cross-check against outputs/h2_spatial_diagnostics.csv (recorded from the original H2 SEM run using this same weights object): n_rectangles = 158, n_isolated = 0. MATCH: TRUE.
Built the raw binary (0/1) adjacency matrix from this SAME reused neighbour list (weights$nb) via spdep::nb2mat(nb, style = "B") — 158 rectangles, 546 undirected edges (1092 matrix entries = 1), symmetric: TRUE. This differs from the row-STANDARDISED "W"-style listw used by the original H2 errorsarlm()/lagsarlm() SEM/SAR models only in MATRIX REPRESENTATION (binary vs. row-normalised), which is the standard convention difference between CAR models (binary adjacency, correlation strength captured by a free rho) and SAR/SEM models (row-standardised weights) — not a different neighbour-list DEFINITION.

## Data and analysis universe (identical to Round 1)
10464 hauls across 158 rectangles — same universe, same fixed-effect construction as Round 1 (see that run log for the full breakdown of the H2-panel/Couce-year join).

## Reusing Round 1's plain non-spatial comparison model
Loaded fit_plain (residual ~ log_hours_total * phase + mean_ln_B_obs + (1 | stat_rec)) directly from /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_model_objects.rds — NOT refit here, per the brief.

## Model fitting
Fitting CAR model: residual ~ log_hours_total * phase + mean_ln_B_obs + adjacency(1 | stat_rec) via spaMM::fitme(..., method = "REML") ...
CAR model fit time: 0.94 sec elapsed (0.71 sec user, 0.04 sec system).
Fit-time comparison (all reused times parsed directly from outputs/h2h3_feasibility_run_log.md, not refit): CAR adjacency model (spaMM, this script) = 0.94 sec elapsed; Round 1 plain (1 | stat_rec) model (glmmTMB) = 0.70 sec elapsed; Round 1 full continuous-spatial exp() model (glmmTMB) = 1.60 sec elapsed. All three are fast (single-digit seconds); the CAR model here is comparable to or faster than both Round 1 models, consistent with spaMM's non-MCMC REML fitting on this modestly-sized (158 x 158) adjacency matrix — fit time is not a practical obstacle to iterating on this model during further development.
Saved fitted model objects (+ analysis data.frame + adjacency matrix): /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_model_objects.rds

## 1. Convergence diagnostics
spaMM's fitme() does not expose an iterative-optimizer status code the way glmmTMB does (no MCMC either — REML/ML via TMB-free internal likelihood maximisation). Reporting what IS available: (a) any warnings/errors captured during the fit, (b) a robustness check refitting from 3 different starting values for rho spanning the admissible range [-0.2634, 0.1317] and confirming they converge to the same estimate (a standard local-optimum check), and (c) whether all fixed-effect standard errors are finite.
  - car_adjacency: warnings during fit = 0 (none); refits from near-lower/zero/near-upper starting rho all converge to rho = 0.130881 (robust to starting value = TRUE); any non-finite fixed-effect SE = FALSE.
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_convergence_diagnostics.csv

## 2. Spatial autocorrelation parameter (rho) identification
rho estimate = 0.130881. Admissible range for this adjacency matrix (from its eigenvalues, 1/lambda_min to 1/lambda_max) = [-0.263426, 0.131705]. 95% parametric bootstrap CI (percentile method, 99 replicates, seed = 123) = [0.082114, 0.131576]. Rho lies 99.8% of the way from the lower bound to the upper bound of the admissible range (i.e. 0.2% of the range's width from the UPPER bound).
CAR variance component (lambda, variance of the rectangle-level CAR random effect) = 0.006476.
FLAG: rho is PINNED NEAR A BOUNDARY of its admissible range (within 1% of the upper or lower bound: estimate 0.130881 vs. upper bound 0.131705 — 99.79% of the way to the boundary). This is a different FORM of non-identifiability than Round 1's (there, the continuous range parameter's CI spanned several orders of magnitude and exceeded the data's spatial extent; here, the discrete adjacency structure DOES converge to a single well-defined point estimate and a comparatively narrow bootstrap CI, but that point estimate sits essentially at the edge of what the model algebraically allows). Reported as observed, not resolved.

## 3. Partial pooling check
Same shrinkage_ratio definition as Round 1: random_intercept / unpooled_intercept, where unpooled_intercept is the per-rectangle mean of (observed residual - the model's OWN population-level [fixed-effects-only, re.form = NA] prediction).
Spearman correlation of log(n_hauls) with shrinkage_ratio across all 158 rectangles: car_adjacency rho = 0.316; simple_nonspatial (reused Round 1 values) rho = 0.338.
Sample of rectangles spanning the haul-count range (CAR adjacency model):
  - 47E6: n_hauls = 16, unpooled_intercept = -0.0596, random_intercept = -0.0062, shrinkage_ratio = 0.103
  - 51F2: n_hauls = 53, unpooled_intercept = 0.0902, random_intercept = 0.0576, shrinkage_ratio = 0.639
  - 49F2: n_hauls = 60, unpooled_intercept = -0.0997, random_intercept = -0.0869, shrinkage_ratio = 0.872
  - 43F5: n_hauls = 62, unpooled_intercept = 0.0112, random_intercept = 0.0184, shrinkage_ratio = 1.647
  - 44F3: n_hauls = 65, unpooled_intercept = 0.0393, random_intercept = 0.0227, shrinkage_ratio = 0.577
  - 45F1: n_hauls = 67, unpooled_intercept = 0.0220, random_intercept = -0.0047, shrinkage_ratio = -0.215
  - 46F0: n_hauls = 69, unpooled_intercept = -0.0122, random_intercept = -0.0282, shrinkage_ratio = 2.307
  - 40F7: n_hauls = 74, unpooled_intercept = -0.2164, random_intercept = -0.1880, shrinkage_ratio = 0.869
  - 39F4: n_hauls = 82, unpooled_intercept = -0.1054, random_intercept = -0.0928, shrinkage_ratio = 0.881
  - 38F7: n_hauls = 111, unpooled_intercept = -0.2638, random_intercept = -0.2548, shrinkage_ratio = 0.966
Same rectangles, simple non-spatial model (reused from Round 1):
  - 47E6: n_hauls = 16, unpooled_intercept = -0.0418, random_intercept = -0.0175, shrinkage_ratio = 0.419
  - 51F2: n_hauls = 53, unpooled_intercept = 0.1353, random_intercept = 0.0926, shrinkage_ratio = 0.684
  - 49F2: n_hauls = 60, unpooled_intercept = -0.1001, random_intercept = -0.0575, shrinkage_ratio = 0.575
  - 43F5: n_hauls = 62, unpooled_intercept = 0.0544, random_intercept = 0.0394, shrinkage_ratio = 0.724
  - 44F3: n_hauls = 65, unpooled_intercept = 0.0718, random_intercept = 0.0648, shrinkage_ratio = 0.903
  - 45F1: n_hauls = 67, unpooled_intercept = 0.0669, random_intercept = 0.0414, shrinkage_ratio = 0.619
  - 46F0: n_hauls = 69, unpooled_intercept = 0.0495, random_intercept = 0.0151, shrinkage_ratio = 0.304
  - 40F7: n_hauls = 74, unpooled_intercept = -0.1408, random_intercept = -0.1399, shrinkage_ratio = 0.994
  - 39F4: n_hauls = 82, unpooled_intercept = -0.0497, random_intercept = -0.0583, shrinkage_ratio = 1.173
  - 38F7: n_hauls = 111, unpooled_intercept = -0.1957, random_intercept = -0.1921, shrinkage_ratio = 0.982
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_partial_pooling.csv
Saved figure: /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_feasibility_round2_partial_pooling.png

## 4. Fixed effects: CAR adjacency model vs plain non-spatial comparison (reused from Round 1)
Descriptive only — NOT a formal model-selection test, and NOT an interpretation of the fishing-pressure or interaction terms as answering H2/H3. Reporting how much fixed-effect estimates and SEs shift between the plain non-spatial model (Round 1, reused) and this CAR adjacency model:
  - (Intercept): estimate 0.6135 -> 0.7431 (+21.1%); SE 0.2164 -> 0.2188 (+1.1%)
  - log_hours_total: estimate -0.0149 -> -0.0125 (+15.8%); SE 0.0139 -> 0.0138 (-1.1%)
  - phase1989-2000: estimate 0.0583 -> 0.0435 (-25.4%); SE 0.1377 -> 0.1373 (-0.3%)
  - phase2001-2007: estimate -0.7240 -> -0.7280 (-0.6%); SE 0.1360 -> 0.1353 (-0.5%)
  - phase2008-2015: estimate -0.3408 -> -0.3382 (+0.8%); SE 0.1330 -> 0.1323 (-0.5%)
  - mean_ln_B_obs: estimate -0.1518 -> -0.1616 (-6.4%); SE 0.0136 -> 0.0149 (+9.6%)
  - log_hours_total:phase1989-2000: estimate -0.0102 -> -0.0089 (+12.2%); SE 0.0146 -> 0.0146 (-0.3%)
  - log_hours_total:phase2001-2007: estimate 0.0872 -> 0.0875 (+0.3%); SE 0.0146 -> 0.0146 (-0.5%)
  - log_hours_total:phase2008-2015: estimate 0.0579 -> 0.0575 (-0.6%); SE 0.0145 -> 0.0144 (-0.5%)
Saved: /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_fixed_effects.csv

## Feasibility verdict
Fit-level behaviour (no warnings/errors during fitting, robust to 3 different starting values for rho, all fixed-effect SEs finite): TRUE.
Rho boundary-pinning flag: TRUE (rho = 0.130881 is 99.79% of the way to the nearer admissible bound out of [-0.263426, 0.131705]).
PLAIN STATEMENT (per the brief — feasible / feasible with caveats / not feasible):
  FEASIBLE WITH CAVEATS. Unlike Round 1's continuous distance-decay model, the CAR adjacency structure fits cleanly by every check available in this framework: no warnings or errors during fitting, the optimum is reached from 3 widely-spaced starting values for rho (near the lower admissible bound, zero, near the upper admissible bound) with no sensitivity to starting value, and all fixed-effect standard errors are finite. This is a materially better-behaved fit than Round 1's exp() structure, which had range/sill CIs spanning orders of magnitude. HOWEVER, the fitted rho (0.130881) sits essentially AT the upper edge of what this adjacency matrix algebraically permits (admissible range [-0.263426, 0.131705]; rho is 99.79% of the way to that upper bound) — a boundary-pinned rho, even with a narrow bootstrap CI, means the model is pushing as much positive spatial correlation into the CAR term as the adjacency structure allows, which is itself a form of identification strain (in the CAR literature this is a well-known feature of rectangle-density adjacency structures with strong, spatially widespread residual correlation — the model wants MORE spatial smoothing than a single free rho against this particular graph can supply). RECOMMENDATION for supervisor discussion: proceed with this CAR structure as the practical choice between the two spatial approaches tried (it is unambiguously better-identified than Round 1's continuous exp() structure), but flag the boundary-adjacent rho as a caveat — e.g. consider whether a Leroux-type CAR (mixing pure CAR with an unstructured component) or the plain (1|stat_rec) model (which has no such boundary issue at all) is preferred if this sensitivity matters for the eventual H2/H3 model.
This verdict is about model FIT BEHAVIOUR only — it does not interpret, and is not based on, the sign, magnitude, or significance of the fishing-pressure or fishing-pressure x phase coefficients.

## Outputs
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_fixed_effects.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_spatial_param.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_partial_pooling.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_convergence_diagnostics.csv
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_model_objects.rds
- /Users/stuartstokeld/north_sea_eeos/outputs/figures/h2h3_feasibility_round2_partial_pooling.png
- /Users/stuartstokeld/north_sea_eeos/outputs/h2h3_feasibility_round2_run_log.md (this file)
