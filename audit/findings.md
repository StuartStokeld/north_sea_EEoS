# Production audit: findings (GATE B)

Branch `test/production-audit` (from `fix/haul-effort-datatype` @ f0589c2). No pipeline code was changed. All
effects below were computed in audit scripts (`audit/scripts/`) against the committed outputs, with the live
pipeline functions swapped or re-parameterised in memory only.

- **Severity.** **critical** changes reported results. **major** could change results, or violates the methods
  or provenance. **minor** is robustness or clarity. **note** is for information.
- **Test tags.** Tests tagged `[F-..]` in `pipeline/tests/testthat/` encode a finding and fail until it is fixed
  (see "Test suite status").
- **Decisions after GATE A** are recorded at the top of `audit/methods_reference.md`.

## 1. Known-lead verdicts (brief §3)

| Lead | Verdict | Evidence | Finding |
|---|---|---|---|
| 1 LngtCode "0" midpoint | **CONFIRMED** | ICES vocab ref 18: "0" = 0.5 cm class in mm. Code adds 0.5 mm, not 2.5 mm. 190,117 rows, 99.999% multiples of 5 mm | F-A1 |
| 2 LngtCode "5" | **CONFIRMED** | ICES: "5" = 5 cm class in cm. All 1,338 rows are 60–145 cm; code reads them as 6–14.5 cm | F-A2 |
| 3 Rows dropped by unit/code filters | **NOT A BUG** (0 rows). Related: F-B8 | — | — |
| 4 `test-h1_metrics.R:11` | **Test wrong, `log_r2()` correct** (decision). By hand R² = 0.456; test corrected | F-A4 | F-A4 |
| 5 Residual sign | **CONFIRMED: code in error** (decision 1). Models use signed ln B_obs − ln B_pred; spec is absolute | F-A5 | F-A5 |
| 6 log base | **NOT A BUG** (natural log throughout). Methods wording "log" → "ln" | F-A6 | F-A6 |
| 7 FP_between averaging | **Methods wording was wrong** (decision 2); the code's haul-weighted mean is kept. FP_within averages to 0 over each rectangle's hauls (tested) | F-A7 | F-A7 |
| 8 B_obs definition | **NOT A BUG in the text; reverses my GATE A verdict.** FishGlob's published `wgt` is length-based (Σ count × a·L^b, l.1110), not `CatCatchWgt`; it matches our recomputation within 1% for 99.4% of haul × taxon pairs. kg → g correct. But see F-B9 and F-B10 | F-B10 |
| 9 FishGlob version | **RESOLVED.** RData SHA-256 `bd98d2ea…dad13` = git blob `8fea3f9b` = FishGlob commit `acde3dc` (2025-12-08), in tags v2.0.2 and v2.1 (Zenodo 10.5281/zenodo.19056965). The 92 P hauls exist because FishGlob extracted DATRAS in July 2021, before they were resubmitted as P | F-A9 |
| 10 Species pool | **NOT A BUG** (S, N, E and m_min share one row set; 31/32/177 reproduce). Gap: F-B8 | — |
| 11 Length type | **CONFIRMED (minor)** | F-A11 |
| 12 m_min units | **NOT A BUG** (g end to end). Wording F-A12 | — |
| 13 REML / queen CAR | **NOT A BUG** (REML; adjacency = queen contiguity, symmetric, connected; tested). Wording F-A13; provenance F-B13 | — |
| 14 Inclusion rules | **NOT A BUG.** ≥ 10 H1 hauls counted before the Couce filter (minimum 13 after); wording F-A14 | — |
| 15 Robustness specs | **NOT A BUG** (1,000 draws at 80/60/40%; rectangle-level permutation; p = (1+#)/(n+1); ×4 Bonferroni; kNN k = 4). Undocumented: Wald-z p-values | F-A15 |
| 16 Provenance (solver, downloads) | `biomass.py` = upstream a96f406 (golden 42/42). **Silent solver failures CONFIRMED** (F-A16). **Mixed DATRAS download dates** (F-A17) | F-A16, F-A17 |
| 17 Methods counts | **REPRODUCED.** Pre-fix 12,117 / 272 / 44 + 4 / 12,069 / 158 / 10,464; post-fix 10,406 / 10,172 / 10,128 / 8,670 / 158 (tested) | — |

## 2. Findings table

Severity order: critical, then major, minor, note. "Effect" is on the current primary results (exactly 30 min;
10,128 H1 hauls; 8,670 H2/H3 hauls, 158 rectangles).

| ID | Severity | Location | What's wrong | Evidence | Effect on results | Suggested fix |
|---|---|---|---|---|---|---|
| **F-A5** | **critical** | `run_h2h3_phase_v2_refit.R:205`, `run_h2h3_phase_v2_reporting.R:179` (and every robustness script using `residual ~`) | H2/H3 response is the signed residual; spec (decision 1) is \|ln B_obs − ln B_pred\| | Formulas use `residual`; `primary_model_v2.rds` response ≠ abs (8,390 / 8,670 differ in sign) | 96.7% of residuals are negative, so absolute slopes ≈ negated signed slopes. **H2 1992–2001 FP_between becomes Bonferroni-significant** (−0.030, p_bonf 0.032 vs signed +0.026, p_bonf 0.18). H3 Bonferroni pattern unchanged (92–01, 08–15 significant; 02–07 drops to raw only, p_bonf 0.064). The absolute model's Gaussian fit is acceptable (model-residual skew 0.31, excess kurtosis 1.5; QQ r = 0.996; floor at 0 visible): `audit/outputs/abs_residual/`, `audit/figures/abs_residual_diagnostics.png` | Use `abs_residual` in the primary and all robustness models; rerun H2/H3. Fix branch after the residual decision is final |
| **F-A2** | **critical** | `pipeline/R/datras_hl_helpers.R:21-22` | LngtCode "5" (5 cm classes, cm) read as mm: 60–145 cm fish become 6–14.5 cm | ICES ref 18; 1,338 rows, mostly cod, NL 1985–2001 | 538 H1 hauls: E ratio 0.05–7.6, B_pred ratio 0.54–15×, 249 hauls change > 10%. **With F-A1 fixed too:** H1 log R² 0.661 → 0.657, cor² 0.971 → 0.973. **H2 abs FP_between 1992–2001 loses Bonferroni significance** (p_bonf 0.032 → 0.098). Signed FP_between 92–01 and 08–15 lose raw significance. H3 unchanged in significance (`audit/outputs/lngtcode/`) | "5": midpoint = class + 2.5 cm |
| **F-A1** | **major** | same | LngtCode "0" (0.5 cm classes, mm): midpoint +0.5 mm instead of +2.5 mm | 190,117 rows; 7,314 hauls have "0" but no "5" | B_pred ratio 0.98–1.13 in those hauls; 4 change > 10% (jointly with F-A2 above) | "0": midpoint = (class + 2.5)/10 cm |
| **F-B10** | **major** (interpretation) | FishGlob `get_datras.R:1110`; `h1_join_helpers.R:33` | B_obs is FishGlob's length-based biomass (FishGlob LW table, lower-bound lengths) from the **same DATRAS length data** as S, N and E. The products are not independent | 99.4% within 1% of our recomputation | H1 agreement is partly built in. With independent measured weight (×100-corrected, complete hauls) cor² is 0.92 vs 0.97 and log R² 0.42 vs 0.62 (`diagnosis_measured_weight.md`) | Decide which B_obs (see F-B16); state the shared length basis in the methods |
| **F-B16** | **major** (decision input) | DATRAS HL `SpeciesCategoryWeight` | The independent measured weight has 65% complete coverage (GB-SCT 10%, NL 16%) and whole country-year blocks recorded about 100× too small (DE 1990–2003, DK 1985–99, GB 1985–2003, NO 1989–2002, SE 1985–2001) | 63/163 country-years entirely at ≈ 0.009 × B_obs | As recorded: unusable (log R² −1.36; spurious large slopes). After a block-level ×100 correction: median ratio to FishGlob 0.87. H2: no FP_between slope significant. H3: 02–07 and 08–15 remain Bonferroni-significant with the same signs | Your decision; see `diagnosis_measured_weight.md`. Confirm the unit blocks with ICES before any use |
| **F-B9** | **major** (minor in effect) | `h1_join_helpers.R:33` (`sum(wgt, na.rm = TRUE)`) | NA FishGlob taxon weights count as 0 g | 1,407 haul × taxon records, 1,149 hauls (11.3%). Cause: FishGlob's LW join on taxa + family + genus + rank fails for 6 taxa (e.g. *Lumpenus*: Lumpenidae vs Stichaeidae), 1 taxon name missing, and 511 weight-only records (no lengths) | 76% of affected hauls lose ≤ 0.1% of biomass; 94 lose > 1%; max 25%. Associated with FP (Spearman 0.42 at rectangle level), but no treatment moves any slope by more than 0.011 (`diagnosis_missing_weights.md`) | Treatment (d), fill from lengths × FishGlob LW by name; report the M1a defect to FishGlob. Your choice of treatment |
| **F-B8** | **major** | `datras_hl_helpers.R:115` (`inner_join` to the lookup); `outputs/fishbase_lw_lookup_v2.csv` | Fish taxa without LW coefficients are silently dropped from S, N, E and m_min | 53 taxa in H1 hauls (29 absent from the lookup, 24 in it with NA a/b), 0.39% of N. 4,902 of 10,128 hauls lose ≥ 1 species (median 8.3% of S, max 33%). Mostly accepted species: *Argentina sphyraena*, *Syngnathus acus*, *Callionymus maculatus*, *Pomatoschistus minutus* | Rebuilt lookup resolves 51/53 (only *Ostrea edulis* and *Todarodes*, both non-fish, remain). Adding them: S +1 in 4,902 hauls; H1 log R² 0.661 → 0.658, median B_pred/B_obs 1.958 → 1.965; 2 extra failed solves. H2/H3 not refitted (deferred to the fix branch) | Adopt the re-runnable `audit/scripts/rebuild_lw_lookup.R` (draft) keyed on DATRAS AphiaID; add the missing invertebrates to NON_FISH (F-B12) |
| **F-A16** | **major** | `equation_of_state/biomass.py:85-91`; `build_eeos_predictions.R:147-193` | Solver failures are only printed, never raised; `biomass.py`'s check uses signed `con > thresh`, so negative deviations aren't even printed | 5 H1 hauls have \|constraint dev\| 15–99% (e.g. S = 2, N = 189,168). One violates the Jensen lower bound B ≥ N(E/N)^{4/3} | 5 of 10,128 hauls carry meaningless B_pred (residuals included in H1 and H2/H3) | Compute \|dev\| after each solve in R; exclude > 5% and report the count. Don't edit upstream `biomass.py` |
| **F-B13** | **major** (provenance) | `run_h2h3_within_between.R:80`, `_reporting.R:69`, `_rectangle_subsampling_car_refit.R`, `_spec_a_car_identifiability.R`, `permutation_bootstrap_FP_between_CAR.R:81`, `run_h2h3_phase_v2_refit.R:59` | Live scripts read 7 files that only exploratory scripts write: the CAR adjacency (`h2h3_feasibility_round2_model_objects.rds`), the kNN lag (F-B15), `primary_model_v2_spec_a.rds` (optional), `h2h3_results_fp_slopes_by_phase.csv`, `spatial_diagnostics_comparison.csv` (optional), `exploratory/outputs/permutation_bootstrap_FP_between_RE_objects.rds`. The LW lookup is built in a non-runnable notebook | `audit/dependency_graph.md` (static parse) | Adjacency verified correct (queen, tested), so no current effect. But rebuilding from the repo is impossible, and a rectangle-set change would break it | Move adjacency, kNN-lag and LW-lookup construction into live `pipeline/` scripts |
| **F-B15** | **major** | `run_h2h3_spec_a_car_identifiability.R:47,148` | Spec A (spatial lag) joins a **stale pre-fix** FP_between_lag (`outputs/fp_between_lag_rectangle.rds`, 4 Aug 2026, untracked, exploratory) | Stored FP_between differs from current by up to 0.42; lag up to 0.30 (r = 0.997) | Lag slopes shift, e.g. 2008–15 −0.044 (p 0.084) → −0.032 (p 0.19); absolute 2008–15 lag p 0.057 → 0.14. No 0.05 crossings (`audit/outputs/spec_a/`) | Recompute the lag inside the live script from current FP_between |
| F-A11 | minor | lookup `fallback_level` | 1 FL, 1 NG/OT, 2 untyped records used as TL; 7 b outside 2.5–3.5 | 0.0001% of N | negligible | Convert FL/SL with FishBase length–length tables or fall back to genus |
| F-A17 | minor | `audit/sources.md` | Three DATRAS snapshots: FishGlob (July 2021), HL (9 Apr 2026), HH (9 Nov 2026) | HL 1,372,057 rows vs Nov CSV 1,376,471 | HL→HH match 100%; P-haul DataType differs between snapshots | Document all three dates in the methods; re-download HL with HH at the next refresh |
| F-B11 | minor | `datras_csv_import.R` (HL import) | The committed HL drops SpeciesSex / SpeciesCategory, so duplicated rows can't be told from sex/category splits | 4,097 exact duplicate rows; the raw CSV has 0 exact duplicates and only the documented 58EJ collision | none now (sums are correct) | Keep Sex and CatIdentifier in the HL RDS |
| F-B12 | minor | `datras_constants.R` (NON_FISH_APHIAIDS) | Duplicate entry (138802); *Ostrea edulis* (140658) and *Todarodes* (138281) missing | Dropped later for lack of LW; would enter S once F-B8 is fixed | none now | Add both; dedupe |
| F-B14 | minor | `build_eeos_predictions.R:155` | `haul_eeos_predictions.rds` has no `haul_key` (`row$haul_key` is NULL) | test | none (joins use haul_id) | Keep `haul_key` |
| F-A10 | minor | `supplementary/fishbase_lookup.Rmd` | Genus/family fallback averages a and b separately (arithmetic) over a pool of North Sea species only | 0.13% of individuals | negligible | Note in methods; the rebuild script adds a FishBase-wide pool tier |
| F-B17 | minor | `build_eeos_predictions.R` | Non-integer N (C counts × 0.5, fractional SubFactor) passed to a solver that uses `np.arange(N)` | 3,770 H1 hauls | not quantified | Document, or round N for the solver only |
| F-B18 | minor | `run_h1_null_model.R` | Draws random numbers without `set.seed` | static check | null-model figures not reproducible | Add a seed (determinism layer deferred) |
| F-A4 | note | `pipeline/tests/testthat/test-h1_metrics.R:11` | Test expected negative R² for a ×4 offset | hand calculation 0.456 | none | Test corrected in this branch |
| F-A6 | note | methods text | "log(hours + 1)" is natural log | code | none | Write "ln" |
| F-A7 | note | methods text (now revised) | FP_between is haul-weighted over years with hauls | decision 2 | none (year-weighted version deferred as sensitivity) | — |
| F-A9 | note | `audit/sources.md` | FishGlob version traced (above) | blob match | none | Cite FishGlob v2.1 / commit acde3dc |
| F-A12–A14 | note | methods text | m_min = midpoint mass of the lightest bin; CAR conditional mean is ρ × **sum** of neighbours (binary W); ≥ 10-haul rule applied before the Couce filter; "valid AphiaID" = non-missing | interpretations accepted | none | Reword |
| F-A15 | note | `h2h3_within_between_helpers.R:229-242` | p-values are Wald-normal (z) | code | none | Document |
| F-B19 | note | `methods_draft.md`, `H1_methods_draft.md` | DataType "C, R or P" (S also allowed); LW "per 1 cm bin" (wrong for "." / "0"); Bonferroni described for H2 only | text | none | Edit drafts |
| F-B20 | note | `lngt_midpoint_cm` | Valid ICES codes "2" and "9" would be dropped silently (none present) | ICES ref 18 | none now | Handle "2" (+1 cm); keep "9" excluded, but warn |

## 3. Interpreted statements (accepted at GATE A)

All seven were accepted. The minor ones are logged above as F-A10, F-B12 and F-A12–A14.

## 4. Scripts out of scope (untested)

- **Stale** (per `outputs/live_pipeline_run_log.md`): `run_h1_presentation_figures.R`,
  `fig_state_variable_distributions.R`, `explore_h1_catchability_scaling.R`, `explore_h1_dominance_partial_r2.R`,
  `run_h2h3_car_reporting_diagnostics.R`, `explore_h1_haul_dominance.R`.
- **Archived after GATE A:** `fig_methods_state_and_fishing_pressure.R`, `fig_study_area_intro.R`,
  `run_methods_study_area_visuals.R`, `run_h2h3_phase_v2_rectangle_maps.R`,
  `run_h2_presentation_primary_figure.R`. They hard-code pre-fix counts. **They were untracked**, so `git mv`
  wasn't possible; I moved them with `mv` to `archive/figure_scripts/` (still untracked, nothing committed).
  `CLAUDE.md` never listed them, so its run order is unchanged.

## 5. Test suite status

`Rscript pipeline/tests/run_all.R` runs 109 tests in about 2.5 minutes. 2 slow tests are skipped by default
(`EEOS_SLOW_TESTS=1`).

| Layer | File(s) | Status |
|---|---|---|
| 4.2 contracts | `test-contracts-datras.R`, `test-contracts-inputs.R` | done |
| 4.3 hand-computed units | `test-units-hand.R`; `test-h1_metrics.R` and `test-datras_hl_helpers.R` corrected | done |
| 4.4 solver | `test-eeos-solver.R` (golden 42/42 passes) | done |
| 4.5 integrity | `test-integrity.R` (accounting, joins, provenance) + `audit/dependency_graph.md/.csv` | done; determinism **deferred** (scope) |
| 4.6 metamorphic | `test-metamorphic.R` (C↔R, SubFactor, row order) | done (scoped) |
| 4.7 simulation | — | **deferred**: run on the corrected (absolute-residual) model after the fixes |
| 4.8 mutation | `audit/scripts/mutation_test.R`, `audit/mutation_report.md` | done (4 mutants) |
| 4.9 docs vs outputs | — | **deferred** until after the fixes and the results rewrite |
| 4.10 static checks | — | **deferred** (also: `lintr` isn't installed) |

The tests that fail on the current code all carry a finding tag:

- F-A1, F-A2: lengths
- F-A5: residual (model object and both script formulas)
- F-A11: LW length type
- F-A16: solver failures, `biomass.py` check, invariant
- F-B8: missing LW
- F-B9: NA weights (contract and unit)
- F-B11: HL fields
- F-B12: non-fish list
- F-B13: provenance (2 tests)
- F-B14: `haul_key`

`run_all.R` reports any untagged failure separately as a regression.

## 6. Follow-ups logged (not done)

- **WoRMS ↔ FishBase synonym resolution** in `rebuild_lw_lookup.R`. The `rfishbase::synonyms()` step was
  dropped by scope decision. Names now match on WoRMS valid or verbatim names. WoRMS records are cached
  (212/212).
- **Rebuilt lookup vs current lookup.** Of the 177 current taxa, 157 are present in the rebuild; 110 have
  identical a and b; 25 change W(30 cm) by more than 10% (FishBase 26.06 vs the May 2026 pull, plus choice of
  record). FishBase 25.04 couldn't be fetched (S3 timeout), so the exact original release is unconfirmed.
- **H2/H3 refits** for F-B8 and the full combined fix set, on the fix branches.
- **Year-weighted FP_between** sensitivity (decision 2).
