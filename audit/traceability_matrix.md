# Traceability matrix — methods text vs live pipeline (GATE A)

Branch `test/production-audit` (from `fix/haul-effort-datatype` @ `f0589c2`). Audit date 2026-11-12.
**No pipeline code changed.** Evidence scripts were run from a scratch directory against the committed
outputs (`outputs/*.rds`, exactly-30-min primary run).

Methods sources:
- **MR** = `audit/methods_reference.md` (supplied, written before the fix)
- **MD** = `display_discussion/methods_draft.md`, **H1D** = `display_discussion/H1_methods_draft.md`,
  **HD** = `display_discussion/H2_H3_methods_detailed.md` (updated drafts on the fix branch)

Status key: **MATCH** · **MISMATCH** (code ≠ text, or two texts contradict each other) · **NOT IMPL** ·
**UNDOC** (code does something no methods text states) · **INTERP** (I had to interpret the text; please
confirm). Test IDs are *planned* (written after GATE A). Prefixes: `C-` contract, `U-` unit, `S-` solver,
`I-` integrity, `M-` metamorphic, `V-` simulation, `D-` docs-vs-outputs.

---

## A. Data selection and survey scope

| # | Methods statement (quoted) | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| A1 | MR: "All analyses used Quarter 1 records … 1985–2015" | `pipeline/R/datras_constants.R:14-16`; `datras_hl_helpers.R:84-86`; `datras_hh_helpers.R:66-70`; `h1_join_helpers.R:27-30` | Q1 and 1985–2015 filters applied to HL, HH and FishGlob | MATCH | C-HL-01, C-FG-01 |
| A2 | MR: "since 1963" / "(January–February)" | — | Descriptive only | n/a | — |
| A3 | MR: "Unaggregated trawl and biological information … downloaded 9 April 2026" | `outputs/datras_hl_raw.rds` (first committed 2026-07-13, mtime 2026-07-30); `outputs/datras_hh_raw.rds` (`hh_source = "ICES CSV"`, written 2026-11-10) | The HL download date is not recorded anywhere in the repo. HH was imported on 10 Nov from an ICES CSV whose download date isn't recorded either | **NEEDS INPUT** (lead 16) | I-PROV-03 |
| A4 | MD/H1D: HL joined to HH; keep `HaulVal == V`, DataType C/R/P/S, `HaulDur == 30` | `datras_hh_helpers.R:161-231`, `datras_constants.R:31-34` | As stated. Hauls are classified into reasons, in order: no_hh 0 · ambiguous_hh_key 9 · invalid 33 · datatype 0 · duration 1,941 · kept 10,406 | MATCH | U-ALIGN-01..05, I-ACC-01 |
| A5 | MD: "restricted to valid hauls of exactly 30 min (DataType C, R or P)"; H1D §1b: "C, R, P or S" | `datras_constants.R:33` | Keeps C/R/P/S (HH has no S hauls) | **MISMATCH (text vs text)**: MD leaves out S. No effect on results | D-01 |
| A6 | *(not in any text)* | `datras_hh_helpers.R:183-201` | 9 HL hauls whose 6-field key matches 2 HH records (1985 NO 58EJ, GOV + HOB gear, haul numbers 40–48) are excluded as `ambiguous_hh_key` | **UNDOC** | I-ACC-02 |

## B. HL length-row processing

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| B1 | MR: rows dropped if "(a) lacked a valid length code or length class" | `datras_hl_helpers.R:91,104` | Drops `LngtClass` NA, and codes outside {".","0","1","5"} (midpoint NA). In the aligned data 0 rows are dropped by either rule | MATCH. Codes "2" and "9" (+ group), both valid ICES codes, would be dropped silently as "invalid"; neither occurs in the data | U-LNGT-05, C-HL-02 |
| B2 | MR: "(b) lacked a valid WoRMS AphiaID" | `datras_hl_helpers.R:87` | Only checks `!is.na(Valid_Aphia)`. Validity isn't checked against WoRMS | **INTERP** ("valid" = non-missing) | C-HL-03 |
| B3 | MR: "(c) had missing or non-positive length counts" | `datras_hl_helpers.R:89-90,106` | `HLNoAtLngt` NA or ≤ 0 dropped; raised count ≤ 0 dropped | MATCH (0 rows dropped) | U-HL-02 |
| B4 | MR: "(d) belonged to non-fish taxa (cephalopods, crustaceans and bivalves)" | `datras_constants.R:4-12` | Hard-coded list of 47 AphiaIDs (138802 appears twice). It drops 31,614 rows. *Ostrea edulis* (140658) and *Todarodes* (138281) are not on the list; they are dropped later only because they have no LW record | MATCH in effect; list incomplete (minor) | C-HL-04, U-HL-03 |
| B5 | MR: "Counts were converted to haul-raised numbers using the DATRAS subsampling factor" | `datras_hh_helpers.R:224-231`, `datras_hl_helpers.R:100` | `n = HLNoAtLngt × SubFactor × f`, with `f = HaulDur/60` for C and 1 otherwise | MATCH (MR doesn't mention the C conversion; MD/H1D do) | U-ALIGN-03, M-CR-01 |
| B6 | MR: "L is the DATRAS length-class midpoint in cm" | `datras_hl_helpers.R:13-25` (`lngt_midpoint_cm`) | "." → (cls+0.5)/10 ✔; "1" → cls+0.5 ✔; **"0" → (cls+0.5)/10 ✘ (should be (cls+2.5)/10)**; **"5" → (cls+0.5)/10 ✘ (should be cls+2.5; class is in cm)** | **MISMATCH**: leads 1–2 CONFIRMED (section E) | U-LNGT-01..04 |
| B7 | H1D §2: coefficients "applied per 1 cm length bin" | `datras_hl_helpers.R:114-131` | Applied per DATRAS bin: 1 mm (".") / 5 mm ("0") / 1 cm ("1") / 5 cm ("5") | **MISMATCH (text)** | D-02 |
| B8 | *(not in any text)* | `datras_hl_helpers.R:114-131` (`inner_join` to lookup) | **30 taxa in HL are missing from the LW lookup and are silently dropped**: 22,161 rows, 0.32% of individuals, at least one species lost in 4,382 of 10,128 H1 hauls (median 1 species, 8% of S). 28 of the 30 are accepted fish species in WoRMS (e.g. *Argentina sphyraena*, *Pomatoschistus minutus*, *Callionymus maculatus*, *Syngnathus acus*, *Alosa fallax*, *Mustelus asterias*, *Sebastes norvegicus*); the other 2 are Squalidae (family) and *Dipturus flossada* (taxon inquirendum). Likely cause (unconfirmed): the lookup's species list came from FishGlob's accepted names, which lump these taxa to genus | **UNDOC** (new finding; probably major because it biases S) | I-ACC-04, C-LW-03 |

## C. Length–weight lookup

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| C1 | MR: "FishBase length–weight coefficients (Froese & Pauly, 2026; accessed May 2026 via rfishbase) using W = aL^b" | `supplementary/fishbase_lookup.Rmd` (not in the live run order) → `outputs/fishbase_lw_lookup_v2.csv` | One record per species: `Type` priority TL > FL > SL > NA > other, then largest `Number`. The Rmd uses an undefined `all_spp` and hard-coded paths, so the lookup can't be rebuilt from the repo. rfishbase/FishBase version not recorded | MATCH in substance; **provenance gap** | I-PROV-02, C-LW-01 |
| C2 | MR: "genus-mean then family-mean a and b were used (n = 31 taxa; 0.13% of individuals)" | Rmd lines 153-203 | 26 genus + 5 family = 31 ✔. a and b are averaged separately and arithmetically, over TL **and NA-type** records. 0.133% reproduces from the lookup's `total_N` (pre-fix counts); on the current aligned data it is 0.108% | MATCH (counts); **INTERP** ("genus-mean" = arithmetic mean of a and of b separately) | U-LW-02 |
| C3 | MR: "Remaining unresolved taxa were dropped (n = 32 names; 0.33% of individuals), leaving 177 taxa" | lookup `fallback_level == "MISSING"` | 32 ✔; 177 finite a,b ✔; 0.326% ✔ (pre-fix `total_N`). This doesn't cover the 30 taxa missing from the lookup altogether (B8) | MATCH as stated; **incomplete** | C-LW-02, D-03 |
| C4 | MR: "W is mass" (units implied g); "L … in cm" | `datras_hl_helpers.R:127` | `mass_g = a·L_cm^b` | MATCH | U-LW-01 |
| C5 | *(not in any text)* (lead 11) | lookup `Type` | 140 TL; 1 FL (*Leptoclinus maculatus*, flagged "FL_needs_conversion" but used unconverted); 2 NG and 1 OT ("other_review"); 2 NA-type. Those present in the current data = 5 taxa, 0.0001% of individuals | **UNDOC** (minor) | C-LW-04 |
| C6 | Brief: LW `b` in ~2.5–3.5 | lookup | 7 records fall outside: *Aphia minuta* 3.55, *Scyliorhinus canicula* 3.615, *Maurolicus muelleri* 3.961, *Lumpenus lampretaeformis* 2.09, *Lesueurigobius friesii* 2.146, *Leptoclinus* 2.16, *Artediellus* 2.29 | contract flag (not necessarily wrong) | C-LW-05 |

## D. State variables, m_min, E

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| D1 | MR: "S and N include only taxa that receive a length–weight coefficient … same species pool as E" | `datras_hl_helpers.R:151-170` | S, N, E_raw, m_min all come from the same `hl_mass` rows (after the LW join and mass > 0) | MATCH (lead 10) | U-STATE-01, M-POOL-01 |
| D2 | MR: "Raw metabolic rate … Σ(abundance × mass^0.75)" | `datras_hl_helpers.R:159` | `sum(NumberAtLength * mass_g^0.75)` | MATCH | U-STATE-02 |
| D3 | MR: "E was normalised by the metabolic rate of the smallest individual in the haul (m_min^0.75)" | `datras_hl_helpers.R:160`, `datras_state_helpers.R:99-106` | `m_min` = smallest **bin-midpoint mass** across all species-bins in the haul; `E = E_raw / m_min^0.75` | **INTERP** ("smallest individual" = midpoint mass of the lightest length bin; individuals below a bin midpoint aren't represented). Units g ✔ (lead 12) | U-STATE-03 |
| D4 | *(not in MR; H1D §3 states it)* | `datras_hl_helpers.R:164-168`; `datras_state_helpers.R:83-96`; `h1_join_helpers.R:60-73` | Keeps hauls with `N ≥ S`, `E_raw` and `m_min` finite and > 0 | MATCH (H1D) / **UNDOC** (MR) | U-STATE-04 |
| D5 | *(not in any text)* | aggregated N | **3,770 of 10,128 H1 hauls have non-integer N** (C counts × 0.5, fractional SubFactor) and are passed as-is to the solver, which uses `np.arange(N)` | **UNDOC** (minor; affects solver semantics) | S-EEOS-05 |

## E. FishGlob join and B_obs

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| E1 | MR: "matched on survey, year, quarter, country, ship and haul number" | `h1_join_helpers.R:40-46` | Key built from haul_id parts 1–5 and 8 (drops gear and StNo). No duplicate keys in FishGlob Q1 1985–2015 | MATCH | U-JOIN-01, I-JOIN-01 |
| E2 | MR: "B_obs … by summing individual mass per haul" (also "Derivation" paragraph) | `h1_join_helpers.R:34`; FishGlob `get_datras.R:364-378,414` | `B_obs = Σ wgt × 1000` (g); FishGlob `wgt = CatCatchWgt/1000` = **reported species catch weight** (C × HaulDur/60), with FishGlob's own ×100 correction for hauls whose LW-estimated/reported ratio is > 20. **Not** a sum of individual masses | **MISMATCH (text)**: lead 8 CONFIRMED. Units kg→g ✔ | U-JOIN-02, D-04 |
| E3 | MR: "44 with missing or non-positive B_obs" | `h1_join_helpers.R:34` (`na.rm = TRUE`) | A haul whose wgt is all NA → B_obs = 0 → excluded ✔. **But 1,149 H1 hauls have some species with NA wgt, which counts as 0: B_obs is underestimated and the haul is kept** | **UNDOC** (probably major for H1) | U-JOIN-03, C-FG-03 |
| E4 | *(not in any text)* | — | B_obs covers every FishGlob fish taxon (including the 30 missing from the lookup and the 32 without coefficients); S/N/E cover only the 177 with coefficients. The two sides use different species pools | **UNDOC** | M-POOL-02 |
| E5 | MR: "All 12,117 FishGlob hauls matched … 272 DATRAS hauls … excluded … 48 … (44 + 4) … 12,069" | `build_eeos_predictions.R` → `h1_dropout_summary.csv`, `h1_filter_exclusions.csv` | Pre-fix archive reproduces exactly: 12,117 / 12,389 DATRAS / 12,117 matched / 272 / 44 bad_B_obs + 4 S<2 / 12,069 | MATCH (lead 17) | D-05 |
| E6 | MD/H1D: 10,128 hauls (C 8,493 / R 1,543 / P 92) | same | Current: 12,117 FishGlob / 10,406 DATRAS / 10,172 matched (1,945 FishGlob hauls unmatched, 234 DATRAS hauls without FishGlob) / 44 dropped (42 bad_B_obs + 2 S<2) / 10,128 | MATCH. The current funnel isn't written up in MR, and "All FishGlob hauls matched" is no longer true | D-05 |
| E7 | MD/CLAUDE: DataType P hauls kept; "local FishGlob build includes these" | `FishGlob_data/` (not a git repo) | `NS-IBTS_clean.RData`: `timestamp = "2021-07"` (DATRAS extraction July 2021), file mtime 19 Mar 2026; local `NEWS.md` = 2.1.0; local `get_datras.R:353` **excludes** P, so the RData was not produced by the local code version | **NEEDS INPUT** (lead 9) | I-PROV-04 |

## F. EEoS prediction

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| F1 | MR: "reference Python EEoS implementation (github.com/micbru/equation_of_state)" | `build_eeos_predictions.R:130-142` | Local `equation_of_state/biomass.py` is **byte-identical to upstream `a96f406`**. Golden: all 42 `data_statevariables.csv` communities reproduce `pBnum` to 8e-16 | MATCH (lead 16, solver part) | S-EEOS-01 |
| F2 | MR: "returned a dimensionless biomass, which was converted to grams by multiplying by the m_min" | `build_eeos_predictions.R:154` | `B_pred = B_pred_norm × m_min` | MATCH | U-EEOS-01 |
| F3 | H1D: hauls failing "S < 2, E ≤ N, missing/non-finite inputs, B_obs ≤ 0" excluded before the solver | `h1_join_helpers.R:111-139` | As stated (E ≤ N: 0 hauls) | MATCH (H1D) / **UNDOC** in MR (E > N) | U-JOIN-04 |
| F4 | *(not in any text)* | `biomass.py:85-91`; `build_eeos_predictions.R:147-193` | `fsolve` non-convergence and the >5% constraint failure only `print()`. They're never raised, so the R `tryCatch` can't catch them. **5 of 10,128 H1 hauls return B_pred from a failed solve** (4 with ~97–99% constraint error) | **UNDOC** (finding) | S-EEOS-03 |
| F5 | MR: "compared … on a strict, unfitted 1:1 line"; "log_r2 = 1 − SS_res/SS_tot" (CLAUDE/H1D) | `h1_common.R:57-70` | As stated | MATCH. Lead 4: `log_r2()` is correct; the test premise is wrong (section G) | U-MET-01 |
| F6 | *(not in any text)* | `build_eeos_predictions.R:155-157` | `haul_key = row$haul_key` is NULL (the joined column is named `join_key`), so `haul_eeos_predictions.rds` **has no `haul_key` column**. Downstream code uses `haul_id`, so nothing breaks today | **UNDOC** (minor) | I-SCHEMA-01 |

## G. Fishing pressure

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| G1 | MR: "summed otter and beam hours (treating missing gear as zero)" | `h2_couce_helpers.R:25-53` | `Data_type == "Reconstructed"`, gear ∈ {Otter, Beam}, pivot with sum, `coalesce(…, 0)` | MATCH (lead 7, part) | U-FP-01 |
| G2 | MR: "mean fishing pressure of log(hours + 1)" | `h2h3_feasibility_helpers.R:65` | `log()` = **natural log** everywhere; no log10 used for FP (the only log10 is in axis labelling) | MATCH; **INTERP** ("log" = ln; suggest writing "ln") (lead 6) | U-FP-02 |
| G3 | MR: "FP_between … each rectangle's mean … over 1985–2015, calculated once per rectangle" | `h2h3_within_between_helpers.R:12-25` on haul rows from `build_feasibility_data()` | **Haul-weighted mean over the rectangle's hauls**, so only years with hauls (and Couce coverage), weighted by haul count. Not a mean over 1985–2015 rectangle-years. MD says "haul-weighted mean", which matches the code | **MISMATCH (MR vs code; MR vs MD)**: lead 7 CONFIRMED. Size: \|Δ FP_between\| vs a mean over all Couce years has median 0.083, max 0.59 (SD of FP_between 1.17); r = 0.991 | U-FP-03, M-FP-01 |
| G4 | MR: "FP_within … that rectangle-year's deviation from the same mean; all hauls in a rectangle-year share one value" | same | Deviation from the haul-weighted mean ✔ shared per rectangle-year ✔. Its **haul-weighted** mean is 0 (1e-15); its rectangle-year mean is not (max 0.47) | MATCH relative to G3's definition | U-FP-04 |
| G5 | MR: rectangle-years without Couce records | `h2h3_feasibility_helpers.R:50-56` | Hauls in such rectangle-years are **dropped** (inner join), not treated as 0. 3 of 4,898 rectangle-years are missing; no zero-hour rows exist | MATCH with "requiring the haul's ICES rectangle-year to have Couce trawling-hours coverage" | U-FP-05 |
| G6 | MR: phases "1985–1991, 1992–2001, 2002–2007, 2008–2015" | `run_h2h3_phase_v2_refit.R:140-148` | `cut(breaks = c(-Inf,1991,2001,2007,Inf), right = TRUE)`; year ranges in data verified inclusive | MATCH | U-PHASE-01 |

## H. H2/H3 model

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| H1 | MR: "EEoS residual is … the magnitude of failure … log(B_obs) compared to log(B_pred), irrespective of sign"; H1D §7: "residual magnitude … as the dependent variable" | `run_h2h3_phase_v2_reporting.R:179`; `run_h2h3_phase_v2_refit.R:205`; `run_h2h3_within_between.R:196,228` | Response is the **signed** `residual = log(B_obs) − log(B_pred)`. MD says "log(B_obs) − log(B_pred) throughout", which matches the code | **MISMATCH (MR + H1D vs code + MD)**: lead 5 CONFIRMED. Not deciding which is right | D-06, M-RESID-01 |
| H2 | MR: "rectangles to contain at least 10 H1 hauls … 158 rectangles" | `h2_common.R:5`; `h2_panel_helpers.R:43-58`; `h2h3_feasibility_helpers.R:40-56` | ≥ 10 counted on **all H1 hauls before** the rectangle-year Couce filter; the rectangle must have some Couce record. After filtering, min hauls/rectangle = 13 | MATCH; **INTERP** (order: count, then rectangle-year filter) (lead 14) | U-INC-01 |
| H3 | MR: "Gaussian … identity link, fitted by restricted maximum likelihood" | `run_h2h3_phase_v2_reporting.R:187-189` | `spaMM::fitme(..., method = "REML")`, gaussian default | MATCH (lead 13) | I-MODEL-01 |
| H4 | MR: "CAR … queen contiguity (shared edge or corner)" | adjacency from `outputs/h2h3_feasibility_round2_model_objects.rds$adjMatrix` | 158×158, binary, symmetric, 0 islands, 1 connected component, 3–8 neighbours; **identical to queen contiguity derived independently from the rectangle codes** | MATCH. **Provenance:** this file is written only by `exploratory/pipeline/run_h2h3_shared_model_feasibility_round2.R` (5 live scripts read it) | C-ADJ-01..03, I-PROV-01 |
| H5 | MR: "each rectangle's intercept, given its neighbours, depends on the **mean** intercept of neighbouring rectangles" | spaMM `adjacency()` with binary W | Proper CAR on **binary** W: conditional mean = ρ·**Σ** of neighbours (precision ∝ I − ρW), not the mean (that would be row-standardised W) | **INTERP / wording mismatch** (please confirm) | — |
| H6 | MR: formula with FP_between × phase + FP_within × phase + CAR | `run_h2h3_phase_v2_reporting.R:179` | Exact match | MATCH | I-MODEL-02 |
| H7 | MR: "report the phase-specific slope" | `h2h3_within_between_helpers.R:197-247` | main + interaction contrast; SE by delta method; **p = 2·pnorm(−\|t\|) (Wald normal approximation)** | MATCH; p-value method **UNDOC** | U-SLOPE-01 |

## I. Robustness

| # | Methods statement | Code location | What the code does | Status | Test |
|---|---|---|---|---|---|
| I1 | MR: "Bonferroni … within each hypothesis family (m = 4; α = 0.0125)"; MD: "across the four phase-specific H2 slopes" | `run_h2h3_h2_multiplicity_subsampling.R:5,201` | `pmin(1, p × 4)` separately for the H2 and H3 families | MATCH (MR); MD mentions H2 only (**text vs text**) | U-BONF-01 |
| I2 | MR: "1,000 random subsamples … 80%, 60% and 40% … without replacement … proportion of converged subsamples … sign … sign and significance (raw p < 0.05)" | `run_h2h3_rectangle_subsampling_car_refit.R:113-117,163,310-356` | As stated; seed 42 + per-iteration seeds; isolates in subsampled adjacency logged but not excluded | MATCH (lead 15) | U-SUB-01, V-* |
| I3 | MR: "1,000 permutations that reassigned the observed FP_between values among rectangles, leaving … FP_within, phase, and the CAR neighbour structure unchanged" | `permutation_bootstrap_FP_between_CAR.R:109-110,162-175,217` | Permutes at rectangle level (stops if FP_between isn't unique per rectangle); N_BOOT 1000; seed 42; p = (1 + #\|null\| ≥ \|obs\|)/(n + 1), two-sided | MATCH; p formula **UNDOC** | U-PERM-01 |
| I4 | MR: "k-nearest-neighbour (k = 4) lag of neighbouring rectangles' FP_between, interacted with phase" | `build_knn_spatial_weights.R`; `h2h3_knn_spatial_helpers.R:170-182`; `run_h2h3_spec_a_car_identifiability.R:174` | Lag = mean of the 4 nearest neighbours (row-standardised, asymmetric); `FP_between_lag × phase_v2`. The Spec A permutation recomputes the lag after each shuffle | MATCH | U-LAG-01 |

## J. Environment

| # | Methods statement | Code | What the code does | Status | Test |
|---|---|---|---|---|---|
| J1 | MR: "R version 4.6.1" | — | R 4.6.1 ✔; Python 3.12.13; spaMM 4.7.17; glmmTMB 1.1.14 | MATCH | I-ENV-01 |

---

## K. Code behaviour no methods text mentions (summary)

1. **Taxa missing from the LW lookup silently dropped** (B8): 30 taxa, 4,382 of 10,128 H1 hauls lose S.
2. **Partial-NA FishGlob weights counted as 0** in B_obs (E3): 1,149 H1 hauls.
3. **Different species pools** for B_obs and S/N/E (E4).
4. **Silent EEoS solver failures** (F4): 5 hauls.
5. Effort-alignment exclusions not in MR: 9 ambiguous HH keys, 33 invalid, 1,941 off-duration (A4, A6).
6. Filters `N ≥ S`, `E > N` (D4, F3); non-integer N passed to the solver (D5).
7. `m_min` defined from bin midpoints (D3).
8. LW record choice and unconverted FL/NG/OT records (C1, C5); separate arithmetic means of a and b (C2).
9. Wald-z p-values and the permutation p formula (H7, I3).
10. CAR adjacency read from an **exploratory** output (H4); LW lookup built outside the live pipeline (C1).
11. Missing `haul_key` column in `haul_eeos_predictions.rds` (F6).
12. **Hard-coded pre-fix counts** that will now `stop()`: `fig_methods_state_and_fishing_pressure.R:169`, `fig_study_area_intro.R:41`, `run_methods_study_area_visuals.R:47` (12,069); `run_h2h3_phase_v2_rectangle_maps.R:156`, `run_h2_presentation_primary_figure.R:145` (10,464). None are in the CLAUDE.md run order or on the stale list. I need to know whether they count as live.
13. Valid ICES LngtCodes "2" and "9" would be dropped as unknown (B1); non-fish list has a duplicate entry and misses 2 present invertebrates (B4).

---

## L. Section 3 lead verdicts (preliminary; full evidence goes into `audit/findings.md` at GATE B)

| Lead | Verdict | One-line evidence |
|---|---|---|
| 1 LngtCode "0" | **CONFIRMED** | ICES ref 18: "0.5 cm length class, reporting units: mm". 99.999% of "0" classes are multiples of 5 mm. The code adds 0.5 mm, not 2.5 mm |
| 2 LngtCode "5" | **CONFIRMED** | ICES ref 18: "5 cm length class, reporting units: cm". All 1,338 rows are 60–145, multiples of 5; in the same species × haul the code-"1" rows reach at most 59 cm, so the "5" rows continue from 60 cm upwards |
| 3 Rows dropped by code/NA filters | **NOT A BUG** (0 rows) · related new finding B8 | All codes in the aligned data are known; no NA classes |
| 4 `test-h1_metrics.R:11` | **NOT A BUG in `log_r2()`; the test is wrong** | By hand: SS_tot = 2(ln 10)² = 10.60, SS_res = 3(ln 4)² = 5.77, so R² = 0.456. A ×4 offset only gives negative R² when it exceeds the spread of log B_obs (testthat isn't installed, so not yet run) |
| 5 Residual sign | **CONFIRMED (methods conflict)** | Code and MD: signed. MR and H1D §7: magnitude |
| 6 log base | **NOT A BUG** (natural log throughout) · INTERP of "log" | — |
| 7 FP_between averaging | **CONFIRMED (MR vs code)** | Haul-weighted, only years with hauls; r = 0.991 with the all-years mean; Couce filters ✔ |
| 8 B_obs definition | **CONFIRMED (text)** + E3 | Reported `CatCatchWgt`, kg → g ✔; partial NA treated as 0 |
| 9 FishGlob P / version | **NEEDS INPUT** | No git metadata; RData timestamp 2021-07; local code excludes P |
| 10 Species pool / counts | **NOT A BUG** (S, N, E, m_min share rows; 31/32/177 reproduced) · % are pre-fix; B8 gap | — |
| 11 Length type | **CONFIRMED (minor)** | 5 non-TL taxa used unconverted, 0.0001% of N |
| 12 m_min units | **NOT A BUG** (g end-to-end) · INTERP D3 | — |
| 13 REML / queen | **NOT A BUG** · provenance flag; H5 wording | — |
| 14 Inclusion rules | **NOT A BUG** · INTERP H2 | min 13 hauls/rectangle |
| 15 Robustness specs | **NOT A BUG** | all specs match |
| 16 Provenance | **PARTIAL**: biomass.py = a96f406 ✔, golden ✔; silent solver failures (F4); **NEEDS INPUT** on download dates | — |
| 17 Counts | **CONFIRMED reproduced** (pre-fix and post-fix) | E5, E6 |

### Leads 1–2 measured on the H1 sample (exactly 30 min, 10,128 hauls)

Recomputed with ICES-correct midpoints ("0" → (cls+2.5)/10, "5" → cls+2.5) and the unchanged solver:

| Group | Hauls | S, N changed | median E_raw ratio | m_min changed | E (normalised) ratio range | B_pred ratio median [5%, 95%] (min–max) | hauls with \|ΔB_pred\| > 10% |
|---|---|---|---|---|---|---|---|
| has code "5" | 538 | 0 | 1.048 (max 7.57) | 429 (median ×1.84) | 0.050 – 7.57 | 1.04 [0.79, 1.61] (0.54 – 15.1) | 254 |
| code "0", no "5" | 7,314 | 0 | 1.004 | 2,150 (median ×1.09) | 0.49 – 1.07 | 1.00 [0.997, 1.04] (0.98 – 1.13) | 4 |
| neither | 2,276 | 0 | 1 | 0 | 1 | 1 | 0 |

H1 headline effect: log_r2 0.6610 → 0.6565; cor² 0.9713 → 0.9733; median B_pred/B_obs 1.958 → 1.968;
correlation between current and corrected residuals 0.980. **The H1 headline barely moves, but 254 hauls
(mostly NL 1985–2001 with large cod) change by >10%, up to 15×.** Those hauls are spatially and temporally
clustered, so the H2/H3 effect needs a refit; that's for GATE B. Note: the log_r2 here (0.661) is the
current exactly-30-min value. CLAUDE.md's −0.223 is the superseded pre-fix value.
