# H1 Pipeline Audit — Checks for Cursor

This document identifies every point in the pipeline where a unit error, column misuse, or assumption mismatch could silently corrupt state variables or the EEoS comparison. Work through each check in order. For each one, confirm the current implementation matches the expected behaviour, or flag and fix.

---

## 1. N — abundance input

**Check:** What column is used as N in the EEoS call?

- **Required:** N must be the total raised individual count per haul from DATRAS HL — specifically `sum(HLNoAtLngt)` across all fish-only species in the haul, after LngtCode unit standardisation and non-fish exclusion.
- **Do not use:** `num_cpue` from FishGlob. This is individuals per hour of trawling (effort-standardised), not a count. It will make N vary with tow duration independently of the actual community, and will be dimensionally inconsistent with E which is built from DATRAS raised counts.
- **Do not use:** `num` from FishGlob. This is total individuals per haul from FishGlob's own pipeline. It is not the same as the DATRAS raised count used to compute E, and mixing sources breaks the internal consistency of S, N, E.

**The N used in the EEoS call must come from the same DATRAS HL rows used to compute E.** S and N should be derived as by-products of the E computation step, not pulled separately from FishGlob.

**Audit action:** Print the source column and a sample of N values alongside `sum(HLNoAtLngt)` per haul from the DATRAS pipeline. Confirm they are identical.

---

## 2. S — species richness input

**Check:** What is S, and is it derived from the same filtered DATRAS HL rows as N and E?

- **Required:** S = number of unique AphiaIDs per haul in the DATRAS HL dataset, after non-fish exclusion, after `HLNoAtLngt > 0` filter, and after LW coverage filter (species with no usable LW parameters cannot contribute to E and should not contribute to S).
- **Do not use:** Species richness from FishGlob (`n_species` or similar). FishGlob's species list is not filtered to the same taxa, the same exclusion criteria, or the same raise structure.

**Audit action:** Confirm S is computed as `n_distinct(AphiaID)` within the same grouped DATRAS HL data used for E, after all filters.

---

## 3. E — total metabolic rate computation

### 3a. LngtCode unit standardisation

**Check:** Are all length bin midpoints converted to centimetres before applying W = aL^b?

Known LngtCode values and their units (confirmed from ICES vocabulary server):
- Code `1`: bins in cm, LngtClass in cm → midpoint = LngtClass + 0.5 (cm)
- Code `0`: bins in 0.5cm, LngtClass in mm → midpoint = (LngtClass + 0.5) / 10 (cm)
- Code `.`: bins in 1mm, LngtClass in mm → midpoint = (LngtClass + 0.5) / 10 (cm)
- Code `5`: bins in cm, LngtClass in cm (confirmed as cm × 10 in the June session, i.e. LngtClass is in mm) → midpoint = (LngtClass + 0.5) / 10 (cm)

**Audit action:** Print unique LngtCode values in the DATRAS HL data. Confirm the conversion applied to each. Print a sample of raw LngtClass vs computed midpoint_cm for each code. Flag any LngtCode not in the list above.

### 3b. FishBase LW parameters — units

**Check:** Are FishBase `a` and `b` values in the correct units for the length midpoints?

FishBase LW parameters use the convention W(g) = a × L(cm)^b. The `a` parameter for a typical demersal fish at TL is in the range 0.001–0.05. Values outside 0.0001–0.5 are suspicious and likely indicate a unit mismatch (e.g., parameters in mm rather than cm, or SL/FL rather than TL).

**Audit action:** From `fishbase_lw_lookup_v2.csv`, print the distribution of `a` values (min, p5, median, p95, max). Flag any species where `a < 0.0001` or `a > 1`. These should be inspected manually — they may have been pulled with incorrect length type or from a non-standard parameterisation.

### 3c. Length type — TL throughout

**Check:** Are FishBase LW parameters filtered to TL (total length) type only?

All NS-IBTS fish species are measured as total length. `LenMeasType` in DATRAS is 86% NA and cannot be used to confirm this per-row; instead the design decision is to filter FishBase LW parameters to TL-type entries as the primary source, with the geographic priority hierarchy (North Sea Q1 → North Sea any season → NE Atlantic → global median) applied within TL-type entries only. If a species has no TL-type entry, it falls to the taxonomic fallback (genus mean → family mean → MISSING), not to FL or SL parameters.

**Audit action:** Confirm the `fishbase_lw_lookup_v2` table has a `type` or `length_type` column indicating TL. Confirm no FL or SL parameters are being used for any species. Print a count of species by parameter type used.

### 3d. Non-fish exclusion

**Check:** Are cephalopods, crustaceans, and bivalves excluded before computing E?

These taxa were explicitly removed by AphiaID in the original DATRAS pipeline session. The implementation should exclude them before the `sum(HLNoAtLngt)` step — they must not appear in N, S, or E.

Known excluded AphiaIDs include cephalopod, crustacean, and bivalve groups. The exclusion removed 7,105 rows from 183,441 → 176,336 fish-only rows.

**Audit action:** Confirm the exclusion filter is applied before N, S, E are computed. Print the number of rows dropped and confirm it matches ~7,105 (proportionally, given the 1985–2015 subset vs 1985–2026 full download).

### 3e. E normalisation — minimum individual = 1

**Check:** Is E normalised so the smallest individual in each haul has metabolic rate = 1?

The micbru `biomass.biomass()` function requires E expressed in units where the minimum individual metabolic rate = 1. This means E_norm = E_raw / m_min^0.75, where m_min is the mass (in grams) of the smallest individual in the haul across all species.

- m_min must be computed from the actual minimum length bin in the haul (the smallest observed individual across all species), not from a species-level minimum.
- The normalisation must be applied per haul, not globally.
- The `haul_state_variables.rds` output should contain both `E_raw` and `E` (normalised), and `m_min`. Confirm these three columns are present and that `E / E_raw` equals `1 / m_min^0.75` for each row.

**Audit action:** For a sample of 10 hauls, verify: `E_raw × (1 / m_min^0.75) == E`. Also confirm `min_epsilon == 1.0` for all hauls (this field in the current output is already uniformly 1, which is correct).

---

## 4. B_obs — observed biomass

**Check:** What are the units of B_obs, and is it consistent with B_pred?

The current `haul_eeos_predictions.rds` stores `B_obs` with values in the range 0.05–18,000 (median ~144). This is described as "already kg×1000, not kg" — i.e., it is stored in grams despite being derived from FishGlob `wgt` which is in kg. Confirm:

- The conversion `B_obs_g = wgt_kg × 1000` was applied in the pipeline before storing `B_obs`.
- `B_pred` (= `B_pred_norm × m_min`) is in grams (mass in grams of smallest individual × dimensionless normalised biomass count).
- Both are on the same scale before computing `ln_B_obs` and `ln_B_pred`.

**Audit action:** Print the FishGlob `wgt` column units from the data documentation or column metadata. Confirm the ×1000 conversion is explicitly present in the pipeline code. Print `B_obs` median alongside `B_pred` median and the ratio — the median B_pred/B_obs ratio is currently ~230, which represents a real scale mismatch likely due to catchability (survey hauls capture a fraction of community biomass), not a unit error. Document this explicitly in the pipeline comments so it is not mistaken for a bug later.

---

## 5. FishGlob join — haul_id format

**Check:** Is the DATRAS-FishGlob join producing a 1:1 match on haul_id?

The FishGlob `haul_id` format is: `Survey Year Quarter Country Ship Gear StNo HaulNo` (space-separated). The DATRAS HH fields that map to this are `Survey`, `Year`, `Quarter`, `Country`, `Ship`, `Gear`, `StNo`, `HaulNo`. This join was validated at 374/374 StatRec matches in the pipeline design session.

**Audit action:** After the join, print:
- Number of DATRAS hauls with E computed: should be ~12,389
- Number successfully joined to FishGlob: shown in `haul_state_variables.rds` as 12,117
- Number with EEoS predictions: 11,548
- Number dropped at each step and the reason (join failure vs EEoS failure vs filter)

The 569 hauls with state variables but no prediction (4.7% drop) need a documented reason — likely EEoS failure due to degenerate state variable combinations (e.g. N < S, or E ≤ N which makes λ₂ undefined).

---

## 6. EEoS call — input validation

**Check:** Are the S, N, E inputs to `biomass.biomass()` valid for every haul before the call?

The micbru implementation computes λ₂ = S / (E − N). This is undefined when E ≤ N. Every haul passed to the EEoS function must satisfy E > N. After normalisation, E is in units where the minimum individual = 1, and N is the count of individuals, so E > N should always hold if the normalisation is correct — but check it empirically.

**Audit action:**
- Print count of hauls where E ≤ N before EEoS call. This should be zero.
- Print count of hauls where S > N (impossible biologically — more species than individuals). Should be zero.
- Print count of hauls where S = 1 (single species). EEoS may behave degenerate here. Flag these.
- The 569 hauls dropped between state variables and predictions are likely failing one of the above. Confirm this.

---

## 7. n_species_with_lw > S anomaly

**Check:** Ten hauls have `n_species_with_lw > S`. This is logically impossible — you cannot have LW parameters for more species than are present.

This likely indicates a join or grouping error where the LW species count is being computed over a broader scope than S. 

**Audit action:** For the 10 affected hauls, print `haul_id`, `S`, and `n_species_with_lw`. Identify whether the extra species come from a join that is not properly restricted to the haul-level species list.

---

## 8. Year and quarter filter

**Check:** Is the analysis restricted to Q1, 1985–2015 throughout the pipeline?

The FishGlob filter is Q1, 1985–2015. The DATRAS HL download covers 1985–2026 and all quarters. The filter must be applied before E computation, not after — otherwise the m_min normalisation could be influenced by non-Q1 hauls.

**Audit action:** Confirm the DATRAS HL data is filtered to `Quarter == 1` and `Year %in% 1985:2015` before any E computation. Print year range and quarter values present in the filtered DATRAS data.

---

## 9. Summary of expected state variable ranges

Use these as sanity checks against the current outputs. Values outside these ranges warrant investigation:

| Variable | Expected range | Current data | Flag if |
|---|---|---|---|
| S | 1–50 | 2–36 ✓ | S > 60 or S = 0 |
| N | 10–10^6 | 8–1.52×10^6 ✓ | N < S |
| E_raw | any positive | 0.08–1.49×10^7 ✓ | E_raw ≤ 0 |
| E (normalised) | > N always | min=2 ✓ | E ≤ N |
| m_min (g) | 0.0001–100 | 0.000002–77.5 — lower tail suspicious | m_min < 0.00001 |
| B_obs (g) | 50–50,000 | 0.05–18,019 — lower tail suspicious | B_obs < 1 |
| B_pred (g) | varies | 0.03–9.27×10^8 | ratio B_pred/B_obs > 10,000 per haul |

The very low m_min values (minimum 0.000002 g) suggest some hauls have implausibly tiny minimum individuals — possibly a length bin at 0.5 cm being treated as a real observation rather than a data entry artefact. Flag any haul where m_min < 0.0001 g (roughly the mass of a 1cm fish) for manual inspection.

---

## 10. Residual sign convention

**Check:** Confirm the residual convention is consistent throughout.

- `residual = ln(B_obs) − ln(B_pred)` — negative means EEoS overestimates (B_pred > B_obs)
- `residual_norm = ln(B_obs) − ln(B_pred_norm)` — this is in mixed units (grams vs normalised metabolic units) and should NOT be used as the dependent variable for H1, H2, or H3
- The correct residual for all downstream analysis is `residual`, not `residual_norm`

**Audit action:** Confirm the column used for H2 and H3 will be `abs_residual` (= |ln(B_obs) − ln(B_pred)|), not `abs_residual_norm`. Add a comment in the code making this explicit.

---

## Priority order

Fix in this order:

1. **N source** (check 1) — most likely to invalidate current results if wrong
2. **S source** (check 2) — same concern
3. **n_species_with_lw > S** (check 7) — indicates a join bug
4. **m_min lower tail** (check 9) — implausibly small values will distort E normalisation
5. **B_obs units confirmation** (check 4) — document the ×1000 conversion explicitly
6. **LngtCode audit** (check 3a) — confirm Code 5 handling
7. **EEoS input validation** (check 6) — understand the 569 dropped hauls
8. **Year/quarter filter** (check 8) — confirm applied before E computation
9. **Residual convention** (check 10) — confirm correct column used downstream
