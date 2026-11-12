# Diagnosis: H1 hauls with FishGlob species weights missing from B_obs

**Status:** diagnosis only, nothing fixed. Script: `audit/scripts/diagnose_missing_weights.R` (re-runnable,
read-only with respect to the pipeline). Tables: `audit/outputs/missing_weights/`. Finding ID: **F-B9** in
`audit/findings.md`.

## What happens

`build_fishglob_haul_table()` (`pipeline/R/h1_join_helpers.R:33`) sets `B_obs = sum(wgt, na.rm = TRUE) * 1000`.
A taxon whose FishGlob `wgt` is NA therefore contributes 0 g, and the haul stays in the analysis with no flag.

1,407 of 145,055 FishGlob haul × taxon records in the 10,128 H1 hauls (0.97%) have NA `wgt`. They fall in
**1,149 hauls (11.3%)**.

### First, what FishGlob `wgt` is (this reverses my GATE A verdict on lead 8)

The last step of FishGlob's `get_datras.R` (l.1110, at the traced commit `acde3dc`) overwrites `wgt` with
`wgtlencpue × area_swept`. That is Σ count × a·L^b / 1000 (kg), using FishGlob's own length–weight table
(`length_weight/length.weight_DATRAS_3August2023.csv`) and the lower bound of each length class. Recomputing it
from our DATRAS HL with FishGlob's conventions matches FishGlob `wgt` within 1% for **99.4%** of 142,909
haul × taxon pairs. So:

- B_obs **is** "summed individual mass", as the methods say. It is not the reported `CatCatchWgt`.
- The ×100 and ×10 corrections in `get_datras.R` apply to the earlier catch-weight column, which this step
  overwrites. They don't reach the published `wgt`.
- B_obs comes from the **same DATRAS length data** as S, N and E (FishGlob's July 2021 extraction), just with a
  different length–weight table and length convention. The two products are not independent measurements of
  biomass. This matters for how the H1 result is read, so I've logged it as its own finding (F-B10).

## 1. Mechanism

Each NA record has exactly one cause:

| Mechanism | Haul × taxon records | Our HL (April 2026) for the same haul × taxon |
|---|---|---|
| **M1a** FishGlob's length–weight join (on taxa + family + genus + rank) fails because the family or genus in the data differs from the table, although the taxon name is there | 760 | complete lengths 732; taxon absent 28 |
| **M1b** taxon name not in FishGlob's length–weight table | 136 | complete lengths 136 |
| **M2** no length counts in FishGlob (`num` NA: taxon recorded with catch weight / total number only) | 511 | taxon absent 510; complete 1 |

None of the NAs come from DATRAS `CatCatchWgt = −9`: the published `wgt` doesn't use `CatCatchWgt`. M1a and M1b
are **FishGlob processing defects**. M2 is a real data gap (no lengths).

**By taxon** (`by_species.csv`). M1 hits every occurrence of the taxon:

| Taxon | Mechanism | Hauls with NA | % of the taxon's hauls |
|---|---|---|---|
| *Lumpenus lampretaeformis* | M1a (data: Lumpenidae; table: Stichaeidae) | 630 | 100% |
| *Phrynorhombus norvegicus* | M1b (table has it as *Zeugopterus norvegicus*) | 136 | 100% |
| *Gaidropsarus vulgaris* | M1a | 98 | 100% |
| *Raniceps raninus*, *Leptoclinus maculatus*, *Gaidropsarus mediterraneus* | M1a | 17, 7, 5 | 100% |
| *Agonus cataphractus*, *Echiichthys vipera*, *Buglossidium*, *Raja*, *Taurulus bubalis*, … | M2 | 80, 65, 51, 34, 33 | 4–80% |

**By DataType:** C 12.7% of hauls (1,079 / 8,493); R 4.1% (63 / 1,543); P 7.6% (7 / 92).
**By country:** SE 40.9% (539 / 1,319); FR 21.4% (332 / 1,555); DE 5.4%; DK 5.9%; GB 5.3%; GB-SCT 3.4%;
NL 2.0%. FR 1988–1994 and SE 1998–2007 have the highest counts. This follows where *Lumpenus* and
*Phrynorhombus* are caught and which countries recorded weight-only taxa.
**By year:** 3–24% of hauls each year, highest 1988–1994 (`by_year.csv`).

## 2. Magnitude per affected haul

I estimated each missing taxon's weight from our DATRAS lengths, using FishGlob's length–weight entry matched
by name and FishGlob's length convention (what FishGlob would have produced without the join defect): 733
records. Where there were no lengths, I used the reported species-category weight from the raw DATRAS CSV
(9 Nov 2026 download; C converted to the haul): 326 records. **348 records (in 225 hauls) stay unresolved**:
33 have lengths but no length–weight entry anywhere, and 315 have neither lengths nor a catch weight. For those
hauls the share below is a lower bound.

Shares are within the affected hauls (n = 1,149):

| Quantile | Share of individuals (%) | Share of biomass (%) |
|---|---|---|
| 25% | 0.000 | 0.003 |
| 50% | 0.04 | 0.02 |
| 75% | 0.16 | 0.09 |
| 90% | 0.53 | 0.72 |
| 95% | 1.06 | 2.17 |
| 99% | 5.0 | 11.3 |
| max | 77.8 | 24.9 |

| Missing biomass share | Hauls |
|---|---|
| ≤ 0.1% | 873 |
| 0.1–1% | 182 |
| 1–5% | 64 |
| 5–10% | 15 |
| 10–25% | 15 |

**Most affected hauls lose a trace.** 76% lose ≤ 0.1% of biomass, and 94 hauls lose more than 1%. The largest
losses (up to 25%) are FR 35TH hauls from 1988–1994. There, a weight-only taxon (no lengths, so it isn't in N
either) carries much of the catch.

## 3. How it shows up in the current results

**H1 residual.** Affected and unaffected hauls have the same residual distribution:

| | Hauls | Signed residual (mean / median) | Absolute residual (mean / median) |
|---|---|---|---|
| Unaffected | 8,979 | −0.652 / −0.673 | 0.669 / 0.675 |
| Affected | 1,149 | −0.656 / −0.666 | 0.664 / 0.666 |

Wilcoxon p = 0.43 (signed) and 0.26 (absolute). Within affected hauls, the missing share is weakly correlated
with the residual: Spearman +0.08 with the signed residual (p = 0.006) and −0.09 with the absolute residual
(p = 0.004). The direction is as expected: missing biomass lowers B_obs, which makes the over-prediction look
larger. But the effect is small.

**H2/H3 panel (8,670 hauls, 158 rectangles).** 600 hauls are affected. 135 rectangles have at least one, and the
worst rectangle has 45.5%.

- **Not random with respect to fishing pressure.**
  - By phase: 11.0% of hauls in 1985–1991, 7.9% in 1992–2001, 3.8% in 2002–2007, 3.9% in 2008–2015.
  - Affected hauls have higher FP_between (mean 9.42 vs 8.96) and higher FP_within (+0.22 vs −0.02).
  - At rectangle level, Spearman(FP_between, % affected) = 0.42 (p = 3e−8).
  - A haul-level logistic regression of "affected" on the model's own fixed effects has a clear FP_between slope
    (+0.58, SE 0.09) and FP_within slope (+0.74, SE 0.12), plus phase interactions. Hauls are not independent
    here, so read this as descriptive.
- **So it is a potential confound in principle.** The missingness is linked to FP through country and region,
  because SE and FR fish the eastern and southern rectangles.
- **But it is too small to matter in practice.** The treatments below show this.

## 4. Effect of each treatment

Primary CAR model only (no permutation or subsampling), fitted on both the signed residual (current code) and
the absolute residual (corrected specification). Under (b), the thresholds use the estimated share, which is a
lower bound for the 225 hauls with an unresolved taxon. Treatment (a) reproduces the live `wb_car_v2` slopes
exactly (max |Δ| 9e−18).

### H1

| Treatment | Hauls | 1:1 log R² | Bias, mean(ln B_obs − ln B_pred) | Bias-removed R² | Median B_pred/B_obs |
|---|---|---|---|---|---|
| (a) current, NA counted as 0 | 10,128 | 0.6610 | −0.652 | 0.9282 | 1.958 |
| (b1) exclude > 1% missing | 10,034 | 0.6599 | −0.654 | 0.9281 | 1.961 |
| (b5) exclude > 5% missing | 10,098 | 0.6604 | −0.653 | 0.9282 | 1.959 |
| (b10) exclude > 10% missing | 10,113 | 0.6608 | −0.653 | 0.9282 | 1.959 |
| (c) exclude all affected | 8,979 | 0.6693 | −0.652 | 0.9297 | 1.960 |
| (d) fill from DATRAS lengths × LW (catch weight where no lengths) | 10,128 | 0.6611 | −0.652 | 0.9280 | 1.957 |

Bias-removed R² = 1 − Σ(r − r̄)² / SS_tot, i.e. the 1:1 line shifted by the mean log bias.

### H2/H3 primary CAR: phase slopes, slope (SE)

`*` raw p < 0.05; `**` Bonferroni p < 0.05 (×4 per family). Full table with p-values:
`audit/outputs/missing_weights/treatments_car_slopes.csv`.

**Signed residual (current code)**

| Treatment | Hauls | FP_between 85–91 | 92–01 | 02–07 | 08–15 | FP_within 85–91 | 92–01 | 02–07 | 08–15 |
|---|---|---|---|---|---|---|---|---|---|
| (a) | 8,670 | +0.004 | +0.026 * | +0.013 | −0.025 * | −0.011 | +0.034 ** | −0.032 ** | +0.032 ** |
| (b1) | 8,589 | +0.004 | +0.025 | +0.012 | −0.026 * | −0.014 | +0.033 ** | −0.031 ** | +0.032 ** |
| (b5) | 8,640 | +0.004 | +0.025 | +0.012 | −0.026 * | −0.012 | +0.035 ** | −0.032 ** | +0.032 ** |
| (b10) | 8,655 | +0.004 | +0.025 * | +0.012 | −0.026 * | −0.011 | +0.035 ** | −0.032 ** | +0.032 ** |
| (c) | 8,070 | +0.001 | +0.025 | +0.011 | −0.025 | −0.020 | +0.039 ** | −0.030 * | +0.032 ** |
| (d) | 8,670 | +0.004 | +0.026 * | +0.013 | −0.026 * | −0.011 | +0.034 ** | −0.032 ** | +0.032 ** |

**Absolute residual (corrected specification)**

| Treatment | Hauls | FP_between 85–91 | 92–01 | 02–07 | 08–15 | FP_within 85–91 | 92–01 | 02–07 | 08–15 |
|---|---|---|---|---|---|---|---|---|---|
| (a) | 8,670 | −0.006 | −0.030 ** | −0.016 | +0.021 | +0.012 | −0.035 ** | +0.027 * | −0.031 ** |
| (b1) | 8,589 | −0.006 | −0.029 ** | −0.016 | +0.021 | +0.015 | −0.034 ** | +0.027 * | −0.031 ** |
| (b5) | 8,640 | −0.006 | −0.029 ** | −0.016 | +0.021 | +0.013 | −0.036 ** | +0.027 * | −0.031 ** |
| (b10) | 8,655 | −0.006 | −0.029 ** | −0.016 | +0.021 | +0.012 | −0.036 ** | +0.027 * | −0.031 ** |
| (c) | 8,070 | −0.001 | −0.029 ** | −0.015 | +0.020 | +0.023 * | −0.040 ** | +0.025 * | −0.031 ** |
| (d) | 8,670 | −0.006 | −0.030 ** | −0.016 | +0.021 | +0.012 | −0.035 ** | +0.027 * | −0.031 ** |

All 158 rectangles are retained under every treatment.

## Reading

- **The treatments barely matter.** No phase slope changes sign under any treatment. The largest change is
  |Δ| ≤ 0.011 under (c), against SEs of 0.009–0.014. The H1 metrics move in the third decimal place. The
  treatment that changes the most is (c), and that is because it drops 11% of hauls (an FR/SE-heavy subset), not
  because it corrects anything.
- **The residual definition matters far more.** 96.7% of hauls have a negative signed residual (over-prediction),
  so |r| ≈ −r and the absolute-residual slopes are close to the negated signed slopes. One consequence: under
  the corrected specification, the FP_between slope in 1992–2001 passes Bonferroni (−0.030, p_bonf < 0.05),
  which it didn't on the signed residual.
- **Treatment choice is yours.** On the evidence, (d) is closest to what FishGlob intended. Fill M1 from DATRAS
  lengths × FishGlob length–weight by name, and M2 from the reported catch weight. 225 hauls would still have
  an unresolved taxon. (a), (b) and (d) are interchangeable for the reported inference.
- **Separately, worth reporting upstream to FishGlob:** the M1a join defect (family/genus mismatches between the
  data and the length–weight table) and the M1b name mismatch (*Phrynorhombus* vs *Zeugopterus norvegicus*).
