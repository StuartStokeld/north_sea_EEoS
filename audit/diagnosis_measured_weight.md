# Diagnosis: DATRAS measured catch weight as an alternative B_obs

**Status:** diagnosis for a decision; nothing changed in the pipeline. Script:
`audit/scripts/diagnose_measured_weight.R`. Tables: `audit/outputs/measured_weight/`. Source: the raw ICES CSV
(HH + HL blocks, downloaded 9 Nov 2026; local only, not in git). That is a different DATRAS snapshot from the HL
behind S/N/E (9 Apr 2026) and from FishGlob's (July 2021). Finding: **F-B16**.

## How to sum it

- **Unit of record.** `SpeciesCategoryWeight` (= `CatCatchWgt`, g) is recorded once per haul × species × sex ×
  category and repeated on every length row. Of 158,615 such keys in the H1 hauls, **none** has more than one
  distinct weight. So the correct sum is: de-duplicate on haul × AphiaID × SpeciesSex × SpeciesCategory, then add
  up the categories and sexes.
- **DataType P subcategories.** 80 species × sex combinations in P hauls have more than one category. These are
  separate subsamples and are summed.
- **Missing values.** −9 → missing. Zero weights (3,591 records) are kept as valid 0 g.
- **Conversion.** DataType C: × HaulDur/60 (per-hour values back to the haul; the DATRAS FAQ). R/P: as recorded.
  Units are g.
- **Species pool.**
  - "Measured haul biomass" = the sum over the taxa that carry length–weight coefficients in the current lookup,
    which is the S/N/E pool.
  - A haul is **complete** when every pool taxon in it has a valid weight.

## Coverage

Of 145,055 haul × species records, 76.4% have a valid weight; for pool taxa it is 76.7%.
**6,603 of 10,128 H1 hauls (65.2%) are complete.**

| DataType | Hauls | % complete | % pool records valid |
|---|---|---|---|
| C | 8,493 | 66.6 | 79.2 |
| R | 1,543 | 55.3 | 60.1 |
| P | 92 | 100 | 100 |

| Country | DE | DK | FR | GB | GB-SCT | NL | NO | SE |
|---|---|---|---|---|---|---|---|---|
| % complete | 78.1 | 86.5 | 72.0 | 100 | **10.3** | **15.7** | 94.5 | 92.6 |

By year, completeness is 34–51% in 1985–1989 and 1993–1995, 63–73% from 1990 to 2008, and 78–96% from 2009
(`coverage_by_year.csv`; country × year in `coverage_by_country_year.csv`). The missingness is strongly
structured: Scottish and Dutch hauls before the late 2000s mostly have no weight.

**Excluded taxa** (non-fish and taxa without length–weight coefficients) are a median 0.03% of measured weight in
complete hauls. The 90th percentile is 3.6% and the maximum 78%.

## Agreement with FishGlob B_obs (complete hauls, n = 6,603)

The ratio measured / FishGlob is **bimodal**:

| Quantile | 1% | 5% | 25% | 50% | 75% | 95% | 99% |
|---|---|---|---|---|---|---|---|
| As recorded | 0.006 | 0.008 | 0.009 | 0.776 | 0.892 | 0.998 | 1.56 |

- **The low cluster (about 0.009, i.e. ≈ 1/100) is a unit problem.** It comes in whole country-year blocks:
  - Of 163 country-years, 63 are entirely low, 96 entirely normal, and 4 mixed. Two of the mixed are 98% low;
    two NO years are 3–6% low.
  - The low blocks are DE 1990–2003, DK 1985–1999, GB 1985–2003, NO 1989–2002 and SE 1985–2001
    (`factor100_by_country_year.csv`).
  - This is the factor-100 problem FishGlob's `get_datras.R` corrects haul by haul (`resc`, l.926–1010) for the
    catch-weight column it later discards.
  - I haven't been able to confirm the cause against ICES documentation. It looks like those submissions were in
    a different unit (e.g. 10 g or hg), or ICES mis-scaled them.
- **After a block-level ×100 correction**, I get the table below. 163 hauls remain outside ×2.
  - The correction is data-informed: blocks are flagged where the median ratio to B_obs is below 0.05.

  | Quantile | 1% | 5% | 25% | 50% | 75% | 95% | 99% |
  |---|---|---|---|---|---|---|---|
  | ×100 blocks | 0.091 | 0.709 | 0.815 | **0.874** | 0.929 | 1.03 | 1.66 |

  - 35% of hauls agree within 10%, and 90% within 25%.
  - The median 0.87 matches the earlier DATRAS-implied ratio (0.89): FishGlob's length-based weight runs about
    13% above the measured weight.

## H1 and primary CAR using measured B_obs (complete hauls only)

The same 6,603 hauls for every row. H2/H3 keeps 5,239 hauls in 151 rectangles, with the adjacency restricted to
those rectangles.

| B_obs | Hauls | 1:1 log R² | Bias (mean ln obs − ln pred) | Bias-removed R² | Median B_pred/B_obs | cor² (ln) |
|---|---|---|---|---|---|---|
| FishGlob, all H1 hauls | 10,128 | 0.661 | −0.652 | 0.928 | 1.96 | 0.971 |
| FishGlob, complete hauls | 6,603 | 0.621 | −0.666 | 0.927 | 1.98 | 0.972 |
| Measured, as recorded | 6,603 | −1.36 | −2.81 | 0.028 | 3.28 | 0.106 |
| Measured, ×100 blocks | 6,603 | 0.419 | −0.831 | 0.872 | 2.25 | 0.916 |

Primary CAR on the absolute residual: phase slopes. `*` raw p < 0.05; `**` Bonferroni p < 0.05.

| B_obs (complete hauls) | FP_between 85–91 | 92–01 | 02–07 | 08–15 | FP_within 85–91 | 92–01 | 02–07 | 08–15 |
|---|---|---|---|---|---|---|---|---|
| FishGlob | −0.005 | −0.049 ** | −0.043 ** | −0.003 | +0.036 * | −0.049 ** | +0.047 ** | −0.042 ** |
| Measured, as recorded | +0.303 ** | −0.388 ** | −0.361 ** | +0.068 | +0.245 ** | −0.118 | +0.950 ** | −0.013 |
| Measured, ×100 blocks | +0.023 | −0.022 | −0.024 | −0.013 | +0.058 * | −0.010 | +0.075 ** | −0.060 ** |

The signed-residual slopes are in `tables.md`.

## Reading

1. **As recorded, measured weight cannot be used.** The unit blocks dominate everything. In the H2/H3 slopes,
   the block structure aliases with phase and FP, producing large spurious "effects".
2. **Even after the ×100 correction:**
   - Coverage is 65% and structured by country and era. GB-SCT and NL are mostly missing, which would remove
     much of the western and southern North Sea from H2.
   - Agreement with FishGlob is only moderate (median 0.87).
   - Measured B_obs fits EEoS less well (cor² 0.92 vs 0.97). That is expected given F-B10: FishGlob B_obs, like
     S, N and E, is built from the length data, while measured weight is independent of it.
3. **The H2/H3 picture shifts with the B_obs definition.**
   - With corrected measured B_obs, no FP_between slope is significant, but FP_within keeps the same 2002–07 and
     2008–15 signs (Bonferroni).
   - Read this as: H2 is sensitive to the B_obs definition, and H3's late-period pattern is more robust.
4. **Decision inputs:**
   - FishGlob B_obs has complete coverage, but it shares its length data with S, N and E (a circularity risk for
     H1).
   - Measured weight is independent, but has 65% coverage and needs an undocumented ×100 unit correction that
     ought to be confirmed with ICES before use.
   - A middle option: keep FishGlob as primary and report the corrected measured weight on complete hauls as a
     sensitivity analysis.
