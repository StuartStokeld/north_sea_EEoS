# FP_between spatial confounding bootstrap — summary (ARCHIVED RE)

> Archived: original spatial-confounding permutation test against the RE
> `(1 | stat_rec)` specification. Current CAR version:
> `pipeline/permutation_bootstrap_FP_between_CAR.R` /
> `outputs/permutation_bootstrap_FP_between_CAR_summary.md`.
> Relocated artifacts: `exploratory/outputs/permutation_bootstrap_FP_between_RE_*`.

Model-appropriate spatial confounding check for H2: randomly reassigns which
rectangle receives which `FP_between` value, then refits the within-between RE model.

- Formula: `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 |      stat_rec)` (REML)
- Phase scheme: **phase_v2** (`phase_v2`)
- Phase levels: 1985-1991 / 1992-2001 / 2002-2007 / 2008-2015

## Settings

- Seed: `42`
- `n_boot`: 1000
- `n_failed` (non-converged / NA): 0
- `n_ok`: 1000
- Cores: 1
- Runtime: 335.3 sec (0.34 sec/replicate)
- Rectangles: 158; hauls: 10464
- Source fit: `/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/primary_model_v2.rds`

### Matched permutation procedure

For each of the 1000 replicates: **one global shuffle** of rectangle-level
`FP_between` → **one refit** of
`residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)` →
**one `fixef()` extract** recording the reference-phase main effect and the three
`FP_between:phase_v2` interaction increments (phase slopes = main + interaction).
Not four independent shuffle-and-refit procedures. The four null distributions
are columns of the same 1000 × K matrix (matched randomizations).

## Observed H2 coefficients

- Reference-phase main effect `FP_between` (1985-1991): **-0.024729** (SE 0.015505; 95% CI [-0.055119, 0.005661]; p = 0.1107)

Phase-specific `FP_between` slopes (presented H2 effects):

- 1985-1991: **-0.024729** (SE 0.015505; 95% CI [-0.055119, 0.005661]; p = 0.1107)
- 1992-2001: **+0.014173** (SE 0.014137; 95% CI [-0.013536, 0.041881]; p = 0.3161)
- 2002-2007: **+0.122331** (SE 0.015465; 95% CI [0.092020, 0.152642]; p = 2.567e-15)
- 2008-2015: **+0.051887** (SE 0.014636; 95% CI [0.023201, 0.080574]; p = 0.0003924)

## Null distribution vs observed

| Target | Observed | Null mean | Null SD | Null 2.5% | Null 97.5% | Empirical p | Inside null 95%? |
|--------|----------|-----------|---------|-----------|------------|-------------|------------------|
| FP_between (reference-phase main effect, 1985-1991) | -0.024729 | +0.000165 | 0.015810 | -0.029975 | 0.031692 | 0.1199 | yes |
| Phase 1985-1991 FP_between slope | -0.024729 | +0.000165 | 0.015810 | -0.029975 | 0.031692 | 0.1199 | yes |
| Phase 1992-2001 FP_between slope | +0.014173 | +0.000284 | 0.015232 | -0.030904 | 0.030465 | 0.3487 | yes |
| Phase 2002-2007 FP_between slope | +0.122331 | +0.000840 | 0.019782 | -0.038040 | 0.040760 | 0.000999 (< 0.001) | no |
| Phase 2008-2015 FP_between slope | +0.051887 | +0.000146 | 0.016808 | -0.031329 | 0.033119 | 0.002997 | no |

## Interpretation (Step 5)

Rule:

- Observed coefficient **inside** the permutation null 95% interval → H2 effect does
  not depend on the spatial arrangement of `FP_between`; no evidence that the
  rectangle intercept and `FP_between` are substituting for one another.
- Observed coefficient **outside** the null 95% interval → rectangle intercept and
  `FP_between` compete for the same spatial signal under the real map; reported H2
  is at least partly attributable to that overlap.

Per-target statements:

### FP_between (reference-phase main effect, 1985-1991)

- Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
- Empirical p = 0.1199; inside null 95% = TRUE
- H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.

### Phase 1985-1991 FP_between slope

- Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
- Empirical p = 0.1199; inside null 95% = TRUE
- H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.

### Phase 1992-2001 FP_between slope

- Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
- Empirical p = 0.3487; inside null 95% = TRUE
- H2 effect does not depend on the spatial arrangement of FP_between. No evidence the rectangle intercept and FP_between are substituting for one another. Reported effect is not an artifact of spatial confounding.

### Phase 2002-2007 FP_between slope

- Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
- Empirical p = 0.000999 (< 0.001); inside null 95% = FALSE
- The rectangle intercept and FP_between are competing for the same spatial signal when FP_between's real spatial arrangement is present. The reported H2 effect is at least partly attributable to this overlap, not purely to the fishing-pressure effect itself.

### Phase 2008-2015 FP_between slope

- Direction: permutation null closer to zero than observed (spatial arrangement strengthens |coef|)
- Empirical p = 0.002997; inside null 95% = FALSE
- The rectangle intercept and FP_between are competing for the same spatial signal when FP_between's real spatial arrangement is present. The reported H2 effect is at least partly attributable to this overlap, not purely to the fishing-pressure effect itself.

## Headline

**Headline (phase-specific H2 slopes, phase_v2):** at least one presented phase-specific `FP_between` slope falls outside its spatial-permutation null 95% interval, indicating that the reported H2 effect depends on the real spatial arrangement of fishing pressure and is at least partly entangled with the rectangle intercept's spatial signal.

## Outputs

- Results CSV: `exploratory/outputs/permutation_bootstrap_FP_between_RE_results.csv`
- Null figure: `exploratory/outputs/figures/permutation_bootstrap_FP_between_RE_null.png`
- Run log: `exploratory/outputs/permutation_bootstrap_FP_between_RE_run_log.md`
- Session info: `exploratory/outputs/permutation_bootstrap_FP_between_RE_sessionInfo.txt`

