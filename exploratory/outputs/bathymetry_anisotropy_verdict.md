# Bathymetry spatial anisotropy — verdict

**Date:** 2026-08-04 10:49 BST
**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`
**Run log:** `outputs/h2h3_bathymetry_anisotropy_run_log.md`

## Confirmatory test (Jammalamadaka–Sarma)

| Test | n | ρ_JS | p | Significant (α=0.05) |
|---|---:|---:|---:|---|
| Primary (signed resid vs along-shelf depth) | 152 | 0.0483 | 0.5523 | FALSE |
| Robustness (BLUP vs along-shelf depth) | 152 | 0.0173 | 0.8292 | FALSE |
| Sensitivity (resid, TID flag at 25%) | 158 | 0.0371 | 0.6411 | FALSE |

## Verdict: **not_confirmed**

**Next step:** Deprioritize further distance-metric refinement; prioritize the spatial-lag covariate test (Decision rule 2 from the Moran/BLUP run).

## Notes

- TID exclusion at 50%: 6 rectangles flagged; retained in tables/maps but excluded from alignment test.
- Directional variograms (descriptive): `outputs/directional_variogram_resid.png`, `outputs/directional_variogram_blup.png`.
- Per-rectangle bathymetry: `outputs/bathymetry_by_rectangle.csv`; bearings: `outputs/bearings_by_rectangle.csv`.
