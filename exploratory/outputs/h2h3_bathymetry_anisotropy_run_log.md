# H2/H3 bathymetry spatial anisotropy check — run log

Diagnostic test of whether primary-model residual spatial correlation is anisotropic and aligned with local shelf (depth-gradient) geometry. Design: `display_discussion/Design_bathymetry_spatial_anisotropy.md`.

## Session
sessionInfo written to: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/h2h3_bathymetry_anisotropy_sessionInfo.txt

## Locked conventions (documented before unblinding confirmatory p)
1. Bearings use gstat plane convention: degrees counterclockwise from +east (0° = east, 90° = north).
2. `bearing_depth` = along-shelf axis = local depth-gradient (cross-shelf) bearing + 90°. Depth = −GEBCO elevation (positive deeper).
3. `bearing_resid` / `bearing_blup` = direction bin (0/45/90/135 ±22.5°) with lowest local semivariance among pairs in lag window [0.5°, 4°] involving that rectangle (strongest residual similarity = along-shelf residual structure).
4. Confirmatory test: fold both bearings to axial [0°, 180), convert to radians, double the angle (axial→circular map), then `circular::cor.circular(..., test = TRUE)` (Jammalamadaka–Sarma).
5. TID: direct soundings = codes {10–17}; land TID=0 excluded from denominator. Flag if pct_sounding < 50% (sensitivity threshold 25%). Flagged rectangles excluded from the circular alignment test only.
6. Primary object: signed rectangle-mean residual; BLUP = robustness only.
7. Alignment confirmed if primary p < 0.05.

## Inputs
Panel rectangles: 158
Loaded residuals: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/residuals_by_rectangle.csv
Loaded BLUPs: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/blups_by_rectangle.csv

## Step 1 — Zonal bathymetry and TID extraction
Running: '/opt/homebrew/bin/python3' '/Users/stuartstokeld/Projects/north_sea_EEoS/pipeline/python/extract_rectangle_bathymetry.py' --shapefile '/Users/stuartstokeld/Projects/north_sea_EEoS/gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp' --panel-ids '/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/_tmp_panel_stat_rec.txt' --bathy '/Users/stuartstokeld/Projects/north_sea_EEoS/data/external/GEBCO_02_Aug_2026_e390ca8d46b0/gebco_2026_n65.0_s50.0_w-5.0_e10.0_geotiff.tif' --tid '/Users/stuartstokeld/Projects/north_sea_EEoS/data/external/GEBCO_02_Aug_2026_e390ca8d46b0/gebco_2026_tid_n65.0_s50.0_w-5.0_e10.0_geotiff.tif' --out '/Users/stuartstokeld/Projects/north_sea_EEoS/outputs/bathymetry_by_rectangle.csv'
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/bathymetry_by_rectangle.csv
TID flags: n_lt50=6 (excluded from alignment); n_lt25=0 (sensitivity)
Depth summary: mean_depth median=63.4 m; depth_range median=60.5 m; grad magnitude median=0.0014176
GEBCO edge distance: min=1.250 deg, median=4.250 deg

## Step 2 — Directional variograms (descriptive)
Residual directional fits:
  dir=0°  nugget=0.0001287  sill=1.779e-08  range=3.598  fit_ok=TRUE
  dir=45°  nugget=0.000183  sill=6.902e-09  range=3.538  fit_ok=TRUE
  dir=90°  nugget=0.0001929  sill=3.476e-08  range=3.86  fit_ok=TRUE
  dir=135°  nugget=0.0001831  sill=6.313e-09  range=4.017  fit_ok=TRUE
BLUP directional fits:
  dir=0°  nugget=0.009232  sill=6.64e-05  range=3.413  fit_ok=TRUE
  dir=45°  nugget=0.01253  sill=4.388e-05  range=3.446  fit_ok=TRUE
  dir=90°  nugget=0.01369  sill=0.0001616  range=3.448  fit_ok=TRUE
  dir=135°  nugget=0.01328  sill=2.223e-05  range=3.976  fit_ok=TRUE
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/directional_variogram_resid.png
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/directional_variogram_blup.png

## Step 3 — Per-rectangle bearings and Jammalamadaka–Sarma test
Alignment sample after TID exclusion: 152 / 158 rectangles
Finite bearing_resid: 152; bearing_blup: 152 (within TID-ok set)

### Confirmatory test (primary residuals)
n=152  rho=0.0483  statistic=0.5944  p.value=0.552272  (axial-doubled)

### Robustness test (BLUPs)
n=152  rho=0.0173  statistic=0.2157  p.value=0.829232  (axial-doubled)
Sensitivity TID<25% exclusion: n=158  rho=0.0371  p=0.64105

## Verdict
Primary p < 0.05? FALSE
Robustness p < 0.05? FALSE
Verdict: not_confirmed
Next step: Deprioritize further distance-metric refinement; prioritize the spatial-lag covariate test (Decision rule 2 from the Moran/BLUP run).
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/bearing_alignment_test.csv
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/bearings_by_rectangle.csv
Wrote: /Users/stuartstokeld/Projects/north_sea_EEoS/outputs/bathymetry_anisotropy_verdict.md
