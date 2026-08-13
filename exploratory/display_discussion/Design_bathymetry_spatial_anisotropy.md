# Design Document: Testing an Ecologically Relevant Spatial Structure via Bathymetry

**Status:** Complete — run 2026-08-04; verdict **not_confirmed**  
**Depends on:** completed BLUP/residual diagnostic (`pipeline/run_h2h3_primary_spatial_autocorr_check.R`)  
**Last updated:** 2026-08-04  
**Outputs:** `outputs/bathymetry_anisotropy_verdict.md`, `outputs/h2h3_bathymetry_anisotropy_run_log.md`  
**Runner:** `pipeline/run_h2h3_bathymetry_anisotropy_check.R`

---

## 1. Aim

The distance-decay spatial model tested earlier assumed an **isotropic** field — that
correlation between two rectangles depends only on the distance between them, not on
direction. That assumption was never tested directly; it was inferred only as one
possible explanation for why the fitted decay range (~17.7°) was unidentifiable and
exceeded the maximum inter-centroid distance in the data (~13°).

This design tests that inference directly: does spatial correlation in our model
residuals run stronger along the shelf than across it? If so, this gives a concrete,
ecologically grounded reason the isotropic assumption failed, and points to a specific
fix — an environment-informed distance metric — rather than defaulting to a more complex
covariance function fit on the same (wrong) distance basis.

**This is a diagnostic test, not a model-fitting exercise.** The goal is to decide between
two competing explanations for the unidentified decay range:

1. Correlation genuinely has no finite spatial range at this resolution (confounded with
   the fixed-effect spatial trend, or a truly borderless field).
2. Correlation has a real, finite range, but geographic (Euclidean lon/lat) distance is
   the wrong metric to detect it — shelf geometry is a better one.

---

## 2. Data

### 2.1 Bathymetry (primary input)

**File:** `data/external/GEBCO_02_Aug_2026_e390ca8d46b0/gebco_2026_n65.0_s50.0_w-5.0_e10.0_geotiff.tif`  
**Source:** GEBCO_2026 Grid, 15 arc-second interval, elevation in metres, GeoTIFF, base
SRTM15+ v2.8 fused with Seabed 2030 Regional Center data.

Used to compute, per rectangle:
- Mean depth
- Depth range (max − min) — flags rectangles straddling the shelf break
- Direction and magnitude of the local depth gradient at the centroid — defines a
  rectangle-specific "cross-shelf axis," independent of a fixed compass bearing

**Resolution note:** 15 arc-seconds ≈ 450 m at these latitudes; each ICES rectangle
(~0.5° lat × 1° lon) will contain thousands of grid cells. Aggregate via **zonal
statistics** against rectangle polygons — not used cell-by-cell.

### 2.2 Type Identifier (TID) grid (QC layer)

**File:** `data/external/GEBCO_02_Aug_2026_e390ca8d46b0/gebco_2026_tid_n65.0_s50.0_w-5.0_e10.0_geotiff.tif`

Flags, per cell, whether depth is from direct soundings vs. satellite-derived/estimated
vs. interpolated values. Used to compute, per rectangle, the proportion of area sourced
from measured soundings. Purpose: rule out the possibility that an apparent directional
signal is an artifact of interpolation rather than real shelf geometry.

### 2.3 SST front data — deferred, not pulled this round

Considered as a second candidate driver of anisotropy (thermal fronts can create
directional structure independent of depth). Held back pending the bathymetry result —
if bathymetry alone explains the observed anisotropy, this may not be needed for this
question. Revisit if bathymetry result is inconclusive or only partially explanatory.

### 2.4 Rectangle geometry (confirmed)

**Shapefile:** `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp`  
**Helper:** `pipeline/R/h2_common.R` → `h2_ices_shapefile_path()`  
**Same geometry** used to rebuild queen contiguity in
`pipeline/run_h2h3_primary_spatial_autocorr_check.R` (validated against archived Moran I
to machine precision). Panel: **158** rectangles from `outputs/h2_rectangle_panel.rds`.

### 2.5 Primary-model residual object (confirmed)

| Role | File | Column |
|------|------|--------|
| **Primary** | `outputs/residuals_by_rectangle.csv` | `resid` (signed rectangle-mean response residual) |
| Robustness | `outputs/blups_by_rectangle.csv` | `blup` |

Source model: `outputs/h2h3_wb_model_objects.rds` (`fit_wb`).  
Residual definition matches the Moran primary row in
`outputs/h2h3_primary_spatial_autocorr_run_log.md`.

---

## 3. Locked pre-hoc decisions

| # | Decision | Lock |
|---|----------|------|
| 1 | **TID confidence threshold** | Flag rectangle if **&lt;50%** of non-land cells are direct-sounding TID classes. Sensitivity report at **25%**. Flagged rectangles retained in maps/tables but **excluded from the gradient-alignment test**. |
| 2 | **Primary residual object** | Signed `resid`. BLUP bearing used only for robustness. |
| 3 | **Along-shelf axis** | Derived **per rectangle** from local depth-gradient direction. A single global compass split is a foil only, not confirmatory. |
| 4 | **Alignment decision rule** | See §3.1. |

### 3.1 Alignment decision rule (confirmatory)

Alignment confirmed if **p &lt; 0.05** on the circular–circular correlation between the
residual-correlation bearing and the local depth-gradient bearing, computed both ways
for robustness:

**Primary**

```r
circular::cor.circular(bearing_resid, bearing_depth, test = TRUE)
```

Returns the Jammalamadaka–Sarma coefficient and its asymptotic p-value directly.

**Robustness (same test, substitute BLUP-correlation bearing)**

```r
circular::cor.circular(bearing_blup, bearing_depth, test = TRUE)
```

- Confirmatory call is the **primary** (residual) test.
- Robustness agreement strengthens the claim; disagreement is reported as partial
  and does not overturn a clear primary result by itself.
- Both circular variables are **per-rectangle** angles (length 158, or fewer after TID
  exclusion), so the test is well-defined. Global directional-variogram bearings are
  descriptive context, not the confirmatory statistic.

**Local residual-correlation bearing (`bearing_resid`):** for each rectangle, among
pairs involving that rectangle in a fixed lag window, take the direction of strongest
similarity (lowest local semivariance / highest correlation). Convert to the
along-contour (along-shelf) axis by rotating the cross-shelf depth-gradient bearing
by 90° when comparing “along-shelf residual structure vs along-shelf depth contours”
— document the convention in the run log before unblinding the p-value.

---

## 4. Method

1. **Zonal extraction.** For each of the 158 rectangles: mean depth, depth range, and
   depth-gradient direction/magnitude at the centroid. Cross-reference against the TID
   grid for percent-sounding-coverage; apply the §3 TID threshold.
2. **Directional variogram on model residuals — not on depth itself.** Depth is the
   candidate *explanation*; the object tested for anisotropy is the primary model's
   rectangle-mean residuals. Split pairwise semivariance by bearing
   (e.g. `gstat::variogram(resid ~ 1, locations = ~lon+lat, data = rect_centroids,
   alpha = c(0,45,90,135))`) and compare range/sill across bearings. Descriptive only.
3. **Per-rectangle bearings + confirmatory test.** Estimate `bearing_resid` (and
   `bearing_blup`) locally; compare to `bearing_depth` via §3.1. A match
   (correlation strongest along-shelf, weakest cross-shelf, aligned with depth
   contours, primary p &lt; 0.05) is the confirmatory result.

### 4.1 Boundary coverage (resolved)

Panel polygon extent: lat **[51.0, 61.5]**, lon **[−4.0, 9.0]**.  
GEBCO extent: lat **[50, 65]**, lon **[−5, 10]**.

All 158 rectangles lie inside the grid with margin. No rectangle approaches the 65°N
edge. Westernmost (WEST = −4) and easternmost (EAST = 9) rectangles retain ≥1° of
bathymetry outside the polygon for gradient estimation. No boundary exclusion required
a priori; still log min distance-to-GEBCO-edge per rectangle as a QC column.

---

## 5. How results feed the next steps

| Result | Interpretation | Next step |
|---|---|---|
| Anisotropy confirmed, aligns with depth geometry (primary p &lt; 0.05) | Isotropic distance was the wrong metric — explains why the original exponential decay was unidentifiable (fitting one range to two different physical decay scales) | Replace raw lon/lat distance with an environment-informed distance (e.g. cross-shelf weighted, or depth-contour distance) in any future spatial-covariance model. Do not re-attempt a plain isotropic model. |
| Anisotropy weak or doesn't align with depth (primary p ≥ 0.05) | Points back to the earlier explanations: confounding with the fixed-effect spatial trend, or a genuinely borderless correlation field | Deprioritize further distance-metric refinement; prioritize the spatial-lag covariate test (testing the source-sink/mechanistic hypothesis directly) — consistent with Decision rule 2 from the Moran/BLUP run |
| Partial (e.g. directional variogram suggests anisotropy but circular test fails, or primary/robustness disagree) | Inconclusive for depth as the metric; possible SST-front or mixed drivers | Do not treat as confirmed alignment; consider SST-front layer before committing to a depth-based distance metric |

Either outcome is decision-relevant — this is not a test that can return a null result
with no bearing on next steps.

---

## 6. Previously open risks — now closed

| Risk | Resolution |
|------|------------|
| Rectangle polygon file unconfirmed | Confirmed: `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp` (same as queen-weights rebuild) |
| Fixed global “along-shelf” bearing | Forbidden for confirmatory test; per-rectangle depth-gradient axis only |
| Boundary coverage at 65°N | Checked; panel max NORTH = 61.5°; safe |
| TID threshold undecided | Locked at 50% (sensitivity 25%) before seeing results |
| SST fronts | Deferred until bathymetry verdict |
| Package availability | Implementation will need `terra`/`sf` (zonal), `gstat` (directional variogram), `circular` (Jammalamadaka–Sarma). Likely outside `renv.lock` — run with `Rscript --vanilla` after install, same pattern as glmmTMB/spaMM scripts |

---

## 7. Deliverables

- Per-rectangle bathymetry table: mean depth, depth range, gradient direction/magnitude,
  percent-sounding-coverage (from TID), distance-to-GEBCO-edge, TID-flag
- Directional variogram plot + fitted range/sill by bearing, on primary model residuals
  (and BLUP robustness panel)
- Per-rectangle bearing table + Jammalamadaka–Sarma coefficient, asymptotic p-value
  (primary and robustness)
- Short written verdict: anisotropy confirmed / not confirmed / partial, with the
  corresponding next-step recommendation from §5

### Proposed output paths

```
outputs/bathymetry_by_rectangle.csv
outputs/directional_variogram_resid.png
outputs/directional_variogram_blup.png
outputs/bearing_alignment_test.csv
outputs/bathymetry_anisotropy_verdict.md
```

### Proposed script

```
pipeline/run_h2h3_bathymetry_anisotropy_check.R
```
(helpers under `pipeline/R/` as needed)
