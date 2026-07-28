# Missing visuals — haul map & composition figures — run log

Fills two gaps flagged in `outputs/exploratory_review/index.md` (Task 3, "core visual suite"): a standalone geographic haul map, and the plain S/D/size_CV composition distributions (as opposed to the residual-ratio overlays already produced for Task 2). H2 residual map deliberately out of scope.

Loaded `h3_pre_rectangle_usability_flags.csv` (197 rectangles), `h1_dominance_haul_table.csv` (12069 hauls), `haul_eeos_predictions.rds` (12069 hauls).
Composition table: joined dominance table to `haul_eeos_predictions.rds` on `haul_id` for `S` only; 12069 hauls, zero-loss join (verified 1:1 match, no recomputation).

## Haul map

Shapefile source: `gis/ICES_rectangles/ICES_Statistical_Rectangles_Eco.shp`, loaded via `load_ices_rectangles_sf()` (`pipeline/R/h2_spatial_helpers.R`), the same ICES StatRec grid already used for the H2/H3 spatial maps (e.g. `h3_pre_C2`/`D2`, `h3_policy_*` zone maps). No new shapefile was needed or added.
haul_map_1: 0 of 197 rectangles in the usability-flags table could not be matched to shapefile geometry and are excluded from the map (total_hauls / coverage counts below use the full 197).
Rectangles without Couce fishing-pressure coverage: 32 of 197 (matches the confirmed universe from Task 1's run log).

## Composition figures

All values (S, D, size_CV) read as-is from `haul_eeos_predictions.rds` / `h1_dominance_haul_table.csv`; none recomputed.

### "Big shoal" cutoff check
Searched `pipeline/*.R`, `pipeline/R/*.R`, and `display_discussion/*.md` for any existing formal cutoff jointly combining D and size_CV into a single "big shoal" classification. None found: the H1 dominance analysis (`explore_h1_haul_dominance.R`, `display_discussion/H1_dominance_results_summary.md`) treats the top decile of D and the bottom decile of size_CV as separate marginal extremes for descriptive/taxonomic breakdown only (e.g. "D, top decile" and "size_CV, bottom decile" reported as distinct rows), never combined into one joint threshold. `display_discussion/Methods_note_for_discussion.md` uses only qualitative language ("high D and low size_CV") with no numeric cutoff.
Per the briefing's fallback, this script marks the top decile of D as a visual reference line only — explicitly not a formal "big shoal" classification. Reused the existing top-decile boundary already computed in `outputs/h1_dominance_by_D_bins.csv` (produced by `explore_h1_haul_dominance.R`) rather than computing a fresh quantile.
D top-decile reference threshold = 0.9006 (from h1_dominance_by_D_bins.csv, n_bins=10 bin=10, n=1206 there vs n=1206 here after the S-join; 10.0% of hauls).

## Outputs
- `/Users/stuartstokeld/north_sea_eeos/outputs/figures/haul_map_1.png` — Total haul count per ICES rectangle, 1985-2015 (choropleth, N=197 rectangles); the 32 rectangles without Couce fishing-pressure coverage are outlined in red.
- `/Users/stuartstokeld/north_sea_eeos/outputs/figures/composition_1.png` — Distribution histograms of S, D, and size_CV per haul (N=12069).
- `/Users/stuartstokeld/north_sea_eeos/outputs/figures/composition_2.png` — Plain scatter of D vs size_CV, N=12069 hauls, density contours only (no ratio/residual colouring).
- `/Users/stuartstokeld/north_sea_eeos/outputs/figures/composition_3.png` — 1206 of 12069 hauls (10.0%) fall at/above the D top-decile reference line (D >= 0.901); reference only, not a formal "big shoal" classification (none exists elsewhere in the H1 dominance analysis).
- /Users/stuartstokeld/north_sea_eeos/outputs/missing_visuals_run_log.md (this file)
