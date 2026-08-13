# H2 rectangle sub-sampling sensitivity

Each iteration drops 25% of rectangles (without replacement) from the 158 ICES rectangles, then refits `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)` on the remaining hauls. Reference signs / significance are from the full-sample primary RE fit (`wb_primary_v2`).

- Iterations requested: 1000
- Converged iterations: 1000 (failed: 0)
- Seed: 42
- Runtime: 247.0 sec (0.25 sec/iter)

## Summary by phase

| phase | full-sample slope | sign | sig (α=0.05) | % sign retained | % significant | % sig & same sign | n_ok |
|---|---|---|---|---|---|---|---|
| 1985-1991 | -0.0247 | − | N | 100.0% | 13.4% | 13.4% | 1000 |
| 1992-2001 | +0.0142 | + | N | 95.4% | 1.8% | 1.8% | 1000 |
| 2002-2007 | +0.1223 | + | Y | 100.0% | 100.0% | 100.0% | 1000 |
| 2008-2015 | +0.0519 | + | Y | 100.0% | 97.5% | 97.5% | 1000 |

## Deliverable sentences

- Phase 1985-1991 (negative full-sample RE slope): original sign held in 100.0% of 1000 rectangle subsamples; significance at α=0.05 held in 13.4% (same-sign significance in 13.4%).
- Phase 1992-2001 (positive full-sample RE slope): original sign held in 95.4% of 1000 rectangle subsamples; significance at α=0.05 held in 1.8% (same-sign significance in 1.8%).
- Phase 2002-2007 (positive full-sample RE slope): original sign held in 100.0% of 1000 rectangle subsamples; significance at α=0.05 held in 100.0% (same-sign significance in 100.0%).
- Phase 2008-2015 (positive full-sample RE slope): original sign held in 100.0% of 1000 rectangle subsamples; significance at α=0.05 held in 97.5% (same-sign significance in 97.5%).

Figure: `outputs/figures/h2_rectangle_subsampling_slope_histograms.png` (per-phase density of subsampled slopes; dashed red = full-sample RE).

### Notes

- Subsampling is at the rectangle level (random-effect / spatial unit), not hauls.
- Part B uses the primary RE model for computational feasibility (~4 min here); Part A multiplicity corrects the presented H2 slopes (`wb_car_v2`).
- Significance in the loop is uncorrected α = 0.05 (as specified in the brief).

