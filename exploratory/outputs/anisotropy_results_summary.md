# Bathymetry / anisotropy results — pre-registered decision rules

**Date:** 2026-08-04 11:02 BST
**Source run:** `outputs/h2h3_bathymetry_anisotropy_run_log.md`; this report: `outputs/anisotropy_preregistered_report_run_log.md`
**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`

Decision rule for “Anisotropy confirmed?”: **asymptotic p < 0.05** on the 
Jammalamadaka–Sarma circular–circular correlation. Coefficient magnitude is 
descriptive only. Permutation p (999 pairing shuffles) is reported for the 
primary row as pre-registered robustness of the p-value, not a substitute rule.

## Table 1 — Headline result

| Test object | Circular correlation coefficient | p (asymptotic) | p (permutation) | Anisotropy confirmed? |
|---|---:|---:|---:|---|
| Signed residuals × local depth-gradient bearing (primary) | 0.0483 | 0.5523 | 0.5616 | No |

Primary pairing: residual-correlation bearing vs local **along-shelf** axis (= per-rectangle depth-gradient direction + 90°). N = 152 after excluding TID-flagged rectangles (<50% direct soundings).

## Table 2 — Robustness / contrast rows

*Not confirmatory. Same p < 0.05 labeling for transparency only.*

| Variant | Coefficient | p (asymptotic) | Confirmed under same rule? | Role |
|---|---:|---:|---|---|
| BLUPs × local depth-gradient bearing | 0.0173 | 0.8292 | No | Robustness check |
| Signed residuals × global compass bearing (foil) | 0.1239 | 0.1256 | No | Foil — not confirmatory |
| TID 25% threshold (sensitivity) | 0.0371 | 0.6411 | No | Sensitivity on exclusion rule |

Foil definition: depth-independent piecewise compass along-shelf — 90° (N–S) for lon < 4°E, 0° (E–W) for lon ≥ 4°E. Permutation p for foil = 0.0751; for BLUP = 0.8398; for TID-25% = 0.6316.

## Table 3 — TID exclusion summary

| Threshold | N flagged | Role |
|---|---:|---|
| <50% direct-sounding (primary exclusion) | 6 | Excluded from Table 1 / confirmatory tests; **kept in bearing plot** |
| <25% direct-sounding (sensitivity) | 0 | Labels the Table 2 sensitivity row only |

Total rectangles = 158. Table 1 N = 152. Deliverable plot `outputs/anisotropy_bearing_comparison_plot.png` marks the 6 TID-flagged rectangles (triangles) and does not remove them.

## Narrative verdict

No evidence that shelf-geometry-based direction explains the earlier unidentified decay range. Consistent with either a genuinely borderless spatial field or confounding with the fixed-effect spatial trend — does not distinguish between these two.

A confirmed alignment is expected to reflect a consistent ~90° rotational offset between residual-correlation bearing and the raw depth-gradient (cross-shelf) direction: Table 1 pairs residual bearing with the along-shelf axis (gradient + 90°), so confirmation means those axes coincide — not that residual correlation points the same compass way as the steepest depth slope.

**Verdict label:** `not_confirmed`

## Deliverables

- `outputs/anisotropy_results_summary.md` (this file)
- `outputs/anisotropy_bearing_comparison_plot.png`
- `outputs/anisotropy_preregistered_report_run_log.md`
