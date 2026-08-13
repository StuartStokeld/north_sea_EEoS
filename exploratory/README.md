# Exploratory / superseded analysis

Research history that informed the final design, but is **not** the live pipeline.

Nothing here has been deleted from git history — it was relocated for visibility.

For the live pipeline and presented results, see:

- [`../outputs/live_pipeline_run_log.md`](../outputs/live_pipeline_run_log.md)
- [`../pipeline/`](../pipeline/)
- [`../display_discussion/`](../display_discussion/)

---

## Layout

| Path | Contents |
|------|----------|
| `pipeline/` | Archived run scripts (ARCHIVED header on each) |
| `pipeline/R/` | Helpers exclusive to those scripts (+ symlinks to live helpers) |
| `pipeline/python/` | Archived bathymetry extraction helper |
| `outputs/` | Tables, run logs, and figures from archived runs |
| `display_discussion/` | Superseded design / results notes |

Re-run only intentionally from `exploratory/pipeline/` (update write paths to
`exploratory/outputs/` before re-running — many scripts still contain historical
top-level `outputs/` paths).

---

## What lives here (including 2026-08-13 archive)

1. **Pre-H3 / zone schemes / blended-term / GAM** — shared-model feasibility, policy zones, temporal robustness on blended term.
2. **Bai–Perron structural-break check** — original data-driven phases (1989 / 2001 / 2008).
3. **H2 dose-response GAM linearity** — development gate only.
4. **Early rectangle SEM** — `run_h2_models.R` + `h2_ols_*` / `h2_sem_*` / Moran map figures.
5. **RE Spec A + RE Spec A permutation** — historical first Spec A path; current claim is CAR Spec A in live `pipeline/`.
6. **Spec B / A+B** — neighbour biomass lag.
7. **RE FP_between permutation** — superseded by live CAR permutation.
8. **OLS residual-proxy rectangle subsampling** — superseded by GLMM/CAR Part B.
9. **Bathymetry anisotropy / BYM smoothness / n_hauls metric / haul-count maps** — supporting diagnostics not in methods model list.
10. **Original-phase standalone proportional / pooled / presentation / raw-correlation scripts** — superseded by `run_h2h3_phase_v2_reporting.R`.
11. **Superseded results note** — `display_discussion/H2_H3_results_interpretation.md` (blended-term era).

A few intermediate files remain under top-level `outputs/` because live scripts
still read them (e.g. `h2h3_wb_model_objects.rds`, CAR `adjMatrix` RDS) — listed
in the live pipeline run log.
