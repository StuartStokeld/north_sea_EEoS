# Exploratory / superseded analysis

Research history that informed the final design, but is **not** the live pipeline
a supervisor should land on first.

Nothing here has been deleted from git history — it was relocated for visibility.

For the live pipeline and presented results, see:

- [`../outputs/live_pipeline_run_log.md`](../outputs/live_pipeline_run_log.md)
- [`../pipeline/`](../pipeline/)
- [`../display_discussion/`](../display_discussion/) (design document, one-pager, results drafts)

---

## Layout (mirrors the main repo)

| Path | Contents |
|------|----------|
| `pipeline/` | Superseded / exploratory run scripts |
| `pipeline/R/` | Helpers exclusive to those scripts (+ symlinks to shared live helpers) |
| `outputs/` | Tables, run logs, and figures from those runs |
| `outputs/exploratory_review/` | Earlier collation of exploratory H1/H3 review materials |

---

## What was moved here

1. **Pre-H3 feasibility** — `run_h3_pre_exploration.R` and outputs (`h3_pre_*`).
2. **Zone schemes** — `run_h3_policy_zones.R` (Scheme A block-merge, Scheme B pressure-tier) and outputs (`h3_policy_*`); related design-support zone figures.
3. **Biomass-included shared-model feasibility** — `run_h2h3_shared_model_feasibility.R`, `…_round2.R`, and most `h2h3_feasibility_*` outputs.
4. **Blended-term (undecomposed) model + GAM** — `run_h2h3_shared_model_results.R`, `run_h2h3_proportional_effects.R`, `run_h2h3_temporal_robustness.R`, and their outputs/figures.
5. **Design-support explorations** — `run_h2h3_design_support.R` (aside from A1 / A4 inputs kept for the live structbreak script).
6. **Step-0 / dominance / biomass H2 runs** — `run_h2_dominance_*`, `run_h2_biomass_*`, `run_h2_robustness_*`, `run_h2_sar_lag_*`, `run_missing_visuals.R` and matching outputs under `exploratory/outputs/`.

A few intermediate files remain under top-level `outputs/` because the live
within-between / structbreak scripts still read them — listed in the live pipeline
run log.
