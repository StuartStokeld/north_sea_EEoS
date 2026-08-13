# North Sea EEoS (`north_sea_EEoS`)

Masters project: Ecological Equation of State analysis in North Sea demersal fish
communities (H1 haul-level biomass; H2/H3 fishing-pressure effects on EEoS
residuals, NS-IBTS Q1 1985–2015).

Public GitHub repo: **`north_sea_EEoS`**. Open **`north_sea_eeos.Rproj`** at this level.

---

## Start here (supervisor review)

| Document | What it is |
|----------|------------|
| [`display_discussion/One page readme`](display_discussion/One%20page%20readme) | One-page project narrative + links |
| [`display_discussion/H2_H3_methods_detailed.md`](display_discussion/H2_H3_methods_detailed.md) | H2/H3 design document (within-between CAR) |
| [`display_discussion/H3_results_draft.md`](display_discussion/H3_results_draft.md) / [`H2_results_draft.md`](display_discussion/H2_results_draft.md) | Hypothesis results notes |
| [`outputs/live_pipeline_run_log.md`](outputs/live_pipeline_run_log.md) | Explicit list of **live** scripts and outputs |

**Live pipeline (summary):** H1 haul-level EEoS → Couce/H2 panel → within-between
prerequisite → **policy-anchored `phase_v2` primary** (RE for H3, CAR for H2) →
Bonferroni / permutation / subsampling / Spec A. Full list:
[`pipeline/README.md`](pipeline/README.md) and the live pipeline run log above.

**Exploratory / superseded work** (SEM, Spec B, RE Spec A, blended-term, zones, anisotropy, etc.): [`exploratory/`](exploratory/).

---

## Folders

| Folder | Purpose |
|--------|---------|
| [`pipeline/`](pipeline/) | **Live** H1 / H2–H3 analysis scripts, helpers, tests |
| [`outputs/`](outputs/) | Live pipeline results (RDS, CSV, figures) + live run log |
| [`display_discussion/`](display_discussion/) | Design docs, one-pager, results drafts/summaries |
| [`exploratory/`](exploratory/) | Superseded / exploratory scripts + outputs (see its README) |
| [`H1_results/`](H1_results/) | Shareable H1 results package |
| [`supplementary/`](supplementary/) | Legacy notebooks and one-off tools |
| [`docs/`](docs/) | Thesis write-up drafts |
| [`gis/`](gis/) | QGIS project, shapefiles |
| `FishGlob_data/` | FishGlob NS-IBTS (local clone; not in this repo) |
| `equation_of_state/` | EEoS Python implementation (local clone; not in this repo) |
| [`local_env_upload_28-7/`](local_env_upload_28-7/) | **Local env archive (28 Jul 2026)** — download scripts, missing-data list, restore notes |
| [`AGENT_ONBOARDING.md`](AGENT_ONBOARDING.md) | **For agents:** clone → status → live pipeline → what to fetch |

**External data not in git?** Run `bash local_env_upload_28-7/download_scripts/download_all.sh` — see [`local_env_upload_28-7/README.md`](local_env_upload_28-7/README.md).

---

## Reproduce the live analysis

Run from the workspace root. Details and data prerequisites:
[`pipeline/README.md`](pipeline/README.md).

**H1:**

```bash
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
Rscript pipeline/run_pipeline_diagnostics.R
```

**H2/H3 primary + methods robustness:**

```bash
Rscript pipeline/import_couce_fishing_effort.R
Rscript pipeline/build_h2_rectangle_panel.R
Rscript --vanilla pipeline/run_h2h3_within_between.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R
Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R
Rscript --vanilla pipeline/run_h2h3_h2_multiplicity_subsampling.R
Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR.R
Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_refit.R
Rscript --vanilla pipeline/run_h2h3_rectangle_subsampling_car_refit.R
Rscript --vanilla pipeline/build_knn_spatial_weights.R
Rscript --vanilla pipeline/run_h2h3_spec_a_car_identifiability.R
Rscript --vanilla pipeline/permutation_bootstrap_FP_between_CAR_spec_a.R
```

Structural-break discovery and the H2 dose-response GAM check are under
[`exploratory/`](exploratory/) (historical; not part of the main pipeline).

---

## Citation

Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). An equation of state
unifies biodiversity, ecosystem functioning, and biomass in ecosystems.
*Communications Biology*, 5, 957.
