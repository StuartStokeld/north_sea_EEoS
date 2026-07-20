# North Sea EEoS (`north_sea_EEoS`)

Masters project: Ecological Equation of State analysis in North Sea demersal fish communities (H1 haul-level biomass; H2 rectangle-level residuals vs fishing pressure, NS-IBTS Q1 1985–2015).

Public GitHub repo: **`north_sea_EEoS`**. Open **`north_sea_eeos.Rproj`** at this level. Data, outputs, Python venv, and renv live here; code is organised in subfolders below.

---

## Folders

| Folder | Purpose |
|--------|---------|
| [`pipeline/`](pipeline/) | H1 and H2 analysis pipelines — run scripts, `R/`, tests, audit checklist |
| [`outputs/`](outputs/) | Pipeline results (RDS, CSV, figures) |
| [`display_discussion/`](display_discussion/) | Reports, briefing HTML/Rmd, review notes |
| [`H1_results/`](H1_results/) | **Shareable H1 results package** — self-contained report + key outputs |
| [`supplementary/`](supplementary/) | Exploratory scripts, legacy notebooks, one-off tools |
| [`docs/`](docs/) | Thesis write-up drafts |
| [`gis/`](gis/) | QGIS project, shapefiles, exported haul layers |
| `FishGlob_data/` | FishGlob NS-IBTS (local clone; not in this repo) |
| `equation_of_state/` | EEoS Python implementation (local clone; not in this repo) |
| `Unaggregated trawl and biological information_*/` | ICES HL CSV source (local; not in this repo) |

---

## Results (to add to as I go)
(H1 results)

Open the pre-knitted report:

- [`H1_results/display_discussion/H1_results_summary.html`](H1_results/display_discussion/H1_results_summary.html)

Or knit from [`H1_results/`](H1_results/) — see [`H1_results/SUPERVISOR_README.txt`](H1_results/SUPERVISOR_README.txt).

Also useful: [`display_discussion/review_supervisor_h1_briefing.html`](display_discussion/review_supervisor_h1_briefing.html) (methods + figures).

---

## Reproduce the analysis

**H1 (haul-level EEoS):**

```bash
cd /path/to/north_sea_eeos
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
Rscript pipeline/run_pipeline_diagnostics.R
Rscript pipeline/explore_h1_haul_dominance.R          # dominance / size-homogeneity (D, size_CV)
Rscript pipeline/explore_h1_dominance_partial_r2.R    # partial R² follow-up
```

**H2 (rectangle-level residuals vs fishing pressure):** run after H1 (`build_eeos_predictions.R`).

```bash
Rscript pipeline/import_couce_fishing_effort.R
Rscript pipeline/build_h2_rectangle_panel.R
Rscript pipeline/run_h2_models.R
```

Required external data and setup: [`pipeline/README.md`](pipeline/README.md).

**H2 results summary:** [`display_discussion/H2_results_summary.md`](display_discussion/H2_results_summary.md)

---

## Citation

Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). An equation of state unifies biodiversity, ecosystem functioning, and biomass in ecosystems. *Communications Biology*, 5, 957.
