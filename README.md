# North Sea EEoS (`north_sea_EEoS`)

PhD project: Hypothesis 1 — EEoS haul-level biomass (NS-IBTS Q1, 1985–2015).

Public GitHub repo: **`north_sea_EEoS`**. Open **`north_sea_eeos.Rproj`** at this level. Data, outputs, Python venv, and renv live here; code is organised in subfolders below.

---

## Folders

| Folder | Purpose |
|--------|---------|
| [`pipeline/`](pipeline/) | H1 analysis pipeline — run scripts, `R/`, tests, audit checklist |
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

## For supervisors (H1 results)

Open the pre-knitted report:

- [`H1_results/display_discussion/H1_results_summary.html`](H1_results/display_discussion/H1_results_summary.html)

Or knit from [`H1_results/`](H1_results/) — see [`H1_results/SUPERVISOR_README.txt`](H1_results/SUPERVISOR_README.txt).

Also useful: [`display_discussion/review_supervisor_h1_briefing.html`](display_discussion/review_supervisor_h1_briefing.html) (methods + figures).

---

## Reproduce the analysis

```bash
cd /path/to/north_sea_eeos
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
Rscript pipeline/run_pipeline_diagnostics.R
```

Required external data and setup: [`pipeline/README.md`](pipeline/README.md).

---

## Citation

Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). An equation of state unifies biodiversity, ecosystem functioning, and biomass in ecosystems. *Communications Biology*, 5, 957.
