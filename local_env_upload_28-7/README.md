# Local environment upload — 28 Jul 2026

Snapshot of what lived on Stuart's Mac before losing local access. **Primary project content stays at the repo root** (`pipeline/`, `outputs/`, `display_discussion/`, `exploratory/`). This folder holds restore scripts, manifests, and items too large or machine-specific for the main tree.

---

## Files NOT in GitHub — must be fetched locally

| Item | Expected path (after restore) | Size (approx.) | How to obtain |
|------|------------------------------|----------------|---------------|
| **ICES DATRAS HL CSV** | `Unaggregated trawl and biological information_*/` … `.csv` | 251 MB | Manual — see [`download_scripts/04_download_datras_hl.md`](download_scripts/04_download_datras_hl.md). **Or skip** if `outputs/datras_hl_raw.rds` suffices. |
| **FishGlob NS-IBTS** | `FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData` | 14 MB | `bash download_scripts/01_clone_fishglob.sh` |
| **EEoS Python** | `equation_of_state/biomass.py` | 17 MB repo | `bash download_scripts/02_clone_equation_of_state.sh` |
| **Couce fishing effort** | `data/external/couce_trawling_effort/NorthSea_trawling_effort_1985to2015_REVIEW_v2.csv` | 2.5 MB | `bash download_scripts/03_download_couce_fishing_effort.sh` (gitignored; script downloads) |
| **Python venv** | `.venv/` | ~217 MB | **Not uploaded.** `bash download_scripts/05_setup_python_venv.sh` + [`python_requirements.txt`](python_requirements.txt) |
| **renv library** | `renv/library/` | ~6 MB | **Not uploaded.** `renv::restore()` from root `renv.lock` |
| **Session junk** | `.RData`, `.Rhistory`, `Rplots.pdf`, `.Rproj.user/` | varies | **Not uploaded** — discard |

### Optional / not required for live pipeline

| Item | Notes |
|------|--------|
| Full `FishGlob_data/` (all surveys) | Clone script gets NS-IBTS only; ~329 MB if you keep entire repo |
| `gis/NS_EC.qgz`, `gis/ns_ibts_hauls_2000_2020.gpkg` | QGIS extras — commit separately if needed |
| `north_sea_eeos_H1_supervisor_package.zip` | Build artifact — not in git |
| Word drafts | Copied to [`docs/write_up/`](docs/write_up/) in this archive |

---

## Quick restore (new machine)

From repo root:

```bash
bash local_env_upload_28-7/download_scripts/download_all.sh
bash local_env_upload_28-7/download_scripts/05_setup_python_venv.sh   # optional
# In R: open north_sea_eeos.Rproj → renv::restore()
# Manual: DATRAS CSV if you need to rebuild datras_hl_raw.rds from source
```

---

## What's in this folder

| Path | Purpose |
|------|---------|
| [`download_scripts/`](download_scripts/) | Clone/download scripts + DATRAS instructions |
| [`MANIFEST_git_status_at_upload.txt`](MANIFEST_git_status_at_upload.txt) | `git status` at upload time |
| [`MANIFEST_files_over_50MB.txt`](MANIFEST_files_over_50MB.txt) | Large local files (why they were excluded) |
| [`python_requirements.txt`](python_requirements.txt) | Pip freeze from `.venv` (28 Jul 2026) |
| [`diagnose_startup.sh`](diagnose_startup.sh) | R/renv startup diagnostic from local machine |
| [`docs/write_up/`](docs/write_up/) | Thesis Word drafts (not in primary tree) |

---

## R / Python versions (local machine, 28 Jul 2026)

- R 4.4.x with `renv` — see root `renv.lock`
- Python 3 venv at `.venv/` — see `python_requirements.txt`
- Ad hoc packages for H2/H3 (`glmmTMB`, `spaMM`, `strucchange`) installed outside renv — run live H2/H3 scripts with `Rscript --vanilla`
