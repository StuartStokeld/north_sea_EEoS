# Agent onboarding — North Sea EEoS repo

**Repo:** https://github.com/StuartStokeld/north_sea_EEoS  
**Purpose:** Masters project — EEoS haul-level biomass (H1), fishing-pressure effects on EEoS residuals spatially (H2) and temporally (H3), NS-IBTS Q1 1985–2015.

Read this after `git clone` to understand layout, what is complete, and what is missing locally.

---

## 1. Start here (human / supervisor)

| Document | Role |
|----------|------|
| [`README.md`](README.md) | Repo map |
| [`display_discussion/One page readme`](display_discussion/One%20page%20readme) | One-page narrative |
| [`display_discussion/H2_H3_methods_justification`](display_discussion/H2_H3_methods_justification) | Methods justification + figures |
| [`outputs/live_pipeline_run_log.md`](outputs/live_pipeline_run_log.md) | **Canonical list of live scripts** |

---

## 2. Repository layout

| Path | Status |
|------|--------|
| `pipeline/` | **Live** H1 + H2/H3 scripts (structbreak, within-between model) |
| `outputs/` | Live pipeline results (RDS, CSV, figures) |
| `display_discussion/` | Design docs, results drafts, presentation figures |
| `exploratory/` | Superseded work (zone schemes, blended-term model, biomass runs, step-0 diagnostics) |
| `H1_results/` | Shareable H1 package |
| `supplementary/` | Legacy notebooks and one-off tools |
| `local_env_upload_28-7/` | **Restore scripts + list of data not in git** |

---

## 3. Live pipeline (presented results)

**H1:** `build_datras_state_variables.R` → `build_eeos_predictions.R` → Harte baseline / ln(E) / null / diagnostics.

**H2/H3:** `import_couce_fishing_effort.R` → `build_h2_rectangle_panel.R` → `run_h2h3_structbreak_check.R` → **`run_h2h3_within_between.R`** (primary, biomass-free, FP_between/FP_within decomposition) → proportional effects / presentation figures.

Full command list: [`outputs/live_pipeline_run_log.md`](outputs/live_pipeline_run_log.md).

---

## 4. External data — not all in GitHub

Before re-running from scratch, fetch dependencies:

```bash
bash local_env_upload_28-7/download_scripts/download_all.sh
```

| Dependency | In git? | Notes |
|------------|---------|-------|
| `outputs/datras_hl_raw.rds` | **Yes** (~3.8 MB) | Processed HL; H1 can start here |
| ICES DATRAS raw CSV | **No** (251 MB) | Manual download — see `local_env_upload_28-7/download_scripts/04_download_datras_hl.md` |
| `FishGlob_data/.../NS-IBTS_clean.RData` | **No** | Clone script |
| `equation_of_state/biomass.py` | **No** | Clone script |
| Couce fishing CSV | **No** (gitignored) | Download script |
| `renv.lock` | **Yes** | Run `renv::restore()` |
| `.venv` | **No** | Recreate from `local_env_upload_28-7/python_requirements.txt` |

Complete table: [`local_env_upload_28-7/README.md`](local_env_upload_28-7/README.md).

---

## 5. Analysis status (Jul 2026)

| Hypothesis | Status |
|------------|--------|
| **H1** | Complete — EEoS fails unfitted 1:1 at haul level; dominance/catchability explored as diagnostics |
| **H2 (rectangle SEM)** | Complete — fishing pressure not robust after spatial correction; biomass confound documented |
| **H2/H3 (shared model)** | **Primary results:** within-between decomposed, biomass-free CAR/RE model; policy-anchored `phase_v2` (1992/2002/2008) |
| Exploratory zone schemes | Superseded — in `exploratory/` |
| Blended-term + GAM model | Superseded — in `exploratory/` |

Results drafts: `display_discussion/H2_results_draft.md`, `H3_results_draft.md`, `One page read me.md`.
(`H2_H3_results_interpretation.md` is superseded — blended-term era.)

---

## 6. Conventions agents should respect

- Run scripts from **repo root**: `Rscript pipeline/...`
- Live H2/H3 glmmTMB/spaMM scripts: **`Rscript --vanilla`** (packages not in renv.lock)
- Primary residual: `log(B_obs) - log(B_pred)`
- Exploratory material stays in `exploratory/` — do not treat as presented results
- Do not delete git history; reorganize visibility only

---

## 7. Known local-only state (28 Jul 2026 upload)

- Uncommitted local edits may exist on `.Rprofile`, H1 summary HTML/Rmd — check `git status`
- Large nested clones (`FishGlob_data/`, `equation_of_state/`) expected beside repo after restore scripts
- Manifests at upload: `local_env_upload_28-7/MANIFEST_*.txt`
