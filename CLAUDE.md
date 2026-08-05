# CLAUDE.md — North Sea EEoS project context

**GitHub:** https://github.com/StuartStokeld/north_sea_EEoS  
**Local path:** `/Users/stuartstokeld/Projects/north_sea_EEoS`  
**Default branch:** `main`  
**RStudio project:** `north_sea_eeos.Rproj`

Use this file as the primary briefing when working in this repo. Prefer it over scattered notes; for deeper detail follow the linked docs below.

---

## What this project is

Masters / PhD research testing the **Ecological Equation of State (EEoS / METE)** of Harte et al. (2022) on North Sea demersal fish from the **NS-IBTS Q1** survey (**1985–2015**).

Three hypotheses:

| ID | Question |
|----|----------|
| **H1** | Can EEoS predict haul-level biomass from state variables S, N, E (unfitted)? |
| **H2** | Do EEoS prediction failures correlate with fishing disturbance **spatially**? |
| **H3** | Do EEoS prediction failures correlate with fishing disturbance **temporally**? |

Core citation: Harte, J., Brush, J. M., Newman, E. A., & Umemura, K. (2022). *Communications Biology*, 5, 957.  
EEoS Python implementation: [micbru/equation_of_state](https://github.com/micbru/equation_of_state) (`biomass.py`), called from R via `reticulate`.

---

## Headline results (presented)

**H1** (n = 12,069 hauls): EEoS does **not** predict absolute haul biomass.  
`log_r² = −0.223` (worse than mean) vs productivity 1:1 baseline `0.736`; `cor² = 0.926` with systematic **~3.7–4.1× overprediction** that grows with biomass magnitude. Relative productivity-ratio structure resembles Harte’s “disturbed systems.” Dominance / “big shoals” is a second, smaller failure source.

**H2/H3** (shared model; n = 10,464 hauls, 158 rectangles): within-between fishing-pressure decomposition, biomass-free CAR/RE model; policy-anchored phases **1992 / 2002 / 2008** (`phase_v2`: 1985–1991 / 1992–2001 / 2002–2007 / 2008–2015).  
Neither H2 nor H3 shows a single uniform disturbance signal. Marginal R² ≈ **0.051**. Fishing pressure × phase effects flip direction across periods.

Narrative one-pager: `display_discussion/One page read me.md`  
Results drafts: `display_discussion/H2_results_draft.md`, `H3_results_draft.md`  
(`H2_H3_results_interpretation.md` is superseded — blended-term era only)

---

## Repo layout

| Path | Role |
|------|------|
| `pipeline/` | **Live** H1 + H2/H3 scripts and helpers |
| `outputs/` | Live results (RDS, CSV, figures) + `live_pipeline_run_log.md` |
| `display_discussion/` | Design docs, results drafts, presentation figures |
| `exploratory/` | **Superseded** work — do not treat as presented results |
| `H1_results/` | Shareable H1 package for supervisors |
| `supplementary/` | Legacy notebooks / one-off tools |
| `gis/` | ICES rectangle shapefiles, QGIS |
| `data/external/` | Couce fishing effort (CSV often gitignored) |
| `FishGlob_data/` | Local clone — **not fully in git** |
| `equation_of_state/` | Local clone of EEoS Python — **not fully in git** |
| `local_env_upload_28-7/` | Restore scripts + list of data missing from git |
| `AGENT_ONBOARDING.md` | Shorter agent clone/restore checklist |
| `renv.lock` | R package lockfile |

Canonical live-script list: `outputs/live_pipeline_run_log.md`.

---

## Live pipeline (run from repo root)

### H1

```bash
Rscript pipeline/build_datras_state_variables.R
Rscript pipeline/build_eeos_predictions.R
Rscript pipeline/run_h1_harte_baseline.R
Rscript pipeline/run_h1_lne_reference.R
Rscript pipeline/run_h1_null_model.R
Rscript pipeline/run_pipeline_diagnostics.R
```

### H2 panel → phases → primary H2/H3 model

```bash
Rscript pipeline/import_couce_fishing_effort.R
Rscript pipeline/build_h2_rectangle_panel.R
Rscript --vanilla pipeline/run_h2h3_structbreak_check.R
Rscript --vanilla pipeline/run_h2h3_within_between.R          # original data-driven phase
Rscript --vanilla pipeline/run_h2h3_phase_v2_refit.R          # PRIMARY (policy-anchored phase_v2)
Rscript --vanilla pipeline/run_h2h3_phase_v2_reporting.R      # slopes / CAR / proportional / figures
Rscript --vanilla pipeline/run_h2h3_wb_proportional_effects.R # original-phase proportional (archive)
```

Primary H2/H3 model: `residual ~ FP_between * phase_v2 + FP_within * phase_v2 + (1 | stat_rec)` 
(with optional CAR sensitivity). **No biomass covariate.** Policy-anchored `phase_v2`
(1992/2002/2008); original data-driven `phase` retained in artifacts for comparison.
Artifact: `outputs/primary_model_v2.rds`.

Use `Rscript --vanilla` for glmmTMB / spaMM / strucchange scripts (packages may be outside `renv.lock`).
Full step table and helpers: `pipeline/README.md`.

---

## Key scientific conventions

| Concept | Convention |
|---------|------------|
| Primary residual | `log(B_obs) - log(B_pred)` |
| Primary R² | `log_r2()` = 1 − SS_res/SS_tot (**can be negative**); not `cor(log B)²` |
| `B_obs`, `B_pred` | grams |
| `E` | normalised metabolic rate `E_raw / m_min^0.75` |
| Headline productivity baseline | `E_calibrated = E × m_min` (not raw `E_raw`) |
| Survey | NS-IBTS Q1, years 1985–2015 |
| Fishing pressure | Couce et al. rectangle-year trawling hours |
| H2 term | `FP_between` — rectangle long-run mean fishing pressure |
| H3 term | `FP_within` — year deviation from that rectangle’s mean |
| Phases (current primary / `phase_v2`) | Policy-anchored: 1992, 2002, 2008 (1985–1991 / 1992–2001 / 2002–2007 / 2008–2015) |
| Phases (original / `phase`) | Data-driven structural breaks at 1989, 2001, 2008 |

Median `B_pred/B_obs` > 1 means systematic scale offset (EEoS over-predicts), not a unit bug.

---

## Conceptual framing agents must respect

Two failure modes for METE/EEoS at haul level:

1. **Disturbance (sought signal):** fishing alters community structure so S, N, E no longer match METE assumptions.
2. **Catchability / unrepresentative hauls (confound):** hauls under-sample standing biomass or capture “big shoals” that violate METE aggregation assumptions.

Do not apply catchability corrections as primary inference without careful justification — they can erase or invent the disturbance signal. Dominance / size-CV explorations are **diagnostics**, not corrections applied to the live H2/H3 model.

See: `display_discussion/Methods_note_for_discussion.md`.

---

## Data not in GitHub

Restore helper:

```bash
bash local_env_upload_28-7/download_scripts/download_all.sh
```

| Dependency | In git? |
|------------|---------|
| `outputs/datras_hl_raw.rds` | Yes (~3.8 MB) — H1 can start here |
| Raw ICES DATRAS HL CSV (~251 MB) | No |
| `FishGlob_data/.../NS-IBTS_clean.RData` | No — clone/restore |
| `equation_of_state/biomass.py` | No — clone/restore |
| Couce fishing CSV | No (often gitignored) |
| `renv.lock` | Yes — `renv::restore()` |
| `.venv` | No — recreate from `local_env_upload_28-7/python_requirements.txt` |

Details: `local_env_upload_28-7/README.md`.

---

## Analysis status (as of Jul 2026)

| Workstream | Status |
|------------|--------|
| H1 haul-level EEoS | Complete |
| H2 rectangle SEM (earlier) | Complete; fishing not robust after spatial correction; biomass confound documented |
| H2/H3 within-between (primary) | Complete — presented results |
| Zone schemes / blended-term + GAM | Superseded → `exploratory/` |

---

## Agent conventions

1. Run scripts from **repo root**: `Rscript pipeline/...`
2. Treat `pipeline/` + `outputs/live_pipeline_run_log.md` as source of truth for presented work.
3. Do **not** present `exploratory/` as current results.
4. Prefer editing live scripts/docs over resurrecting superseded paths.
5. Do not delete git history; reorganize visibility only.
6. Python EEoS lives in `equation_of_state/biomass.py`; R orchestration in `pipeline/`.
7. When changing H2/H3 inference, update `display_discussion/` drafts and the live run log if script roles change.

---

## Best reading order for a new session

1. This file (`CLAUDE.md`)
2. `display_discussion/One page read me.md`
3. `outputs/live_pipeline_run_log.md`
4. `pipeline/README.md` (if running or editing code)
5. `display_discussion/H2_H3_methods_detailed.md` + `H2_H3_results_interpretation.md`
6. `AGENT_ONBOARDING.md` (clone/restore specifics)

---

## Related external repos

| Repo | URL |
|------|-----|
| **This project** | https://github.com/StuartStokeld/north_sea_EEoS |
| EEoS implementation | https://github.com/micbru/equation_of_state |
