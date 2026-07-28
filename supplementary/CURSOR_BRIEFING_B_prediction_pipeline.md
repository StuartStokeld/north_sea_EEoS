# Cursor Briefing: RFEC B Prediction Pipeline
## Task: Build the EEoS haul-level B prediction pipeline in R

---

## 1. Project context

This is a PhD research project (RFEC — Reference Free Ecosystem Condition) testing whether the Ecological Equation of State (EEoS), derived from Maximum Entropy Theory of Ecology (METE), predicts community biomass in North Sea demersal fish communities.

EEoS predicts total community biomass B from three state variables:
- **S** — species richness per haul
- **N** — total abundance per haul
- **E** — total community metabolic rate per haul (normalised, dimensionless)

The pipeline you are building assembles S, N, E, and observed B per haul, calls the Python EEoS function to get predicted B, then computes residuals. This is Hypothesis 1 (H1).

---

## 2. Project file locations

All files are under `/Users/stuartstokeld/north_sea_eeos/`

### Input files already built

| File | Description |
|------|-------------|
| `outputs/datras_haul_E.rds` | Haul-level E (computed, normalisation still needed — see Section 5) |
| `outputs/datras_haul_mean_length.rds` | Species × haul mean length in mm |
| `outputs/fishbase_lw_lookup_v2.csv` | FishBase LW parameters per species (a, b, fallback_level) |
| `FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData` | FishGlob NS-IBTS haul data — loads object named `data` |

### Output files to create

| File | Description |
|------|-------------|
| `outputs/haul_state_variables.rds` | S, N, E (normalised), B per haul — the assembled input table |
| `outputs/haul_eeos_predictions.rds` | S, N, E, B_obs, B_pred, residual per haul |

---

## 3. FishGlob data structure

Load with:
```r
load("/Users/stuartstokeld/north_sea_eeos/FishGlob_data/outputs/Cleaned_data/NS-IBTS_clean.RData")
# Object is named 'data', 364,196 rows × 42 columns
```

Key columns:
- `haul_id` — composite haul identifier, space-separated string in format:
  `"NS-IBTS 2010 1 DE 06NI GOV 107 5"` (Survey Year Quarter Country Ship Gear StNo HaulNo)
- `accepted_name` — species name
- `num_cpue` — abundance (individuals per km²) — **use this for N, not raw `num`**
- `wgt` — biomass in grams per haul per species
- `stat_rec` — ICES rectangle code (e.g. "44G1")
- `year`, `quarter`

Derive per-haul S, N, B:
```r
fishglob_haul <- data %>%
  filter(quarter == 1, year >= 1985, year <= 2015) %>%
  group_by(haul_id) %>%
  summarise(
    S     = n_distinct(accepted_name),
    N     = sum(num_cpue, na.rm = TRUE),
    B_obs = sum(wgt, na.rm = TRUE),
    stat_rec = first(stat_rec),
    year     = first(year),
    .groups = "drop"
  )
```

---

## 4. haul_E data structure

`datras_haul_E.rds` has these columns:
- `haul_key` — composite key in format: `"NS-IBTS_1985_1_DE_06DA_1"` (underscore-separated: Survey_Year_Quarter_Country_Platform_HaulNumber)
- `E` — **unnormalised** haul-level metabolic rate (see Section 5)
- `Survey`, `Year`, `Quarter`, `Country`, `Platform`, `HaulNumber`
- `n_species_with_lw` — number of species contributing to E

---

## 5. Join key between FishGlob and haul_E

The FishGlob `haul_id` uses space separators and 8 fields:
`"NS-IBTS 2010 1 DE 06NI GOV 107 5"`
→ Survey, Year, Quarter, Country, Ship, Gear, StNo, HaulNo

The DATRAS `haul_key` uses underscore separators and 6 fields:
`"NS-IBTS_2010_1_DE_06NI_5"`
→ Survey, Year, Quarter, Country, Platform, HaulNumber

**The join drops Gear and StNo** — these are not in the DATRAS HL download.
Build a 6-field join key from FishGlob by parsing haul_id:

```r
library(stringr)

fishglob_haul <- fishglob_haul %>%
  mutate(
    haul_parts = str_split(haul_id, " "),
    fg_survey   = map_chr(haul_parts, 1),
    fg_year     = map_chr(haul_parts, 2),
    fg_quarter  = map_chr(haul_parts, 3),
    fg_country  = map_chr(haul_parts, 4),
    fg_ship     = map_chr(haul_parts, 5),
    # field 6 is Gear — skip
    # field 7 is StNo — skip
    fg_haulno   = map_chr(haul_parts, 8),
    join_key    = paste(fg_survey, fg_year, fg_quarter, 
                        fg_country, fg_ship, fg_haulno, sep = "_")
  )
```

Then join on `join_key == haul_key`.

---

## 6. E normalisation — CRITICAL

The current `haul_E$E` is **not yet normalised**. The EEoS implementation (`micbru/equation_of_state`) requires E to be dimensionless — normalised so that the smallest individual metabolic rate in the haul equals 1.

The normalisation formula is:
```
E_normalised = E_raw / m_min^(3/4)
```

Where `m_min` is the mass in grams of the smallest individual across all species in that haul.

To apply this you need to:
1. Load `datras_haul_mean_length.rds` — species × haul mean length in mm
2. Join to `fishbase_lw_lookup_v2.csv` on AphiaID → aphia_id to get a and b
3. Compute `mass_g = a * (length_mm/10)^b` per species per haul
4. Find `m_min` = minimum mass_g across all species in each haul
5. Compute `m_min_epsilon = m_min^0.75`
6. `E_normalised = E_raw / m_min_epsilon`

**Validation check**: after normalisation, the minimum epsilon (m^0.75) in any haul should equal exactly 1.

---

## 7. EEoS Python function

The reference implementation is at: https://github.com/micbru/equation_of_state

The key function is in `equation_of_state.py`. It takes S, N, E (normalised) and returns predicted B.

Call it from R using `reticulate`:
```r
library(reticulate)
# Point to your Python environment if needed
# use_python("/usr/local/bin/python3")

source_python("/path/to/equation_of_state.py")

# Call per haul — the function signature is approximately:
# B_pred = equation_of_state(S, N, E)
# Confirm exact function name and signature by reading equation_of_state.py
```

You will need to clone the repo or have it available locally:
```bash
git clone https://github.com/micbru/equation_of_state.git
```

The repo is likely already present at:
`/Users/stuartstokeld/north_sea_eeos/BioticHomogenization-main/` — check if equation_of_state.py is in there, or in a separate clone.

---

## 8. Expected output schema

`haul_eeos_predictions.rds` should have one row per haul with:

| Column | Description |
|--------|-------------|
| `haul_id` | FishGlob haul identifier |
| `haul_key` | DATRAS haul key |
| `year` | Survey year |
| `stat_rec` | ICES rectangle |
| `S` | Species richness |
| `N` | Total abundance (num_cpue) |
| `E` | Normalised metabolic rate |
| `B_obs` | Observed biomass (grams, from FishGlob wgt) |
| `B_pred` | EEoS predicted biomass |
| `ln_B_obs` | log(B_obs) |
| `ln_B_pred` | log(B_pred) |
| `residual` | ln_B_obs − ln_B_pred |
| `abs_residual` | abs(residual) |
| `n_species_with_lw` | Species count contributing to E |

---

## 9. Filters to apply before EEoS evaluation

- Quarter == 1 only
- Year 1985–2015 only
- Exclude hauls where S < 2 (EEoS undefined for single-species communities)
- Exclude hauls where E <= N after normalisation (mathematically undefined: lambda_2 = S/(E-N) requires E > N)
- Exclude hauls where B_obs == 0 or is NA

---

## 10. Implementation notes

- Work in R throughout; use `reticulate` only for the EEoS function call
- The EEoS function may need to be called row-by-row (use `purrr::map` or `rowwise()`) — check whether it is vectorised
- Wrap each EEoS call in `tryCatch` — some hauls may fail to converge; log these separately rather than erroring out
- Save failed hauls to `outputs/eeos_failed_hauls.rds` with the error message
- Report convergence failure rate in the output summary

---

## 11. Suggested R script structure

Create: `/Users/stuartstokeld/north_sea_eeos/build_eeos_predictions.R`

Sections:
1. Load libraries and data
2. Derive FishGlob haul-level S, N, B
3. Build join key from haul_id
4. Join FishGlob to haul_E
5. Compute E normalisation (m_min per haul)
6. Apply filters
7. Call EEoS function per haul
8. Compute residuals
9. Save outputs and print summary statistics

---

## 12. Summary statistics to print at end

```r
cat("Total hauls assembled:", nrow(haul_joined), "\n")
cat("Hauls passing all filters:", nrow(haul_filtered), "\n")
cat("Hauls with successful EEoS prediction:", nrow(haul_predictions), "\n")
cat("Hauls failed to converge:", nrow(failed_hauls), "\n")
cat("E normalisation check — min epsilon == 1:", 
    all(abs(haul_predictions$min_epsilon - 1) < 1e-10), "\n")
cat("R² (log scale):", ..., "\n")
cat("Median absolute residual:", median(haul_predictions$abs_residual), "\n")
```
