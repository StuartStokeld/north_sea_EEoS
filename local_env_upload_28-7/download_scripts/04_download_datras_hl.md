# Download ICES DATRAS unaggregated HL data (~251 MB)

**Not stored in GitHub** — exceeds the 100 MB file limit and is excluded by `.gitignore`.

## What the pipeline needs

After download, place the export so this path resolves:

```
north_sea_eeos/
  Unaggregated trawl and biological information_<DATE>/
    Unaggregated trawl and biological information_<DATE>.csv
```

The pipeline auto-detects any folder matching `Unaggregated trawl and biological information_*` at the repo root (`pipeline/R/datras_csv_import.R::find_ices_hl_csv()`).

Alternatively, if you already have `outputs/datras_hl_raw.rds` (3.8 MB, **is** in the GitHub repo), you can skip the CSV entirely — the H1 pipeline reads the RDS directly.

## How to obtain the CSV

1. Go to [ICES DATRAS](https://datras.ices.dk/) → **Data Products** → **Unaggregated trawl and biological information**.
2. Request/export for **NS-IBTS**, **Quarter 1**, **years 1985–2015** (match `pipeline/R/datras_constants.R`).
3. Download the zip/csv export to your machine.
4. Unzip into a folder at the repo root named like the ICES export (e.g. `Unaggregated trawl and biological information_2026-06-01 18_45_33/`).

## Import into the pipeline

```bash
cd /path/to/north_sea_eeos
Rscript pipeline/import_datras_hl_from_csv.R
# or rely on build_datras_state_variables.R to import automatically
```

## Local reference (this machine, 28 Jul 2026)

- Original folder: `Unaggregated trawl and biological information_2026-06-01 18_45_33/`
- CSV size: ~251 MB
- Processed RDS in repo: `outputs/datras_hl_raw.rds`
