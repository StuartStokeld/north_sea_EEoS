H1 — EEoS haul-level biomass prediction: H1_results package
============================================================

Contents
--------
  display_discussion/H1_results_summary.Rmd   Report source (knit to HTML)
  display_discussion/H1_results_summary.html  Pre-knitted copy (read-only backup)
  display_discussion/H1_results_summary.md    Markdown companion (email-friendly)
  display_discussion/H1_metrics_review_H2_H3.md  Which H1 stats carry to H2/H3 (DV = mean absolute residual)
  outputs/haul_eeos_predictions.rds           Raw EEoS predictions (S, N, E, B_obs, B_pred)
  outputs/haul_h1_benchmarks.rds              Haul data used by the report (+ ln(E) benchmark)
  outputs/h1_harte_baseline_metrics.csv   Harte unfitted baseline metrics (Tier 1)
  outputs/h1_model_comparison.rds             Unified model comparison (Tier 1 + Tier 2)
  outputs/h1_lne_coefficients.csv             ln(E) OLS regression coefficients (Tier 2)
  outputs/null_distribution_summary.csv       Null permutation summary (referenced in report)
  outputs/h1_model_comparison.csv             Human-readable comparison table (optional)

Folder layout
-------------
Keep this structure when unzipping. The Rmd resolves paths relative to
display_discussion/../outputs/ — do not flatten into a single folder.

How to knit
-----------
1. Unzip this archive (or open the H1_results/ folder in the repo).
2. Open display_discussion/H1_results_summary.Rmd in RStudio (or run from R).
3. Install packages if needed:

     install.packages(c("dplyr", "ggplot2", "kableExtra", "patchwork", "rmarkdown"))

4. Knit to HTML (Knit button, or):

     rmarkdown::render("display_discussion/H1_results_summary.Rmd")

   Run from the directory that contains both display_discussion/ and outputs/.

Data notes
----------
  haul_eeos_predictions.rds   n = 12,069 hauls; NS-IBTS Q1 1985–2015
  haul_h1_benchmarks.rds    Same hauls with ln(E) benchmark columns added
  B_obs, B_pred               Both in grams

No Python or pipeline scripts are required to knit this report.

Contact
-------
Stuart Stokeld — North Sea EEoS PhD project
