# ARCHIVED (exploratory) — not part of the live methods model list.
# Re-run only intentionally from exploratory/pipeline/.
# Historical outputs live under exploratory/outputs/ (update write paths before re-running).
#
# Structure bathymetry/anisotropy outputs under pre-registered decision rules.
#
# Does not re-fit models or change thresholds. Reads completed diagnostic CSVs,
# adds the pre-registered permutation p-value and global-compass foil, and writes
# the locked report tables + bearing plot + run log.
#
# Run: Rscript --vanilla pipeline/run_h2h3_anisotropy_preregistered_report.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run: Rscript pipeline/run_h2h3_anisotropy_preregistered_report.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2h3_bathymetry_anisotropy_helpers.R"))
ensure_bathymetry_libs(project_root)

path_bearings <- file.path(project_root, "outputs", "bearings_by_rectangle.csv")
path_bathy <- file.path(project_root, "outputs", "bathymetry_by_rectangle.csv")
stopifnot(file.exists(path_bearings), file.exists(path_bathy))

path_summary <- file.path(project_root, "outputs", "anisotropy_results_summary.md")
path_plot <- file.path(project_root, "outputs", "anisotropy_bearing_comparison_plot.png")
path_run_log <- file.path(project_root, "outputs", "anisotropy_preregistered_report_run_log.md")

ALPHA <- 0.05
N_PERM <- 999L
SEED <- 20260804L

run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

yesno <- function(p, alpha = ALPHA) {
  if (!is.finite(p)) return("No (p non-finite)")
  if (p < alpha) "Yes" else "No"
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.0001) return(sprintf("%.2e", p))
  sprintf("%.4f", p)
}

fmt_rho <- function(r) {
  if (!is.finite(r)) return("NA")
  sprintf("%.4f", r)
}

logmsg("# Anisotropy pre-registered report — run log")
logmsg("")
logmsg(paste0("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")))
logmsg("**Inputs:** `outputs/bearings_by_rectangle.csv`, `outputs/bathymetry_by_rectangle.csv`")
logmsg("**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`")
logmsg("")
logmsg("## Pre-registered rules applied (exact)")
logmsg("")
logmsg(
  "1. **TID:** exclude rectangles with pct_sounding < 50% from confirmatory tests; ",
  "keep them in the bearing plot. Report 25% as labeled sensitivity only."
)
logmsg(
  "2. **Primary object:** signed residual-correlation bearing. BLUP = Table 2 ",
  "robustness only."
)
logmsg(
  "3. **Along-shelf axis (Table 1):** per-rectangle local depth-gradient direction ",
  "rotated +90° (along-shelf). Global compass split = Table 2 foil only."
)
logmsg(
  "4. **Decision:** Jammalamadaka–Sarma via `circular::cor.circular(..., test=TRUE)`. ",
  "Table 1 verdict uses **asymptotic p < 0.05**. Permutation p (999 shuffles of the ",
  "bearing pairing) reported alongside; coefficient magnitude is descriptive only."
)
logmsg("")
logmsg("### Foil operationalization (Table 2)")
logmsg(
  "Global compass along-shelf (depth-independent): 90° (N–S) for lon < 4°E; ",
  "0° (E–W) for lon ≥ 4°E. Same JS test and p < 0.05 labeling as primary, but ",
  "**not confirmatory**."
)

bear <- utils::read.csv(path_bearings, stringsAsFactors = FALSE)
bathy <- utils::read.csv(path_bathy, stringsAsFactors = FALSE)
# prefer tid flags from bearings (already joined); cross-check bathy
n_flag_50 <- sum(bear$tid_flag_lt50)
n_flag_25 <- sum(bear$tid_flag_lt25)
logmsg("")
logmsg(sprintf("N rectangles: %d", nrow(bear)))
logmsg(sprintf("TID flagged <50%%: %d", n_flag_50))
logmsg(sprintf("TID flagged <25%%: %d", n_flag_25))

use50 <- !bear$tid_flag_lt50
use25 <- !bear$tid_flag_lt25

# --- Table 1 primary: resid × local along-shelf depth ---
logmsg("")
logmsg("## Table 1 computation")
logmsg(
  "Pairing: bearing_resid × bearing_depth (= grad_bearing_along = depth-gradient + 90°). ",
  "Sample: TID-OK at 50% threshold."
)
primary <- jammalamadaka_sarma_permute(
  bear$bearing_resid[use50],
  bear$bearing_depth[use50],
  n_perm = N_PERM,
  seed = SEED
)
logmsg(sprintf(
  "Primary: n=%d  rho=%s  p_asymp=%s  p_perm=%s  confirmed=%s",
  primary$n, fmt_rho(primary$rho), fmt_p(primary$p_asymptotic),
  fmt_p(primary$p_permutation), yesno(primary$p_asymptotic)
))

# --- Table 2 rows ---
logmsg("")
logmsg("## Table 2 computations")

blup <- jammalamadaka_sarma_permute(
  bear$bearing_blup[use50],
  bear$bearing_depth[use50],
  n_perm = N_PERM,
  seed = SEED + 1L
)
logmsg(sprintf(
  "BLUP robustness: n=%d  rho=%s  p_asymp=%s  p_perm=%s  label=%s",
  blup$n, fmt_rho(blup$rho), fmt_p(blup$p_asymptotic),
  fmt_p(blup$p_permutation), yesno(blup$p_asymptotic)
))

bear$bearing_global_compass <- global_compass_along_shelf_deg(bear$lon)
foil <- jammalamadaka_sarma_permute(
  bear$bearing_resid[use50],
  bear$bearing_global_compass[use50],
  n_perm = N_PERM,
  seed = SEED + 2L
)
logmsg(sprintf(
  "Global compass foil: n=%d  rho=%s  p_asymp=%s  p_perm=%s  label=%s",
  foil$n, fmt_rho(foil$rho), fmt_p(foil$p_asymptotic),
  fmt_p(foil$p_permutation), yesno(foil$p_asymptotic)
))

sens <- jammalamadaka_sarma_permute(
  bear$bearing_resid[use25],
  bear$bearing_depth[use25],
  n_perm = N_PERM,
  seed = SEED + 3L
)
logmsg(sprintf(
  "TID 25%% sensitivity: n=%d  rho=%s  p_asymp=%s  p_perm=%s  label=%s",
  sens$n, fmt_rho(sens$rho), fmt_p(sens$p_asymptotic),
  fmt_p(sens$p_permutation), yesno(sens$p_asymptotic)
))

# --- Verdict narrative (locked branches) ---
prim_conf <- is.finite(primary$p_asymptotic) && primary$p_asymptotic < ALPHA
foil_conf <- is.finite(foil$p_asymptotic) && foil$p_asymptotic < ALPHA

if (prim_conf && !foil_conf) {
  narrative <- paste0(
    "Correlation has a real, finite range, but geographic distance was the wrong ",
    "metric — local shelf geometry explains the directional pattern the isotropic ",
    "distance-decay model couldn't detect."
  )
  verdict_label <- "confirmed_shelf_geometry"
} else if (!prim_conf) {
  narrative <- paste0(
    "No evidence that shelf-geometry-based direction explains the earlier ",
    "unidentified decay range. Consistent with either a genuinely borderless ",
    "spatial field or confounding with the fixed-effect spatial trend — does not ",
    "distinguish between these two."
  )
  verdict_label <- "not_confirmed"
} else {
  # Table 1 confirms AND foil also confirms
  narrative <- paste0(
    "Ambiguous result requiring follow-up: the primary local depth-gradient ",
    "alignment test meets p < 0.05, but the global compass foil also meets ",
    "p < 0.05 under the same rule — do not force into either verdict category."
  )
  verdict_label <- "ambiguous_primary_and_foil"
}

offset_note <- paste0(
  "A confirmed alignment is expected to reflect a consistent ~90° rotational ",
  "offset between residual-correlation bearing and the raw depth-gradient ",
  "(cross-shelf) direction: Table 1 pairs residual bearing with the along-shelf ",
  "axis (gradient + 90°), so confirmation means those axes coincide — not that ",
  "residual correlation points the same compass way as the steepest depth slope."
)

logmsg("")
logmsg("## Verdict branch")
logmsg("verdict_label: ", verdict_label)
logmsg("narrative: ", narrative)

# --- Plot (all rectangles; TID flagged marked, not removed) ---
save_bearing_comparison_plot(bear, path_plot)
logmsg("")
logmsg("Wrote plot: ", path_plot)
logmsg(
  "Plot includes all ", nrow(bear), " rectangles; ", n_flag_50,
  " TID-flagged (<50%) shown as triangles."
)

# --- Summary markdown ---
summary_md <- c(
  "# Bathymetry / anisotropy results — pre-registered decision rules",
  "",
  paste0("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste0(
    "**Source run:** `outputs/h2h3_bathymetry_anisotropy_run_log.md`; ",
    "this report: `outputs/anisotropy_preregistered_report_run_log.md`"
  ),
  paste0("**Design:** `display_discussion/Design_bathymetry_spatial_anisotropy.md`"),
  "",
  "Decision rule for “Anisotropy confirmed?”: **asymptotic p < 0.05** on the ",
  "Jammalamadaka–Sarma circular–circular correlation. Coefficient magnitude is ",
  "descriptive only. Permutation p (999 pairing shuffles) is reported for the ",
  "primary row as pre-registered robustness of the p-value, not a substitute rule.",
  "",
  "## Table 1 — Headline result",
  "",
  "| Test object | Circular correlation coefficient | p (asymptotic) | p (permutation) | Anisotropy confirmed? |",
  "|---|---:|---:|---:|---|",
  sprintf(
    "| Signed residuals × local depth-gradient bearing (primary) | %s | %s | %s | %s |",
    fmt_rho(primary$rho),
    fmt_p(primary$p_asymptotic),
    fmt_p(primary$p_permutation),
    yesno(primary$p_asymptotic)
  ),
  "",
  paste0(
    "Primary pairing: residual-correlation bearing vs local **along-shelf** axis ",
    "(= per-rectangle depth-gradient direction + 90°). N = ", primary$n,
    " after excluding TID-flagged rectangles (<50% direct soundings)."
  ),
  "",
  "## Table 2 — Robustness / contrast rows",
  "",
  "*Not confirmatory. Same p < 0.05 labeling for transparency only.*",
  "",
  "| Variant | Coefficient | p (asymptotic) | Confirmed under same rule? | Role |",
  "|---|---:|---:|---|---|",
  sprintf(
    "| BLUPs × local depth-gradient bearing | %s | %s | %s | Robustness check |",
    fmt_rho(blup$rho), fmt_p(blup$p_asymptotic), yesno(blup$p_asymptotic)
  ),
  sprintf(
    "| Signed residuals × global compass bearing (foil) | %s | %s | %s | Foil — not confirmatory |",
    fmt_rho(foil$rho), fmt_p(foil$p_asymptotic), yesno(foil$p_asymptotic)
  ),
  sprintf(
    "| TID 25%% threshold (sensitivity) | %s | %s | %s | Sensitivity on exclusion rule |",
    fmt_rho(sens$rho), fmt_p(sens$p_asymptotic), yesno(sens$p_asymptotic)
  ),
  "",
  paste0(
    "Foil definition: depth-independent piecewise compass along-shelf — 90° (N–S) ",
    "for lon < 4°E, 0° (E–W) for lon ≥ 4°E. Permutation p for foil = ",
    fmt_p(foil$p_permutation), "; for BLUP = ", fmt_p(blup$p_permutation),
    "; for TID-25% = ", fmt_p(sens$p_permutation), "."
  ),
  "",
  "## Table 3 — TID exclusion summary",
  "",
  sprintf("| Threshold | N flagged | Role |"),
  "|---|---:|---|",
  sprintf(
    "| <50%% direct-sounding (primary exclusion) | %d | Excluded from Table 1 / confirmatory tests; **kept in bearing plot** |",
    n_flag_50
  ),
  sprintf(
    "| <25%% direct-sounding (sensitivity) | %d | Labels the Table 2 sensitivity row only |",
    n_flag_25
  ),
  "",
  paste0(
    "Total rectangles = ", nrow(bear), ". Table 1 N = ", primary$n,
    ". Deliverable plot `outputs/anisotropy_bearing_comparison_plot.png` marks the ",
    n_flag_50, " TID-flagged rectangles (triangles) and does not remove them."
  ),
  "",
  "## Narrative verdict",
  "",
  narrative,
  "",
  offset_note,
  "",
  paste0("**Verdict label:** `", verdict_label, "`"),
  "",
  "## Deliverables",
  "",
  "- `outputs/anisotropy_results_summary.md` (this file)",
  "- `outputs/anisotropy_bearing_comparison_plot.png`",
  "- `outputs/anisotropy_preregistered_report_run_log.md`"
)
writeLines(summary_md, path_summary)
logmsg("Wrote summary: ", path_summary)

writeLines(run_log, path_run_log)
logmsg("Wrote run log: ", path_run_log)

invisible(NULL)
