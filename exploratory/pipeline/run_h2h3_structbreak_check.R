# ARCHIVED (exploratory): Bai–Perron / BIC structural-break discovery for the
# original data-driven phases. Not part of current primary (phase_v2) or
# sensitivity reporting. Re-run only intentionally:
#   Rscript --vanilla exploratory/pipeline/run_h2h3_structbreak_check.R
#
# Structural-break check — fishing-pressure time series.
# See CURSOR_BRIEFING "Structural Break Check — Fishing Pressure Time Series"
# (chat-supplied, not a repo file) for the full spec this script implements.
#
# PURPOSE: Section A of the H2/H3 design-support task (run_h2h3_design_support.R)
# described a four-phase pattern in annual mean fishing hours (rise to 1990,
# gradual decline through the 2000s, levelling off 2010-2015) with boundaries
# chosen BY EYE. This script checks whether that pattern is statistically
# supported by a formal structural-break method, using the method's own
# model-selection criterion (BIC) to pick the number of breaks — NOT assumed
# in advance to match the 4 visual phases.
#
# THIS SCRIPT DOES NOT: force a match to the four-phase pattern, adjust the
# visual boundaries to fit the statistical result, re-fit the within-phase
# trends (A3) under any new break structure, or recommend which set of
# boundaries (visual or statistical) to use going forward. Read-only on
# outputs/h2h3_designA1_year_fishing_summary.csv (not recomputed from raw
# Couce data).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = TRUE)
} else if (dir.exists("R")) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  stop("Run: Rscript --vanilla exploratory/pipeline/run_h2h3_structbreak_check.R")
}
source(file.path(script_dir, "R", "h1_common.R"))
ctx <- load_pipeline_context()
project_root <- ctx$project_root
script_dir <- ctx$script_dir
source(file.path(script_dir, "R", "h2h3_structbreak_helpers.R"))

# ---------------------------------------------------------------------------
# Provisional constants / method choices — named once, surfaced in the run log.
# ---------------------------------------------------------------------------
VISUAL_BOUNDARIES <- c(1990L, 2000L, 2010L) # from the briefing: eyeballed boundaries behind Section A's 4-phase description
BREAK_FORMULA <- mean_hours ~ year          # breaks allow BOTH intercept and slope to change (matches A3's per-phase linear-trend framing); NOT a pure mean-shift ("~1") model
STRUCCHANGE_H <- 0.15                       # strucchange::breakpoints() package default minimum-segment fraction; NOT tuned for this task
CANDIDATE_BREAKS_MIN_REQUIRED <- 5L         # brief requires reporting >=0-5 breaks; actual max reported is whatever STRUCCHANGE_H allows (see run log)

# ---------------------------------------------------------------------------
# Paths (archived under exploratory/; live input A1 stays at top-level outputs/)
# ---------------------------------------------------------------------------
out_root <- file.path(project_root, "exploratory", "outputs")
path_in_A1 <- file.path(project_root, "outputs", "h2h3_designA1_year_fishing_summary.csv")
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_out_criterion <- file.path(out_root, "h2h3_designA4_structbreak_criterion.csv")
path_out_years <- file.path(out_root, "h2h3_designA4_structbreak_years.csv")
path_out_fig <- file.path(fig_dir, "h2h3_designA4_structbreak_series.png")
path_out_run_log <- file.path(out_root, "h2h3_designA4_structbreak_run_log.md")

stopifnot(file.exists(path_in_A1))

# ---------------------------------------------------------------------------
# Run log accumulator
# ---------------------------------------------------------------------------
run_log <- character(0)
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  run_log <<- c(run_log, msg)
}

logmsg("# Structural-break check — fishing-pressure time series — run log")
logmsg("")
logmsg("Checks whether Section A's eyeballed 4-phase pattern is statistically supported. Reports the comparison as numbers only — does not force a match, does not adjust the visual boundaries, does not recommend which set of boundaries to use.")

# ---------------------------------------------------------------------------
# Package availability
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Method / package")
if (!requireNamespace("strucchange", quietly = TRUE)) {
  stop(
    "Package 'strucchange' is not installed and no fallback was implemented in this script. ",
    "Per the briefing, an equivalent well-established structural-break method could be substituted, ",
    "but strucchange::breakpoints() is the one actually used below — install it (e.g. ",
    "install.packages('strucchange')) before re-running."
  )
}
library(strucchange)
logmsg(
  "Package: strucchange (version ", as.character(utils::packageVersion("strucchange")), ", Bai & Perron dynamic-",
  "programming least-squares breakpoint estimation, `breakpoints()`). NOTE: strucchange was NOT part of this ",
  "project's renv-managed dependency set at the time this script was written — it was installed ad hoc ",
  "(`install.packages('strucchange', repos = 'https://cloud.r-project.org')`) into the ambient R library, not ",
  "added to renv.lock. Flagged as an environment note, not resolved here (out of scope for this task to restructure ",
  "the project's dependency management)."
)
logmsg(
  "Model: `", deparse(BREAK_FORMULA), "` — breaks allow BOTH the intercept and the linear slope (hours/year) to ",
  "change at each breakpoint, matching how Section A's within-phase trends (A3) are themselves framed (separate ",
  "linear fits per phase), rather than a pure mean-shift ('~1') model. This is a method/formula choice, not a ",
  "value dictated by the briefing, so it is surfaced explicitly here."
)

# ---------------------------------------------------------------------------
# Data (same series as Section A — read directly, not recomputed)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Data")
year_summary <- read_csv(path_in_A1, show_col_types = FALSE) %>% arrange(year)
logmsg(
  "Read ", path_in_A1, " directly (", nrow(year_summary), " year rows, ", min(year_summary$year), "-",
  max(year_summary$year), "), NOT recomputed from raw Couce data — same full ", "215-rectangle Couce universe as ",
  "Section A. Series used: `mean_hours` (mean annual Couce fishing hours across all rectangles with coverage that year)."
)

# ---------------------------------------------------------------------------
# Fit breakpoints() and build the criterion table (item 2)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Break-count model-selection criterion (item 2)")

bp_full <- breakpoints(BREAK_FORMULA, data = year_summary, h = STRUCCHANGE_H)
criterion_table <- build_breakpoint_criterion_table(bp_full)
n_breaks_max <- max(criterion_table$n_breaks)
n_breaks_opt <- criterion_table$n_breaks[criterion_table$is_bic_optimal][1]

logmsg(
  "`h` = ", STRUCCHANGE_H, " (package default; minimum segment length = h*n = ", round(STRUCCHANGE_H * nrow(year_summary), 2),
  " ~= ", ceiling(STRUCCHANGE_H * nrow(year_summary)), " observations) determines the maximum number of breaks ",
  "the dynamic program can evaluate for n=", nrow(year_summary), " observations: max feasible = ", n_breaks_max,
  " breaks (>= the ", CANDIDATE_BREAKS_MIN_REQUIRED, " the briefing asked for; all ", n_breaks_max + 1L,
  " candidate counts from 0 to ", n_breaks_max, " are reported below, not just 0-5)."
)
write_csv(criterion_table, path_out_criterion)
logmsg("Saved: ", path_out_criterion)
for (i in seq_len(nrow(criterion_table))) {
  r <- criterion_table[i, ]
  logmsg(sprintf("  - m=%d breaks: RSS=%.3e, BIC=%.4f%s", r$n_breaks, r$rss, r$bic, ifelse(r$is_bic_optimal, "  <- BIC-optimal", "")))
}
logmsg(
  "BIC-optimal break count = ", n_breaks_opt, ". NOTE: BIC at m=", n_breaks_opt, " (",
  sprintf("%.4f", criterion_table$bic[criterion_table$n_breaks == n_breaks_opt]), ") is only marginally lower than ",
  "at m=", n_breaks_opt - 1L, " (", sprintf("%.4f", criterion_table$bic[criterion_table$n_breaks == n_breaks_opt - 1L]),
  ") — reported as a numeric fact about how decisive the BIC minimum is, not interpreted further."
)

# ---------------------------------------------------------------------------
# Break years + CIs at the optimal count (item 3)
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Selected break year(s) at the BIC-optimal count, with 95% CIs (item 3)")

year_table <- build_breakpoint_year_table(bp_full, n_breaks_opt, year_summary$year)
write_csv(year_table, path_out_years)
logmsg("Saved: ", path_out_years)
for (i in seq_len(nrow(year_table))) {
  r <- year_table[i, ]
  ci_str <- if (r$ci_computed) sprintf("[%d, %d]", r$ci_lower_year, r$ci_upper_year) else "NA (confint() could not be computed)"
  logmsg(sprintf("  - break %d: year = %d, 95%% CI = %s", r$break_number, r$break_year, ci_str))
}

# ===========================================================================
# Comparison against the visual boundaries (item 4) — numbers only
# ===========================================================================
logmsg("")
logmsg("## Comparison: statistically located break years vs visual boundaries (item 4)")
logmsg(
  "Visual boundaries (from the briefing, eyeballed from h3_pre_D1/D3 and Section A): ",
  paste(VISUAL_BOUNDARIES, collapse = ", "), ". Statistically located break years (BIC-optimal, m=", n_breaks_opt,
  "): ", paste(year_table$break_year, collapse = ", "), "."
)
comparison_table <- compare_breaks_to_visual(year_table$break_year, VISUAL_BOUNDARIES)
for (i in seq_len(nrow(comparison_table))) {
  r <- comparison_table[i, ]
  logmsg(sprintf("  - [%s] %d -> nearest %d (|diff| = %d years)", r$direction, r$reference_year, r$nearest_other_year, r$abs_diff_years))
}

# ===========================================================================
# Figure: series with statistical breaks + visual boundaries marked
# ===========================================================================
logmsg("")
logmsg("## Figure")

break_lines <- tibble(year = year_table$break_year, kind = "Statistical break (BIC-optimal)")
visual_lines <- tibble(year = VISUAL_BOUNDARIES, kind = "Visual boundary (eyeballed)")
ci_bands <- year_table %>%
  filter(ci_computed) %>%
  transmute(xmin = ci_lower_year, xmax = ci_upper_year)

p <- ggplot(year_summary, aes(x = year, y = mean_hours))
if (nrow(ci_bands) > 0L) {
  p <- p + geom_rect(
    data = ci_bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#b2182b", alpha = 0.12
  )
}
p <- p +
  geom_line(colour = "#08306b", linewidth = 1) +
  geom_point(colour = "#08306b", size = 1.6) +
  geom_vline(data = visual_lines, aes(xintercept = year, colour = kind), linetype = "dashed", linewidth = 0.9) +
  geom_vline(data = break_lines, aes(xintercept = year, colour = kind), linetype = "solid", linewidth = 0.9) +
  scale_colour_manual(values = c("Visual boundary (eyeballed)" = "grey40", "Statistical break (BIC-optimal)" = "#b2182b"), name = NULL) +
  labs(
    x = "Year", y = "Mean Couce fishing hours (215-rectangle universe)",
    title = "Structural-break check: statistical breaks vs visual phase boundaries",
    subtitle = sprintf("strucchange::breakpoints(mean_hours ~ year), BIC-optimal m=%d breaks; shaded = 95%% CI on break year", n_breaks_opt),
    caption = "Dashed grey = visual boundaries (1990/2000/2010, eyeballed from h3_pre_D1/D3). Solid red = statistically located breaks. Numbers only, no recommendation."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(path_out_fig, p, width = 9.5, height = 6, dpi = 150)
logmsg("Saved: ", path_out_fig)

# ---------------------------------------------------------------------------
# Outputs index
# ---------------------------------------------------------------------------
logmsg("")
logmsg("## Outputs")
logmsg("- ", path_out_criterion)
logmsg("- ", path_out_years)
logmsg("- ", path_out_fig)
logmsg("- ", path_out_run_log, " (this file)")

writeLines(run_log, path_out_run_log)
cat("\nSaved run log:", path_out_run_log, "\n")
cat("=== Structural-break check complete — numbers only, no boundary recommendation made. ===\n")
