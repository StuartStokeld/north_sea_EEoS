# Helpers for the H2/H3 within-between fishing-pressure decomposition
# (run_h2h3_within_between.R). Separates persistent between-rectangle fishing
# pressure (H2) from within-rectangle year-to-year deviations (H3).

suppressPackageStartupMessages({
  library(dplyr)
})

#' Add FP_between (rectangle mean of log_hours_total) and FP_within
#' (deviation from that mean). FP_within has mean ~0 within each rectangle
#' by construction.
add_fp_within_between <- function(dat) {
  required <- c("stat_rec", "log_hours_total")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("add_fp_within_between: missing columns: ", paste(missing, collapse = ", "))
  }
  dat %>%
    group_by(stat_rec) %>%
    mutate(
      FP_between = mean(log_hours_total, na.rm = TRUE),
      FP_within = log_hours_total - FP_between
    ) %>%
    ungroup()
}

#' Sanity checks on the decomposition: within-rectangle mean of FP_within
#' near zero; correlation of FP_between with FP_within near zero.
fp_within_between_checks <- function(dat) {
  within_means <- dat %>%
    group_by(stat_rec) %>%
    summarise(mean_within = mean(FP_within, na.rm = TRUE), .groups = "drop")
  cor_bw <- suppressWarnings(stats::cor(dat$FP_between, dat$FP_within, use = "complete.obs"))
  list(
    cor_between_within = cor_bw,
    max_abs_within_rect_mean = max(abs(within_means$mean_within), na.rm = TRUE),
    mean_abs_within_rect_mean = mean(abs(within_means$mean_within), na.rm = TRUE),
    n_rectangles = nrow(within_means),
    var_between = stats::var(dat$FP_between, na.rm = TRUE),
    var_within = stats::var(dat$FP_within, na.rm = TRUE)
  )
}

#' True if `term` is an interaction of `component` with a phase level
#' (either `component:phaseX` or `phaseX:component` — glmmTMB may emit either).
is_fp_phase_interaction <- function(term, component) {
  grepl(paste0("^", component, ":phase"), term) |
    grepl(paste0("^phase[^:]+:", component, "$"), term)
}

#' Label within-between fixed-effect terms by hypothesis.
label_wb_terms <- function(term) {
  dplyr::case_when(
    term == "FP_between" ~
      "H2 (between-rectangle FP, reference phase 1985-1988)",
    is_fp_phase_interaction(term, "FP_between") ~
      "H2 (between-rectangle FP x phase interaction vs reference)",
    term == "FP_within" ~
      "H3 (within-rectangle FP deviation, reference phase 1985-1988)",
    is_fp_phase_interaction(term, "FP_within") ~
      "H3 (within-rectangle FP x phase interaction vs reference)",
    grepl("^phase", term) ~ "control (phase main effect)",
    term == "(Intercept)" ~ "control (intercept)",
    TRUE ~ "other"
  )
}

#' Plain-language term labels for the within-between FE table.
annotate_wb_fixed_effects <- function(fe_table) {
  fe_table %>%
    mutate(
      hypothesis = label_wb_terms(term),
      hypothesis_group = dplyr::case_when(
        grepl("^H2", hypothesis) ~ "H2_spatial_between",
        grepl("^H3", hypothesis) ~ "H3_temporal_within",
        TRUE ~ "control"
      ),
      term_plain = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "FP_between" ~
          "Between-rectangle FP (rectangle mean log(hours+1)), phase 1985-1988",
        term == "FP_within" ~
          "Within-rectangle FP deviation (year - rectangle mean), phase 1985-1988",
        term %in% c("phase1989-2000") ~ "Phase 1989-2000 (vs 1985-1988)",
        term %in% c("phase2001-2007") ~ "Phase 2001-2007 (vs 1985-1988)",
        term %in% c("phase2008-2015") ~ "Phase 2008-2015 (vs 1985-1988)",
        grepl("FP_between", term) & grepl("phase1989-2000", term) ~
          "Between-FP x phase 1989-2000 (change in between slope vs 1985-1988)",
        grepl("FP_between", term) & grepl("phase2001-2007", term) ~
          "Between-FP x phase 2001-2007 (change in between slope vs 1985-1988)",
        grepl("FP_between", term) & grepl("phase2008-2015", term) ~
          "Between-FP x phase 2008-2015 (change in between slope vs 1985-1988)",
        grepl("FP_within", term) & grepl("phase1989-2000", term) ~
          "Within-FP x phase 1989-2000 (change in within slope vs 1985-1988)",
        grepl("FP_within", term) & grepl("phase2001-2007", term) ~
          "Within-FP x phase 2001-2007 (change in within slope vs 1985-1988)",
        grepl("FP_within", term) & grepl("phase2008-2015", term) ~
          "Within-FP x phase 2008-2015 (change in within slope vs 1985-1988)",
        TRUE ~ term
      )
    )
}

#' Resolve the interaction coefficient name for `term_name` × phase level,
#' accepting either `term:phaseX` or `phaseX:term`.
find_fp_phase_interaction_name <- function(beta_names, term_name, phase) {
  candidates <- c(
    paste0(term_name, ":phase", phase),
    paste0("phase", phase, ":", term_name)
  )
  hit <- intersect(candidates, beta_names)
  if (length(hit) == 0L) {
    stop(
      "Expected interaction term not found for ", term_name, " × ", phase,
      ". Looked for: ", paste(candidates, collapse = " / "),
      ". Available: ", paste(beta_names, collapse = ", ")
    )
  }
  hit[[1]]
}

#' Per-phase slopes for a named FP component (FP_between or FP_within)
#' from a glmmTMB fit, with delta-method SE / 95% CI.
extract_wb_phase_slopes <- function(fit, term_name, model_id, hypothesis_group) {
  b <- glmmTMB::fixef(fit)$cond
  V <- stats::vcov(fit)$cond
  bounds <- data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    year_start = c(1985L, 1989L, 2001L, 2008L),
    year_end = c(1988L, 2000L, 2007L, 2015L),
    stringsAsFactors = FALSE
  )
  ref_phase <- bounds$phase[1]
  if (!term_name %in% names(b)) {
    stop("Term not found in fixed effects: ", term_name)
  }
  rows <- lapply(seq_len(nrow(bounds)), function(i) {
    ph <- bounds$phase[i]
    if (identical(ph, ref_phase)) {
      est <- unname(b[[term_name]])
      se <- sqrt(V[[term_name, term_name]])
    } else {
      int_term <- find_fp_phase_interaction_name(names(b), term_name, ph)
      cvec <- stats::setNames(rep(0, length(b)), names(b))
      cvec[[term_name]] <- 1
      cvec[[int_term]] <- 1
      est <- as.numeric(sum(cvec * b))
      se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    }
    z <- est / se
    data.frame(
      model_id = model_id,
      hypothesis_group = hypothesis_group,
      component = term_name,
      phase = ph,
      year_start = bounds$year_start[i],
      year_end = bounds$year_end[i],
      fp_slope = est,
      fp_slope_se = se,
      fp_slope_lo = est - 1.96 * se,
      fp_slope_hi = est + 1.96 * se,
      statistic = z,
      p_value = 2 * stats::pnorm(-abs(z)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Same as extract_wb_phase_slopes but for a spaMM CAR fit.
extract_wb_phase_slopes_spamm <- function(fit, term_name, model_id, hypothesis_group) {
  b <- spaMM::fixef(fit)
  V <- as.matrix(stats::vcov(fit))
  bounds <- data.frame(
    phase = c("1985-1988", "1989-2000", "2001-2007", "2008-2015"),
    year_start = c(1985L, 1989L, 2001L, 2008L),
    year_end = c(1988L, 2000L, 2007L, 2015L),
    stringsAsFactors = FALSE
  )
  ref_phase <- bounds$phase[1]
  if (!term_name %in% names(b)) {
    stop("Term not found in spaMM fixed effects: ", term_name)
  }
  rows <- lapply(seq_len(nrow(bounds)), function(i) {
    ph <- bounds$phase[i]
    if (identical(ph, ref_phase)) {
      est <- unname(b[[term_name]])
      se <- sqrt(V[[term_name, term_name]])
    } else {
      int_term <- find_fp_phase_interaction_name(names(b), term_name, ph)
      cvec <- stats::setNames(rep(0, length(b)), names(b))
      cvec[[term_name]] <- 1
      cvec[[int_term]] <- 1
      est <- as.numeric(sum(cvec * b))
      se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    }
    z <- est / se
    data.frame(
      model_id = model_id,
      hypothesis_group = hypothesis_group,
      component = term_name,
      phase = ph,
      year_start = bounds$year_start[i],
      year_end = bounds$year_end[i],
      fp_slope = est,
      fp_slope_se = se,
      fp_slope_lo = est - 1.96 * se,
      fp_slope_hi = est + 1.96 * se,
      statistic = z,
      p_value = 2 * stats::pnorm(-abs(z)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Names of coefficients that jointly define all phase-specific slopes for
#' a component (main effect + its phase interactions; either name order).
wb_all_slope_terms <- function(beta_names, term_name) {
  ints <- beta_names[is_fp_phase_interaction(beta_names, term_name)]
  c(term_name, ints)
}

#' Names of interaction-only terms for a component (H3-style change-across-phases).
wb_interaction_terms <- function(beta_names, term_name) {
  beta_names[is_fp_phase_interaction(beta_names, term_name)]
}
