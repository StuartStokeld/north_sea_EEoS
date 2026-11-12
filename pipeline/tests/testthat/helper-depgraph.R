# Static script -> file dependency graph (audit layer 4.5, provenance / staleness).
#
# Parses R scripts (no execution) and records which project files each script reads and
# writes. Paths are resolved from string literals combined through file.path() / paste0(),
# simple `var <- path` assignments, the h2_output_paths() list and a few path helpers.
# Unresolvable paths are recorded as NA so that gaps are visible rather than silent.

DG_READ_FUNS <- c("readRDS", "read_csv", "read.csv", "read_rds", "fread", "load", "read.table",
                  "read_tsv", "readLines", "st_read", "read_sf", "read.dbf", "read_excel")
DG_WRITE_FUNS <- c("saveRDS", "write_csv", "write.csv", "write_rds", "save", "ggsave", "png",
                   "pdf", "writeLines", "sink", "fwrite", "write.table", "write_tsv", "st_write")

# helper functions that read fixed files (pipeline/R)
DG_HELPER_READS <- list(
  load_aligned_hl = c("outputs/datras_hl_raw.rds", "outputs/datras_hh_raw.rds"),
  ensure_haul_E_raw = "outputs/haul_state_variables.rds",
  augment_haul_harte_baseline = "outputs/haul_state_variables.rds"
)

dg_h2_paths <- c(
  couce_year = "outputs/h2_couce_year_effort.rds",
  couce_rectangle = "outputs/h2_couce_rectangle_effort.rds",
  couce_diagnostics = "outputs/h2_couce_import_diagnostics.csv",
  panel = "outputs/h2_rectangle_panel.rds", panel_csv = "outputs/h2_rectangle_panel.csv"
)

dg_fun_name <- function(e) {
  f <- e[[1]]
  if (is.name(f)) return(as.character(f))
  if (is.call(f) && as.character(f[[1]]) %in% c("::", ":::")) return(as.character(f[[3]]))
  NA_character_
}

dg_resolve <- function(e, vars) {
  if (identical(e, quote(expr = ))) return(NA_character_)
  if (is.character(e) && length(e) == 1L) return(e)
  if (is.name(e)) {
    nm <- as.character(e)
    if (nm %in% c("project_root", "root", "repo_root", "workspace_root", "proj_root")) return("")
    if (!is.null(vars[[nm]])) return(vars[[nm]])
    return(NA_character_)
  }
  if (is.call(e)) {
    fn <- dg_fun_name(e)
    if (identical(fn, "$") && is.name(e[[2]]) && as.character(e[[2]]) == "paths") {
      p <- dg_h2_paths[as.character(e[[3]])]
      return(if (is.na(p)) NA_character_ else unname(p))
    }
    if (fn %in% c("file.path", "paste0", "paste")) {
      args <- as.list(e)[-1]
      if (!is.null(names(args))) args <- args[!names(args) %in% c("sep", "collapse")]
      parts <- vapply(args, dg_resolve, "", vars = vars)
      if (anyNA(parts)) return(NA_character_)
      parts <- parts[nzchar(parts)]
      sep <- if (fn == "file.path") "/" else ""
      return(gsub("/+", "/", paste(parts, collapse = sep)))
    }
    if (fn %in% c("normalizePath", "path.expand")) return(dg_resolve(e[[2]], vars))
    if (fn == "test_path") return(NA_character_)
  }
  NA_character_
}

dg_target_arg <- function(e, fn) {
  a <- as.list(e)[-1]
  nms <- names(a)
  pick <- function(n) if (!is.null(nms) && n %in% nms) a[[n]] else NULL
  for (n in c("file", "filename", "path", "con", "dsn")) {
    if (!is.null(pick(n))) return(pick(n))
  }
  pos <- if (is.null(nms)) a else a[!nzchar(nms)]
  if (fn %in% c("saveRDS", "write_csv", "write.csv", "write_rds", "fwrite", "write.table",
                "write_tsv", "st_write", "writeLines")) {
    return(if (length(pos) >= 2L) pos[[2]] else NULL)
  }
  if (fn == "save") return(NULL) # save(...) without file=
  if (length(pos) >= 1L) pos[[1]] else NULL
}

#' Reads/writes of one script: data.frame(script, direction, path, call).
dg_parse_script <- function(script, rel_to) {
  exprs <- tryCatch(as.list(parse(script, keep.source = FALSE)), error = function(e) list())
  vars <- list()
  out <- list()
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    fn <- dg_fun_name(e)
    if (!is.na(fn) && fn %in% c("<-", "=") && length(e) == 3L && is.name(e[[2]])) {
      p <- dg_resolve(e[[3]], vars)
      if (!is.na(p)) vars[[as.character(e[[2]])]] <<- p
    }
    if (!is.na(fn) && (fn %in% DG_READ_FUNS || fn %in% DG_WRITE_FUNS)) {
      tgt <- dg_target_arg(e, fn)
      p <- if (is.null(tgt)) NA_character_ else dg_resolve(tgt, vars)
      out[[length(out) + 1L]] <<- data.frame(
        direction = if (fn %in% DG_READ_FUNS) "read" else "write",
        path = p, call = fn, stringsAsFactors = FALSE)
    }
    if (!is.na(fn) && fn %in% names(DG_HELPER_READS)) {
      for (p in DG_HELPER_READS[[fn]]) {
        out[[length(out) + 1L]] <<- data.frame(direction = "read", path = p, call = fn,
                                               stringsAsFactors = FALSE)
      }
    }
    args <- as.list(e)[-1]
    for (i in seq_along(args)) {
      if (!identical(args[[i]], quote(expr = )) && is.call(args[[i]])) walk(args[[i]])
    }
    invisible()
  }
  for (e in exprs) walk(e)
  res <- do.call(rbind, c(out, list(data.frame(direction = character(), path = character(),
                                                 call = character()))))
  if (nrow(res)) {
    res$path <- sub("^/", "", res$path)
    res$path <- sub("^\\./", "", res$path)
  }
  res$script <- rep(sub(paste0("^", rel_to, "/?"), "", script), nrow(res))
  res[, c("script", "direction", "path", "call")]
}

#' Live scripts = those named in the CLAUDE.md run order (H1 and H2/H3 blocks).
dg_live_scripts <- function(root) {
  txt <- readLines(file.path(root, "CLAUDE.md"))
  hits <- regmatches(txt, regexpr("pipeline/[A-Za-z0-9_]+\\.R", txt))
  unique(hits)
}

#' Full graph over pipeline/*.R (and exploratory/**/*.R as potential producers).
dg_build <- function(root) {
  scripts <- c(list.files(file.path(root, "pipeline"), pattern = "\\.R$", full.names = TRUE),
               list.files(file.path(root, "exploratory"), pattern = "\\.R$", full.names = TRUE,
                          recursive = TRUE))
  g <- do.call(rbind, lapply(scripts, dg_parse_script, rel_to = root))
  live <- dg_live_scripts(root)
  g$script_class <- ifelse(g$script %in% live, "live",
                           ifelse(grepl("^exploratory/", g$script), "exploratory", "pipeline_other"))
  unique(g)
}

#' For each file read by a live script: who produces it.
dg_live_inputs <- function(g, raw_inputs) {
  reads <- unique(g[g$script_class == "live" & g$direction == "read" & !is.na(g$path), c("script", "path")])
  writes <- g[g$direction == "write" & !is.na(g$path), ]
  reads$producers_live <- vapply(reads$path, function(p) paste(sort(unique(writes$script[writes$path == p & writes$script_class == "live"])), collapse = "; "), "")
  reads$producers_other <- vapply(reads$path, function(p) paste(sort(unique(writes$script[writes$path == p & writes$script_class != "live"])), collapse = "; "), "")
  reads$raw_input <- vapply(reads$path, function(p) any(startsWith(p, raw_inputs)), TRUE)
  reads$status <- with(reads, ifelse(nzchar(producers_live), "ok: live producer",
                              ifelse(raw_input, "ok: documented raw input",
                              ifelse(grepl("exploratory/", producers_other), "FLAG: produced only by exploratory script",
                              ifelse(nzchar(producers_other), "FLAG: produced only by a non-run-order script",
                                     "FLAG: no producer found")))))
  reads
}

DG_RAW_INPUTS <- c(
  "FishGlob_data/", "equation_of_state/", "data/external/", "gis/",
  "outputs/datras_hl_raw.rds", "outputs/datras_hh_raw.rds", # committed raw inputs (import scripts optional)
  "outputs/fishbase_lw_lookup_v2.csv" # built outside the pipeline (supplementary notebook) - see F-B13
)
