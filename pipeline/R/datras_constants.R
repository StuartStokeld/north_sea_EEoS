# Shared DATRAS pipeline constants (north_sea_eeos)

#' Non-fish AphiaIDs excluded before S, N, E computation (cephalopods, crustaceans, bivalves).
NON_FISH_APHIAIDS <- c(
  153131, 138138, 140270, 138139, 140271, 11734, 11760,
  141452, 140625, 11716, 140621, 141454, 11718, 141449,
  140624, 140600, 138481, 138482, 138483, 141448, 141450,
  141443, 141444, 138477, 11725, 494117, 138802, 153083,
  11707, 11709, 138802, 1437787, 1437789, 140605,
  107254, 107276, 107205, 107253, 107347, 107350, 107398, 106789,
  140712, 140687, 658802, 167039, 167569
)

ANALYSIS_YEAR_MIN <- 1985L
ANALYSIS_YEAR_MAX <- 2015L
ANALYSIS_QUARTER <- 1L

KNOWN_LNGT_CODES <- c("0", "1", "5", ".")

# --- Haul effort alignment (HH record) ------------------------------------------
# DATRAS HL counts are per haul for DataType R/P/S but per hour for DataType C
# (DATRAS FAQ; for C, SubFactor = 1). Hauls are aligned to HH before S, N, E are
# built: valid hauls only, standard 30-min tows (NS-IBTS manual, ICES 2020 SISP 10
# Rev. 11, s2.5.4 and Annex 2), C converted back to raw counts with HaulDur / 60.
# No rescaling to a common duration (S cannot be rescaled).
#
# Duration window: PRIMARY analysis = hauls of EXACTLY 30 min (HaulDur == 30;
# supervisor decision, 10 Nov 2026). SENSITIVITY = 26-34 min (30 +/- 4), run with
# EEOS_HAUL_DUR_MIN=26 EEOS_HAUL_DUR_MAX=34 (outputs in outputs/sensitivity_window_26_34/).

HAUL_DUR_MIN <- 30
HAUL_DUR_MAX <- 30
HH_DATATYPES_KEEP <- c("C", "R", "P", "S")
HH_HAULVAL_KEEP <- "V"

#' Minimum share of HL hauls that must match an HH record (code-mismatch guard).
HH_HL_MIN_MATCH_RATE <- 0.95

#' Haul-duration window in minutes, overridable for sensitivity runs via
#' environment variables EEOS_HAUL_DUR_MIN / EEOS_HAUL_DUR_MAX
#' (e.g. the 26-34 sensitivity: EEOS_HAUL_DUR_MIN=26 EEOS_HAUL_DUR_MAX=34 Rscript ...).
#' Defaults: HAUL_DUR_MIN / HAUL_DUR_MAX (exactly 30 min, primary).
haul_dur_window <- function() {
  read_env <- function(var, default) {
    val <- Sys.getenv(var, unset = "")
    if (!nzchar(val)) {
      return(default)
    }
    num <- suppressWarnings(as.numeric(val))
    if (length(num) != 1L || is.na(num) || num < 0) {
      stop(var, " must be a non-negative number, got: '", val, "'")
    }
    num
  }
  dur_min <- read_env("EEOS_HAUL_DUR_MIN", HAUL_DUR_MIN)
  dur_max <- read_env("EEOS_HAUL_DUR_MAX", HAUL_DUR_MAX)
  if (dur_min > dur_max) {
    stop("Haul-duration window is empty: min ", dur_min, " > max ", dur_max)
  }
  c(min = dur_min, max = dur_max)
}
