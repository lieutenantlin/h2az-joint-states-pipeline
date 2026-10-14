# Project-wide settings. The runner changes to the project root before sourcing this file.
PROJECT_ROOT <- normalizePath(getwd(), mustWork = TRUE)
NORMALIZATION <- "background" # Change to "spikein" only with validated calibration BAMs in every row.
THREADS <- suppressWarnings(as.integer(Sys.getenv("THREADS", "")))
if (length(THREADS) != 1L || is.na(THREADS) || THREADS < 1L) {
  THREADS <- parallel::detectCores(logical = FALSE)
}
# Base R may return NA when it cannot determine the physical core count.
if (length(THREADS) != 1L || is.na(THREADS) || THREADS < 1L) THREADS <- 1L
THREADS <- as.integer(THREADS)
BLACKLIST <- file.path(PROJECT_ROOT, "resources", "mm10-blacklist.v2.bed.gz")

if (!NORMALIZATION %in% c("background", "spikein")) {
  stop("NORMALIZATION must be 'background' or 'spikein'.", call. = FALSE)
}
if (!file.exists(BLACKLIST)) stop("Missing mm10 blacklist: ", BLACKLIST, call. = FALSE)
