#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
source("config/settings.R")

inventory <- "results/01_inventory/samples_validated.csv"
if (!file.exists(inventory)) stop("Run scripts/01_validate_inventory.R first", call. = FALSE)
samples <- utils::read.csv(inventory, stringsAsFactors = FALSE, check.names = FALSE)
if (!nrow(samples)) stop("Validated inventory is empty", call. = FALSE)
if (!file.exists(BLACKLIST)) stop(sprintf("Blacklist is missing: %s", BLACKLIST), call. = FALSE)
blacklist <- read_bed_gr(BLACKLIST)
peak_path <- function(row) {
  if (nzchar(row[["peaks"]])) return(row[["peaks"]])
  if (identical(row[["layout"]], "PE")) return(file.path("results/02_peaks", paste0(row[["sample_id"]], ".stringent.bed")))
  file.path("results/02_peaks", paste0(row[["sample_id"]], "_peaks.narrowPeak"))
}
strata <- list()
support <- data.frame(mark = character(), group = character(), replicate_peaks = integer(), supported_regions = integer(), stringsAsFactors = FALSE)
for (mark in MARKS) for (group in GROUPS) {
  x <- samples[samples$mark == mark & samples$group == group, , drop = FALSE]
  peaks <- lapply(seq_len(nrow(x)), function(i) {
    path <- peak_path(x[i, , drop = FALSE])
    if (!file.exists(path)) stop(sprintf("Peak file is missing for %s: %s", x$sample_id[[i]], path), call. = FALSE)
    read_bed_gr(path)
  })
  regions <- replicate_supported(peaks, min_reps = 2L)
  strata[[paste(mark, group, sep = "_")]] <- regions
  support <- rbind(support, data.frame(mark = mark, group = group, replicate_peaks = sum(vapply(peaks, length, integer(1))), supported_regions = length(regions), stringsAsFactors = FALSE))
  if (!length(regions)) warning(sprintf("No supported regions for %s x %s", mark, group), call. = FALSE)
}
master <- build_master(strata, blacklist)
if (!length(master)) stop("Master region set is empty after chromosome and blacklist filtering", call. = FALSE)
out <- "results/03_master"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(master, file.path(out, "master_regions.rds"))
master_ids <- master_locus_ids(master)
master_bed <- data.frame(as.character(GenomicRanges::seqnames(master)), GenomicRanges::start(master) - 1L,
                         GenomicRanges::end(master), master_ids, stringsAsFactors = FALSE)
utils::write.table(master_bed, file.path(out, "master_regions.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
utils::write.table(master_bed[, 1:3], file.path(out, "master_regions.counting.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
utils::write.table(support, file.path(out, "stratum_support_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
widths <- GenomicRanges::width(master)
width_summary <- data.frame(statistic = c("n", "min", "q1", "median", "mean", "q3", "max"),
                            width_bp = c(length(widths), min(widths), stats::quantile(widths, .25, names = FALSE), stats::median(widths), mean(widths), stats::quantile(widths, .75, names = FALSE), max(widths)))
utils::write.table(width_summary, file.path(out, "width_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
widest <- master_ids[order(widths, decreasing = TRUE)[seq_len(min(10L, length(master)))]]
utils::write.table(master_bed[match(widest, master_ids), , drop = FALSE], file.path(out, "widest10.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d fixed master regions.\n", length(master)))
