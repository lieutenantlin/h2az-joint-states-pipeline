#!/usr/bin/env Rscript
# Development-anchored figures: how Bap1 cKO affects loci that change during normal (control) development.
# Reads completed Step 04 outputs for all marks; runs before Step 05 and does not need joint eligibility.
source("R/pipeline_functions.R")
source("R/developmental_figures.R")
source("config/settings.R")
for (pkg in c("ggplot2", "DESeq2", "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("%s is required for developmental figures", pkg), call. = FALSE)
}

master_path <- "results/03_master/master_regions.rds"
inventory_path <- "results/01_inventory/samples_validated.csv"
mark_dir <- function(mark) file.path("results/04_diffbind", mark)
result_paths <- stats::setNames(file.path(vapply(MARKS, mark_dir, ""), paste0(MARKS, "_all_loci_results.tsv")), MARKS)
dds_paths <- stats::setNames(file.path(vapply(MARKS, mark_dir, ""), paste0(MARKS, "_dds.rds")), MARKS)
map_paths <- stats::setNames(file.path(vapply(MARKS, mark_dir, ""), "dba_row_to_locus_id.tsv"), MARKS)
norm_paths <- stats::setNames(file.path(vapply(MARKS, mark_dir, ""), "normalization.tsv"), MARKS)
inputs <- c(master_path, inventory_path, result_paths, dds_paths, map_paths, norm_paths)
missing <- inputs[!file.exists(inputs)]
if (length(missing)) stop(paste(c("Missing inputs; run Steps 01-04 first:", missing), collapse = "\n  "), call. = FALSE)

master <- readRDS(master_path)
master_ids <- master_locus_ids(master)
samples <- utils::read.csv(inventory_path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
samples_by_mark <- lapply(stats::setNames(MARKS, MARKS), function(mark) samples[samples$mark == mark, , drop = FALSE])

# The saved Step 04 normalization must match the current setting so captions describe the analysis actually run.
for (mark in MARKS) {
  norm <- utils::read.delim(norm_paths[[mark]], stringsAsFactors = FALSE, colClasses = "character")
  if (!identical(norm$sample_id, samples_by_mark[[mark]]$sample_id)) {
    stop(sprintf("%s sample IDs do not match the validated inventory; rerun Step 04", norm_paths[[mark]]), call. = FALSE)
  }
  modes <- unique(norm$normalization_mode)
  if (length(modes) != 1L || !identical(modes, NORMALIZATION)) {
    stop(sprintf("%s records normalization '%s' but config/settings.R selects '%s'; rerun Step 04 or restore the matching setting",
                 norm_paths[[mark]], paste(modes, collapse = ","), NORMALIZATION), call. = FALSE)
  }
}

calls <- read_developmental_calls(result_paths, master_ids)
signals <- lapply(stats::setNames(MARKS, MARKS), function(mark) {
  row_map <- utils::read.delim(map_paths[[mark]], stringsAsFactors = FALSE, colClasses = "character")
  tryCatch(developmental_signal(readRDS(dds_paths[[mark]]), row_map, samples_by_mark[[mark]], master),
           error = function(e) stop(sprintf("%s: %s", mark, conditionMessage(e)), call. = FALSE))
})
cohorts <- developmental_cohorts(calls, master_ids)
responses <- developmental_responses(calls, cohorts)
counts <- developmental_response_counts(responses)
manifest <- heatmap_row_manifest(responses)
scatter_summary <- developmental_scatter_summary(responses)

out <- "results/developmental_responses"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(cohorts, "developmental_cohorts.tsv")
write_tsv(responses, "developmental_locus_responses.tsv")
write_tsv(counts, "response_frequencies.tsv")
write_tsv(manifest, "heatmap_rows.tsv")
write_tsv(scatter_summary, "scatter_summary.tsv")
for (mark in MARKS) {
  for (kind in c("log2_counts", "row_z")) {
    m <- signals[[mark]][[kind]]
    x <- data.frame(locus_id = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
    name <- if (kind == "log2_counts") "log2_normalized_counts" else "row_z"
    write_tsv(x, sprintf("%s_%s.tsv", mark, name))
  }
}

render_pdf <- function(path, plots, width, height) {
  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off())
  for (p in plots) print(p)
}
render_png <- function(path, plot, width, height) {
  grDevices::png(path, width = width, height = height, res = 150, type = "cairo")
  on.exit(grDevices::dev.off())
  print(plot)
}

# Remove only this stage's paged heatmap PNGs so a smaller cohort cannot leave obsolete pages.
page_pattern <- sprintf("^heatmap_(%s)_(gain|loss)_(%s)_page_[0-9]{3,}\\.png$", paste(MARKS, collapse = "|"), paste(MARKS, collapse = "|"))
stale <- list.files(out, pattern = page_pattern, full.names = TRUE)
if (length(stale)) unlink(stale)

multi_batch <- "batch" %in% names(samples) &&
  any(vapply(samples_by_mark, function(x) length(unique(x$batch)) > 1L, logical(1)))
for (anchor in MARKS) {
  p <- plot_developmental_frequencies(counts, anchor, developmental_caption(NORMALIZATION, "frequency"))
  render_pdf(file.path(out, sprintf("response_frequency_%s.pdf", anchor)), list(p), 14, 8)
  render_png(file.path(out, sprintf("response_frequency_%s.png", anchor)), p, 2100, 1200)
  p <- plot_developmental_scatter(responses, scatter_summary, anchor, developmental_caption(NORMALIZATION, "scatter"))
  render_pdf(file.path(out, sprintf("developmental_vs_late_cKO_%s.pdf", anchor)), list(p), 14, 8)
  render_png(file.path(out, sprintf("developmental_vs_late_cKO_%s.png", anchor)), p, 2100, 1200)
  for (cohort in DEVELOPMENTAL_COHORTS) {
    rows <- manifest[manifest$anchor_mark == anchor & manifest$cohort == cohort, , drop = FALSE]
    pages <- if (nrow(rows)) sort(unique(rows$page)) else 1L
    for (response in MARKS) {
      stem <- sprintf("heatmap_%s_%s_%s", anchor, cohort, response)
      plots <- lapply(pages, function(page) {
        plot_developmental_heatmap(
          signals[[response]]$row_z, rows[rows$page == page, , drop = FALSE], samples_by_mark[[response]],
          sprintf("%s signal at %s developmental-%s loci (page %d of %d)", response, anchor, cohort, page, length(pages)),
          developmental_caption(NORMALIZATION, "heatmap", multi_batch))
      })
      render_pdf(file.path(out, paste0(stem, ".pdf")), plots, 14, 10)
      for (i in seq_along(pages)) render_png(file.path(out, sprintf("%s_page_%03d.png", stem, pages[[i]])), plots[[i]], 2100, 1500)
    }
  }
}
cohort_sizes <- with(cohorts[cohorts$cohort %in% DEVELOPMENTAL_COHORTS, ], table(anchor_mark, cohort))
print(cohort_sizes)
cat("Wrote development-anchored response tables and figures to", out, "\n")
