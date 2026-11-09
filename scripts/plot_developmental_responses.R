#!/usr/bin/env Rscript
# Development-anchored figures: how Bap1 cKO affects loci that change during normal (control) development.
# Reads completed Step 04 outputs for all marks; runs before Step 05 and does not need joint eligibility.
source("R/pipeline_functions.R")
source("R/developmental_figures.R")
source("config/settings.R")
for (pkg in c("ggplot2", "DESeq2", "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("%s is required for developmental figures", pkg), call. = FALSE)
}
if (!isTRUE(capabilities("cairo"))) {
  stop("Cairo support is required for developmental PNG and SVG output", call. = FALSE)
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

out <- "results/developmental_responses"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name, directory = out) {
  utils::write.table(x, file.path(directory, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
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

multi_batch <- "batch" %in% names(samples) &&
  any(vapply(samples_by_mark, function(x) length(unique(x$batch)) > 1L, logical(1)))

render_view <- function(all_peaks = FALSE) {
  view_out <- if (all_peaks) file.path(out, "all_peaks") else out
  cohort_levels <- if (all_peaks) "all" else DEVELOPMENTAL_COHORTS
  responses <- developmental_responses(calls, cohorts, all_peaks = all_peaks)
  counts <- developmental_response_counts(responses, cohort_levels = cohort_levels)
  manifest <- heatmap_row_manifest(responses, cohort_levels = cohort_levels)
  scatter_summary <- developmental_scatter_summary(responses, cohort_levels = cohort_levels)
  dir.create(view_out, recursive = TRUE, showWarnings = FALSE)
  write_tsv(cohorts, "developmental_cohorts.tsv", view_out)
  write_tsv(responses, "developmental_locus_responses.tsv", view_out)
  write_tsv(counts, "response_frequencies.tsv", view_out)
  write_tsv(manifest, "heatmap_rows.tsv", view_out)
  write_tsv(scatter_summary, "scatter_summary.tsv", view_out)

  # Heatmaps are exported only as compiled multipage PDFs; remove legacy per-page PNG/SVG files from this view.
  page_pattern <- sprintf("^heatmap_(%s)_(%s)_(%s)_page_[0-9]{3,}\\.(png|svg)$",
                          paste(MARKS, collapse = "|"), paste(cohort_levels, collapse = "|"),
                          paste(MARKS, collapse = "|"))
  stale <- list.files(view_out, pattern = page_pattern, full.names = TRUE)
  if (length(stale)) unlink(stale)

  for (anchor in MARKS) {
    p <- plot_developmental_frequencies(
      counts, anchor, developmental_caption(NORMALIZATION, "frequency", all_peaks = all_peaks),
      cohort_levels = cohort_levels)
    stem <- file.path(view_out, sprintf("response_frequency_%s", anchor))
    render_pdf(paste0(stem, ".pdf"), list(p), 14, 8)
    render_png(paste0(stem, ".png"), p, 2100, 1200)
    render_developmental_svg(paste0(stem, ".svg"), p, 14, 8)
    p <- plot_developmental_scatter(
      responses, scatter_summary, anchor, developmental_caption(NORMALIZATION, "scatter", all_peaks = all_peaks),
      cohort_levels = cohort_levels)
    stem <- file.path(view_out, sprintf("developmental_vs_late_cKO_%s", anchor))
    render_pdf(paste0(stem, ".pdf"), list(p), 14, 8)
    render_png(paste0(stem, ".png"), p, 2100, 1200)
    render_developmental_svg(paste0(stem, ".svg"), p, 14, 8)
    for (cohort in cohort_levels) {
      rows <- manifest[manifest$anchor_mark == anchor & manifest$cohort == cohort, , drop = FALSE]
      pages <- if (nrow(rows)) sort(unique(rows$page)) else 1L
      for (response in MARKS) {
        stem <- sprintf("heatmap_%s_%s_%s", anchor, cohort, response)
        plots <- lapply(pages, function(page) {
          title <- if (all_peaks) {
            sprintf("%s signal at %s all developmental-comparison peaks (page %d of %d)",
                    response, anchor, page, length(pages))
          } else {
            sprintf("%s signal at %s developmental-%s loci (page %d of %d)",
                    response, anchor, cohort, page, length(pages))
          }
          plot_developmental_heatmap(
            signals[[response]]$row_z, rows[rows$page == page, , drop = FALSE], samples_by_mark[[response]],
            title, developmental_caption(NORMALIZATION, "heatmap", multi_batch, all_peaks = all_peaks))
        })
        render_pdf(file.path(view_out, paste0(stem, ".pdf")), plots, 14, 10)
      }
    }
  }
}
render_view()
if (DEVELOPMENTAL_ALL_PEAKS) render_view(all_peaks = TRUE)
cohort_sizes <- with(cohorts[cohorts$cohort %in% DEVELOPMENTAL_COHORTS, ], table(anchor_mark, cohort))
print(cohort_sizes)
cat("Wrote development-anchored response tables and figures to", out, "\n")
