#!/usr/bin/env Rscript
# A separate presentation view of completed developmental results: one mark/cohort per figure.
# No differential refit or new cohort selection; existing cross-mark outputs remain untouched.
source("R/pipeline_functions.R")
source("R/developmental_figures.R")
source("config/settings.R")
if (!requireNamespace("ggplot2", quietly = TRUE)) stopf("ggplot2 is required for same-mark developmental figures")

input <- "results/developmental_responses"
inventory_path <- "results/01_inventory/samples_validated.csv"
cohort_path <- file.path(input, "developmental_cohorts.tsv")
response_path <- file.path(input, "developmental_locus_responses.tsv")
z_paths <- stats::setNames(file.path(input, paste0(MARKS, "_row_z.tsv")), MARKS)
norm_paths <- stats::setNames(file.path("results/04_diffbind", MARKS, "normalization.tsv"), MARKS)
missing <- c(inventory_path, cohort_path, response_path, z_paths, norm_paths)
missing <- missing[!file.exists(missing)]
if (length(missing)) {
  stopf("Missing inputs; first run scripts/plot_developmental_responses.R after Steps 01-04:\n  %s", paste(missing, collapse = "\n  "))
}
read_table <- function(path, required) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA")
  if (!all(required %in% names(x))) stopf("%s lacks required columns: %s", path, paste(setdiff(required, names(x)), collapse = ", "))
  x
}
as_numeric <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (any(is.na(value) & !is.na(x))) stopf("Non-numeric values in %s", label)
  value
}
samples <- utils::read.csv(inventory_path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
if (!all(c("sample_id", "mark", "group") %in% names(samples)) || anyDuplicated(samples$sample_id)) {
  stopf("%s must contain unique sample_id values plus mark and group", inventory_path)
}
cohorts <- read_table(cohort_path, c("anchor_mark", "locus_id", "cohort", "developmental_lfc", "developmental_padj"))
if (!all(cohorts$anchor_mark %in% MARKS) || !all(cohorts$cohort %in% unname(COHORT_BY_STATE)) ||
    anyNA(cohorts$locus_id) || anyDuplicated(cohorts[c("anchor_mark", "locus_id")])) {
  stopf("%s contains invalid or duplicate developmental cohort keys", cohort_path)
}
for (column in c("developmental_lfc", "developmental_padj")) cohorts[[column]] <- as_numeric(cohorts[[column]], paste(cohort_path, column))
required <- c("anchor_mark", "cohort", "locus_id", "response_mark", "developmental_lfc", "developmental_padj",
              unlist(lapply(names(RESPONSE_PREFIXES), function(p) paste0(p, c("_lfc", "_padj", "_state", "_reason")))))
responses <- read_table(response_path, required)
responses <- responses[responses$anchor_mark == responses$response_mark, , drop = FALSE]
if (anyNA(responses[c("anchor_mark", "response_mark", "cohort", "locus_id")]) ||
    !all(responses$anchor_mark %in% MARKS) || !all(responses$cohort %in% DEVELOPMENTAL_COHORTS) ||
    anyDuplicated(responses[c("anchor_mark", "locus_id")])) stopf("%s contains invalid or duplicate same-mark response keys", response_path)
for (column in c("developmental_lfc", "developmental_padj", unlist(lapply(names(RESPONSE_PREFIXES), function(p) paste0(p, c("_lfc", "_padj")))))) {
  responses[[column]] <- as_numeric(responses[[column]], paste(response_path, column))
}
for (prefix in names(RESPONSE_PREFIXES)) {
  if (!all(responses[[paste0(prefix, "_state")]] %in% RESPONSE_STATE_LEVELS)) stopf("%s has invalid %s states", response_path, prefix)
}

# Preflight all marks before creating outputs. Matrix columns and normalization must match the inventory.
samples_by_mark <- lapply(stats::setNames(MARKS, MARKS), function(mark) samples[samples$mark == mark, , drop = FALSE])
signals <- lapply(stats::setNames(MARKS, MARKS), function(mark) {
  s <- samples_by_mark[[mark]]
  if (!nrow(s) || !setequal(s$group, GROUPS)) stopf("Validated inventory lacks required groups for %s", mark)
  norm <- read_table(norm_paths[[mark]], c("sample_id", "normalization_mode"))
  if (!identical(norm$sample_id, s$sample_id)) stopf("%s sample IDs do not match the validated inventory", norm_paths[[mark]])
  if (!identical(unique(norm$normalization_mode), NORMALIZATION)) {
    stopf("%s normalization disagrees with config/settings.R; restore matching settings or rerun Step 04 and the developmental figure stage", norm_paths[[mark]])
  }
  expected <- cohorts[cohorts$anchor_mark == mark, , drop = FALSE]
  first_ids <- cohorts$locus_id[cohorts$anchor_mark == MARKS[[1L]]]
  if (!nrow(expected) || !setequal(expected$locus_id, first_ids)) stopf("%s does not cover the same loci for every mark", cohort_path)
  selected <- expected[expected$cohort %in% DEVELOPMENTAL_COHORTS, , drop = FALSE]
  observed <- responses[responses$anchor_mark == mark, , drop = FALSE]
  ix <- match(selected$locus_id, observed$locus_id)
  if (nrow(observed) != nrow(selected) || anyNA(ix) || !identical(selected$cohort, observed$cohort[ix]) ||
      !isTRUE(all.equal(selected$developmental_lfc, observed$developmental_lfc[ix], check.attributes = FALSE)) ||
      !isTRUE(all.equal(selected$developmental_padj, observed$developmental_padj[ix], check.attributes = FALSE))) {
    stopf("%s same-mark membership/effects disagree with %s for %s", response_path, cohort_path, mark)
  }
  z <- read_table(z_paths[[mark]], c("locus_id", s$sample_id))
  if (!identical(names(z), c("locus_id", s$sample_id)) || anyDuplicated(z$locus_id) ||
      !setequal(z$locus_id, expected$locus_id)) stopf("%s locus/sample identities disagree with developmental inputs", z_paths[[mark]])
  values <- lapply(s$sample_id, function(id) as_numeric(z[[id]], paste(z_paths[[mark]], id)))
  m <- do.call(cbind, values)
  dimnames(m) <- list(z$locus_id, s$sample_id)
  if (!all(is.finite(m))) stopf("%s contains nonfinite row-z signal", z_paths[[mark]])
  m
})
counts <- developmental_response_counts(responses)
counts <- counts[counts$anchor_mark == counts$response_mark, , drop = FALSE]
manifest <- heatmap_row_manifest(responses)
summary <- developmental_scatter_summary(responses)
summary <- summary[summary$anchor_mark == summary$response_mark, , drop = FALSE]

out <- "results/developmental_same_mark"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(cohorts, "developmental_cohorts.tsv")
write_tsv(responses, "developmental_locus_responses.tsv")
write_tsv(counts, "response_frequencies.tsv")
write_tsv(manifest, "heatmap_rows.tsv")
write_tsv(summary, "scatter_summary.tsv")
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
# Only remove obsolete paged PNGs owned by this stage, never source/cross-mark output files.
pattern <- sprintf("^heatmap_(%s)_(gain|loss)_page_[0-9]{3,}\\.png$", paste(MARKS, collapse = "|"))
stale <- list.files(out, pattern = pattern, full.names = TRUE)
if (length(stale)) unlink(stale)
for (mark in MARKS) for (cohort in DEVELOPMENTAL_COHORTS) {
  caption <- function(family) developmental_caption(NORMALIZATION, family,
    multi_batch = "batch" %in% names(samples) && length(unique(samples_by_mark[[mark]]$batch)) > 1L, same_mark = TRUE)
  p <- plot_developmental_frequencies(counts, mark, caption("frequency"), response_marks = mark, cohort_levels = cohort) +
    ggplot2::ggtitle(sprintf("%s developmental-%s loci: same-mark cKO responses", mark, cohort))
  stem <- sprintf("response_frequency_%s_%s", mark, cohort)
  render_pdf(file.path(out, paste0(stem, ".pdf")), list(p), 10, 8)
  render_png(file.path(out, paste0(stem, ".png")), p, 1500, 1200)
  p <- plot_developmental_scatter(responses, summary, mark, caption("scatter"), response_marks = mark, cohort_levels = cohort) +
    ggplot2::labs(title = sprintf("%s developmental-%s loci: developmental versus late cKO effect", mark, cohort),
                  y = sprintf("%s late cKO vs control log2 fold change", mark)) +
    ggplot2::theme(legend.box = "vertical")
  stem <- sprintf("developmental_vs_late_cKO_%s_%s", mark, cohort)
  render_pdf(file.path(out, paste0(stem, ".pdf")), list(p), 10, 8)
  render_png(file.path(out, paste0(stem, ".png")), p, 1500, 1200)
  rows <- manifest[manifest$anchor_mark == mark & manifest$cohort == cohort, , drop = FALSE]
  pages <- if (nrow(rows)) sort(unique(rows$page)) else 1L
  plots <- lapply(pages, function(page) plot_developmental_heatmap(
    signals[[mark]], rows[rows$page == page, , drop = FALSE], samples_by_mark[[mark]],
    sprintf("%s signal at its developmental-%s loci (page %d of %d)", mark, cohort, page, length(pages)), caption("heatmap")))
  stem <- sprintf("heatmap_%s_%s", mark, cohort)
  render_pdf(file.path(out, paste0(stem, ".pdf")), plots, 14, 10)
  for (i in seq_along(pages)) render_png(file.path(out, sprintf("%s_page_%03d.png", stem, pages[[i]])), plots[[i]], 2100, 1500)
}
cat("Wrote separate same-mark developmental cohort figures to", out, "\n")
