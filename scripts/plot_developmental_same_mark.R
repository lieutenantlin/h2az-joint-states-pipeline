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
all_response_path <- file.path(input, "all_peaks", "developmental_locus_responses.tsv")
if (DEVELOPMENTAL_ALL_PEAKS && !file.exists(all_response_path)) {
  stopf("Missing all-peak developmental responses; rerun scripts/plot_developmental_responses.R with H2AZ_DEVELOPMENTAL_ALL_PEAKS=1")
}
if (!isTRUE(capabilities("cairo"))) stopf("Cairo support is required for developmental PNG and SVG output")
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
cohort_required <- c("anchor_mark", "locus_id", "cohort", "developmental_lfc", "developmental_padj")
if (DEVELOPMENTAL_ALL_PEAKS) cohort_required <- c(cohort_required, "developmental_state", "developmental_reason")
cohorts <- read_table(cohort_path, cohort_required)
if (!all(cohorts$anchor_mark %in% MARKS) || !all(cohorts$cohort %in% unname(COHORT_BY_STATE)) ||
    anyNA(cohorts$locus_id) || anyDuplicated(cohorts[c("anchor_mark", "locus_id")])) {
  stopf("%s contains invalid or duplicate developmental cohort keys", cohort_path)
}
for (column in c("developmental_lfc", "developmental_padj")) cohorts[[column]] <- as_numeric(cohorts[[column]], paste(cohort_path, column))
required <- c("anchor_mark", "cohort", "locus_id", "response_mark", "developmental_lfc", "developmental_padj",
              unlist(lapply(names(RESPONSE_PREFIXES), function(p) paste0(p, c("_lfc", "_padj", "_state", "_reason")))))
read_responses <- function(path, all_peaks = FALSE) {
  columns <- required
  if (all_peaks) columns <- c(columns, "developmental_state", "developmental_reason")
  x <- read_table(path, columns)
  x <- x[x$anchor_mark == x$response_mark, , drop = FALSE]
  levels <- if (all_peaks) "all" else DEVELOPMENTAL_COHORTS
  if (anyNA(x[c("anchor_mark", "response_mark", "cohort", "locus_id")]) ||
      !all(x$anchor_mark %in% MARKS) || !all(x$cohort %in% levels) ||
      anyDuplicated(x[c("anchor_mark", "locus_id")])) stopf("%s contains invalid or duplicate same-mark response keys", path)
  for (column in c("developmental_lfc", "developmental_padj", unlist(lapply(names(RESPONSE_PREFIXES), function(p) paste0(p, c("_lfc", "_padj")))))) {
    x[[column]] <- as_numeric(x[[column]], paste(path, column))
  }
  for (prefix in names(RESPONSE_PREFIXES)) {
    if (!all(x[[paste0(prefix, "_state")]] %in% RESPONSE_STATE_LEVELS)) stopf("%s has invalid %s states", path, prefix)
  }
  if (all_peaks && !all(x$developmental_state %in% RESPONSE_STATE_LEVELS)) stopf("%s has invalid developmental states", path)
  x
}
views <- list(default = list(path = response_path, responses = read_responses(response_path), levels = DEVELOPMENTAL_COHORTS))
if (DEVELOPMENTAL_ALL_PEAKS) {
  if (!all(cohorts$developmental_state %in% RESPONSE_STATE_LEVELS) ||
      !identical(unname(COHORT_BY_STATE[cohorts$developmental_state]), cohorts$cohort)) {
    stopf("%s contains invalid developmental states or cohort assignments", cohort_path)
  }
  views$all <- list(path = all_response_path, responses = read_responses(all_response_path, TRUE), levels = "all")
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
  for (view in views) {
    all_peaks <- identical(view$levels, "all")
    selected <- if (all_peaks) expected else expected[expected$cohort %in% DEVELOPMENTAL_COHORTS, , drop = FALSE]
    observed <- view$responses[view$responses$anchor_mark == mark, , drop = FALSE]
    ix <- match(selected$locus_id, observed$locus_id)
    expected_cohort <- if (all_peaks) rep("all", nrow(selected)) else selected$cohort
    mismatch <- nrow(observed) != nrow(selected) || anyNA(ix) ||
      !identical(expected_cohort, observed$cohort[ix]) ||
      !isTRUE(all.equal(selected$developmental_lfc, observed$developmental_lfc[ix], check.attributes = FALSE)) ||
      !isTRUE(all.equal(selected$developmental_padj, observed$developmental_padj[ix], check.attributes = FALSE))
    if (all_peaks && !mismatch) {
      mismatch <- !identical(selected$developmental_state, observed$developmental_state[ix]) ||
        !identical(selected$developmental_reason, observed$developmental_reason[ix])
    }
    if (mismatch) stopf("%s same-mark membership/effects disagree with %s for %s", view$path, cohort_path, mark)
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
render_view <- function(view, out) {
  responses <- view$responses
  levels <- view$levels
  all_peaks <- identical(levels, "all")
  counts <- developmental_response_counts(responses, cohort_levels = levels)
  counts <- counts[counts$anchor_mark == counts$response_mark, , drop = FALSE]
  manifest <- heatmap_row_manifest(responses, cohort_levels = levels)
  summary <- developmental_scatter_summary(responses, cohort_levels = levels)
  summary <- summary[summary$anchor_mark == summary$response_mark, , drop = FALSE]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  write_tsv(cohorts, "developmental_cohorts.tsv")
  write_tsv(responses, "developmental_locus_responses.tsv")
  write_tsv(counts, "response_frequencies.tsv")
  write_tsv(manifest, "heatmap_rows.tsv")
  write_tsv(summary, "scatter_summary.tsv")
  # Heatmaps are exported only as compiled multipage PDFs; remove legacy per-page PNG/SVG files from this view.
  cohort_pattern <- if (all_peaks) "all" else "gain|loss"
  pattern <- sprintf("^heatmap_(%s)_(%s)_page_[0-9]{3,}\\.(png|svg)$", paste(MARKS, collapse = "|"), cohort_pattern)
  stale <- list.files(out, pattern = pattern, full.names = TRUE)
  if (length(stale)) unlink(stale)
  for (mark in MARKS) for (cohort in levels) {
    caption <- function(family) developmental_caption(NORMALIZATION, family,
      multi_batch = "batch" %in% names(samples) && length(unique(samples_by_mark[[mark]]$batch)) > 1L,
      same_mark = TRUE, all_peaks = all_peaks)
    label <- if (all_peaks) "all developmental-comparison peaks" else sprintf("developmental-%s loci", cohort)
    p <- plot_developmental_frequencies(counts, mark, caption("frequency"), response_marks = mark, cohort_levels = cohort) +
      ggplot2::ggtitle(sprintf("%s %s: same-mark cKO responses", mark, label))
    stem <- sprintf("response_frequency_%s_%s", mark, cohort)
    render_pdf(file.path(out, paste0(stem, ".pdf")), list(p), 10, 8)
    render_png(file.path(out, paste0(stem, ".png")), p, 1500, 1200)
    render_developmental_svg(file.path(out, paste0(stem, ".svg")), p, 10, 8)
    p <- plot_developmental_scatter(responses, summary, mark, caption("scatter"), response_marks = mark, cohort_levels = cohort) +
      ggplot2::labs(title = sprintf("%s %s: developmental versus late cKO effect", mark, label),
                    y = sprintf("%s late cKO vs control log2 fold change", mark)) +
      ggplot2::theme(legend.box = "vertical")
    stem <- sprintf("developmental_vs_late_cKO_%s_%s", mark, cohort)
    render_pdf(file.path(out, paste0(stem, ".pdf")), list(p), 10, 8)
    render_png(file.path(out, paste0(stem, ".png")), p, 1500, 1200)
    render_developmental_svg(file.path(out, paste0(stem, ".svg")), p, 10, 8)
    rows <- manifest[manifest$anchor_mark == mark & manifest$cohort == cohort, , drop = FALSE]
    pages <- if (nrow(rows)) sort(unique(rows$page)) else 1L
    plots <- lapply(pages, function(page) plot_developmental_heatmap(
      signals[[mark]], rows[rows$page == page, , drop = FALSE], samples_by_mark[[mark]],
      sprintf("%s signal at its %s (page %d of %d)", mark, label, page, length(pages)), caption("heatmap")))
    stem <- sprintf("heatmap_%s_%s", mark, cohort)
    render_pdf(file.path(out, paste0(stem, ".pdf")), plots, 14, 10)
  }
}
out <- "results/developmental_same_mark"
render_view(views$default, out)
if (DEVELOPMENTAL_ALL_PEAKS) render_view(views$all, file.path(out, "all_peaks"))
cat("Wrote separate same-mark developmental cohort figures to", out, "\n")
