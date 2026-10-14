#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
# Step 07: select a small, ranked set of loci for manual IGV review.
wide_path <- "results/05_joint/joint_locus_states.tsv"
if (!file.exists(wide_path)) stop("Run Step 05 before Step 07", call. = FALSE)
wide <- utils::read.delim(wide_path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA")
result_tables <- lapply(MARKS, function(mark) {
  path <- file.path("results/04_diffbind", mark, paste0(mark, "_all_loci_results.tsv"))
  if (!file.exists(path)) stop(sprintf("Missing differential-results table: %s", path), call. = FALSE)
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
  if (!all(c("locus_id", "contrast", "mark", "baseMean", "padj") %in% names(x))) stop(sprintf("Malformed table: %s", path), call. = FALSE)
  x
})
calls <- do.call(rbind, result_tables)
# Review representative concordant, discordant, and no-detected-change joint-state patterns.
review_states <- c("+++", "---", "+0-", "000")
selected <- list()
for (contrast in names(CONTRASTS)) for (state in review_states) {
  loci <- wide$locus_id[wide$contrast == contrast & wide$joint == state]
  if (!length(loci)) next
  metrics <- lapply(loci, function(locus) {
    x <- calls[calls$locus_id == locus & calls$contrast == contrast, , drop = FALSE]
    if (nrow(x) != length(MARKS) || !setequal(x$mark, MARKS)) stop("Review input is incomplete", call. = FALSE)
    x <- x[match(MARKS, x$mark), , drop = FALSE]
    # Changed states rank by their weakest changed-mark FDR; 000 ranks by signal strength.
    if (identical(state, "000")) return(min(as.numeric(x$baseMean), na.rm = TRUE))
    changed <- strsplit(state, "", fixed = TRUE)[[1L]] != "0"
    max(as.numeric(x$padj)[changed], na.rm = TRUE)
  })
  metric <- unlist(metrics, use.names = FALSE)
  if (identical(state, "000")) order_index <- order(-metric, loci) else order_index <- order(metric, loci)
  keep <- order_index[seq_len(min(3L, length(order_index)))]
  selected[[length(selected) + 1L]] <- data.frame(contrast = contrast, state = state, locus_id = loci[keep],
                                                   rank = seq_along(keep), ranking_metric = metric[keep],
                                                   stringsAsFactors = FALSE)
}
out <- "results/07_review"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
review <- if (length(selected)) do.call(rbind, selected) else data.frame(contrast = character(), state = character(), locus_id = character(), rank = integer(), ranking_metric = numeric())
utils::write.table(review, file.path(out, "review_loci.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat(sprintf("Wrote %d loci for manual IGV review.\n", nrow(review)))
