#!/usr/bin/env Rscript

# Step 04: count one mark on the fixed master loci and test all planned four-group contrasts.
main <- function(mark) {
  source("R/pipeline_functions.R")
  source("config/settings.R")
  if (!requireNamespace("DiffBind", quietly = TRUE) || !requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DiffBind and DESeq2 are required for Step 04", call. = FALSE)
  }
  if (!mark %in% MARKS) stop(sprintf("Unknown mark: %s", mark), call. = FALSE)
  inventory <- "results/01_inventory/samples_validated.csv"
  master_path <- "results/03_master/master_regions.rds"
  if (!file.exists(inventory) || !file.exists(master_path)) stop("Run Steps 01 and 03 before Step 04", call. = FALSE)
  samples <- utils::read.csv(inventory, stringsAsFactors = FALSE, check.names = FALSE)
  samples <- samples[samples$mark == mark, , drop = FALSE]
  if (!nrow(samples)) stop(sprintf("No validated samples for %s", mark), call. = FALSE)
  if (!identical(sort(unique(samples$group)), sort(GROUPS))) stop("All four required groups are required; reduced-comparison designs are not selected automatically", call. = FALSE)
  # The planned contrast set needs every EC/EK/LC/LK group; no reduced design is inferred.
  master <- readRDS(master_path)
  master_ids <- master_locus_ids(master)
  if (!length(master)) stop("master_regions.rds is empty", call. = FALSE)
  if (length(unique(samples$layout)) != 1L) stop("Layout must be identical within a mark", call. = FALSE)
  if (!identical(NORMALIZATION, "background") && !identical(NORMALIZATION, "spikein")) stop("NORMALIZATION must be background or spikein", call. = FALSE)

  # DiffBind raw peak files require a fourth score column. This disposable file only initializes
  # dba(); dba.count() receives the exact master GRanges below, while the authoritative on-disk
  # master_regions.counting.bed remains three-column BED as specified.
  tmp_peaks <- tempfile(pattern = paste0(mark, "_master_"), fileext = ".raw.bed")
  on.exit(unlink(tmp_peaks), add = TRUE)
  utils::write.table(data.frame(as.character(GenomicRanges::seqnames(master)), GenomicRanges::start(master) - 1L,
                                GenomicRanges::end(master), rep.int(1, length(master))),
                     tmp_peaks, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  ss <- data.frame(SampleID = samples$sample_id, Factor = samples$mark,
                   Condition = samples$group, Replicate = samples$rep_index, bamReads = samples$bam,
                   Peaks = rep.int(tmp_peaks, nrow(samples)), PeakCaller = rep.int("raw", nrow(samples)),
                   ScoreCol = rep.int(4L, nrow(samples)), stringsAsFactors = FALSE, check.names = FALSE)
  if ("batch" %in% names(samples)) ss$Tissue <- samples$batch
  if (identical(NORMALIZATION, "spikein")) ss$Spikein <- samples$spikein_bam

  dbo <- DiffBind::dba(sampleSheet = ss)
  dbo$config$singleEnd <- identical(samples$layout[[1L]], "SE")
  dbo$config$cores <- as.integer(THREADS)
  dbo$config$RunParallel <- as.integer(THREADS) > 1L
  master_for_count <- master
  names(master_for_count) <- NULL
  S4Vectors::mcols(master_for_count) <- NULL
  dbo <- DiffBind::dba.count(dbo, peaks = master_for_count, score = DiffBind::DBA_SCORE_READS,
                             summits = FALSE, filter = 0, minCount = 0, bSubControl = FALSE,
                             bUseSummarizeOverlaps = TRUE)
  # DiffBind 3.20's dba.count() uses bSubControl while counting but does not retain an
  # explicit FALSE in the DBA object; pv.normalize otherwise defaults a missing field to TRUE.
  # Persist the documented requested setting so normalization cannot silently choose subtraction.
  dbo$bSubControl <- FALSE
  # Record and verify the native DiffBind row order before relabeling rows with master IDs.
  raw_binding <- DiffBind::dba.peakset(dbo, bRetrieve = TRUE, DataType = DiffBind::DBA_DATA_GRANGES)
  if (!inherits(raw_binding, "GRanges")) stop("DiffBind did not return a GRanges raw binding matrix", call. = FALSE)
  raw_index_native <- map_to_master(raw_binding, master)
  dba_to_master <- order(raw_index_native)
  raw_binding <- raw_binding[dba_to_master]
  if (!identical(raw_index_native[dba_to_master], seq_along(master))) {
    stop("Raw binding rows cannot be aligned to master order", call. = FALSE)
  }
  raw_counts <- as.matrix(S4Vectors::mcols(raw_binding))
  expected_ids <- as.character(samples$sample_id)
  if (!identical(colnames(raw_counts), expected_ids)) {
    stop(sprintf("DiffBind raw count columns do not exactly match sample IDs/order (got: %s)", paste(colnames(raw_counts), collapse = ", ")), call. = FALSE)
  }
  rownames(raw_counts) <- master_ids
  row_mapping <- data.frame(dba_row_index = seq_along(master),
                            locus_id = master_ids[raw_index_native], stringsAsFactors = FALSE)

  # Use the configured normalization mode and save the returned factors for interpretation.
  normalization <- if (identical(NORMALIZATION, "spikein")) {
    DiffBind::dba.normalize(dbo, method = DiffBind::DBA_DESEQ2, normalize = DiffBind::DBA_NORM_NATIVE, spikein = TRUE)
  } else {
    DiffBind::dba.normalize(dbo, method = DiffBind::DBA_DESEQ2, normalize = DiffBind::DBA_NORM_NATIVE, background = TRUE)
  }
  norm_info <- DiffBind::dba.normalize(normalization, method = DiffBind::DBA_DESEQ2, bRetrieve = TRUE)
  needed_norm <- c("norm.method", "norm.factors", "lib.method", "lib.sizes")
  if (!is.list(norm_info) || !all(needed_norm %in% names(norm_info)) ||
      length(norm_info$norm.factors) != nrow(samples) || length(norm_info$lib.sizes) != nrow(samples)) {
    stop("DiffBind returned incomplete DESeq2 normalization metadata", call. = FALSE)
  }
  norm_table <- data.frame(sample_id = expected_ids, normalization_mode = NORMALIZATION,
                           norm_method = norm_info$norm.method, norm_factor = as.numeric(norm_info$norm.factors),
                           lib_method = norm_info$lib.method, lib_size = as.numeric(norm_info$lib.sizes),
                           background = isTRUE(norm_info$background),
                           control_subtract = isTRUE(norm_info$control.subtract), stringsAsFactors = FALSE)
  out <- file.path("results/04_diffbind", mark)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(norm_table, file.path(out, "normalization.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(row_mapping, file.path(out, "dba_row_to_locus_id.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  # Include batch only when validation established a multi-batch, estimable design.
  design <- if ("batch" %in% names(samples) && length(unique(samples$batch)) > 1L) "~ Tissue + Condition" else "~ Condition"
  dbo <- DiffBind::dba.contrast(normalization, design = design, contrast = c("Condition", "EK", "EC"),
                                reorderMeta = list(Condition = GROUPS))
  dbo <- DiffBind::dba.analyze(dbo, method = DiffBind::DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
  dds <- DiffBind::dba.analyze(dbo, bRetrieveAnalysis = DiffBind::DBA_DESEQ2)
  saveRDS(dds, file.path(out, paste0(mark, "_dds.rds")))

  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  group_candidates <- names(cd)[vapply(cd, function(x) setequal(as.character(x), GROUPS), logical(1))]
  if (length(group_candidates) != 1L) stop("Could not identify exactly one DESeq2 group column containing EC/EK/LC/LK", call. = FALSE)
  group_col <- group_candidates[[1L]]
  sample_groups <- as.character(cd[[group_col]])
  dds_counts_native <- DESeq2::counts(dds, normalized = FALSE)
  if (!identical(colnames(dds_counts_native), expected_ids) || !identical(sample_groups, as.character(samples$group))) {
    stop("DESeq2 count columns/sample metadata do not exactly preserve validated sample order", call. = FALSE)
  }
  dds_ranges <- tryCatch(SummarizedExperiment::rowRanges(dds), error = function(e) NULL)
  if (inherits(dds_ranges, "GRanges") && length(dds_ranges) == length(master)) {
    dds_index <- map_to_master(dds_ranges, master)
    if (!identical(dds_index, raw_index_native)) stop("DESeq2 row ranges do not preserve DiffBind binding-matrix order", call. = FALSE)
  } else if (!identical(rownames(dds_counts_native), as.character(seq_along(master)))) {
    stop("DiffBind DESeq2 binding-matrix row indices are not the expected 1..N native order", call. = FALSE)
  }
  dds_counts <- dds_counts_native[dba_to_master, , drop = FALSE]
  if (!identical(dim(dds_counts), dim(raw_counts)) || !isTRUE(all.equal(unname(dds_counts), unname(raw_counts), tolerance = 0, check.attributes = FALSE))) {
    stop("DESeq2 counts do not equal the raw DiffBind binding matrix after recorded master reordering", call. = FALSE)
  }

  ref <- if (is.factor(cd[[group_col]])) levels(cd[[group_col]])[[1L]] else {
    group_names <- DESeq2::resultsNames(dds)
    candidate <- GROUPS[vapply(GROUPS, function(g) !any(grepl(paste0("^", group_col, "_", g, "_vs_"), group_names)), logical(1))]
    if (length(candidate) != 1L) stop("Cannot determine the DESeq2 reference group", call. = FALSE)
    candidate
  }
  # A state of 0 is no detected change; N records an untestable locus/contrast.
  all_results <- lapply(names(CONTRASTS), function(contrast_name) {
    weights <- CONTRASTS[[contrast_name]]
    res <- DESeq2::results(dds, contrast = contrast_vector(DESeq2::resultsNames(dds), group_col, ref, weights), alpha = ALPHA)
    result <- as.data.frame(res)[dba_to_master, , drop = FALSE]
    zero <- zero_in_groups(dds_counts, sample_groups, weights)
    counted <- rep(TRUE, length(master))
    data.frame(locus_id = master_ids, contrast = contrast_name, mark = mark,
               baseMean = result$baseMean, log2FoldChange = result$log2FoldChange, lfcSE = result$lfcSE,
               pvalue = result$pvalue, padj = result$padj, zero_counts = zero, counted = counted,
               state = classify_state(result$log2FoldChange, result$pvalue, result$padj, zero),
               reason = state_reason(result$log2FoldChange, result$pvalue, result$padj, zero, counted),
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  all_results <- do.call(rbind, all_results)
  utils::write.table(all_results, file.path(out, paste0(mark, "_all_loci_results.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  grDevices::pdf(file.path(out, "MA_plots.pdf"), width = 7, height = 6, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (contrast_name in names(CONTRASTS)) {
    res <- DESeq2::results(dds, contrast = contrast_vector(DESeq2::resultsNames(dds), group_col, ref, CONTRASTS[[contrast_name]]), alpha = ALPHA)
    DESeq2::plotMA(res, alpha = ALPHA, main = paste(mark, contrast_name))
  }
  cat(sprintf("Completed DiffBind/DESeq2 analysis for %s across %d master loci.\n", mark, length(master)))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: scripts/04_diffbind_per_mark.R <H2AZ|H2AZac|H2AZub>", call. = FALSE)
main(args[[1L]])
