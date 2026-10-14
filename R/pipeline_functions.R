# Shared constants and helpers for the ordered H2A.Z joint-state pipeline.
MARKS <- c("H2AZ", "H2AZac", "H2AZub")
GROUPS <- c("EC", "EK", "LC", "LK")
KEEP_CHR <- paste0("chr", c(1:19, "X"))
ALPHA <- 0.05
CONTRASTS <- list(
  early_KO_vs_control = c(EC = -1, EK = 1, LC = 0, LK = 0),
  late_KO_vs_control = c(EC = 0, EK = 0, LC = -1, LK = 1),
  control_late_vs_early = c(EC = -1, EK = 0, LC = 1, LK = 0),
  KO_late_vs_early = c(EC = 0, EK = -1, LC = 0, LK = 1),
  stage_genotype_interaction = c(EC = 1, EK = -1, LC = -1, LK = 1)
)
PRIMARY_CONTRASTS <- names(CONTRASTS)[1:4]
.state_symbols <- c("+", "0", "-")
STATE_LEVELS <- with(
  expand.grid(H2AZub = .state_symbols, H2AZac = .state_symbols,
              H2AZ = .state_symbols, stringsAsFactors = FALSE),
  paste0(H2AZ, H2AZac, H2AZub)
)

# BED/GRanges and fixed-master-region integrity helpers.

stopf <- function(...) stop(sprintf(...), call. = FALSE)

read_bed_gr <- function(path) {
  if (!file.exists(path)) stopf("BED file does not exist: %s", path)
  x <- tryCatch(utils::read.table(path, header = FALSE, sep = "\t", quote = "",
                                  comment.char = "", stringsAsFactors = FALSE,
                                  fill = FALSE),
                error = function(e) stopf("Cannot read BED file %s: %s", path, e$message))
  if (!nrow(x)) stopf("BED file is empty: %s", path)
  if (ncol(x) < 3L) stopf("BED file needs at least three columns: %s", path)
  start_numeric <- suppressWarnings(as.numeric(x[[2L]]))
  end_numeric <- suppressWarnings(as.numeric(x[[3L]]))
  starts <- suppressWarnings(as.integer(x[[2L]]))
  ends <- suppressWarnings(as.integer(x[[3L]]))
  if (anyNA(start_numeric) || anyNA(end_numeric) || any(start_numeric != starts) ||
      any(end_numeric != ends) || anyNA(starts) || anyNA(ends) || any(starts < 0L) ||
      any(ends <= starts)) {
    stopf("BED coordinates are invalid: %s", path)
  }
  # BED starts are zero-based; GRanges starts are one-based.
  GenomicRanges::GRanges(seqnames = x[[1L]], ranges = IRanges::IRanges(starts + 1L, ends))
}

replicate_supported <- function(list_of_replicate_GRanges, min_reps = 2L) {
  if (!length(list_of_replicate_GRanges)) return(GenomicRanges::GRanges())
  if (min_reps < 1L) stopf("min_reps must be at least one")
  if (any(!vapply(list_of_replicate_GRanges, inherits, logical(1), what = "GRanges"))) {
    stopf("replicate_supported requires a list of GRanges")
  }
  nonempty <- Filter(length, list_of_replicate_GRanges)
  if (!length(nonempty)) return(GenomicRanges::GRanges())
  # Count each candidate at most once per replicate, so duplicate peaks cannot inflate support.
  candidates <- GenomicRanges::reduce(Reduce(function(left, right) c(left, right), nonempty), min.gapwidth = 0L)
  supported <- integer(length(candidates))
  for (replicate in nonempty) {
    hits <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(candidates, replicate)))
    supported[hits] <- supported[hits] + 1L
  }
  candidates[supported >= min_reps]
}

build_master <- function(list_of_supported, blacklist) {
  if (!inherits(blacklist, "GRanges")) stopf("blacklist must be a GRanges")
  if (!length(list_of_supported)) stopf("No supported regions were supplied")
  nonempty <- Filter(length, list_of_supported)
  if (!length(nonempty)) stopf("No supported regions were supplied")
  # Keep directly adjacent regions distinct while merging true overlaps.
  master <- GenomicRanges::reduce(Reduce(function(left, right) c(left, right), nonempty), min.gapwidth = 0L)
  master <- master[as.character(GenomicRanges::seqnames(master)) %in% KEEP_CHR]
  if (length(blacklist) && length(master)) {
    drop <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(master, blacklist)))
    if (length(drop)) master <- master[-drop]
  }
  if (!length(master)) return(master)
  chr_rank <- match(as.character(GenomicRanges::seqnames(master)), KEEP_CHR)
  master <- master[order(chr_rank, GenomicRanges::start(master), GenomicRanges::end(master))]
  ids <- sprintf("%s:%d-%d", as.character(GenomicRanges::seqnames(master)),
                 GenomicRanges::start(master) - 1L, GenomicRanges::end(master))
  names(master) <- ids
  S4Vectors::mcols(master)$locus_id <- ids
  master
}

master_locus_ids <- function(master) {
  if (!inherits(master, "GRanges")) stopf("master must be a GRanges")
  ids <- if ("locus_id" %in% names(S4Vectors::mcols(master))) as.character(S4Vectors::mcols(master)$locus_id) else names(master)
  if (is.null(ids) || length(ids) != length(master) || anyNA(ids) || anyDuplicated(ids)) {
    stopf("Master regions lack unique persisted locus_id values")
  }
  ids
}

map_to_master <- function(counted, master) {
  if (!inherits(counted, "GRanges") || !inherits(master, "GRanges")) {
    stopf("map_to_master requires two GRanges objects")
  }
  if (length(counted) != length(master)) {
    stopf("Counted intervals (%d) do not match master intervals (%d)", length(counted), length(master))
  }
  # Permit only a one-base coordinate shift while requiring a unique bidirectional mapping.
  hits <- GenomicRanges::findOverlaps(counted, master, type = "equal", maxgap = 1L)
  qn <- tabulate(S4Vectors::queryHits(hits), nbins = length(counted))
  sn <- tabulate(S4Vectors::subjectHits(hits), nbins = length(master))
  if (any(qn != 1L) || any(sn != 1L)) {
    stopf("Counted and master intervals are not in strict one-to-one correspondence (only 1 bp coordinate tolerance is allowed)")
  }
  master_index <- integer(length(counted))
  master_index[S4Vectors::queryHits(hits)] <- S4Vectors::subjectHits(hits)
  master_index
}

# DESeq2 contrast and per-mark state helpers.

contrast_vector <- function(results_names, group_col, ref, weights) {
  if (is.null(names(weights)) || !setequal(names(weights), GROUPS)) {
    stopf("Contrast weights must be named EC, EK, LC, LK")
  }
  if (!ref %in% GROUPS) stopf("Unknown reference group: %s", ref)
  if (abs(sum(weights)) > sqrt(.Machine$double.eps)) {
    stopf("Cell-mean contrast weights must sum to zero")
  }
  out <- stats::setNames(numeric(length(results_names)), results_names)
  expected <- setdiff(GROUPS, ref)
  for (level in expected) {
    coefficient <- paste0(group_col, "_", level, "_vs_", ref)
    if (!coefficient %in% results_names) {
      stopf("DESeq2 coefficient is absent: %s (available: %s)", coefficient,
            paste(results_names, collapse = ", "))
    }
    out[[coefficient]] <- weights[[level]]
  }
  out
}

zero_in_groups <- function(counts, sample_groups, weights) {
  counts <- as.matrix(counts)
  if (length(sample_groups) != ncol(counts)) stopf("sample_groups must match count columns")
  groups <- names(weights)[weights != 0]
  if (!length(groups) || !all(groups %in% sample_groups)) stopf("Contrast groups are absent from counts")
  rowSums(counts[, sample_groups %in% groups, drop = FALSE] != 0L) == 0L
}

# 0 means no detected change; N means the locus is untestable for this contrast.
classify_state <- function(lfc, pvalue, padj, zero_counts = FALSE) {
  missing <- is.na(lfc) | is.na(pvalue) | is.na(padj) | zero_counts
  out <- rep("0", length(lfc))
  out[!missing & padj < ALPHA & lfc > 0] <- "+"
  out[!missing & padj < ALPHA & lfc < 0] <- "-"
  out[missing] <- "N"
  out
}

state_reason <- function(lfc, pvalue, padj, zero_counts = FALSE, counted = TRUE) {
  n <- max(length(lfc), length(pvalue), length(padj), length(zero_counts), length(counted))
  lfc <- rep(lfc, length.out = n); pvalue <- rep(pvalue, length.out = n)
  padj <- rep(padj, length.out = n); zero_counts <- rep(zero_counts, length.out = n)
  counted <- rep(counted, length.out = n)
  out <- rep("tested", n)
  out[!counted] <- "not_counted"
  out[counted & zero_counts] <- "zero_counts_in_contrast_groups"
  out[counted & !zero_counts & is.na(pvalue)] <- "pvalue_NA_outlier_or_allzero"
  out[counted & !zero_counts & !is.na(pvalue) & is.na(padj)] <- "padj_NA_independent_filtering"
  out[counted & !zero_counts & !is.na(pvalue) & !is.na(padj) & is.na(lfc)] <- "lfc_NA"
  out
}

# Joint-state construction and eligibility/counting helpers.

joint_states <- function(calls) {
  needed <- c("locus_id", "contrast", "mark", "state")
  if (!all(needed %in% names(calls))) stopf("calls must contain: %s", paste(needed, collapse = ", "))
  if (anyDuplicated(calls[c("locus_id", "contrast", "mark")])) stopf("Duplicate locus/contrast/mark calls")
  if (!all(calls$mark %in% MARKS) || !all(calls$contrast %in% names(CONTRASTS))) {
    stopf("Calls contain unknown marks or contrasts")
  }
  loci <- unique(as.character(calls$locus_id))
  complete <- expand.grid(locus_id = loci, contrast = names(CONTRASTS), mark = MARKS,
                          stringsAsFactors = FALSE)
  observed <- merge(complete, calls[needed], by = needed[1:3], all.x = TRUE, sort = FALSE)
  if (nrow(calls) != nrow(complete) || anyNA(observed$state)) {
    stopf("Joint input must contain exactly every locus x five contrasts x three marks")
  }
  keys <- unique(observed[c("locus_id", "contrast")])
  keys <- keys[order(keys$locus_id, match(keys$contrast, names(CONTRASTS))), , drop = FALSE]
  result <- keys
  for (mark in MARKS) {
    x <- observed[observed$mark == mark, c("locus_id", "contrast", "state"), drop = FALSE]
    names(x)[3L] <- mark
    result <- merge(result, x, by = c("locus_id", "contrast"), all.x = TRUE, sort = FALSE)
  }
  result <- result[match(paste(keys$locus_id, keys$contrast), paste(result$locus_id, result$contrast)), , drop = FALSE]
  result$joint <- do.call(paste0, result[MARKS])
  result
}

eligible_loci <- function(wide, contrasts) {
  if (!all(c("locus_id", "contrast", MARKS) %in% names(wide))) {
    stopf("wide must contain locus_id, contrast, and all mark state columns")
  }
  contrast_names <- if (is.list(contrasts)) names(contrasts) else as.character(contrasts)
  if (!length(contrast_names) || !all(contrast_names %in% wide$contrast)) {
    stopf("Requested contrasts are absent")
  }
  state_ok <- apply(wide[MARKS], 1L, function(x) all(!is.na(x) & x != "N"))
  complete_by_locus <- vapply(split(seq_len(nrow(wide)), wide$locus_id), function(ix) {
    rows <- wide[ix, , drop = FALSE]
    all(contrast_names %in% rows$contrast) && all(state_ok[ix][match(contrast_names, rows$contrast)])
  }, logical(1))
  as.character(names(complete_by_locus)[complete_by_locus])
}

# Joint-state counting, plotting, and IGV-oriented BED export helpers.

count_states <- function(wide, loci, contrasts) {
  if (!all(c("locus_id", "contrast", "joint") %in% names(wide))) {
    stopf("wide must contain locus_id, contrast, and joint")
  }
  contrast_names <- if (is.list(contrasts)) names(contrasts) else as.character(contrasts)
  pieces <- lapply(contrast_names, function(contrast) {
    x <- wide[wide$contrast == contrast & wide$locus_id %in% loci, , drop = FALSE]
    counts <- table(factor(x$joint, levels = STATE_LEVELS))
    data.frame(contrast = contrast, state = STATE_LEVELS, n = as.integer(counts),
               eligible = nrow(x), percent = if (nrow(x)) 100 * as.integer(counts) / nrow(x) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

plot_joint_states <- function(tab, caption = "+ increased, - decreased, 0 no detected change (BH FDR < 0.05, per mark)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stopf("ggplot2 is required for plotting")
  needed <- c("contrast", "state", "n", "eligible")
  if (!all(needed %in% names(tab))) stopf("tab must contain: %s", paste(needed, collapse = ", "))
  tab$state <- factor(as.character(tab$state), levels = STATE_LEVELS)
  tab$state_label <- vapply(strsplit(as.character(tab$state), "", fixed = TRUE), paste, collapse = "\n", FUN.VALUE = character(1))
  labels <- unique(tab[c("contrast", "eligible")])
  labels$contrast_label <- sprintf("%s (n = %d)", labels$contrast, labels$eligible)
  tab$contrast_label <- factor(labels$contrast_label[match(tab$contrast, labels$contrast)], levels = labels$contrast_label)
  ggplot2::ggplot(tab, ggplot2::aes(x = state, y = n)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(labels = stats::setNames(tab$state_label, tab$state), drop = FALSE) +
    ggplot2::facet_wrap(~ contrast_label, ncol = 2, scales = "fixed", axes = "all") +
    ggplot2::labs(x = "Rows, top to bottom: H2A.Z / H2A.Zac / H2A.Zub", y = "# of loci",
                  caption = caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(family = "mono"))
}

write_state_beds <- function(states, master, output_dir) {
  if (!inherits(master, "GRanges")) stopf("master must be a GRanges")
  if (!all(c("locus_id", "contrast", "joint") %in% names(states))) {
    stopf("states must contain locus_id, contrast, and joint")
  }
  ids <- master_locus_ids(master)
  if (anyNA(match(states$locus_id, ids))) stopf("State table contains locus IDs absent from master regions")
  for (contrast in unique(states$contrast)) {
    x <- states[states$contrast == contrast, , drop = FALSE]
    dir <- file.path(output_dir, contrast)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    stale <- Sys.glob(file.path(dir, sprintf("%s_???.bed", contrast)))
    if (length(stale)) unlink(stale)
    for (state in unique(x$joint)) {
      y <- x[x$joint == state, , drop = FALSE]
      if (!nrow(y)) next
      dir <- file.path(output_dir, contrast)
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      # Use filesystem-safe P/Z/M codes for +, 0, and - joint-state symbols.
      code <- chartr("+0-", "PZM", state)
      ranges <- master[match(y$locus_id, ids)]
      bed <- data.frame(as.character(GenomicRanges::seqnames(ranges)),
                        GenomicRanges::start(ranges) - 1L,
                        GenomicRanges::end(ranges),
                        y$locus_id, stringsAsFactors = FALSE)
      utils::write.table(bed, file.path(dir, sprintf("%s_%s.bed", contrast, code)), sep = "\t",
                         row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  invisible(NULL)
}
