# Development-anchored Bap1 response helpers: cohorts, response frequencies, replicate signal,
# and figure builders. Source after R/pipeline_functions.R.

DEVELOPMENTAL_CONTRAST <- "control_late_vs_early"
DEVELOPMENTAL_COHORTS <- c("gain", "loss")
COHORT_BY_STATE <- c("+" = "gain", "-" = "loss", "0" = "no_detected_change", "N" = "untestable")
COHORT_LABELS <- c(gain = "Developmental gain (control late > early)", loss = "Developmental loss (control late < early)")
RESPONSE_PREFIXES <- c(early = "early_KO_vs_control", late = "late_KO_vs_control",
                       ko_developmental = "KO_late_vs_early", interaction = "stage_genotype_interaction")
RESPONSE_STAGES <- c("early", "late")
STAGE_LABELS <- c(early = "Early cKO vs control", late = "Late cKO vs control")
RESPONSE_STATE_LEVELS <- c("+", "-", "0", "N")
RESPONSE_STATE_LABELS <- c("+" = "Increased", "-" = "Decreased", "0" = "No detected change", "N" = "Untestable")
RESPONSE_STATE_COLORS <- c(Increased = "#D55E00", Decreased = "#0072B2", `No detected change` = "#999999", Untestable = "#CC79A7")
INTERACTION_SHAPES <- c(Increased = 2, Decreased = 6, `No detected change` = 1, Untestable = 4)
HEATMAP_GROUP_ORDER <- c("EC", "LC", "EK", "LK")
HEATMAP_GROUP_LABELS <- c(EC = "Early control", LC = "Late control", EK = "Early cKO", LK = "Late cKO")
HEATMAP_PAGE_SIZE <- 500L
HEATMAP_DISPLAY_LIMIT <- 3
CALL_COLUMNS <- c("locus_id", "contrast", "mark", "log2FoldChange", "pvalue", "padj", "zero_counts", "state", "reason")

# Step 04 result loading and validation.

read_developmental_calls <- function(paths, master_ids) {
  if (is.null(names(paths)) || !setequal(names(paths), MARKS)) stopf("Result paths must be named by MARKS")
  required <- c("locus_id", "contrast", "mark", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
                "zero_counts", "counted", "state", "reason")
  pieces <- lapply(MARKS, function(mark) {
    path <- paths[[mark]]
    if (!file.exists(path)) stopf("Missing Step 04 result table: %s", path)
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA")
    missing <- setdiff(required, names(x))
    if (length(missing)) stopf("%s lacks required columns: %s", path, paste(missing, collapse = ", "))
    if (!all(x$mark == mark)) stopf("%s contains rows for a mark other than %s", path, mark)
    if (anyDuplicated(x[c("locus_id", "contrast")])) stopf("%s has duplicate locus_id/contrast keys", path)
    if (nrow(x) != length(master_ids) * length(CONTRASTS) || !setequal(x$locus_id, master_ids) ||
        !setequal(x$contrast, names(CONTRASTS))) {
      stopf("%s must contain exactly master loci x planned contrasts", path)
    }
    for (column in c("log2FoldChange", "pvalue", "padj")) {
      value <- suppressWarnings(as.numeric(x[[column]]))
      if (any(is.na(value) & !is.na(x[[column]]))) stopf("%s has non-numeric %s values", path, column)
      x[[column]] <- value
    }
    if (!all(x$zero_counts %in% c("TRUE", "FALSE"))) stopf("%s has invalid zero_counts values", path)
    x$zero_counts <- x$zero_counts == "TRUE"
    if (!all(x$state %in% RESPONSE_STATE_LEVELS)) stopf("%s has invalid states", path)
    expected <- classify_state(x$log2FoldChange, x$pvalue, x$padj, x$zero_counts)
    if (!identical(expected, x$state)) stopf("%s has states inconsistent with classify_state()", path)
    x[CALL_COLUMNS]
  })
  do.call(rbind, pieces)
}

# Cohort definition and cKO response joins.

developmental_cohorts <- function(calls, master_ids) {
  pieces <- lapply(MARKS, function(mark) {
    x <- calls[calls$mark == mark & calls$contrast == DEVELOPMENTAL_CONTRAST, , drop = FALSE]
    ix <- match(master_ids, x$locus_id)
    if (anyNA(ix) || nrow(x) != length(master_ids)) stopf("Developmental calls for %s do not cover master loci exactly", mark)
    x <- x[ix, , drop = FALSE]
    data.frame(anchor_mark = mark, locus_id = master_ids, cohort = unname(COHORT_BY_STATE[x$state]),
               developmental_state = x$state, developmental_lfc = x$log2FoldChange,
               developmental_padj = x$padj, developmental_reason = x$reason, stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

developmental_responses <- function(calls, cohorts) {
  members <- cohorts[cohorts$cohort %in% DEVELOPMENTAL_COHORTS, , drop = FALSE]
  members <- members[order(match(members$anchor_mark, MARKS), match(members$cohort, DEVELOPMENTAL_COHORTS),
                           seq_len(nrow(members))), , drop = FALSE]
  pieces <- lapply(MARKS, function(response_mark) {
    out <- data.frame(anchor_mark = members$anchor_mark, cohort = members$cohort, locus_id = members$locus_id,
                      response_mark = rep(response_mark, nrow(members)),
                      developmental_lfc = members$developmental_lfc, developmental_padj = members$developmental_padj,
                      stringsAsFactors = FALSE)
    for (prefix in names(RESPONSE_PREFIXES)) {
      x <- calls[calls$mark == response_mark & calls$contrast == RESPONSE_PREFIXES[[prefix]], , drop = FALSE]
      ix <- match(members$locus_id, x$locus_id)
      if (anyNA(ix)) stopf("%s calls for %s lack developmental cohort loci", RESPONSE_PREFIXES[[prefix]], response_mark)
      out[[paste0(prefix, "_lfc")]] <- x$log2FoldChange[ix]
      out[[paste0(prefix, "_padj")]] <- x$padj[ix]
      out[[paste0(prefix, "_state")]] <- x$state[ix]
      out[[paste0(prefix, "_reason")]] <- x$reason[ix]
    }
    out$member_index <- seq_len(nrow(members))
    out
  })
  out <- do.call(rbind, pieces)
  out <- out[order(match(out$anchor_mark, MARKS), match(out$cohort, DEVELOPMENTAL_COHORTS),
                   match(out$response_mark, MARKS), out$member_index), , drop = FALSE]
  out$member_index <- NULL
  rownames(out) <- NULL
  out
}

# Denominators are the selected response mark/stage's testable cohort loci; 0 stays counted.
developmental_response_counts <- function(responses) {
  rows <- list()
  for (anchor in MARKS) for (cohort in DEVELOPMENTAL_COHORTS) for (response in MARKS) for (stage in RESPONSE_STAGES) {
    states <- responses[[paste0(stage, "_state")]][responses$anchor_mark == anchor & responses$cohort == cohort &
                                                      responses$response_mark == response]
    total <- length(states)
    untestable <- sum(states == "N")
    testable <- total - untestable
    n <- vapply(c("+", "-", "0"), function(s) sum(states == s), integer(1))
    if (sum(n) != testable) stopf("Response counts disagree with testable denominator")
    rows[[length(rows) + 1L]] <- data.frame(anchor_mark = anchor, cohort = cohort, response_mark = response,
                                            stage = stage, state = c("+", "-", "0"), n = unname(n),
                                            cohort_total = total, testable = testable, untestable = untestable,
                                            percent = if (testable) 100 * unname(n) / testable else NA_real_,
                                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Replicate signal from the saved Step 04 DESeq2 object, keyed by the persisted row map.

developmental_signal <- function(dds, row_mapping, samples, master) {
  if (!inherits(dds, "DESeqDataSet")) stopf("Expected a DESeqDataSet")
  master_ids <- master_locus_ids(master)
  n <- nrow(dds)
  if (n != length(master_ids)) stopf("DESeq2 rows (%d) do not match master loci (%d)", n, length(master_ids))
  nf <- DESeq2::normalizationFactors(dds)
  sf <- DESeq2::sizeFactors(dds)
  if (!is.null(nf)) {
    if (!all(is.finite(nf)) || any(nf <= 0)) stopf("normalizationFactors must be positive and finite")
  } else if (is.null(sf) || !all(is.finite(sf)) || any(sf <= 0)) {
    stopf("DESeq2 object lacks positive finite size or normalization factors")
  }
  if (!identical(colnames(dds), as.character(samples$sample_id))) {
    stopf("DESeq2 sample IDs do not exactly match the validated inventory order")
  }
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  group_cols <- names(cd)[vapply(cd, function(x) setequal(as.character(x), GROUPS), logical(1))]
  if (length(group_cols) != 1L) stopf("Could not identify exactly one DESeq2 group column containing EC/EK/LC/LK")
  if (!identical(as.character(cd[[group_cols]]), as.character(samples$group))) {
    stopf("DESeq2 group labels do not match the validated inventory")
  }
  if (!all(c("dba_row_index", "locus_id") %in% names(row_mapping))) stopf("Row map lacks dba_row_index/locus_id")
  index <- suppressWarnings(as.integer(row_mapping$dba_row_index))
  if (nrow(row_mapping) != n || anyNA(index) || !identical(sort(index), seq_len(n)) ||
      anyDuplicated(row_mapping$locus_id) || !setequal(row_mapping$locus_id, master_ids)) {
    stopf("Row map must bijectively map DESeq2 rows 1..%d to master locus IDs", n)
  }
  native_ids <- character(n)
  native_ids[index] <- as.character(row_mapping$locus_id)
  ranges <- tryCatch(SummarizedExperiment::rowRanges(dds), error = function(e) NULL)
  if (inherits(ranges, "GRanges")) {
    if (!identical(master_ids[map_to_master(ranges, master)], native_ids)) {
      stopf("DESeq2 row ranges disagree with the persisted row map")
    }
  } else if (!identical(rownames(dds), as.character(seq_len(n)))) {
    stopf("DESeq2 row names are not the expected native 1..N order")
  }
  normalized <- DESeq2::counts(dds, normalized = TRUE)
  if (!all(is.finite(normalized)) || any(normalized < 0)) stopf("Normalized counts must be finite and nonnegative")
  normalized <- normalized[match(master_ids, native_ids), , drop = FALSE]
  dimnames(normalized) <- list(master_ids, as.character(samples$sample_id))
  log2_counts <- log2(normalized + 1)
  centers <- rowMeans(log2_counts)
  spreads <- apply(log2_counts, 1L, stats::sd)
  row_z <- (log2_counts - centers) / spreads
  row_z[is.na(spreads) | spreads == 0, ] <- 0
  list(log2_counts = log2_counts, row_z = row_z)
}

# Deterministic heatmap row manifest shared by every response mark.

heatmap_row_manifest <- function(responses, page_size = HEATMAP_PAGE_SIZE) {
  pieces <- list()
  for (anchor in MARKS) for (cohort in DEVELOPMENTAL_COHORTS) {
    x <- responses[responses$anchor_mark == anchor & responses$cohort == cohort & responses$response_mark == anchor, , drop = FALSE]
    if (!nrow(x)) next
    x <- x[order(match(x$late_state, RESPONSE_STATE_LEVELS), -abs(x$developmental_lfc), x$locus_id, method = "radix"), , drop = FALSE]
    position <- seq_len(nrow(x)) - 1L
    pieces[[length(pieces) + 1L]] <- data.frame(anchor_mark = anchor, cohort = cohort,
                                                page = position %/% page_size + 1L,
                                                row_in_page = position %% page_size + 1L,
                                                locus_id = x$locus_id, anchor_late_state = x$late_state,
                                                stringsAsFactors = FALSE)
  }
  if (!length(pieces)) {
    return(data.frame(anchor_mark = character(), cohort = character(), page = integer(), row_in_page = integer(),
                      locus_id = character(), anchor_late_state = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, pieces)
}

developmental_scatter_summary <- function(responses) {
  rows <- list()
  for (anchor in MARKS) for (cohort in DEVELOPMENTAL_COHORTS) for (response in MARKS) {
    x <- responses[responses$anchor_mark == anchor & responses$cohort == cohort & responses$response_mark == response, , drop = FALSE]
    finite <- is.finite(x$developmental_lfc) & is.finite(x$late_lfc)
    rows[[length(rows) + 1L]] <- data.frame(anchor_mark = anchor, cohort = cohort, response_mark = response,
                                            cohort_total = nrow(x), plotted = sum(finite),
                                            excluded_nonfinite = sum(!finite),
                                            interaction_detected = sum(finite & x$interaction_state %in% c("+", "-")),
                                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Captions shared by the three figure families.

developmental_caption <- function(normalization, family, multi_batch = FALSE) {
  normalization_text <- if (identical(normalization, "background")) {
    "Background normalization: relative enrichment; no absolute genome-wide gain/loss claim."
  } else {
    "Calibrated spike-in normalization."
  }
  parts <- c(
    "+/- = per-mark BH FDR < 0.05 (no fold-change threshold); 0 = no detected change, not demonstrated stability; N = untestable.",
    normalization_text,
    "H2AZac and H2AZub are separate enrichment signals, not stoichiometric modification fractions.",
    "Descriptive frequencies conditional on control-based cohort selection (control late vs early), not an independently validated estimate of Bap1 dependence.",
    "Same-mark early/late responses and interactions reuse control samples used for selection, so noise can favor early increases, late decreases and negative interactions in gain cohorts (opposite for loss cohorts).",
    "Cross-mark panels avoid shared libraries but may share animals and biological variation."
  )
  if (identical(family, "heatmap")) {
    parts <- c(parts, "Fill: row z-score of log2(DESeq2-normalized count + 1) across all samples of the displayed mark; colors clipped to [-3, 3]. Control trajectories are selected by construction; normalized counts are descriptive, not model estimates.")
    if (multi_batch) parts <- c(parts, "Multiple batches present: displayed normalized counts are not batch-corrected.")
  } else if (identical(family, "scatter")) {
    parts <- c(parts, "Same-mark x and y share late-control samples, so an apparent anticorrelation is not independent evidence of reversal. Shape = response-mark interaction (LK-EK)-(LC-EC): an altered trajectory, not evidence of direct BAP1 binding nor a failure-to-acquire classification.")
  }
  paste(strwrap(paste(parts, collapse = " "), width = 190), collapse = "\n")
}

# Figure builders; callers own graphics devices.

plot_developmental_frequencies <- function(counts, anchor_mark, caption) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stopf("ggplot2 is required for plotting")
  x <- counts[counts$anchor_mark == anchor_mark, , drop = FALSE]
  if (!nrow(x)) stopf("No response counts for anchor mark %s", anchor_mark)
  x$cohort_f <- factor(unname(COHORT_LABELS[x$cohort]), levels = unname(COHORT_LABELS))
  x$response_f <- factor(x$response_mark, levels = MARKS)
  x$stage_f <- factor(unname(STAGE_LABELS[x$stage]), levels = unname(STAGE_LABELS))
  fill_levels <- unname(RESPONSE_STATE_LABELS[c("+", "-", "0")])
  x$state_f <- factor(unname(RESPONSE_STATE_LABELS[x$state]), levels = fill_levels)
  bars <- x[x$testable > 0 & x$n > 0, , drop = FALSE]
  bars$label <- sprintf("%d (%.1f%%)", bars$n, bars$percent)
  summary <- unique(x[c("cohort_f", "response_f", "stage_f", "cohort_total", "testable", "untestable")])
  summary$label <- sprintf("%d/%d testable\nuntestable = %d", summary$testable, summary$cohort_total, summary$untestable)
  empty <- summary[summary$testable == 0, , drop = FALSE]
  empty$label <- ifelse(empty$cohort_total == 0, "No developmental loci", "No testable response loci")
  ggplot2::ggplot() +
    ggplot2::geom_col(data = bars, ggplot2::aes(x = stage_f, y = percent, fill = state_f), width = 0.7) +
    ggplot2::geom_text(data = bars, ggplot2::aes(x = stage_f, y = percent, label = label, group = state_f),
                       position = ggplot2::position_stack(vjust = 0.5), size = 2.6) +
    ggplot2::geom_text(data = summary, ggplot2::aes(x = stage_f, y = 113, label = label), size = 2.6, lineheight = 0.9) +
    ggplot2::geom_text(data = empty, ggplot2::aes(x = stage_f, y = 50, label = label), size = 2.8, fontface = "italic") +
    ggplot2::facet_grid(cohort_f ~ response_f, drop = FALSE,
                        labeller = ggplot2::labeller(response_f = function(v) paste("Response mark:", v),
                                                     cohort_f = ggplot2::label_wrap_gen(24))) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 124), breaks = seq(0, 100, 25)) +
    ggplot2::scale_fill_manual(values = RESPONSE_STATE_COLORS[fill_levels], limits = fill_levels, drop = FALSE,
                               name = "cKO response") +
    ggplot2::labs(title = sprintf("Bap1 cKO responses at %s developmental loci", anchor_mark),
                  x = NULL, y = "% of testable cohort loci", caption = caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7), legend.position = "bottom")
}

plot_developmental_heatmap <- function(signal_z, rows, samples, title, caption) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stopf("ggplot2 is required for plotting")
  theme <- ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7))
  if (!nrow(rows)) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0, label = "No developmental loci in this cohort", size = 6) +
             ggplot2::labs(title = title, caption = caption) + ggplot2::theme_void() + theme)
  }
  samples <- samples[order(match(samples$group, HEATMAP_GROUP_ORDER), seq_len(nrow(samples))), , drop = FALSE]
  if (anyNA(match(rows$locus_id, rownames(signal_z)))) stopf("Heatmap rows are absent from the signal matrix")
  if (anyNA(match(samples$sample_id, colnames(signal_z)))) stopf("Heatmap samples are absent from the signal matrix")
  z <- signal_z[rows$locus_id, samples$sample_id, drop = FALSE]
  long <- data.frame(locus_id = rep(rows$locus_id, times = ncol(z)),
                     sample_id = rep(samples$sample_id, each = nrow(z)),
                     group = rep(samples$group, each = nrow(z)),
                     state = rep(rows$anchor_late_state, times = ncol(z)),
                     display = pmin(pmax(as.vector(z), -HEATMAP_DISPLAY_LIMIT), HEATMAP_DISPLAY_LIMIT),
                     stringsAsFactors = FALSE)
  long$locus_f <- factor(long$locus_id, levels = rev(rows$locus_id))
  long$sample_f <- factor(long$sample_id, levels = samples$sample_id)
  long$group_f <- factor(unname(HEATMAP_GROUP_LABELS[long$group]), levels = unname(HEATMAP_GROUP_LABELS[HEATMAP_GROUP_ORDER]))
  state_labels <- paste("Anchor late cKO:", RESPONSE_STATE_LABELS[RESPONSE_STATE_LEVELS])
  long$state_f <- factor(paste("Anchor late cKO:", RESPONSE_STATE_LABELS[long$state]), levels = state_labels)
  y_text <- if (nrow(rows) > 50L) ggplot2::element_blank() else ggplot2::element_text(size = 6)
  ggplot2::ggplot(long, ggplot2::aes(x = sample_f, y = locus_f, fill = display)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(state_f ~ group_f, scales = "free", space = "free",
                        labeller = ggplot2::labeller(state_f = ggplot2::label_wrap_gen(16))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                                  limits = c(-HEATMAP_DISPLAY_LIMIT, HEATMAP_DISPLAY_LIMIT),
                                  name = "Row z\nlog2(norm+1)") +
    ggplot2::labs(title = title, x = "Replicate sample", y = sprintf("Developmental loci (n = %d on page)", nrow(rows)),
                  caption = caption) +
    ggplot2::theme_bw() + theme +
    ggplot2::theme(axis.text.y = y_text, axis.ticks.y = if (nrow(rows) > 50L) ggplot2::element_blank() else ggplot2::element_line(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                   panel.spacing = ggplot2::unit(2, "pt"), strip.text.y = ggplot2::element_text(angle = 0, size = 7),
                   panel.grid = ggplot2::element_blank())
}

plot_developmental_scatter <- function(responses, summary, anchor_mark, caption) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stopf("ggplot2 is required for plotting")
  x <- responses[responses$anchor_mark == anchor_mark, , drop = FALSE]
  x <- x[is.finite(x$developmental_lfc) & is.finite(x$late_lfc), , drop = FALSE]
  s <- summary[summary$anchor_mark == anchor_mark, , drop = FALSE]
  extent <- function(v) { m <- if (length(v)) max(abs(v)) else 0; if (m > 0) m * 1.05 else 1 }
  ex <- extent(x$developmental_lfc)
  ey <- extent(x$late_lfc)
  factorize <- function(d) {
    d$cohort_f <- factor(unname(COHORT_LABELS[d$cohort]), levels = unname(COHORT_LABELS))
    d$response_f <- factor(d$response_mark, levels = MARKS)
    d
  }
  x <- factorize(x)
  s <- factorize(s)
  state_levels <- unname(RESPONSE_STATE_LABELS[RESPONSE_STATE_LEVELS])
  x$late_f <- factor(unname(RESPONSE_STATE_LABELS[x$late_state]), levels = state_levels)
  x$interaction_f <- factor(unname(RESPONSE_STATE_LABELS[x$interaction_state]), levels = state_levels)
  s$label <- ifelse(s$cohort_total == 0, "No developmental loci",
                    ifelse(s$plotted == 0, "No finite effect estimates",
                           sprintf("cohort n = %d; plotted = %d\nnonfinite excluded = %d; interaction calls = %d",
                                   s$cohort_total, s$plotted, s$excluded_nonfinite, s$interaction_detected)))
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_point(data = x, ggplot2::aes(x = developmental_lfc, y = late_lfc, colour = late_f, shape = interaction_f),
                        alpha = 0.5, size = 1.2, show.legend = TRUE) +
    ggplot2::geom_text(data = s, ggplot2::aes(x = -ex, y = ey, label = label), hjust = 0, vjust = 1, size = 2.6, lineheight = 0.9) +
    ggplot2::facet_grid(cohort_f ~ response_f, drop = FALSE,
                        labeller = ggplot2::labeller(response_f = function(v) paste("Response mark:", v),
                                                     cohort_f = ggplot2::label_wrap_gen(24))) +
    ggplot2::coord_cartesian(xlim = c(-ex, ex), ylim = c(-ey, ey)) +
    ggplot2::scale_colour_manual(values = RESPONSE_STATE_COLORS, limits = state_levels, drop = FALSE,
                                 name = "Late cKO response") +
    ggplot2::scale_shape_manual(values = INTERACTION_SHAPES, limits = state_levels, drop = FALSE,
                                name = "Interaction") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 16, alpha = 1, size = 2.5)),
                    shape = ggplot2::guide_legend(override.aes = list(colour = "black", alpha = 1, size = 2.5))) +
    ggplot2::labs(title = sprintf("%s control developmental change versus late cKO effect", anchor_mark),
                  x = sprintf("%s control late vs early log2 fold change", anchor_mark),
                  y = "Response mark late cKO vs control log2 fold change", caption = caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7), legend.position = "bottom")
}
