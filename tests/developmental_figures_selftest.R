#!/usr/bin/env Rscript
# Focused checks for development-anchored cohorts, denominators and signal mapping; not biological validation.
source("R/pipeline_functions.R")
source("R/developmental_figures.R")
expect_equal <- function(actual, expected, message = "values differ") {
  if (!isTRUE(all.equal(actual, expected, check.attributes = FALSE))) stop(message, call. = FALSE)
}
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) e$message)
  if (is.null(message) || !grepl(pattern, message, fixed = TRUE)) stop(sprintf("Expected error containing '%s', got '%s'", pattern, message), call. = FALSE)
}
for (pkg in c("DESeq2", "GenomicRanges")) if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required", call. = FALSE)

# Calls consistent with classify_state(): + / - significant, 0 not significant, N untestable.
call_values <- function(state) {
  switch(state,
         "+" = list(lfc = 1.5, p = 0.001, padj = 0.01),
         "-" = list(lfc = -1.5, p = 0.001, padj = 0.01),
         "0" = list(lfc = 0.1, p = 0.5, padj = 0.6),
         "N" = list(lfc = NA_real_, p = NA_real_, padj = NA_real_))
}
make_calls <- function(ids, states = list()) {
  grid <- expand.grid(locus_id = ids, contrast = names(CONTRASTS), mark = MARKS, stringsAsFactors = FALSE)
  grid$state <- "0"
  for (key in names(states)) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
    hit <- grid$mark == parts[[1L]] & grid$contrast == parts[[2L]]
    grid$state[hit] <- states[[key]][match(grid$locus_id[hit], ids)]
  }
  values <- lapply(grid$state, call_values)
  grid$log2FoldChange <- vapply(values, `[[`, numeric(1), "lfc")
  grid$pvalue <- vapply(values, `[[`, numeric(1), "p")
  grid$padj <- vapply(values, `[[`, numeric(1), "padj")
  grid$zero_counts <- FALSE
  grid$reason <- ifelse(grid$state == "N", "pvalue_NA_outlier_or_allzero", "tested")
  grid[CALL_COLUMNS]
}

ids <- paste0("chr1:", seq(0, 700, 100), "-", seq(50, 750, 100))
# Loci 1-5 gain, 6 loses, 7 unchanged, 8 untestable in H2AZub development.
calls <- make_calls(ids, list(
  "H2AZub|control_late_vs_early" = c("+", "+", "+", "+", "+", "-", "0", "N"),
  "H2AZub|late_KO_vs_control" = c("+", "-", "0", "0", "N", "-", "+", "+"),
  "H2AZub|early_KO_vs_control" = c("N", "N", "+", "0", "0", "0", "0", "0"),
  "H2AZac|late_KO_vs_control" = c("-", "-", "-", "N", "0", "0", "0", "0")
))
cohorts <- developmental_cohorts(calls, ids)
ub <- cohorts[cohorts$anchor_mark == "H2AZub", ]
expect_equal(ub$cohort, c(rep("gain", 5), "loss", "no_detected_change", "untestable"), "Control states must define cohorts")
responses <- developmental_responses(calls, cohorts)
counts <- developmental_response_counts(responses)
pick <- function(anchor, cohort, response, stage) counts[counts$anchor_mark == anchor & counts$cohort == cohort &
                                                          counts$response_mark == response & counts$stage == stage, ]
late <- pick("H2AZub", "gain", "H2AZub", "late")
expect_equal(late$n, c(1L, 1L, 2L), "Late same-mark counts must keep no-detected-change loci")
expect_equal(unique(late$cohort_total), 5L)
expect_equal(unique(late$testable), 4L)
expect_equal(unique(late$untestable), 1L)
expect_equal(late$percent, c(25, 25, 50))
early <- pick("H2AZub", "gain", "H2AZub", "early")
expect_equal(c(unique(early$testable), unique(early$untestable)), c(3L, 2L), "Early denominator must be independent of late N calls")
expect_equal(early$n, c(1L, 0L, 2L))
cross <- pick("H2AZub", "gain", "H2AZac", "late")
expect_equal(c(unique(cross$testable), cross$n), c(4L, 0L, 3L, 1L), "Cross-mark denominators use only the response mark's N calls")
loss <- pick("H2AZub", "loss", "H2AZub", "late")
expect_equal(c(unique(loss$cohort_total), loss$n), c(1L, 0L, 1L, 0L), "Loss cohort must contain only the control loss locus")
if (any(c(ids[7], ids[8]) %in% responses$locus_id[responses$anchor_mark == "H2AZub"])) {
  stop("Unchanged and untestable control loci must not enter developmental cohorts", call. = FALSE)
}
empty <- pick("H2AZ", "gain", "H2AZ", "late")
expect_equal(c(unique(empty$cohort_total), empty$n), c(0L, 0L, 0L, 0L), "Empty cohorts keep complete zero rows")
if (!all(is.na(empty$percent))) stop("Empty cohorts must have NA percentages", call. = FALSE)
expect_equal(nrow(counts), length(MARKS) * 2L * length(MARKS) * 2L * 3L, "Response grid must be complete")

all_n <- make_calls(ids, list("H2AZ|control_late_vs_early" = c("+", "+", rep("0", 6)),
                              "H2AZ|late_KO_vs_control" = rep("N", 8)))
all_n_counts <- developmental_response_counts(developmental_responses(all_n, developmental_cohorts(all_n, ids)))
x <- all_n_counts[all_n_counts$anchor_mark == "H2AZ" & all_n_counts$cohort == "gain" & all_n_counts$response_mark == "H2AZ" & all_n_counts$stage == "late", ]
expect_equal(c(unique(x$cohort_total), unique(x$testable), x$n), c(2L, 0L, 0L, 0L, 0L), "All-N responses are not 100% no change")
if (!all(is.na(x$percent))) stop("All-N responses must have NA percentages", call. = FALSE)

# cKO and interaction outcomes cannot alter control-defined membership.
changed <- calls
hit <- changed$contrast != "control_late_vs_early"
changed[hit, ] <- make_calls(ids, list(
  "H2AZub|control_late_vs_early" = c("+", "+", "+", "+", "+", "-", "0", "N"),
  "H2AZ|late_KO_vs_control" = rep("+", 8), "H2AZub|stage_genotype_interaction" = rep("-", 8)))[hit, ]
changed_responses <- developmental_responses(changed, developmental_cohorts(changed, ids))
expect_equal(changed_responses[c("anchor_mark", "cohort", "locus_id", "response_mark")],
             responses[c("anchor_mark", "cohort", "locus_id", "response_mark")], "cKO outcomes must not change membership")
sel <- responses$anchor_mark == "H2AZub" & responses$cohort == "gain" & responses$response_mark == "H2AZac"
expect_equal(responses$late_state[sel], c("-", "-", "-", "N", "0"), "Cross-mark responses must stay attached to anchor locus IDs")

# Result-table validation from disposable TSVs.
tmp <- tempfile("h2az-dev-calls-")
dir.create(tmp)
paths <- stats::setNames(file.path(tmp, paste0(MARKS, ".tsv")), MARKS)
for (mark in MARKS) {
  x <- calls[calls$mark == mark, ]
  x$baseMean <- 10; x$lfcSE <- 0.2; x$counted <- TRUE
  utils::write.table(x, paths[[mark]], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
expect_equal(nrow(read_developmental_calls(paths, ids)), nrow(calls))
bad <- utils::read.delim(paths[["H2AZ"]], colClasses = "character")
bad$state[bad$state == "0"][1L] <- "+"
utils::write.table(bad, paths[["H2AZ"]], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
expect_error(read_developmental_calls(paths, ids), "inconsistent with classify_state")
bad <- rbind(bad[1L, ], bad[-nrow(bad), ])
utils::write.table(bad, paths[["H2AZ"]], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
expect_error(read_developmental_calls(paths, ids), "duplicate locus_id/contrast keys")
unlink(tmp, recursive = TRUE)

# Replicate signal: native row permutation, size factors, normalization factors and identity checks.
master <- build_master(list(GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1L, 201L, 401L), c(100L, 300L, 500L)))),
                       GenomicRanges::GRanges())
master_ids <- master_locus_ids(master)
samples <- data.frame(sample_id = c("ec1", "ec2", "ek1", "ek2", "lc1", "lc2", "lk1", "lk2"),
                      group = rep(GROUPS, each = 2L), stringsAsFactors = FALSE)
by_master <- rbind(c(10, 20, 30, 40, 50, 60, 70, 80), rep(c(5, 10), 4), c(0, 2, 4, 6, 8, 10, 12, 14))
native_order <- c(3L, 1L, 2L)  # native row i holds master locus native_order[i]
make_dds <- function() {
  counts <- by_master[native_order, ]
  storage.mode(counts) <- "integer"
  dimnames(counts) <- list(as.character(seq_len(3L)), samples$sample_id)
  cd <- data.frame(Condition = factor(samples$group, levels = GROUPS), row.names = samples$sample_id)
  suppressMessages(DESeq2::DESeqDataSetFromMatrix(counts, cd, ~ Condition))
}
row_map <- data.frame(dba_row_index = 1:3, locus_id = master_ids[native_order], stringsAsFactors = FALSE)
size_factors <- c(1, 2, 1, 2, 1, 2, 1, 2)
dds <- make_dds()
DESeq2::sizeFactors(dds) <- size_factors
signal <- developmental_signal(dds, row_map, samples, master)
expect_equal(rownames(signal$log2_counts), master_ids, "Signal rows must be in master order")
expect_equal(signal$log2_counts[master_ids[1], "ec2"], log2(20 / 2 + 1), "Size factors must normalize the named locus/sample")
expect_equal(signal$log2_counts[master_ids[3], "lk2"], log2(14 / 2 + 1))
expect_equal(unname(signal$row_z[master_ids[2], ]), rep(0, 8), "Constant rows must have zero z-scores")
z1 <- log2(by_master[1, ] / size_factors + 1)
expect_equal(unname(signal$row_z[master_ids[1], ]), (z1 - mean(z1)) / stats::sd(z1))
dds_nf <- make_dds()
factors <- matrix(1, 3, 8, dimnames = dimnames(dds_nf))
factors[which(native_order == 1L), "ec1"] <- 4
DESeq2::normalizationFactors(dds_nf) <- factors
signal_nf <- developmental_signal(dds_nf, row_map, samples, master)
expect_equal(signal_nf$log2_counts[master_ids[1], "ec1"], log2(10 / 4 + 1), "Per-locus normalization factors must be honored")
expect_equal(signal_nf$log2_counts[master_ids[3], "ec1"], log2(0 + 1))
expect_error(developmental_signal(dds, transform(row_map, locus_id = master_ids[c(1, 1, 2)]), samples, master), "bijectively")
expect_error(developmental_signal(dds, row_map[1:2, ], samples, master), "bijectively")
expect_error(developmental_signal(dds, row_map, samples[c(2, 1, 3:8), ], master), "sample IDs")
swapped <- samples; swapped$group[1:2] <- "LC"; swapped$group[5:6] <- "EC"
expect_error(developmental_signal(dds, row_map, swapped, master), "group labels")
dds_ranges <- dds
SummarizedExperiment::rowRanges(dds_ranges) <- unname(master[c(1L, 2L, 3L)])
expect_error(developmental_signal(dds_ranges, row_map, samples, master), "row ranges disagree")
SummarizedExperiment::rowRanges(dds_ranges) <- unname(master[native_order])
expect_equal(developmental_signal(dds_ranges, row_map, samples, master)$log2_counts, signal$log2_counts,
             "Matching row ranges must be accepted")

# Heatmap pagination and scatter accounting.
many <- sprintf("chr2:%06d-%06d", seq_len(501L) * 1000L, seq_len(501L) * 1000L + 50L)
many_calls <- make_calls(many, list("H2AZ|control_late_vs_early" = rep("+", 501L),
                                    "H2AZ|late_KO_vs_control" = c("N", rep("0", 499L), "+")))
many_calls$log2FoldChange[many_calls$mark == "H2AZ" & many_calls$contrast == "control_late_vs_early"] <- seq(1, 2, length.out = 501L)
many_responses <- developmental_responses(many_calls, developmental_cohorts(many_calls, many))
manifest <- heatmap_row_manifest(many_responses)
h <- manifest[manifest$anchor_mark == "H2AZ" & manifest$cohort == "gain", ]
expect_equal(as.integer(table(h$page)), c(500L, 1L), "501 loci must paginate as 500 + 1")
expect_equal(anyDuplicated(h$locus_id), 0L)
expect_equal(h$locus_id[1:2], c(many[501], many[500]), "Rows order by late state, then descending developmental LFC")
expect_equal(h$locus_id[501], many[1], "Untestable anchor responses sort last")
if (nrow(manifest[manifest$anchor_mark != "H2AZ", ])) stop("Other anchors have no cohorts here", call. = FALSE)

scatter_responses <- responses
row <- which(scatter_responses$anchor_mark == "H2AZub" & scatter_responses$cohort == "gain" & scatter_responses$response_mark == "H2AZub")
scatter_responses$late_lfc[row[1]] <- NA_real_
scatter_responses$late_lfc[row[5]] <- 0.4  # finite estimate with an N state stays plotted
scatter_responses$interaction_state[row[2:3]] <- c("+", "-")
scatter_responses$interaction_state[row[1]] <- "-"
summary <- developmental_scatter_summary(scatter_responses)
s <- summary[summary$anchor_mark == "H2AZub" & summary$cohort == "gain" & summary$response_mark == "H2AZub", ]
expect_equal(c(s$cohort_total, s$plotted, s$excluded_nonfinite, s$interaction_detected), c(5L, 4L, 1L, 2L),
             "Scatter excludes only nonfinite estimates and counts interactions among plotted loci")
cat("developmental_figures_selftest.R: all developmental cohort invariants passed\n")
