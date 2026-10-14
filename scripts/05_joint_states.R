#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
# Step 05: combine complete per-mark contrast calls into joint states and shared denominators.

master_path <- "results/03_master/master_regions.rds"
if (!file.exists(master_path)) stop("Run Step 03 before Step 05", call. = FALSE)
master <- readRDS(master_path)
expected_loci <- master_locus_ids(master)
required <- c("locus_id", "contrast", "mark", "baseMean", "padj", "state", "reason")
# Reject incomplete mark tables so missing calls cannot be interpreted as no detected change.
tables <- lapply(MARKS, function(mark) {
  path <- file.path("results/04_diffbind", mark, paste0(mark, "_all_loci_results.tsv"))
  if (!file.exists(path)) stop(sprintf("Missing differential-results table: %s", path), call. = FALSE)
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA")
  if (!all(required %in% names(x))) stop(sprintf("Malformed differential-results table: %s", path), call. = FALSE)
  if (!all(x$mark == mark)) stop(sprintf("Result table %s contains an unexpected mark value", path), call. = FALSE)
  if (anyDuplicated(x[c("locus_id", "contrast")]) || nrow(x) != length(expected_loci) * length(CONTRASTS) ||
      !setequal(x$locus_id, expected_loci) || !setequal(x$contrast, names(CONTRASTS))) {
    stop(sprintf("Result table %s must contain exactly master loci x planned contrasts", path), call. = FALSE)
  }
  x$state <- as.character(x$state)
  if (!all(x$state %in% c("+", "0", "-", "N"))) stop(sprintf("Result table %s has invalid states", path), call. = FALSE)
  x
})
calls <- do.call(rbind, tables)
wide <- joint_states(calls)
if (!setequal(unique(wide$locus_id), expected_loci)) stop("Differential result locus IDs do not exactly equal master loci", call. = FALSE)
wide <- wide[order(match(wide$locus_id, expected_loci), match(wide$contrast, names(CONTRASTS))), , drop = FALSE]
primary_loci <- eligible_loci(wide, PRIMARY_CONTRASTS)
interaction_loci <- eligible_loci(wide, "stage_genotype_interaction")
if (!length(primary_loci)) stop("No loci are jointly eligible across all four primary contrasts", call. = FALSE)
if (!length(interaction_loci)) stop("No loci are jointly eligible for the interaction contrast", call. = FALSE)
# Primary panels share loci testable across all four primary contrasts; interaction is separate.
primary_counts <- count_states(wide, primary_loci, PRIMARY_CONTRASTS)
interaction_counts <- count_states(wide, interaction_loci, "stage_genotype_interaction")
reason_table <- as.data.frame(with(calls, table(mark, contrast, reason, useNA = "ifany")), stringsAsFactors = FALSE)
names(reason_table) <- c("mark", "contrast", "reason", "n")
reason_table <- reason_table[reason_table$n > 0L, , drop = FALSE]
# Preserve per-mark testability reasons alongside the state-count summaries.
eligibility <- data.frame(contrast = c(PRIMARY_CONTRASTS, "stage_genotype_interaction"),
                          eligible = c(rep.int(length(primary_loci), length(PRIMARY_CONTRASTS)), length(interaction_loci)),
                          eligible_primary = length(primary_loci),
                          eligible_interaction = length(interaction_loci),
                          denominator_rule = c(rep.int("shared across all four primary contrasts", length(PRIMARY_CONTRASTS)), "interaction-specific"),
                          stringsAsFactors = FALSE)
out <- "results/05_joint"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.table(wide, file.path(out, "joint_locus_states.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(primary_counts, file.path(out, "state_counts_primary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(interaction_counts, file.path(out, "state_counts_interaction.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(reason_table, file.path(out, "testability_by_mark_contrast_reason.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(eligibility, file.path(out, "eligibility_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Joined %d master loci; %d are eligible for every primary panel.\n", length(expected_loci), length(primary_loci)))
