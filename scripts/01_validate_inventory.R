#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
source("config/settings.R")
# Step 01: validate analysis inputs and write a deterministic inventory for downstream stages.

required_columns <- c("sample_id", "include", "mark", "stage", "genotype", "animal_id", "batch",
                      "sex", "layout", "bam", "peaks", "igg_bam", "spikein_bam", "bigwig")
path <- "config/samples.csv"
if (!file.exists(path)) stop("Missing config/samples.csv", call. = FALSE)
samples <- utils::read.csv(path, colClasses = "character", na.strings = character(), check.names = FALSE,
                           stringsAsFactors = FALSE)
errors <- character()
if (!identical(names(samples), required_columns)) {
  errors <- c(errors, sprintf("samples.csv columns must be exactly: %s", paste(required_columns, collapse = ", ")))
}
if (!nrow(samples)) errors <- c(errors, "samples.csv has no data rows")
if (nrow(samples) && identical(names(samples), required_columns)) {
  samples[is.na(samples)] <- ""
  if (any(!samples$include %in% c("TRUE", "FALSE"))) {
    errors <- c(errors, "include must be TRUE or FALSE")
  }
  if (anyDuplicated(samples$sample_id)) errors <- c(errors, "sample_id values must be unique")
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", samples$sample_id))) {
    errors <- c(errors, "sample_id must match ^[A-Za-z][A-Za-z0-9_]*$")
  }
  # Excluded rows remain in the source sheet but do not enter analysis or peak calling.
  included <- samples[samples$include == "TRUE", , drop = FALSE]
  if (!nrow(included)) errors <- c(errors, "No included samples")
  if (nrow(included)) {
    valid_mark <- included$mark %in% MARKS
    valid_stage <- included$stage %in% c("early", "late")
    valid_genotype <- included$genotype %in% c("control", "cKO")
    valid_layout <- included$layout %in% c("PE", "SE")
    if (any(!valid_mark)) errors <- c(errors, "mark must be one of H2AZ, H2AZac, H2AZub")
    if (any(!valid_stage)) errors <- c(errors, "stage must be early or late")
    if (any(!valid_genotype)) errors <- c(errors, "genotype must be control or cKO")
    if (any(!valid_layout)) errors <- c(errors, "layout must be PE or SE")
    if (any(!nzchar(included$animal_id))) errors <- c(errors, "animal_id is required for every included row")
    if (any(!nzchar(included$batch))) errors <- c(errors, "batch is required for every included row")
    if (anyDuplicated(included[c("mark", "animal_id")])) errors <- c(errors, "Only one included row per animal_id per mark is allowed")
    for (mark in intersect(MARKS, unique(included$mark))) {
      layouts <- unique(included$layout[included$mark == mark])
      if (length(layouts) > 1L) errors <- c(errors, sprintf("layout must be identical within mark %s", mark))
    }
    # Encode the planned stage x genotype groups used by every later contrast.
    included$group <- NA_character_
    valid_group <- valid_stage & valid_genotype
    included$group[valid_group] <- ifelse(included$stage[valid_group] == "early" & included$genotype[valid_group] == "control", "EC",
                                    ifelse(included$stage[valid_group] == "early" & included$genotype[valid_group] == "cKO", "EK",
                                    ifelse(included$stage[valid_group] == "late" & included$genotype[valid_group] == "control", "LC", "LK")))
    for (i in seq_len(nrow(included))) {
      bam <- included$bam[[i]]
      if (!nzchar(bam) || !file.exists(bam)) {
        errors <- c(errors, sprintf("BAM is missing for %s: %s", included$sample_id[[i]], bam))
      } else if (!file.exists(paste0(bam, ".bai")) && !file.exists(sub("\\.bam$", ".bai", bam, ignore.case = TRUE))) {
        errors <- c(errors, sprintf("BAM index is missing for %s: %s", included$sample_id[[i]], bam))
      }
      for (field in c("peaks", "igg_bam")) {
        value <- included[[field]][[i]]
        if (nzchar(value) && !file.exists(value)) errors <- c(errors, sprintf("%s is missing for %s: %s", field, included$sample_id[[i]], value))
      }
      if (identical(NORMALIZATION, "spikein")) {
        value <- included$spikein_bam[[i]]
        if (!nzchar(value) || !file.exists(value)) errors <- c(errors, sprintf("spikein_bam is required and missing for %s", included$sample_id[[i]]))
      }
    }
    if (!identical(NORMALIZATION, "background") && !identical(NORMALIZATION, "spikein")) {
      errors <- c(errors, "NORMALIZATION must be 'background' or 'spikein'")
    }
    if (all(valid_mark & valid_group)) {
      replication <- with(included, table(factor(mark, levels = MARKS), factor(group, levels = GROUPS)))
      cat("Biological replicate table (mark x group):\n")
      print(replication)
      short <- which(replication < 2L, arr.ind = TRUE)
      if (nrow(short)) errors <- c(errors, apply(short, 1L, function(z) sprintf("Need >=2 biological replicates for %s x %s", rownames(replication)[z[[1L]]], colnames(replication)[z[[2L]]])))
      # A batch term is used later only when it is estimable with the four groups.
      for (mark in MARKS) {
        x <- included[included$mark == mark, , drop = FALSE]
        if (length(unique(x$batch)) > 1L) {
          group <- factor(x$group, levels = GROUPS)
          batch <- factor(x$batch)
          model <- stats::model.matrix(~ batch + group)
          if (qr(model)$rank < ncol(model)) {
            cat(sprintf("Confounded batch/group table for %s:\n", mark)); print(table(batch, group))
            errors <- c(errors, sprintf("batch + group design is rank deficient for mark %s", mark))
          }
        }
        if (length(unique(x$sex[nzchar(x$sex)])) > 1L) {
          cat(sprintf("Sex by group for %s (informational):\n", mark)); print(table(x$group, x$sex))
        }
      }
    }
    # Stable order and within-group indices make the downstream sample sheet reproducible.
    included <- included[order(match(included$mark, MARKS), match(included$group, GROUPS), included$sample_id), , drop = FALSE]
    included$rep_index <- ave(included$sample_id, interaction(included$mark, included$group, drop = TRUE), FUN = seq_along)
  }
}
if (length(errors)) stop(paste(c("Sample inventory validation failed:", paste0("- ", unique(errors))), collapse = "\n"), call. = FALSE)
dir.create("results/01_inventory", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(included, "results/01_inventory/samples_validated.csv", row.names = FALSE, quote = TRUE)
# Only included samples without supplied peaks become Step 02 jobs.
jobs <- included[!nzchar(included$peaks), c("sample_id", "layout", "bam", "igg_bam"), drop = FALSE]
job_path <- "results/01_inventory/peak_jobs.tsv"
if (nrow(jobs)) {
  utils::write.table(jobs, job_path, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (file.exists(job_path)) {
  unlink(job_path)
}
cat(sprintf("Validated %d included samples.\n", nrow(included)))
