#!/usr/bin/env Rscript
# Executable Step 01 validation fixtures; these do not establish alignment or analysis validity.
expect_true <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
expect_equal <- function(actual, expected, message = "values differ") {
  if (!identical(actual, expected)) stop(message, call. = FALSE)
}
expect_failure <- function(result, pattern) {
  status <- attr(result, "status")
  output <- paste(result, collapse = "\n")
  if (is.null(status) || status == 0L || !grepl(pattern, output, fixed = TRUE)) {
    stop(sprintf("Expected failure containing '%s', got status %s and output:\n%s",
                 pattern, if (is.null(status)) "0" else status, output), call. = FALSE)
  }
}

project_root <- normalizePath(".", mustWork = TRUE)
canonical_columns <- c("sample_id", "include", "mark", "stage", "genotype", "animal_id", "batch",
                       "sex", "layout", "bam", "peaks", "igg_bam", "spikein_bam", "bigwig")
optional_columns <- c("animal_id", "batch")

groups <- data.frame(
  group = c("EC", "EK", "LC", "LK"),
  stage = c("early", "early", "late", "late"),
  genotype = c("control", "cKO", "control", "cKO"),
  stringsAsFactors = FALSE
)

make_rows <- function(include_excluded = FALSE) {
  rows <- do.call(rbind, lapply(c("H2AZ", "H2AZac", "H2AZub"), function(mark) {
    do.call(rbind, lapply(seq_len(nrow(groups)), function(group_index) {
      do.call(rbind, lapply(1:2, function(replicate_index) {
        group <- groups[group_index, ]
        sample_id <- sprintf("%s_%s_r%d", mark, group$group, replicate_index)
        data.frame(
          sample_id = sample_id,
          include = "TRUE",
          mark = mark,
          stage = group$stage,
          genotype = group$genotype,
          animal_id = sprintf("animal_%s_%d", group$group, replicate_index),
          batch = c("A", "B")[replicate_index],
          sex = c("F", "M")[replicate_index],
          layout = "SE",
          bam = file.path("fixtures", paste0(sample_id, ".bam")),
          peaks = file.path("fixtures", paste0(sample_id, ".bed")),
          igg_bam = "",
          spikein_bam = file.path("fixtures", paste0(sample_id, "_spikein.bam")),
          bigwig = "",
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  if (include_excluded) {
    rows <- rbind(rows, data.frame(
      sample_id = "excluded_metadata_free", include = "FALSE", mark = "not_a_mark", stage = "not_a_stage",
      genotype = "not_a_genotype", animal_id = "", batch = "", sex = "", layout = "not_a_layout",
      bam = "", peaks = "", igg_bam = "", spikein_bam = "", bigwig = "", stringsAsFactors = FALSE
    ))
  }
  rows
}

copy_project_files <- function(root) {
  for (directory in c("R", "config", "scripts", "resources", "fixtures")) dir.create(file.path(root, directory), recursive = TRUE)
  copies <- c("R/pipeline_functions.R", "config/settings.R", "scripts/01_validate_inventory.R", "resources/mm10-blacklist.v2.bed.gz")
  for (path in copies) {
    expect_true(file.copy(file.path(project_root, path), file.path(root, path), overwrite = TRUE),
                sprintf("Could not copy %s", path))
  }
}

write_fixture <- function(root, rows, headers, header_override = NULL) {
  dir.create(file.path(root, "fixtures"), showWarnings = FALSE)
  for (bam in unique(rows$bam[nzchar(rows$bam)])) {
    file.create(file.path(root, bam))
    file.create(file.path(root, paste0(bam, ".bai")))
  }
  for (bam in unique(rows$spikein_bam[nzchar(rows$spikein_bam)])) file.create(file.path(root, bam))
  for (peak in unique(rows$peaks[nzchar(rows$peaks)])) file.create(file.path(root, peak))
  path <- file.path(root, "config", "samples.csv")
  utils::write.table(rows[, headers, drop = FALSE], path, sep = ",", row.names = FALSE, col.names = TRUE,
                     quote = TRUE, na = "")
  if (!is.null(header_override)) {
    contents <- readLines(path, warn = FALSE)
    contents[[1L]] <- paste(header_override, collapse = ",")
    writeLines(contents, path)
  }
}

run_fixture <- function(rows, headers, header_override = NULL) {
  root <- tempfile("h2az-inventory-")
  dir.create(root)
  previous_directory <- getwd()
  on.exit({
    setwd(previous_directory)
    unlink(root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  copy_project_files(root)
  write_fixture(root, rows, headers, header_override)
  setwd(root)
  output <- system2(file.path(R.home("bin"), "Rscript"), "scripts/01_validate_inventory.R",
                    stdout = TRUE, stderr = TRUE)
  validated <- NULL
  if (is.null(attr(output, "status"))) {
    validated <- utils::read.csv(file.path("results", "01_inventory", "samples_validated.csv"),
                                 check.names = FALSE, stringsAsFactors = FALSE)
  }
  list(output = output, validated = validated)
}

base_rows <- make_rows()
expect_true(all(table(base_rows$animal_id) == length(unique(base_rows$mark))),
            "Fixture animal IDs must repeat only across marks")
group_labels <- with(base_rows, ifelse(stage == "early" & genotype == "control", "EC",
                                ifelse(stage == "early" & genotype == "cKO", "EK",
                                ifelse(stage == "late" & genotype == "control", "LC", "LK"))))
expected_sample_ids <- with(base_rows, sample_id[order(match(mark, c("H2AZ", "H2AZac", "H2AZub")),
                                                       match(group_labels, c("EC", "EK", "LC", "LK")), sample_id)])
expected_groups <- group_labels[order(match(base_rows$mark, c("H2AZ", "H2AZac", "H2AZub")),
                                      match(group_labels, c("EC", "EK", "LC", "LK")), base_rows$sample_id)]

for (absent in list(character(), "animal_id", "batch", optional_columns)) {
  headers <- canonical_columns[!canonical_columns %in% absent]
  fixture <- run_fixture(base_rows, headers)
  expect_true(is.null(attr(fixture$output, "status")), sprintf("Valid header failed: %s", paste(headers, collapse = ",")))
  validated <- fixture$validated
  expect_equal(validated$sample_id, expected_sample_ids, "Validated sample identities must have stable ordering")
  expect_equal(validated$group, expected_groups, "Validated group assignments must be retained")
  expect_equal(names(validated)[names(validated) %in% optional_columns], optional_columns[optional_columns %in% headers],
               "Validated optional metadata columns must match supplied input")
}

excluded_fixture <- run_fixture(make_rows(include_excluded = TRUE), canonical_columns)
expect_true(is.null(attr(excluded_fixture$output, "status")), "Excluded rows must not require optional metadata")
expect_equal(nrow(excluded_fixture$validated), 24L, "Excluded rows must not enter the validated inventory")

single_batch <- base_rows
single_batch$batch <- "A"
single_batch_fixture <- run_fixture(single_batch, canonical_columns)
expect_true(is.null(attr(single_batch_fixture$output, "status")), "A supplied single batch must remain valid")

duplicate_animal <- base_rows
duplicate_animal$animal_id[duplicate_animal$mark == "H2AZ" & group_labels == "EC"] <- "animal_1"
expect_failure(run_fixture(duplicate_animal, canonical_columns)$output, "Only one included row per animal_id per mark is allowed")

blank_animal <- base_rows
blank_animal$animal_id[[1L]] <- ""
expect_failure(run_fixture(blank_animal, canonical_columns)$output, "animal_id is required for every included row")

empty_animal <- base_rows
empty_animal$animal_id <- ""
expect_failure(run_fixture(empty_animal, canonical_columns)$output, "animal_id is required for every included row")

blank_batch <- base_rows
blank_batch$batch[[1L]] <- ""
expect_failure(run_fixture(blank_batch, canonical_columns)$output, "batch is required for every included row")

empty_batch <- base_rows
empty_batch$batch <- ""
expect_failure(run_fixture(empty_batch, canonical_columns)$output, "batch is required for every included row")

confounded_batch <- base_rows
confounded_batch$batch <- ifelse(confounded_batch$genotype == "control", "A", "B")
expect_failure(run_fixture(confounded_batch, canonical_columns)$output, "batch + group design is rank deficient")

missing_required <- canonical_columns[canonical_columns != "sex"]
expect_failure(run_fixture(base_rows, missing_required)$output, "columns")
duplicate_header <- canonical_columns
duplicate_header[[length(duplicate_header)]] <- "peaks"
expect_failure(run_fixture(base_rows, canonical_columns, duplicate_header)$output, "columns")
unexpected_column <- base_rows
unexpected_column$unexpected <- ""
expect_failure(run_fixture(unexpected_column, c(canonical_columns, "unexpected"))$output, "columns")
reordered <- canonical_columns
reordered[c(1L, 2L)] <- reordered[c(2L, 1L)]
expect_failure(run_fixture(base_rows, reordered)$output, "columns")
reduced_headers <- canonical_columns[canonical_columns != "animal_id"]
misplaced_batch <- append(reduced_headers[reduced_headers != "batch"], "batch", after = match("sex", reduced_headers[reduced_headers != "batch"]))
expect_failure(run_fixture(base_rows, misplaced_batch)$output, "columns")
reversed_optional <- canonical_columns
reversed_optional[c(6L, 7L)] <- reversed_optional[c(7L, 6L)]
expect_failure(run_fixture(base_rows, reversed_optional)$output, "columns")

insufficient_replicates <- base_rows[-1L, , drop = FALSE]
expect_failure(run_fixture(insufficient_replicates, canonical_columns)$output, "Need >=2 biological replicates")

cat("inventory_selftest.R: all inventory validation fixtures passed\n")
