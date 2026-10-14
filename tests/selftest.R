#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
expect_equal <- function(actual, expected, message = "values differ") if (!identical(actual, expected)) stop(message, call. = FALSE)
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = function(e) e$message)
  if (is.null(message) || !grepl(pattern, message, fixed = TRUE)) stop(sprintf("Expected error containing '%s', got '%s'", pattern, message), call. = FALSE)
}

bed <- tempfile(fileext = ".bed")
writeLines("chr1\t0\t100", bed)
gr <- read_bed_gr(bed)
expect_equal(GenomicRanges::start(gr), 1L, "BED start must convert from zero to one based")
expect_equal(GenomicRanges::end(gr), 100L, "BED end must remain unchanged")
writeLines("chr1\t0.5\t100", bed)
expect_error(read_bed_gr(bed), "BED coordinates are invalid")

r1 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(1L, 1001L), c(100L, 1100L)))
r2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50L, 150L))
r3 <- GenomicRanges::GRanges("chr2", IRanges::IRanges(5000L, 5100L))
supported <- replicate_supported(list(r1, r2, r3))
expect_equal(length(supported), 1L, "Only the two-replicate locus should remain")
expect_equal(as.integer(GenomicRanges::start(supported)), 1L)
expect_equal(as.integer(GenomicRanges::end(supported)), 150L)

duplicate_one_replicate <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(1L, 1L), c(100L, 100L)))
expect_equal(length(replicate_supported(list(duplicate_one_replicate, GenomicRanges::GRanges()))), 0L,
             "Duplicate peaks from one replicate must not count as two replicates")

master_input <- GenomicRanges::GRanges(c("chr1", "chr1", "chr1", "chrM", "chrUn", "chr2"),
  IRanges::IRanges(c(1L, 151L, 301L, 1L, 1L, 1L), c(150L, 300L, 400L, 100L, 100L, 100L)))
blacklist <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1L, 100L))
master <- build_master(list(master_input), blacklist)
expect_equal(master_locus_ids(master), c("chr1:0-150", "chr1:150-300", "chr1:300-400"), "Master must retain adjacent ranges separately")
master_rds <- tempfile(fileext = ".rds"); saveRDS(master, master_rds)
expect_equal(master_locus_ids(readRDS(master_rds)), master_locus_ids(master), "Master locus IDs must persist in RDS")
nonoverlap_blacklist <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1000L, 1100L))
expect_equal(length(build_master(list(master_input[1:3]), nonoverlap_blacklist)), 3L, "Non-overlapping blacklist must not remove master ranges")
shift_one <- master; GenomicRanges::start(shift_one)[1L] <- 2L; GenomicRanges::end(shift_one)[1L] <- 151L
expect_equal(map_to_master(shift_one, master), seq_along(master), "One-base coordinate shift must map")
shift_four <- master; GenomicRanges::start(shift_four)[1L] <- 5L; GenomicRanges::end(shift_four)[1L] <- 154L
expect_error(map_to_master(shift_four, master), "strict one-to-one")

ambiguous_master <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(1L, 3L), c(100L, 102L)))
ambiguous_counted <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(2L, 3L), c(101L, 102L)))
expect_error(map_to_master(ambiguous_counted, ambiguous_master), "strict one-to-one")

means <- c(EC = 2, EK = 3.5, LC = 1, LK = 4)
for (ref in c("EC", "LC")) {
  names <- c("Intercept", paste0("batch_B_vs_A"), paste0("Condition_", setdiff(GROUPS, ref), "_vs_", ref))
  for (contrast in CONTRASTS) {
    v <- contrast_vector(names, "Condition", ref, contrast)
    estimated <- sum(contrast * means)
    reconstructed <- sum(v[paste0("Condition_", setdiff(GROUPS, ref), "_vs_", ref)] *
                         (means[setdiff(GROUPS, ref)] - means[[ref]]))
    if (!isTRUE(all.equal(unname(estimated), unname(reconstructed)))) stop("Numeric contrast does not recover requested group mean contrast", call. = FALSE)
  }
}
expect_equal(classify_state(c(1, -1, .2, 1, NA, 1, 1, 1),
                            c(.01, .01, .6, .01, .01, .01, NA, .01),
                            c(.01, .01, .6, .05, .01, .01, .01, NA),
                            c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)),
             c("+", "-", "0", "0", "N", "N", "N", "N"))
partial_zero_counts <- rbind(c(EC = 0, EK = 0, LC = 5, LK = 5), c(EC = 1, EK = 0, LC = 0, LK = 0))
expect_equal(zero_in_groups(partial_zero_counts, colnames(partial_zero_counts), CONTRASTS$early_KO_vs_control), c(TRUE, FALSE),
             "zero_in_groups must consider only groups with nonzero contrast weights")
expect_equal(STATE_LEVELS[1:4], c("+++", "++0", "++-", "+0+"))
expect_equal(tail(STATE_LEVELS, 1L), "---")

calls <- expand.grid(locus_id = c("locus1", "locus2"), contrast = names(CONTRASTS), mark = MARKS, stringsAsFactors = FALSE)
calls$state <- "0"
calls$state[calls$locus_id == "locus1" & calls$contrast == "early_KO_vs_control" & calls$mark == "H2AZ"] <- "+"
calls$state[calls$locus_id == "locus1" & calls$contrast == "early_KO_vs_control" & calls$mark == "H2AZub"] <- "-"
calls$state[calls$locus_id == "locus2" & calls$contrast == "late_KO_vs_control" & calls$mark == "H2AZac"] <- "N"
wide <- joint_states(calls)
expect_equal(wide$joint[wide$locus_id == "locus1" & wide$contrast == "early_KO_vs_control"], "+0-")
primary_loci <- eligible_loci(wide, PRIMARY_CONTRASTS)
expect_equal(primary_loci, "locus1", "An N state must exclude a locus")
tab <- count_states(wide, primary_loci, PRIMARY_CONTRASTS)
if (nrow(tab) != 27L * 4L || any(tapply(tab$n, tab$contrast, sum) != unique(tab$eligible))) stop("State panels must contain all 27 states and sum to eligible", call. = FALSE)
expect_equal(levels(factor(STATE_LEVELS, levels = STATE_LEVELS)), STATE_LEVELS, "Fixed state ordering must be retained")
expect_error(joint_states(calls[-1L, ]), "exactly every locus x five contrasts x three marks")
cat("selftest.R: all helper invariants passed\n")
