#!/usr/bin/env Rscript
source("R/pipeline_functions.R")
source("config/settings.R")
# Step 06: render joint-state summaries and export eligible state-specific BEDs for review.
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for Step 06", call. = FALSE)
read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing Step 05 output: %s", path), call. = FALSE)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = "NA")
}
wide <- read_tsv("results/05_joint/joint_locus_states.tsv")
primary <- read_tsv("results/05_joint/state_counts_primary.tsv")
interaction <- read_tsv("results/05_joint/state_counts_interaction.tsv")
summary <- read_tsv("results/05_joint/eligibility_summary.tsv")
for (x in list(primary, interaction)) {
  if (!all(c("contrast", "state", "n", "eligible") %in% names(x))) stop("Malformed state-count table", call. = FALSE)
}
primary$n <- as.numeric(primary$n); primary$eligible <- as.numeric(primary$eligible)
interaction$n <- as.numeric(interaction$n); interaction$eligible <- as.numeric(interaction$eligible)
# Guard against figures whose state bars disagree with their reported eligibility denominators.
if (!all(vapply(PRIMARY_CONTRASTS, function(ct) sum(primary$n[primary$contrast == ct]) == unique(primary$eligible[primary$contrast == ct]), logical(1)))) {
  stop("Primary state bars do not sum to their eligible denominators", call. = FALSE)
}
if (sum(interaction$n[interaction$contrast == "stage_genotype_interaction"]) !=
    unique(interaction$eligible[interaction$contrast == "stage_genotype_interaction"])) {
  stop("Interaction state bars do not sum to their eligible denominator", call. = FALSE)
}
out <- "results/06_figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
# Carry the selected normalization's interpretive limitation into each figure caption.
normalization_caption <- if (identical(NORMALIZATION, "background")) {
  "+ increased, - decreased, 0 no detected change (BH FDR < 0.05, per mark). Background normalization: relative enrichment; no genome-wide gain/loss claim."
} else {
  "+ increased, - decreased, 0 no detected change (BH FDR < 0.05, per mark). Calibrated spike-in normalization."
}
primary_plot <- plot_joint_states(primary, caption = normalization_caption)
grDevices::pdf(file.path(out, "joint_states_primary.pdf"), width = 13, height = 7)
print(primary_plot)
grDevices::dev.off()
grDevices::png(file.path(out, "joint_states_primary.png"), width = 1950, height = 1050, res = 150, type = "cairo")
print(primary_plot)
grDevices::dev.off()
interaction_plot <- plot_joint_states(interaction, caption = normalization_caption) +
  ggplot2::ggtitle("Stage x genotype interaction: + = cKO effect more positive at late stage")
grDevices::pdf(file.path(out, "joint_states_interaction.pdf"), width = 13, height = 7)
print(interaction_plot)
grDevices::dev.off()
grDevices::png(file.path(out, "joint_states_interaction.png"), width = 1950, height = 1050, res = 150, type = "cairo")
print(interaction_plot)
grDevices::dev.off()
primary_loci <- eligible_loci(wide, PRIMARY_CONTRASTS)
# Check the recomputed primary BED membership count against the exported denominator.
primary_eligible <- as.integer(summary$eligible[summary$contrast == PRIMARY_CONTRASTS[[1L]]])
if (length(primary_loci) != primary_eligible) stop("Primary BED denominator does not match eligibility summary", call. = FALSE)
# The interaction uses its own eligibility denominator, unlike the shared primary panels.
interaction_loci <- eligible_loci(wide, "stage_genotype_interaction")
state_rows <- rbind(
  wide[wide$contrast %in% PRIMARY_CONTRASTS & wide$locus_id %in% unique(primary_loci), c("locus_id", "contrast", "joint"), drop = FALSE],
  wide[wide$contrast == "stage_genotype_interaction" & wide$locus_id %in% unique(interaction_loci), c("locus_id", "contrast", "joint"), drop = FALSE]
)
write_state_beds(state_rows, readRDS("results/03_master/master_regions.rds"), file.path(out, "state_beds"))
cat("Wrote primary and interaction figures plus populated state BEDs.\n")
