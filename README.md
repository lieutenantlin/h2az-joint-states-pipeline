# H2A.Z joint-state pipeline

This directory is the runnable R pipeline described in `../H2A_Z_JOINT_STATES_PLAN.md`. It takes existing mm10 CUT&RUN BAMs for H2AZ, H2AZac and H2AZub and makes joint-state summaries and figures. The input table is intentionally header-only: add one row per biological antibody library to `config/samples.csv` before running.

## Requirements and setup

Use a Linux/macOS shell and conda for the fully isolated, user-local setup. System R may instead be used with the official Bioconductor/CRAN route below. Peak-calling binaries and SEACR are stored outside this source tree.

Preferred reproducible setup (isolated Conda environments):

```sh
bash setup_r_dependencies.sh
bash setup_peak_tools.sh
```

`setup/r-dependencies.yml` documents the R packages; `setup/peak-tools.yml` documents the external peak tools. The verified R environment contains DiffBind 3.20.0, DESeq2 1.50.2, GenomicRanges 1.62.1, ggplot2 4.0.3, Rsamtools 2.26.0 and csaw 1.44.0. To select that R explicitly for a run, use `H2AZ_USE_CONDA_R=1 bash run_pipeline.sh`; the runner sets both `R_LIBS` and `R_LIBS_USER` to that environment's R library to avoid mixing inherited user packages. Otherwise the runner uses the `Rscript` on the caller's `PATH`; the existence of the optional conda R environment does not silently shadow system R. To override its location, set `H2AZ_R_ENV` for both setup and run.

To use an existing system R instead, install the Bioconductor packages via `BiocManager` and ggplot2 from CRAN (do not install ggplot2 with `BiocManager`):

```sh
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("DiffBind", "DESeq2", "GenomicRanges", "Rsamtools", "csaw"), ask=FALSE, update=FALSE); install.packages("ggplot2", repos="https://cloud.r-project.org")'
Rscript -e 'for (x in c("DiffBind","DESeq2","GenomicRanges","ggplot2")) {stopifnot(requireNamespace(x, quietly=TRUE)); cat(x, as.character(packageVersion(x)), "\n")}'
```

On hosts where an R-level Bioconductor update would upgrade many existing packages, use the isolated conda environment above or an isolated `R_LIBS_USER`; do not update the system library as a workaround.

Peak tools default to `$HOME/.local/share/conda/envs/h2az-joint-states/bin`; SEACR's official 1.3 script/helper default to `$HOME/.local/share/h2az-joint-states/`. `config/tools.sh` adds the peak-tools bin directory when present. Override locations with `H2AZ_TOOLS_ENV` and `H2AZ_TOOL_HOME` during setup, then use `H2AZ_TOOLS_ENV`, `SEACR` and `MACS3` when running. The `setup_peak_tools.sh` command prints installed versions and paths.

Check `config/settings.R`: `NORMALIZATION` defaults to `"background"` (relative enrichment, not absolute genome-wide gain/loss); `THREADS` is detected from physical cores; `BLACKLIST` points at the downloaded ENCODE mm10 blacklist v2 at `resources/mm10-blacklist.v2.bed.gz`. The bundled file came from [Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz). Set normalization to `"spikein"` only when all included libraries have validated, comparable calibration BAMs.

## Input files and sample sheet

The default CSV omits the independently optional `animal_id` and `batch` columns:

```text
sample_id,include,mark,stage,genotype,sex,layout,bam,peaks,igg_bam,spikein_bam,bigwig
```

Only `animal_id` and `batch` may be omitted; all remaining columns must appear in the order above. Either optional column may be supplied independently in its original position, after `genotype` and before `sex` (`animal_id` before `batch` when both are present). The complete supported header is:

```text
sample_id,include,mark,stage,genotype,animal_id,batch,sex,layout,bam,peaks,igg_bam,spikein_bam,bigwig
```

Omission means removing the column, not leaving its cells empty. Every included row must have a nonempty value in each supplied optional column; all-empty and partially empty columns are errors. Excluded rows do not participate in these metadata checks. Supplied animal IDs must be unique within each mark; the same animal may appear across different marks. Omitted columns stay absent in the validated inventory: no animal IDs or batch labels are invented.

Independent biological replicates are still required when `animal_id` is omitted; the user is responsible for their independence. Technical replicates are not biological replicates. Supply `batch` when batch adjustment is needed: multiple batch levels select `~ Tissue + Condition`, while an omitted or single-level batch selects `~ Condition`. Removing batch metadata does not correct confounding.

Use `TRUE`/`FALSE` in `include`; only included rows enter the analysis. Values are `mark` = `H2AZ`, `H2AZac` or `H2AZub`; `stage` = `early` or `late`; `genotype` = `control` or `cKO`; `layout` = `PE` or `SE`. Use one row per animal per mark: merge technical runs first. Use valid DiffBind-safe IDs matching `[A-Za-z][A-Za-z0-9_]*`. The BAM must be indexed, mm10-aligned, use `chr` chromosome names, and have consistent alignment/MAPQ/duplicate filtering across samples. Fill `peaks` with a precomputed peak file path to reuse it; otherwise Step 2 calls peaks. `igg_bam` and `spikein_bam` are optional only under the plan's conditions; `bigwig` is for visualization, not statistical inference. Supply file paths readable from the project root; paths with commas cannot be represented safely in this CSV format.

For a missing PE peak set, the pipeline uses SEACR 1.3 and requires samtools and bedtools. It makes one bedGraph per target and reuses each IgG bedGraph for libraries referencing the same IgG BAM. Without IgG, SEACR calls the top 1% of signal blocks; this is exploratory, not an FDR. For a missing SE peak set, it uses MACS3 with `-q 0.01`; without IgG it uses `--nolambda`. Peak calls are replicate-specific. Existing peaks avoid the peak-calling dependencies at run time.

The retained real-tool smoke fixture uses 200 loci across chr1, strand-balanced 50 bp reads (plus starts centered at peak center −100 bp; minus starts centered at +50 bp) and dispersed background. MACS3 built a 197 bp model and produced 200 narrowPeak records both with and without IgG using the planned options; its log, BAMs and results are under `/tmp/h2az-peak-smoke/`. This validates execution, not biological signal.

At least two biological replicates are needed per mark/group. **All four groups (EC, EK, LC and LK) must be represented for each analyzed mark.** The differential scripts fit the planned four-group contrasts and the joint analysis expects the four primary comparisons; setting `include=FALSE` to remove an entire group does not automatically produce a valid reduced-comparison pipeline. If an entire group is unavailable, the contrast/design and requested panels must be re-planned rather than interpreting missing groups as supported results. Multi-level batch effects must be estimable and not confounded with group.

## Run

Fill the sample sheet and run from this directory:

```sh
bash run_pipeline.sh
```

The runner resolves its own project directory, sources the settings for tools, runs `tests/selftest.R` and `tests/inventory_selftest.R` before validating the user inventory, calls peaks only when `results/01_inventory/peak_jobs.tsv` is created, then runs Steps 3–7 in order and stops on the first failure. Results are written under `results/<step>/`; reruns overwrite the same outputs. Do not treat a successful synthetic self-test as biological validation.

## Pipeline script guide

| Order | Path | Role |
|---|---|---|
| Setup | `setup_r_dependencies.sh` | Provisions the isolated R and Bioconductor environment. |
| Setup | `setup_peak_tools.sh` | Provisions samtools, bedtools, MACS3, and external SEACR 1.3 files. |
| Support | `config/tools.sh` | Resolves tool paths, threads, and optional Conda R. |
| Support | `config/settings.R` | Defines normalization, threads, and the mm10 blacklist. |
| Check | `tests/selftest.R` | Checks shared helper invariants before processing samples; it does not validate biological conclusions. |
| Check | `tests/inventory_selftest.R` | Runs Step 01 in disposable projects to check all four supported headers and invalid metadata/header cases; placeholder files prove validation only, not alignment or analysis. |
| 01 | `scripts/01_validate_inventory.R` | Validates the sample sheet and writes the included, ordered inventory and any missing-peak jobs. |
| 02 | `scripts/02_call_peaks.sh` | Calls missing PE peaks with SEACR or missing SE peaks with MACS3. |
| 03 | `scripts/03_build_master_regions.R` | Builds the fixed blacklist-filtered master regions from replicate-supported peaks. |
| 04 | `scripts/04_diffbind_per_mark.R` | Counts and tests each mark against the master regions for all planned contrasts. |
| 05 | `scripts/05_joint_states.R` | Joins mark-level calls into joint states and computes primary/interaction eligibility. |
| 06 | `scripts/06_plot_and_export.R` | Produces figures and state-specific BED files for eligible loci. |
| 07 | `scripts/07_review_loci.R` | Ranks representative loci for manual IGV review. |
| Support | `R/pipeline_functions.R` | Supplies shared region, contrast, state, plotting, and BED-export helpers. |

## Output manifest

- `results/01_inventory/`: validated sample CSV and, only when peaks are missing, `peak_jobs.tsv`.
- `results/02_peaks/`: per-library SEACR `.stringent.bed` or MACS3 `_peaks.narrowPeak` outputs plus generated bedGraphs.
- `results/03_master/`: fixed master RDS, four-column BED with locus IDs, three-column DiffBind counting BED, support/width summaries and `widest10.bed`.
- `results/04_diffbind/<MARK>/`: normalized library metadata, DESeq2 object, all-locus results and MA plots for each mark.
- `results/05_joint/`: per-locus joint calls, primary/interaction state-count tables, testability reasons and eligibility denominators.
- `results/06_figures/`: primary and interaction PDF/PNG figures and populated per-state BEDs under `state_beds/<contrast>/`.
- `results/07_review/`: ranked loci for manual IGV review.

Implementation details that differ from the original plan's prototype notes: chromosome filtering uses `match()` against `KEEP_CHR` and coordinate sorting rather than unavailable `keepSeqlevels`; PNGs are written with a Cairo device directly rather than `ggsave()` to avoid the reported graphics API mismatch. DiffBind initialization uses a temporary four-column raw peak file with a score column; the persisted `master_regions.counting.bed` remains the planned three-column BED. DiffBind 3.20.0 does not retain the `bSubControl=FALSE` argument from `dba.count()` as an explicit object field, so the script records `dbo$bSubControl <- FALSE` before normalization. It also writes `dba_row_to_locus_id.tsv` and verifies native row ordering against DESeq2 counts before reordering results to master loci.

Integration evidence: the synthetic 36-BAM / 200-locus full run completed with background normalization; spike-in normalization and PE count-path checks also passed. The synthetic check validates pipeline behavior only, not biological concordance.

Optional-metadata verification (2026-10-26): `tests/selftest.R` and `tests/inventory_selftest.R` passed; the inventory suite also passed with spike-in settings in a disposable project. A disposable runner check passed both self-tests before rejecting the header-only user inventory. The real analysis smoke used the configured Conda R environment (`H2AZ_USE_CONDA_R=1` via `config/tools.sh`): R 4.5.3, DiffBind 3.20.0 and DESeq2 1.50.2. System R lacked DiffBind and DESeq2; invoking Conda R without the configured library isolation encountered an incompatible user `rlang`.

Both omitted-metadata and populated-metadata smoke projects validated 24 SE libraries and retained 200 master loci with the bundled blacklist. The exact deterministic counts (`20 + ((window_index + replicate_index) %% 7)`, plus eight cKO reads in even windows) completed counting/background normalization but failed DESeq2's dispersion fit in both modes because estimates were too close to the minimum. A supplemental fixture used `set.seed(261026)` and negative-binomial counts with those means (`size = 20`), retaining the same 200 windows and 2,000 background reads per library. It completed real Steps 01, 03 and 04 for H2AZ in both modes: saved DESeq2 sample IDs matched the validated H2AZ inventory, the designs were `~ Condition` without batch and `~ Tissue + Condition` with balanced A/B batches, and every retained locus had all five contrasts. Temporary projects were removed; production results were untouched. These synthetic checks establish execution, not biological validity.

## Interpretation and biological review

`+` and `-` mean statistically detected increase/decrease (per-mark BH FDR < 0.05). `0` means **no detected change**, not demonstrated stability; `N` means untestable and is excluded from the shared eligibility denominator. No fold-change threshold is applied. H2AZac and H2AZub are analyzed as their own signals, not as stoichiometric modification fractions. With background normalization, interpretation is relative enrichment only; absolute genome-wide gain/loss claims are not supported.

Load state BEDs, `master_regions.bed`, and replicate bigWigs in IGV as specified in the plan. Inspect representative loci and the widest regions manually. This repository/task has not established biological agreement in IGV: that assessment requires the lab's actual data and remains the user's validation step.

## Sources

See the plan's Sources section for the DiffBind/DESeq2 references, [SEACR](https://github.com/FredHutch/SEACR), [MACS3](https://macs3-project.github.io/MACS/docs/callpeak.html), and the blacklist source.
