# Changelog

All notable changes to this project are documented here.

---

## 2026-06-16 — minorfixes

**Changed**
- `RNAseq_shiny_v2.5.R`: updated BAM file path from `/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/` to `/project/def-rallard/COMMUN/raredisease_rnaseq/results_nextflow_rnasplice_09_05_2026/star_salmon/`.
- `RNAseq_shiny_v2.5.R`: merged the gene model slider and plot into a single card so the slider is aligned with the plot.
- `scripts/Shiny/rnaseq_shinyhelper_functions.R`, `RNAseq_shiny_v2.5.R`, `RNAseq_shiny_v2.4.R`: renamed `plotting_coverage()` to `genemodel_plot()`.
- `scripts/Shiny/global.R`, `scripts/Shiny/reactive_module.R`, `RNAseq_shiny_v2.5.R`: the app no longer throws an error when `gwASE.tsv` or `gwImprinted.tsv` are absent — missing files are handled gracefully with NULL checks, empty table fallbacks, and a logged warning.
- `scripts/FRASER/fraser.R`: wrapped the FRASER pipeline in a `tryCatch()` block to handle the intermittent featureCounts SAF error ("no features were loaded"), which occurs when a candidate gene region has no split reads. On failure, an empty CSV is written to prevent the sample from being re-run indefinitely.
- `scripts/cp_and_cleanup.R`: section 7 cp commands now guarded by `Sys.glob()` checks — files are only copied if they exist.
- `.gitignore`: added `scripts/**/.*` to exclude hidden files (e.g. `.DS_Store`) at any depth under `scripts/`.
- `RNAseq_shiny_v2.5.R` show alll hpo terms for a proband.
---

## 2026-06-15 — code-standardisation

**Added**
- Unit tests for `plot_expression_cohort()` and `plot_expression_family()` in `scripts/tests/testthat/test_helper_functions.R`.
- `plot_expression_cohort()` and `plot_expression_family()` helper functions extracted from `RNAseq_shiny_v2.5.R` into `scripts/Shiny/rnaseq_shinyhelper_functions.R`.

**Changed**
- Full code standardisation pass across all R scripts: consistent `=` assignments, `TRUE`/`FALSE` literals, section headers, indentation, and removal of dead commented-out code.
- Hard-coded paths updated for `/project/` HPC mount.
- FRASER and consensus scripts now skip already-processed samples.

---

## 2026-05-19 — ase

**Changed**
- Refactored `ASE.R`.
- Refactored `RNAseq_shiny` (now `RNAseq_shiny_v2.5.R`); modularised reactive inputs into `reactive_module.R`.

---

## 2026-04-24 — ase

**Changed**
- Allele-specific expression analysis v1.
- Minor UI improvements to the Shiny dashboard.
- Consensus: added exon-skipping event detection.

---

## 2026-03-27 — lr_quant

**Changed**
- Long-read RNA-seq analysis added under `scripts/lr_quant/`, based on the [PacBio Iso-Seq reference pipeline](https://isoseq.how/getting-started.html).
- Added download button for the gene prioritisation table.

**Fixed**
- Bug fixes for various edge cases.

---

## 2026-03-27 — bugfixes

**Fixed**
- MANE transcript deduplication (duplicate entries can exist).
- Gene prioritisation table and FRASER summary statistics.

---

## 2026-03-26

**Fixed**
- Bug fixes related to mitochondrial (MT) genes.

---

## 2026-03-24

**Changed**
- Allele-specific expression module added to the Shiny dashboard.

---

## 2026-03-16 — gene-prioritisation

**Fixed**
- Gene prioritisation table logic updated and corrected.

---

## 2026-02-26 — gene-prioritisation

**Fixed**
- Added HPO term names to the gene prioritisation table.

---

## 2026-02-24

**Fixed**
- Minor improvements to gene prioritisation table content.
- Pass per-exon min p-value and max z-score to the genome-wide OUTRIDER view.

---

## 2026-02-20 — gene-prioritisation

**Added**
- Gene prioritisation table combining HPO terms, FRASER, and OUTRIDER gene lists.

**Fixed**
- Minor dashboard improvements.

---

## 2026-02-19 — manhattan-plot

**Fixed**
- OUTRIDER Manhattan plot: point shape now reflects fold-change direction (over/under-expression).

---

## 2026-02-18

**Changed**
- Split the proband selector into a gene-level and sample-level selector.

---

## 2026-02-17 — gwFRASER

**Added**
- Genome-wide OUTRIDER Manhattan plot (manhattan plot function extended accordingly).

**Changed**
- Standardised OUTRIDER output column names.
- BAM subsets now limited to ±5 kb around candidate genes to reduce storage.

**Fixed**
- Updated README.
- Corrected number of genes shown on the genome-wide FRASER plot.

---

## 2026-01-29 — OUTRIDER-per-exon

**Added**
- OUTRIDER results for all candidate exons.
- Genome-wide FRASER results (Manhattan plot and data table).

---

## 2026-01-28 — modern-UI

**Changed**
- Migrated Shiny UI to `bslib` with Bootstrap 5 theme.

---

## 2026-01-22 — error-handling

**Changed**
- Produce placeholder plots when no gene or expression data is available.
- Removed automatic unzipping of the data archive (too cumbersome).

**Fixed**
- Several bugs triggered when multiple candidate genes were present.

---

## 2026-01-20 — new-genes

**Added**
- OUTRIDER p-value threshold selector in the dashboard.

**Changed**
- Parameters now read from a single `configs.json` file.
- `RNAseq_shiny_v2.4.R`: refactored proband index `i` to be defined reactively.

**Fixed**
- Minor refactor to accommodate new gene types (including lncRNAs).
- Bug: not all BAM alignments shown when multiple candidate genes present.
- Bug: incorrect `gsub('03$')` pattern for sample `HSJ_032`.

---

## 2026-01-09

**Changed**
- Minor refactor to accommodate new samples.
- Moved as many parameters as possible into the `.slurm` scripts.

---

## 2026-01-06 — coverage-slider

**Added**
- Slider widget for the coverage plot genomic window.

---

## 2025-12-18 — v2.4

**Changed**
- Updated README, LICENSE, and VERSION.json.

---

## 2025-12-17 — v2.4

**Added**
- `rnaseq_helper_functions.R`: coverage plotting and gene model preparation extracted from `fraser.R`.
- Option to display ±1 kb or ±5 kb coverage window.

**Changed**
- Reference transcript now selected by MANE (Matched Annotation from NCBI and EMBL-EBI).
- Major refactor of `fraser.R`.
- Gene model preparation moved to `featureCounts.R`.

**Fixed**
- Minor Shiny app simplifications.

---

## 2025-12-05 — v2.3

**Added**
- OUTRIDER per-exon analysis and corresponding table in the dashboard.

**Changed**
- Separate coverage plots when multiple mutations are present.
- FASTA sequence displayed as HTML with mutations highlighted in bold.
- Minor refactor of `outrider.R`.

**Fixed**
- Coverage plot hidden when no mutation is defined.
- Incorrect mutation location for sample `HSJ_003_03`.

---

## 2025-11-27 — v2.2

**Added**
- Consensus FASTA sequence generation (`consensus/`).
- This changelog.

**Changed**
- Migrated and updated `RNAseq_version.json` to `scripts/`.
- Standardised script naming to `*.sh` / `*.slurm` / `*.R`.

**Fixed**
- FRASER: removed unused `*closeup.pdf` plot generation.
