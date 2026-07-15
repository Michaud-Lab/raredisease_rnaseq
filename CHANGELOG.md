# Changelog

All notable changes to this project are documented here.

---

## 2026-07-15 — refactor-fraser-consensus

**Changed**
- `scripts/FRASER/fraser.R`: per-candidate logic wrapped into a `fraser_pipeline(candidates, i)` function, called from a loop over all candidates; added a progress print statement before the loop.
- `scripts/consensus/consensus.R`: per-candidate logic wrapped into a `consensus_pipeline(candidates, i)` function, called from a loop over all candidates. The script now reads `candidate_genes.csv`, `candidate_genes_extra.csv`, and `candidate_genes_automated.csv` itself (matching `fraser.R`) instead of receiving `candidates`/`candidate_genes_extra` paths and a proband index `i` as command-line arguments.
- `scripts/consensus/consensus.slurm`: updated to match the new `consensus.R` argument signature (`workdir`, `fc_exons`, `genome`, `gtf`); removed the per-candidate `Rscript` loop and `candidate_genes`/`candidate_genes_extra` line-count logic, since looping is now handled inside `consensus.R`.
- `scripts/rnaseq_helper_functions.R`: `candidate_genes_automated()` haemoglobin gene exclusion regex broadened from `HBA|HBB|HLA` to anchored prefixes covering `HBA, HBB, HLA, HBG, HBD, HBQ, HBE, HBZ`.
- `scripts/Shiny/rnaseq_shinyhelper_functions.R`: `plot_hb_fraction()` legend moved to a horizontal orientation anchored to the top-right corner of the plot (was default right-side vertical legend).
- `RNAseq_shiny_v2.5.R`: bolded the "Information", "Software Version", and MANE reference card headers.

---

## 2026-07-13 — centralized-package-loading

**Changed**
- Pipeline and analysis scripts now load packages through the shared `load_install_library()` helper instead of individual `library()` calls, so missing packages are installed on demand. Consecutive `library()` calls were collapsed into single grouped `load_install_library(c(...))` calls, and each standalone script gained a `source()` of the helper file:
  - `scripts/ASE/ASE.R`
  - `scripts/OUTRIDER/outrider.R`
  - `scripts/featureCounts/featureCounts.R`
  - `scripts/featureCounts/rnaseq_helper_functions.R` (inside `gene_annotation()`; relies on the helper being sourced by its caller)
  - `scripts/lr_quant/lr_quant.R`
  - `scripts/FRASER/fraser.R`
  - `scripts/FRASER/fraser_gw.R`
  - `scripts/consensus/consensus.R`
- The recursive `library()` calls inside `load_install_library()` itself, plus `scripts/Shiny/create_credentials_db.R` and the `testthat` suite, were intentionally left using `library()`.

---

## 2026-07-08 — password_protection

**Added**
- `scripts/Shiny/create_credentials_db.R`: one-time setup script that creates `data/credentials.sqlite` (bcrypt-hashed passwords via `shinymanager::create_db()`). Passwords are replaced with `"CHANGE_ME"` placeholders — fill them in locally before running. The script aborts with a clear error if any placeholder is left unchanged. The database file is gitignored.
- `RNAseq_shiny_v2.5.R`: new `--use_password` command-line flag; when passed, the UI is wrapped with `shinymanager::secure_app()` and the server calls `shinymanager::secure_server()` to enforce login. Without the flag the app launches as before with no authentication. Flags can be combined: `Rscript RNAseq_shiny_v2.5.R --data_minimal --use_password`.
- `scripts/Shiny/global.R`: added `shinymanager` to the package list.
- `data/credentials.sqlite` added to `.gitignore`.

---

## 2026-07-06 — data-minimal-and-bugfixes

**Added**
- `data_minimal/`: minimal synthetic dataset (2 probands: RNA_101_03_PAX / MCM5, RNA_106_03_PAX / TTI1; 5 cohort samples; 25 chromosomes represented in OUTRIDER and FRASER results) for fast local development and testing of the Shiny app.

**Changed**
- `RNAseq_shiny_v2.5.R`: added `--data_minimal` command-line argument; when passed, the app loads from `data_minimal/` instead of `data/`.
- `scripts/Shiny/global.R`: `params$datadir` now resolves to `data_minimal/` when `use_data_minimal` is `TRUE`; added `header = TRUE` to all four `read.table()` calls (`fc_exons_raw.tsv`, `fc_exons_tpm.tsv`, `fc_genes_tpm.tsv`, `fc_genes_raw_ALL.tsv`) to support files without an implicit row-name column.
- `scripts/Shiny/rnaseq_shinyhelper_functions.R`: `genemodel_plot()` — `candidate$position` is now coerced with `as.character()` before `strsplit()` and resulting NAs are filtered out, preventing a crash when position is empty; `manhattan_plot()` — returns `plot(0, main = 'no data available')` early when the input data frame is NULL, empty, or contains no rows for the requested sample.

---

## 2026-07-03 — outrider-hb-filter

**Changed**
- `scripts/OUTRIDER/outrider.R`: haemoglobin genes (`HBB`, `HBA-1`, `HBA-2`) are now excluded from `genes_counts` before both `OutriderDataSet()` calls (gene-level and exon-level). Gene-level filtering uses `%in%` on rownames; exon-level filtering uses `grepl('^(HBB|HBA-1|HBA-2)_', ...)` to match the `geneID_ensemblID_transcriptID_exonID` rowname format.

---

## 2026-07-02 — candidates-extra-column

**Changed**
- `scripts/featureCounts/featureCounts.R`, `scripts/FRASER/fraser.R`, `scripts/OUTRIDER/outrider.R`, `scripts/consensus/consensus.R`, `scripts/cp_and_cleanup.R`, `scripts/Shiny/global.R`: added `candidates_extra = candidates_extra[, colnames(candidates)]` before every `rbind(candidates, candidates_extra)` to drop any extra columns in `candidate_genes_extra.csv` before merging. This prevents column-mismatch errors when the extra CSV has more columns than `candidate_genes.csv`.

---

## 2026-06-23 — new-plots

**Added**
- `scripts/Shiny/rnaseq_shinyhelper_functions.R`: new `plot_hb_fraction()` function — stacked bar chart showing the fraction of total reads attributed to a user-specified set of genes (default: HBA1, HBA2, HBB, HBG1, HBG2) per sample. Hover tooltip shows sample name, gene ID, percentage of total reads, and raw read count.
- `scripts/tests/testthat/test_helper_functions.R`: 7 unit tests for `plot_hb_fraction()` covering return type, bar chart trace type, fraction range (0–100%), correct trace count, single-gene edge case, custom gene list, and exclusion of the `Length` column.

**Changed**
- `RNAseq_shiny_v2.5.R`: replaced the "Search expression" datatable with a per-gene fraction-of-total-reads bar plot using `plot_hb_fraction()`; gene selection uses `selectizeInput` with server-side loading to handle the full genome-wide gene list efficiently. A fixed haemoglobin multi-gene stacked bar plot is shown below as a second card in the same tab.

---

## 2026-06-18 — candidates-extra

**Changed**
- `scripts/featureCounts/featureCounts.R`, `scripts/FRASER/fraser.R`, `scripts/OUTRIDER/outrider.R`, `scripts/consensus/consensus.R`, `scripts/cp_and_cleanup.R`, `scripts/Shiny/global.R`: all scripts now merge `candidates_extra` from a local CSV (path stored in `configs.json` under `general.candidate_genes_extra`) into `candidates` after loading `candidate_genes.csv`. The CSV is downloaded from Google Sheets by running `update_genes.sh` before the pipeline. The URL is passed as a command-line argument in slurm scripts; `cp_and_cleanup.R` and `global.R` read it directly from `configs.json`.
- `scripts/FRASER/fraser.slurm`, `scripts/OUTRIDER/outrider.slurm`, `scripts/consensus/consensus.slurm`: updated to read `candidate_genes_extra` from `configs.json` and pass it as a new trailing argument to their respective R scripts.
- `README.md`: added `update_genes.sh` as Step 2 in the pipeline overview diagram and usage section; renumbered subsequent steps accordingly.

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
- `RNAseq_shiny_v2.5.R`: added a `candidates_table` to the Proband Selection tab listing all probands with their candidate gene IDs (comma-separated per proband); table is placed to the right of the `selectInput` using a `layout_columns(col_widths = c(3, 9))` layout; displays 10 rows per page.

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
