### 2026-01-09 
**Changed** 

* Minor refactor to accomodate new samples
* Moved parameters as much as possible to the .slurm scripts.

### 2026-01-06  - feat: add slider for coverage

**Added**
* Add a slider widget to the coverage plots

### 2025-12-18 - v2.4: update_readme
		
**Changed**		
* Update to the README.md, LICENCE, VERSION.json		

### 2025-12-17 - v2.4

**Changed**
* Major refactor to fraser.r
   * Moved the coverage plotting to a new rnaseq_helper_functions.R , instead simply create files (.csvs + .rda) to generate plots later as part of RNAseq_shiny_v2.4.R
      * Added a option to plot a +/-1kb or +/-5kb window of coverage and refactoring of syntax. 
   * Moved the gene model preparation to a new rnaseq_helper_functions.R function, then run it as part of featurecounts.R
* Reference transcript according to: Matched Annotation from NCBI and EMBL-EBI (MANE)

**Fixed**
* Minor refactor to simplify Shiny App

### 2025-12-05 - v2.3

**Added**
*  New analysis (outrider per exon) and new table for dataviz in the Shiny App. 

**Changed**
*   Show  a seperate coverage plot in cased there are several mutations.
*   Show fasta sequence as HTML with mutations in bold.
*   Minor refactor to outrider.r to streamline the code.

**Fixed**
*   Remove coverage plot when no mutation is present.
*   Wrong mutation location for HSJ_003_03 fixed.

### 2025-11-27 - v2.2

**Added**
*   New feature to generate a consensus fasta sequence (`consensus/`).
*   Created this changelog.

**Changed**
*   migrated and updated RNAseq_version.json to `scripts/`.
*   standardised script naming to `*.sh` / `*.slurm` / `*.R`.

**Fixed**
*   FRASER:Don't generate the *closeup.pdf plots because we don't use it anymore.
