### 2026-02-18 - main

**Changed**
* Split the selector into a gene and a sample specific selector.


### 2026-02-17 - gwFRASER_minorfix

**Changed**
* Add a plot of gwOUTRIDER (and modify the manhattan plot function accordingly)
* Modify the output of outrider to standardise column names
* Subset only +/5Kb of candidate genes to save space.

**fixed**
* readme.md
* how many genes to show on gwFRASER


### 2026-01-29 - Outrider perexon

**Changed**
* Add results of OUTRIDER for all candidate exons.
* Add results of FRASER genome-wide (as a Manhattan plot & data.table)


### 2026-01-28 -feat: a modern UI to Shiny App

**Changed**
* Use bslib to provide a modern UI to Shiny using bootstrap 5.0.  


### 2026-01-22 - error messages

**Changed**
* Produce default plots when no gene/expression is available. 
* Removed the unzipping into a temp folder of data, its just cumbersome for now

**fixed** 
* Several Bug fixes when multiple genes were present.

### 2026-01-20 -  new_genes

**Changed**
* Parameters are now pulled from a single configs.json file.
* Added a selector for the OUTRIDER pvalue threshold
* RNAseq_shiny_v2.4.R: refactor to define an 'i' variable reactively.

**Fixed**
* Minor refactor to accomodate new genes (even for lncRNA)
* Bug fix that did not show all the .bam alignments when multiple genes present.
* Bug fix for HSJ_032 gsub('03$') pattern.

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
