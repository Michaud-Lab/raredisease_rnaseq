# raredisease_rnaseq

RNA-seq analysis pipeline for rare disease candidate gene investigation at CHU Sainte-Justine. Starting from raw sequencing reads, the pipeline produces per-gene and per-exon expression counts, aberrant expression and splicing calls, allele-specific expression results, and an interactive Shiny dashboard for clinical review.

**Maintainer:** sebastien.renaut.hsj@ssss.gouv.qc.ca

---

## Pipeline overview

```
Raw reads (.fastq.gz)
      │
      ▼
nextflow rnasplice         →  BAM files, QC (MultiQC)
      │
      ▼
update_genes.sh            →  Run locally to download candidate_genes_extra.csv and push to serve (Fir)
      │
      ▼
featureCounts.slurm        →  Per-gene and per-exon raw counts, TPM, candidate genes generation.
      │
      ├──▶ outrider.slurm  →  Aberrant gene / exon expression per candidate gene and genome-wide (OUTRIDER)
      │
      ├──▶ fraser.slurm    →  Aberrant splicing per candidate gene and genome-wide (FRASER)
      │
      ├──▶ ase.slurm       →  Allele-specific expression (GATK ASEReadCounter)
      │
      └──▶ consensus.slurm →  Consensus/alternate FASTA sequences per candidate gene
                    │
                    ▼
           cp_and_cleanup.R    →  Consolidate outputs into data/  →  data.zip
                    │
                    ▼
        RNAseq_shiny_v2.5.R   →  Interactive Shiny dashboard
```

> **Cohort size:** A minimum of ~10 samples is required for OUTRIDER and FRASER to produce statistically meaningful outlier calls; >30 is recommended.

> **All-in-one:** Steps 3–8 (featureCounts through cp_and_cleanup.R) can also be run as a single Slurm job via `scripts/run_pipeline.slurm`, instead of submitting each step separately.

---

## Requirements

### HPC modules
```bash
module load nextflow
module load apptainer
module load r        # R >= 4.5
```

### R packages
Each script loads its own dependencies at runtime. Key packages: `OUTRIDER`, `FRASER`, `dplyr`, `tidyr`, `data.table`, `shiny`, `bslib`, `plotly`, `DT`, `igvShiny`, `ggtranscript`, `patchwork`, `biomaRt`.

---

## Input files
```bash
datadir="/project/def-rallard/COMMUN/raredisease_rnaseq/data/input"
scriptsdir="/project/def-rallard/COMMUN/raredisease_rnaseq/scripts"
```

| Files | Description |
|------|-------------|
| `$datadir/nextflow_config.json` | Nextflow configuration |
| `$datadir/nextflow_params.json` | Nextflow rnasplice parameters |
| `$datadir/configs.json` | Paths to all necessary files |

### Input files (nextflow)

| Nextflow | Description | format |
|------|-------------|-------------|
| `$datadir/nextflow_contrast.csv` | Nextflow contrast file (currently unused) | .csv
| `$datadir/nextflow_samples.csv` | .fastq.gz files |.csv
| `sequences/` | Raw paired-end FASTQ files | fastq.gz

### Input files (configs.json)
| configs.json | Description | format |
|------|-------------|-------------|
| `workdir` | Candidate genes and mutations per proband | directory
| `candidate_genes` | Candidate genes and mutations per proband | `candidate_genes.csv`
| `candidate_genes_extra` | More candidate genes | `candidate_genes_extra.csv`
| `rnasplice_bamdir` | aligned files | `*.bam` and `*bam.bai`
| `genome_in` | Reference genome annotation (GTF) | `Homo_sapiens.GRCh38.114.gtf`
| `masterlog` | Anonymized clinical metadata (age, sex, HPO, etc.) | `.xlsx`
| `MANE` | MANE transcript reference | `MANE.GRCh38.v1.5.refseq_genomic.gtf`
| `ens_gene` | Ensembl gene ID mapping file | `ensembl_geneid.tsv`
| `ref_file` | Reference genome (GRCh38) | `Homo_sapiensChr.GRCh38.dna.primary_assembly.fa`
| `ref_annot` | Reference genome annotation (GTF) | `Homo_sapiens.GRCh38.114.gtf`

> **Notes:** Verify all inputs are in place before running.

## Usage

Run each step in order. Steps 3–7 submit Slurm jobs and can run in parallel once Steps 1–2 are complete. Alternatively, `sbatch scripts/run_pipeline.slurm` runs Steps 3–8 sequentially as a single Slurm job.

**Step 1 — Alignment and QC**. 

Aligns raw paired-end FASTQ reads to GRCh38 using the [nf-core rnasplice](https://nf-co.re/rnasplice) Nextflow pipeline, producing sorted/indexed BAM files and a MultiQC report. This is the only step that touches raw sequencing data — every downstream step reads from its BAM output.

```bash
nextflow run rnasplice \
  -params-file data/input/nextflow_params.json \
  -c data/input/nextflow_config.json \
  -resume -bg \
  -w /home/renaut/scratch/nextflow_rnasplice/work
```

> **Notes:** Run this everytime you have a new batch of samples sequenced. You will have to modify `$datadir/nextflow_params.json` to specify a new `outdir` and `input` to list samples in `nextflow_samples.csv`.

**Step 2 — Update candidate genes from Google Sheets**. 

Downloads the clinical team's shared Google Sheet as `candidate_genes_extra.csv` and pushes it to the compute server, so any candidate gene added since the last run is picked up by the rest of the pipeline.

```bash
bash $scriptsdir/update_genes.sh
```

> **Notes:** Run this if/when you have new genes to add to the report. If you have new genes, you will need to re-run Steps 3–9 to update the report. This will be relatively quick since most analyses are skipped if already present.

**Step 3 — Feature counts**

`featureCounts.slurm` reads paths out of `configs.json` and calls `featureCounts.sh`, which:
1. Builds an exon-tagged GTF from `genome_in` (once; skipped if it already exists), and a MANE transcript lookup (`MANE.tsv`) from the MANE reference.
2. Runs `subread featureCounts` twice on every BAM in `rnasplice_bamdir`: once per-gene against the plain GTF, once per-exon against the exon-tagged GTF — both `-p -B -C` (paired, both ends aligned, no chimeras) with multi-overlaps counted (`-O`; fractional for the per-gene run). Both are cached: skipped if `feature_counts_pergene.txt` already exists in `workdir/featureCounts/`.
3. Hands off to `featureCounts.R`, which picks the MANE-selected transcript per gene (falling back to the highest-coverage transcript for genes absent from MANE, e.g. lncRNAs), computes normalized TPM at the gene and exon level (excluding haemoglobin genes `HBB`/`HBA1`/`HBA2`/`HBD` from the exon-level TPM), restricts output to candidate genes/probands, and (re)builds the master candidate list `candidate_genes_ALL.csv`.

   `candidate_genes_ALL.csv` is assembled from three sources, deduplicated on (`geneID`, `ensembl`, `proband`), then merged with the anonymized clinical metadata (`masterlog`):
   - **Manually curated** — every row of `candidate_genes.csv`.
   - **Extra candidates** — every row of `candidate_genes_extra.csv` (the Google Sheet pulled in Step 2).
   - **Automated, from the *previous* run's genome-wide results** (`gwFRASER.tsv`, `gw_genes_OUTRIDER.tsv`, `gwASE.tsv` — skipped on a first run, before those files exist), via `candidate_genes_automated()` in [rnaseq_helper_functions.R](scripts/rnaseq_helper_functions.R):
     - *FRASER:* per sample, splicing outliers with `padjust < 1e-4`, top 5 by p-adjust.
     - *OUTRIDER:* per sample, expression outliers with `pValue < 1e-6`, top 2 by p-value.
     - *ASE:* per sample+gene, allele-specific expression with `pvalue < 1e-49` at ≥2 heterozygous markers.
     - Each of the three excludes haemoglobin/HLA/`SELPLG` genes (`HBA*`, `HBB`, `HLA*`, `HBG`, `HBD`, `HBQ`, `HBE`, `HBZ`, `HBM`, `SELPLG`) and is annotated back to genomic coordinates via `gene_annotation()`; genes that can't be resolved to a locus are dropped.
     - *HPO-driven* — for each proband, `gene_prioritization()` cross-references their HPO terms (from `masterlog`, against `genes_to_phenotype.txt`) with that same previous run's genome-wide FRASER/OUTRIDER results, keeping hits with adjusted gene-score p-value `< 0.01` (top 2 per proband).

   This file drives every subsequent step and the dashboard, and because the automated/HPO sources read the *prior* run's genome-wide outputs, a gene surfaced this way only appears in the report starting on the run after it was first flagged.

```bash
sbatch $scriptsdir/featureCounts/featureCounts.slurm
```

**Step 4 — Aberrant expression (OUTRIDER)**

Runs [OUTRIDER](https://github.com/gagneurlab/OUTRIDER) on the gene- and exon-level count matrices from Step 3 to flag statistically aberrant expression, both genome-wide and specifically for each candidate gene/proband pair. Haemoglobin genes are excluded before fitting, since their extreme expression would otherwise dominate the model.

```bash
sbatch $scriptsdir/OUTRIDER/outrider.slurm
```

**Step 5 — Aberrant splicing (FRASER)**

For each candidate gene/proband pair, subsets the BAM to a ±5 kb window around the gene and runs [FRASER](https://github.com/gagneurlab/FRASER) to detect aberrant splicing (intron retention, exon skipping) against the reference cohort. Also generates the per-gene sashimi coverage plots shown in the dashboard's Gene model/FRASER tabs. Genes with no splice junctions in a proband are skipped rather than failing the run.

```bash
sbatch $scriptsdir/FRASER/fraser.slurm
```

**Step 6 — Allele-specific expression**

Runs GATK `ASEReadCounter` on each proband's RNA-seq BAM at heterozygous SNV positions from their WGS genotypes, then tests for allele-specific expression (binomial test) genome-wide, flagging known imprinted and X-linked genes separately.

```bash
sbatch $scriptsdir/ASE/ase.slurm
```

**Step 7 — Consensus sequences**

Builds a consensus/alternate FASTA sequence for each candidate gene by incorporating the proband's variants into the reference sequence, for side-by-side comparison in the dashboard's "fasta" tab.

```bash
sbatch $scriptsdir/consensus/consensus.slurm
```

**Step 8 — Consolidate outputs**

Copies every step's outputs (FRASER, OUTRIDER, ASE, consensus, featureCounts tables, sashimi plots, MultiQC report) into a single `data/` directory. This is also when `candidate_genes_ALL.csv` gets its final FRASER/OUTRIDER/ASE/HPO annotation. Finishes by zipping `data/` for transfer.

```bash
Rscript $scriptsdir/cp_and_cleanup.R
```

**Step 9a — Launch the Shiny dashboard** *(requires R ≥ 4.5, run in RStudio or locally)*

Launches the interactive dashboard (`RNAseq_shiny_v2.5.R`) for browsing per-proband candidate genes: expression, splicing, ASE results, gene prioritization, and IGV alignments — reading everything from the `data/` directory produced by Step 8.

```bash
Rscript RNAseq_shiny_v2.5.R
```
**Step 9b — Launch the Shiny dashboard (synthetic dataset)**
To launch with the minimal synthetic dataset (`data_minimal/`) for local development and testing without real patient data, modify `use_data_minimal = TRUE`

**Step 9c — Launch the Shiny dashboard with password protection**
To require a username and password before accessing the dashboard, pass `use_password = TRUE`. A credentials database must exist at `data/credentials.sqlite`; create it once by running.  
```

> Long-read FRASER analysis (`lr_quant/lr_quant.slurm`) is available but not part of the standard run.

---

## Shiny dashboard tabs

| Tab | Content |
|-----|---------|
| Proband Selection | Proband/gene selector, clinical summary |
| Expression | Per-exon TPM table (MANE transcript) |
| Plot | Cohort and family expression plots (plotly) |
| OUTRIDER | Genome-wide and candidate gene aberrant expression |
| Gene model | FRASER splicing significance + coverage plot |
| IGV | In-browser BAM alignment viewer |
| FRASER | Genome-wide splicing Manhattan plot and table |
| Gene Prioritization | HPO × OUTRIDER × FRASER gene ranking |
| ASE | Allele-specific expression Manhattan plot and tables |
| fasta | Reference/alternate consensus sequences |
| Search expression | Cohort-wide gene expression lookup |
| MultiQC | Alignment quality control report |

---

## Testing

Unit tests are located in `scripts/tests/testthat/test_helper_functions.R` and cover `gwFRASER_table()`, `manhattan_plot()`, `gene_prioritization()`, `plot_expression_cohort()`, and `plot_expression_family()`.

Run from the project root:
```r
testthat::test_file("scripts/tests/testthat/test_helper_functions.R")
```

---

## Citations

- **OUTRIDER:** Brechtmann F, et al. *OUTRIDER: a statistical method for detecting aberrantly expressed genes in RNA sequencing data.* Am J Hum Genet. 2018;103(6):907–917.
- **FRASER:** Mertes C, et al. *Detection of aberrant splicing events in RNA-seq data using FRASER.* Nat Commun. 2021;12(1):529.
- **featureCounts:** Liao Y, Smyth GK, Shi W. *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.* Bioinformatics. 2014;30(7):923–930.
- **Nextflow rnasplice:** [zenodo.org/records/15194198](https://zenodo.org/records/15194198)
- **PacBio Iso-Seq:** [isoseq.how/getting-started.html](https://isoseq.how/getting-started.html)
