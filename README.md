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
featureCounts.slurm        →  Per-gene and per-exon raw counts, TPM
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

Run each step in order. Steps 3–7 submit Slurm jobs and can run in parallel once Steps 1–2 are complete.

**Step 1 — Alignment and QC**. 
```bash
nextflow run rnasplice \
  -params-file data/input/nextflow_params.json \
  -c data/input/nextflow_config.json \
  -resume -bg \
  -w /home/renaut/scratch/nextflow_rnasplice/work
```

> **Notes:** Run this everytime you have a new batch of samples sequenced. You will have to modify `$datadir/nextflow_params.json` to specify a new `outdir` and `input` to list samples in `nextflow_samples.csv`.

**Step 2 — Update candidate genes from Google Sheets**. 
```bash
bash $scriptsdir/update_genes.sh
```

> **Notes:** Run this if/when you have new genes to add to the report. If you have new genes, you will need to re-run Steps 3–9 to update the report. This will be relatively quick since most analyses are skipped if already present.

**Step 3 — Feature counts**
```bash
sbatch $scriptsdir/featureCounts/featureCounts.slurm
```

**Step 4 — Aberrant expression (OUTRIDER)**
```bash
sbatch $scriptsdir/OUTRIDER/outrider.slurm
```

**Step 5 — Aberrant splicing (FRASER)**
```bash
sbatch $scriptsdir/FRASER/fraser.slurm
```

**Step 6 — Allele-specific expression**
```bash
sbatch $scriptsdir/ASE/ase.slurm
```

**Step 7 — Consensus sequences**
```bash
sbatch $scriptsdir/consensus/consensus.slurm
```

**Step 8 — Consolidate outputs**
```bash
Rscript $scriptsdir/cp_and_cleanup.R
```

**Step 9a — Launch the Shiny dashboard** *(requires R ≥ 4.5, run in RStudio or locally)*
```bash
Rscript RNAseq_shiny_v2.5.R
```
**Step 9b — Launch the Shiny dashboard (synthetic dataset)**
To launch with the minimal synthetic dataset (`data_minimal/`) for local development and testing without real patient data:
```bash
Rscript RNAseq_shiny_v2.5.R --data_minimal
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
