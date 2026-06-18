# raredisease_rnaseq

RNA-seq analysis pipeline for rare disease candidate gene investigation at CHU Sainte-Justine. Starting from raw sequencing reads, the pipeline produces per-gene and per-exon expression counts, aberrant expression and splicing calls, allele-specific expression results, and an interactive Shiny dashboard for clinical review.

**Maintainer:** sebastien.renaut.hsj@ssss.gouv.qc.ca

---

## Pipeline overview

```
Raw reads (.fastq)
      │
      ▼
nextflow rnasplice          →  BAM files, QC (MultiQC)
      │
      ▼
update_genes.sh             →  Download candidate_genes_extra.csv from Google Sheets
      │
      ▼
featureCounts.slurm         →  Per-gene and per-exon raw counts, TPM
      │
      ├──▶ outrider.slurm   →  Aberrant gene / exon expression (OUTRIDER)
      │
      ├──▶ fraser.slurm     →  Aberrant splicing per candidate gene (FRASER)
      │
      ├──▶ ase.slurm        →  Allele-specific expression (GATK ASEReadCounter)
      │
      └──▶ consensus.slurm  →  Consensus/alternate FASTA sequences
                    │
                    ▼
             cp_and_cleanup.R   →  Consolidate outputs → data.zip
                    │
                    ▼
          RNAseq_shiny_v2.5.R  →  Interactive Shiny dashboard
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

Install `ggtranscript` from GitHub if not available via CRAN:
```r
remotes::install_github("dzhang32/ggtranscript")
```

---

## Input files

| File | Description |
|------|-------------|
| `data/input/candidate_genes.csv` | Candidate genes and mutations per proband |
| `data/input/CHUSJ_Master_Linking_Log_modif.xlsx` | Clinical metadata (age, sex, HPO terms, etc.) |
| `data/input/nextflow_config.json` | Nextflow configuration |
| `data/input/nextflow_params.json` | Nextflow rnasplice parameters |
| `data/input/nextflow_contrast.csv` | Nextflow contrast file (currently unused) |
| `sequences/` | Raw paired-end FASTQ files |
| `reference/` | Reference genome and annotation files (GRCh38) |

Verify all inputs are in place before running:
```bash
datadir="${HOME}/scratch/raredisease_rnaseq/data"
scriptsdir="${HOME}/scratch/raredisease_rnaseq/scripts"

ls $datadir/input/candidate_genes.csv
ls $datadir/input/CHUSJ_Master_Linking_Log_modif.xlsx
ls $datadir/input/nextflow_config.json
ls $datadir/input/nextflow_params.json
ls sequences/
ls reference/
```

---

## Usage

Run each step in order. Steps 2–6 submit Slurm jobs and can run in parallel once Step 1 is complete.

**Step 1 — Alignment and QC**
```bash
nextflow run rnasplice \
  -params-file data/input/nextflow_params.json \
  -c data/input/nextflow_config.json \
  -resume -bg \
  -w /home/renaut/scratch/nextflow_rnasplice/work
```

**Step 2 — Update candidate genes from Google Sheets**
```bash
bash $scriptsdir/update_genes.sh
```

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

**Step 9 — Launch the Shiny dashboard** *(requires R ≥ 4.5, run in RStudio or locally)*
```bash
Rscript RNAseq_shiny_v2.5.R
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
