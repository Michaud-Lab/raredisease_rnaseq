
### 0. workdir path  
workdir='/home/renaut/scratch/nextflow_rnasplice'
scriptsdir='/home/renaut/scratch/scripts'

### 1. Setup  
module load nextflow
module load apptainer
module load r 

### 2. mkdir to store results  
mkdir -p $workdir

### 3. Check the required files  
ls $workdir/data/candidate_genes_3.txt   # list of candidate genes to perform analyses on. Typically one / proband. Also contains the candidate mutation(s).
ls $workdir/data/config.json   # Nextflow configuration file
ls $workdir/data/params.json   # rnasplice parameters 
ls $workdir/data/contrast.csv  # Contrast file for rnasplice, altough this is not actually used/meaningfull for us.
ls $workdir/data/RNAseq_version.json #versions 
ls sequences/*  #The raw sequencing files (.fastq paired-end data)

### 4. Check the required scripts  
ls $scriptsdir/*

### 5. Launch the Nextflow pipeline and then the post-processing
nextflow run rnasplice -params-file data/params.json -c data/config.json -resume -bg 
sbatch $scriptsdir/featureCounts/featureCounts.slurm
sbatch $scriptsdir/OUTRIDER/outrider.slurm
sbatch $scriptsdir/FRASER/fraser.slurm
sbatch $scriptsdir/consensus/consensus.slurm
Rscript $scriptsdir/cp_and_cleanup2.R

### 7. Prepare the Shiny app (Run in Rstudio, or Rscript locally, requires R>=4.5)  
Rscript RNAseq_shiny_v2.2.R
