### Code owner      
sebastien.renaut.hsj@ssss.gouv.qc.ca    

### Summary        
 Typically, one would run analyses one step at a time, starting with `nextflow rnasplice` to produce the QC, trimming and aligment (`.bam files`). Then, run `featureCounts` to quantify gene expression per gene and per exon.  Then, run `Outrider` (aberrant gene expression), `fraser` (aberrant splicing pattterns) and `consensus` (consensus sequences). One these analyses are performed, run `cp_and_cleanup.R` to gather all the necessary results into a single zipped `data` folder. With this zipped `data` folder, we can then produce a Shiny web App with `RNAseq_shiny_v2.4.R` script. Note that:    

* You need a resonnably sized cohort to produce results that are statistically sound. Ten is a probably a bare minimum for Outrider and Fraser.      
* Paths are either defined as relative or in a few instances at the beginning of a script.       
* Typically, each analyse is run from a subdirectory by calling a `.slurm` file. The parameters in these files may need to be adjusted according to your HPC set-up. The `.slurm` file will call a  `.R` and `.sh` script. Each analysis requires a number of R libraries to be installed, defined at the beginning of the script.            
* There is currently no unit test, or example dataset.       

### Setup         
`datadir='${HOME}/scratch/raredisease_rnaseq/data'`     
`scriptsdir='${HOME}/scratch/raredisease_rnaseq/scripts'`       

`module load nextflow`      
`module load apptainer`            
`module load r`         

### Check the required files  
`ls $datadir/candidate_genes_3.txt`   # list of candidate genes & mutations      
`ls $datadir/*.xlsx`   # An excel file that contains more info about each sample.       
`ls $datadir/config.json`   # Nextflow configuration file      
`ls $datadir/params.json`   # rnasplice parameters     
`ls $datadir/contrast.csv`  # Contrast file for nextflow. Currently not meaningfull     
`ls sequences/*`  # The raw sequencing files (.fastq paired-end data)       
`ls reference/*`  # The reference genome and annotations        
`ls $scriptsdir/*`   # The required scripts     

### Usage
* Nextflow pipeline & post-processing:    
`nextflow run rnasplice -params-file data/params.json -c data/config.json -resume -bg`      
`sbatch $scriptsdir/featureCounts/featureCounts.slurm`      
`sbatch $scriptsdir/OUTRIDER/outrider.slurm`        
`sbatch $scriptsdir/FRASER/fraser.slurm`        
`sbatch $scriptsdir/consensus/consensus.slurm`      
`Rscript $scriptsdir/cp_and_cleanup.R`      
    
* Shiny app (Run in Rstudio, or Rscript locally, requires R>=4.5)       
`Rscript RNAseq_shiny_v2.4.R`   

### Citations       
* OUTRIDER: Brechtmann, Felix, et al. "OUTRIDER: a statistical method for detecting aberrantly expressed genes in RNA sequencing data." The American Journal of Human Genetics 103.6 (2018): 907-917.       
* Mertes, Christian, et al. "Detection of aberrant splicing events in RNA-seq data using FRASER." Nature communications 12.1 (2021): 529.       
* Liao, Yang, Gordon K. Smyth, and Wei Shi. "featureCounts: an efficient general purpose program for assigning sequence reads to genomic features." Bioinformatics 30.7 (2014): 923-930.        
* [Nextflow rnasplice:](https://zenodo.org/records/15194198)