#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64


module load samtools
module load r

#read configs.json
read workdir candidate_genes < <(
  jq -r '.general | [.workdir, .candidate_genes] | @tsv' ../configs.json
)

temp=$(wc -l < $candidate_genes)
count=$((temp - 1))

#R script launches the bash script to prepare data (i.e. subset the big bam files) and continues with FRASER analysis and plots.
#Can add new genes specified in candidate_genes_3.txt to the for loop below:
for i in $(seq 1 $count); do
  Rscript fraser.R $i $workdir
  echo $i
done


#fraser genome wide
cpu=64
Rscript fraser_gw.R /home/renaut/scratch/raredisease_rnaseq/FRASER/bams_chr_subset/ 16 $cpu


