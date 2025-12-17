#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1


module load bedtools
module load bcftools
module load samtools
module load r


#Rscript consensus.R 1 $HOME/scratch/nextflow_rnasplice/

workdir="$HOME/scratch/raredisease_rnaseq/"

#R script launches the bash script to annotated variant file and continues with consensus.R to produce the reference sequence.
for i in {1..32}; do
  Rscript consensus.R $i $workdir
  echo $i
done







