#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1


module load samtools
module load r

workdir="$HOME/scratch/raredisease_rnaseq/"

#R script launches the bash script to prepare data (i.e. subset the big bam files) and continues with FRASER analysis and plots.
for i in {1..32}; do
  Rscript fraser.R $i $workdir
  echo $i
done
