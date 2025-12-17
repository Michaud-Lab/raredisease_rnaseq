#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load r

workdir="$HOME/scratch/raredisease_rnaseq/"

#R script launches the bash script to prepare data and continues with analysis.
Rscript outrider.R $workdir
Rscript outrider_perexon.R $workdir

