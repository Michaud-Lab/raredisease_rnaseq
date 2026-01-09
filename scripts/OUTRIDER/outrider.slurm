#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load r

#params
workdir="$HOME/scratch/raredisease_rnaseq/"
candidate_genes=$workdir"/data/candidate_genes_3.txt"
fc_pergene=$workdir"featureCounts/feature_counts_pergene.txt"
fc_perexon=$workdir"featureCounts/fc_exons_raw.tsv_ALL"

#R script launches the bash script to prepare data and continues with analysis.
Rscript outrider.R $workdir $candidate_genes $fc_pergene $fc_perexon


