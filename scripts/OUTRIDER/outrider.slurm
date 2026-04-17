#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18

module load r

#read configs.json
read workdir candidate_genes candidate_genes_LR < <(
  jq -r '.general | [.workdir, .candidate_genes, .candidate_genes_LR] | @tsv' ../configs.json
)

read fc_pergene fc_perexon < <(
  jq -r '.OUTRIDER | [.fc_pergene, .fc_perexon] | @tsv' ../configs.json
)

#R script launches the bash script to prepare data and continues with analysis.
Rscript outrider.R $workdir $candidate_genes $fc_pergene $fc_perexon $candidate_genes_LR

#echo  $workdir $candidate_genes $fc_pergene $fc_perexon $candidate_genes_LR
