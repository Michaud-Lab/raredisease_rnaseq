#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load r

#read configs.json
read workdir candidate_genes < <(
  jq -r '.general | [.workdir, .candidate_genes] | @tsv' ../configs.json
)

read fc_pergene fc_perexon < <(
  jq -r '.OUTRIDER | [.fc_pergene, .fc_perexon] | @tsv' ../configs.json
)

#R script launches the bash script to prepare data and continues with analysis.
Rscript outrider.R $workdir $candidate_genes $fc_pergene $fc_perexon


