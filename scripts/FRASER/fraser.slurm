#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1


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
