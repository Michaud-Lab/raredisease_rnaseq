#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

#
module load StdEnv/2023
module load smrtlink/25.3
module load samtools
module load r


#PacBio alignment
if [ ! -f "/home/renaut/scratch/raredisease_rnaseq/lr_quant/bc01/mapped.bam" ]; then
  ./lr_quant.sh 64 
else
  echo "~~~ Not running lr_quant.sh because I already have the bam files~~~"
fi

#read configs.json
read workdir candidate_genes_LR workdir_LR< <(
  jq -r '.general | [.workdir, .candidate_genes_LR, .workdir_LR] | @tsv' ../configs.json
  )

temp=$(wc -l < $candidate_genes_LR)
count=$((temp - 1)) 


#
for i in $(seq 1 $count); do
  Rscript lr_quant.R $workdir"FRASER/" $workdir_LR $i
  echo $i
done


