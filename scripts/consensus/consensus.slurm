#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1


module load bedtools
module load bcftools
module load samtools
module load r

#read configs.json
read workdir gtf candidate_genes candidate_genes_extra < <(
  jq -r '.general | [.workdir, .gtf, .candidate_genes, .candidate_genes_extra] | @tsv' ../../data/input/configs.json
)

read fc_exons genome < <(
  jq -r '.consensus | [.fc_exons, .genome] | @tsv' ../../data/input/configs.json
)

#
temp=$(wc -l < $candidate_genes)
temp_extra=$(wc -l < $candidate_genes_extra)
count=$((temp - 1 + temp_extra - 2))


#R script launches the bash script to annotated variant file and continues with consensus.R to produce the reference sequence.
for i in $(seq 1 $count); do
  Rscript consensus.R $i $workdir $candidate_genes $fc_exons $genome $gtf $candidate_genes_extra
done


echo 'Done Rscript Rscript consensus.R' $i $workdir $candidate_genes $fc_exons $genome $gtf $candidate_genes_extra



