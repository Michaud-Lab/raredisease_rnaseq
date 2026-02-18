#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1


module load bedtools
module load bcftools
module load samtools
module load r

#read configs.json
read workdir candidate_genes < <(
  jq -r '.general | [.workdir, .candidate_genes] | @tsv' ../configs.json
)

read fc_exons ref_file ref_annot < <(
  jq -r '.consensus | [.fc_exons, .ref_file, .ref_annot] | @tsv' ../configs.json
)

#
temp=$(wc -l < $candidate_genes)
count=$((temp - 1))

#
fc_exons=$HOME/$fc_exons
ref_file=$HOME/$ref_file
ref_annot=$HOME/$ref_annot

#R script launches the bash script to annotated variant file and continues with consensus.R to produce the reference sequence.
#Can add new genes here too.
for i in $(seq 85 $count); do
  Rscript consensus.R $i $workdir $candidate_genes $fc_exons $ref_file $ref_annot
done



#echo $i $workdir $candidate_genes $fc_exons $ref_file $ref_annot



