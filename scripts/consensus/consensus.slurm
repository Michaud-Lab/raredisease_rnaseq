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

#params
workdir="$HOME/scratch/raredisease_rnaseq/"
candidate_genes=$workdir"/data/candidate_genes_3.txt"
fc_exons=$workdir"featureCounts/fc_exons_tpm.tsv"
ref_file=$workdir"../reference/Homo_sapiens/Homo_sapiensChr.GRCh38.dna.primary_assembly.fa"
ref_annot=$workdir"../reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"

temp=$(wc -l < $candidate_genes)
count=$((temp - 1))

#R script launches the bash script to annotated variant file and continues with consensus.R to produce the reference sequence.
#Can add new genes here too.
for i in $(seq 9 44); do
  Rscript consensus.R $i $workdir $candidate_genes $fc_exons $ref_file $ref_annot
  echo $i
done







