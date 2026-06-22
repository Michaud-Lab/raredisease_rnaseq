#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module load samtools
module load r

#read configs.json
read workdir rnasplice_bamdir candidate_genes candidate_genes_extra < <(
  jq -r '.general | [.workdir, .rnasplice_bamdir, .candidate_genes, .candidate_genes_extra] | @tsv' ../../data/input/configs.json
)


read fraser_temp_bamdir fraser_bamdir fraser_chr_bamdir  < <(
  jq -r '.FRASER | [.fraser_temp_bamdir, .fraser_bamdir, .fraser_chr_bamdir] | @tsv' ../../data/input/configs.json
)

temp=$(wc -l < $candidate_genes)
temp_extra=$(wc -l < $candidate_genes_extra)
count=$((temp - 1 + temp_extra - 2))

#R script launches the bash script to prepare data (i.e. subset the big bam files) and continues with FRASER analysis and plots.
#Can add new genes specified in candidate_genes_3.txt to the for loop below:
for i in $(seq 1 $count); do
  Rscript fraser.R $i $workdir $rnasplice_bamdir $fraser_temp_bamdir $fraser_bamdir $candidate_genes_extra
done


#fraser genome wide
#cpu=32
#for i in $(seq 1 22) X Y MT; do
# Rscript fraser_gw.R $fraser_chr_bamdir $i $cpu $rnasplice_bamdir
# echo 'Done chr:'$i
#done




echo 'All done FRASER'
