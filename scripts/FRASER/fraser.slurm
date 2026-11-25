#!/bin/bash

#SBATCH --time=03:59:00
#SBATCH --mem=100G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2


module load samtools
module load r

#read configs.json
read workdir rnasplice_bamdir reference_cohorte < <(
  jq -r '.general | [.workdir, .rnasplice_bamdir, .reference_cohorte] | @tsv' ../../data/input/configs.json
)


read fraser_temp_bamdir fraser_bamdir fraser_chr_bamdir  < <(
  jq -r '.FRASER | [.fraser_temp_bamdir, .fraser_bamdir, .fraser_chr_bamdir] | @tsv' ../../data/input/configs.json
)

#Rscript fraser.R $workdir "$rnasplice_bamdir" $fraser_temp_bamdir $fraser_bamdir $reference_cohorte

#fraser genome wide
cpu=4

#all chromosomes are defined and looped over in fraser_gw.R
Rscript fraser_gw.R $fraser_chr_bamdir $cpu "$rnasplice_bamdir" $reference_cohorte

echo 'All done FRASER'

