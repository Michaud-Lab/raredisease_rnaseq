#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mem=340G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

module load bcftools
module load bedtools
module load gatk
module load r


#
./ASE.sh
