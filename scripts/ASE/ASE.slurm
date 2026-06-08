#!/bin/bash

#SBATCH --time=40:00:00
#SBATCH --mem=34G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

module load bcftools
module load bedtools
module load gatk
module load r


./ASE.sh
#./ASE_commonvariants.sh

