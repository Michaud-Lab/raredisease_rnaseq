#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4


module load bedtools
module load gatk
module load r



./ASE.sh
