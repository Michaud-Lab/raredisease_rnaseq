#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module load subread/2.0.6
module load r


./featureCounts.sh

echo 'Done featureCounts.sh'

