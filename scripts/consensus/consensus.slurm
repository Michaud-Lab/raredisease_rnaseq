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
read workdir gtf < <(
  jq -r '.general | [.workdir, .gtf] | @tsv' ../../data/input/configs.json
)

read fc_exons genome < <(
  jq -r '.consensus | [.fc_exons, .genome] | @tsv' ../../data/input/configs.json
)


Rscript consensus.R $workdir $fc_exons $genome $gtf

echo 'Done Rscript consensus.R' $workdir $fc_exons $genome $gtf



