#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --mem=240G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18

module load bcftools
module load bedtools
module load gatk
module load r


#read configs.json
read workdir gtf rnasplice_bamdir < <(
  jq -r '.general | [.workdir, .gtf, .rnasplice_bamdir] | @tsv' ../../data/input/configs.json
)

read genome2 < <(
  jq -r '.ASE | [.genome2] | @tsv' ../../data/input/configs.json
)

#echo $genome2 $gtf $rnasplice_bamdir $workdir
./ASE.sh $genome2 $gtf "$rnasplice_bamdir" $workdir
