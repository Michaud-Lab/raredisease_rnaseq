#!/bin/bash

#SBATCH --time=03:59:00
#SBATCH --mem=100G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2


module load samtools
module load r

#read configs.json
read workdir rnasplice_bamdir < <(
  jq -r '.general | [.workdir, .rnasplice_bamdir] | @tsv' ../../data/input/configs.json
)


read fraser_temp_bamdir fraser_bamdir fraser_chr_bamdir  < <(
  jq -r '.FRASER | [.fraser_temp_bamdir, .fraser_bamdir, .fraser_chr_bamdir] | @tsv' ../../data/input/configs.json
)

Rscript fraser.R $workdir "$rnasplice_bamdir" $fraser_temp_bamdir $fraser_bamdir

#fraser genome wide
#cpu=4

#for i in $(seq 19 22) X Y MT; do
#for i in $(seq 4 9); do
#i=16
#Rscript fraser_gw.R $fraser_chr_bamdir $i $cpu "$rnasplice_bamdir"
  #echo  $fraser_chr_bamdir $i $cpu "$rnasplice_bamdir"
  echo 'Done chr:'$i
#done

echo 'All done FRASER'

