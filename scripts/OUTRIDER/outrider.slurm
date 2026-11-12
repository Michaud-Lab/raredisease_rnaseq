#!/bin/bash

#SBATCH --time=00:59:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18

module load r
cpu=16

# Set to TRUE to (re-)run the expensive per-gene OUTRIDER fit; when FALSE (default),
# outrider.R loads the cached table_genes.rds from a previous run instead. Ignored
# once table_genes.rds already exists (it's always loaded from cache in that case).
force_outrider=TRUE

#read configs.json
read workdir gtf reference_cohorte< <(
  jq -r '.general | [.workdir, .gtf, .reference_cohorte] | @tsv' ../../data/input/configs.json
)

read fc_pergene fc_perexon < <(
  jq -r '.OUTRIDER | [.fc_pergene, .fc_perexon] | @tsv' ../../data/input/configs.json
)

#R script launches the bash script to prepare data and continues with analysis.
Rscript outrider.R $workdir $fc_pergene $fc_perexon $cpu $gtf $force_outrider $reference_cohorte
#echo $workdir $fc_pergene $fc_perexon $cpu  $gtf $force_outrider

Rscript outrider_one_vs_ref.R $workdir $fc_pergene $fc_perexon $cpu $gtf $force_outrider
