#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module load subread/2.0.6
module load r


#read configs.json
read workdir candidate_genes < <(
  jq -r '.general | [.workdir, .candidate_genes] | @tsv' ../configs.json
)

read fc_exons fc_genes cpu genome_in genome_out bamdir MANE ens_gene masterlog  < <(
  jq -r '.featureCounts | [.fc_exons, .fc_genes, .cpu, .genome_in, .genome_out, .bamdir, .MANE, .ens_gene, .masterlog] | @tsv' ../configs.json
)

#run
./featureCounts.sh $cpu $genome_in $genome_out $bamdir $workdir $MANE $candidate_genes $ens_gene $masterlog $fc_exons $fc_genes
echo 'Done featureCounts.sh'

