#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module load subread/2.0.6
module load r


#read configs.json
read workdir candidate_genes candidate_genes_LR candidate_genes_extra < <(
  jq -r '.general | [.workdir, .candidate_genes, .candidate_genes_LR, .candidate_genes_extra] | @tsv' ../../data/input/configs.json
)

read fc_exons fc_genes cpu genome_in genome_out bamdir MANE ens_gene masterlog  < <(
  jq -r '.featureCounts | [.fc_exons, .fc_genes, .cpu, .genome_in, .genome_out, .bamdir, .MANE, .ens_gene, .masterlog] | @tsv' ../../data/input/configs.json
)

echo $candidate_genes_extra

#run
./featureCounts.sh $cpu $genome_in $genome_out $bamdir $workdir $MANE $candidate_genes $ens_gene $masterlog $fc_exons $fc_genes $candidate_genes_LR $candidate_genes_extra
#echo 'Done featureCounts.sh'

