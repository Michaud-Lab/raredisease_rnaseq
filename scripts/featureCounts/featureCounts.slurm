#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16

cpu=16
module load subread/2.0.6
module load r


#read configs.json
read workdir rnasplice_bamdir gtf candidate_genes candidate_genes_LR candidate_genes_extra < <(
  jq -r '.general | [.workdir, .rnasplice_bamdir, .gtf, .candidate_genes, .candidate_genes_LR, .candidate_genes_extra] | @tsv' ../../data/input/configs.json
)

read fc_exons fc_genes gtf_out gtf_MANE ens_gene masterlog  < <(
  jq -r '.featureCounts | [.fc_exons, .fc_genes, .gtf_out, .gtf_MANE, .ens_gene, .masterlog] | @tsv' ../../data/input/configs.json
)

echo $candidate_genes_extra

#run
./featureCounts.sh $cpu $gtf $gtf_out $rnasplice_bamdir $workdir $gtf_MANE $candidate_genes $ens_gene $masterlog $fc_exons $fc_genes $candidate_genes_LR $candidate_genes_extra
#echo 'Done featureCounts.sh'

