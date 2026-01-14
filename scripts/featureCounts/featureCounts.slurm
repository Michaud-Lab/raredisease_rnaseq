#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module load subread/2.0.6
module load r

cpu=16
workdir=$HOME"/scratch/raredisease_rnaseq/"
genome_in=$HOME"/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
genome_out=$HOME"/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.exonid_transcriptid.gtf"
bamdir=$workdir"results_06_01_2026/star_salmon/"
MANE=$HOME"/scratch/reference/MANE/MANE.GRCh38.v1.5.refseq_genomic.gtf"
candidate_genes=$workdir"data/candidate_genes_3.txt"

ens_gene=$workdir"data/ensembl_geneid.tsv"
masterlog=$workdir"data/CHUSJ_Master_Linking_Log_modif.xlsx"
fc_exons=$workdir"featureCounts/feature_counts_perexon_pertranscript.txt"
fc_genes=$workdir"featureCounts/feature_counts_pergene.txt"

./featureCounts.sh $cpu $genome_in $genome_out $bamdir $workdir $MANE $candidate_genes $ens_gene $masterlog $fc_exons $fc_genes
echo 'Done featureCounts.sh'

