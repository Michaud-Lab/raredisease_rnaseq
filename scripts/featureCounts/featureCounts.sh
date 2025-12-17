#!/bin/bash

cpu=16
genome_in="$HOME/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
genome_out="$HOME/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.exonid_transcriptid.gtf"
bamdir="$HOME/scratch/raredisease_rnaseq/results_07_10_2025/star_salmon/"
workdir="$HOME/scratch/raredisease_rnaseq/"

mkdir -p $workdir/featureCounts  

#Prepare genome per exon if it does not exist.
if [ ! -f "$output" ]; then
   echo "Preparing the genome reference"
   awk 'BEGIN{OFS="\t"}
   $3=="exon" {
       # extract existing attributes
       match($0, /gene_id "([^"]+)"/, g)
       match($0, /transcript_id "([^"]+)"/, t)
       match($0, /exon_number "([^"]+)"/, e)
       if(g[1]!="" && t[1]!="" && e[1]!="") {
           # create new attributes string, keeping column 9 as a single field
           $9 = "gene_id \"" g[1] "\"; exon_number \"" e[1] "\"; exon_id \"" g[1] "_" t[1] "_" e[1] "\";"
       }
       print
   }' $genome_in >$genome_out
else
   echo "File $genome_out already exists. No action taken."
fi

#list all the bams
ls -1 $bamdir*bam >bamlist

# Read BAM files into a single string
bams=$(tr '\n' ' ' < "bamlist")

#per exon (we could fraction reads in case they are multi-mapped, instead I choose to keep the single longest transcript with featureCounts.R for quantification per exon)
featureCounts -a $genome_out -o $workdir/featureCounts/feature_counts_perexon_pertranscript.txt -T $cpu  -p -B -C -g exon_id -t exon -O $bams

echo 'Done per exon expression'

#per gene
featureCounts -a $genome_in -o $workdir/featureCounts/feature_counts_pergene.txt -T $cpu -p -B -C -g gene_id -t exon -O --fraction $bams

#
echo 'Done featureCount'
rm bamlist

#get the MANE
awk '{print $14}' $workdir/../reference/MANE/MANE.GRCh38.v1.5.refseq_genomic.gtf | grep 'Ensembl:ENST' | awk '{gsub(/Ensembl:|;|\.[0-9]+/,""); print}'| uniq -c >$workdir/featureCounts/MANE.tsv

#Rscript for data clean-up
Rscript featureCounts.R $workdir

echo 'Done featureCounts.R Rscript'

