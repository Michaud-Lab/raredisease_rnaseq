#!/bin/bash

cpu=$1
genome_in=$2
genome_out=$3
bamdir=$4
workdir=$5
MANE=$6 
candidate_genes=$7
ens_gene=$8
masterlog$9
fc_exons=${10} 
fc_genes${11}

$ens_gene $masterlog $fc_exons $fc_genes 



mkdir -p $workdir/featureCounts  

#Prepare genome per exon if it does not exist.
if [ ! -f "$genome_out" ]; then
   echo "~~~ Preparing the genome reference ~~~"
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
   echo "~~~ File $genome_out already exists. No action taken ~~~"
fi

#get the MANE
if [ ! -f "$workdir/featureCounts/MANE.tsv" ]; then
  awk '{print $14}' $MANE | grep 'Ensembl:ENST' | awk '{gsub(/Ensembl:|;|\.[0-9]+/,""); print}'| uniq -c >$workdir/featureCounts/MANE.tsv
else
  echo "~~~ File $workdir/featureCounts/MANE.tsv already exists. No action taken ~~~"
fi

#list all the bams
ls -1 $bamdir*bam >bamlist
bams=$(tr '\n' ' ' < "bamlist")
rm bamlist

#per gene and exon
if [ ! -f "$workdir/featureCounts/feature_counts_pergene.txt" ]; then
  featureCounts -a $genome_in -o $workdir/featureCounts/feature_counts_pergene.txt -T $cpu -p -B -C -g gene_id -t exon -O --fraction $bams
  featureCounts -a $genome_out -o $workdir/featureCounts/feature_counts_perexon_pertranscript.txt -T $cpu  -p -B -C -g exon_id -t exon -O $bams
else
  echo "~~~ File $workdir/featureCounts/feature_counts_pergene.txt already exists. No action taken ~~~"
fi 

#Rscript for data clean-up
Rscript featureCounts.R $workdir $candidate_genes $ens_gene $masterlog $fc_exons $fc_genes


echo 'Done featureCounts.R Rscript'

