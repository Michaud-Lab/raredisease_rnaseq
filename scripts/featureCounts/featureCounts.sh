#!/bin/bash

cpu=$1
genome_in=$2
genome_out=$3
bamdir=$4
workdir=$5
MANE=$6
ens_gene=$7
masterlog=$8
fc_exons=$9
fc_genes="${10}"
candidate_genes="${11}"   # e.g. data/input/candidate_genes.csv (column "proband" lists the probands to keep)

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
   awk '{print $10,$14}' $MANE | grep 'Ensembl:ENST' | awk '{gsub(/Ensembl:|;|\.[0-9]+/,""); print}'| uniq -c >$workdir/featureCounts/MANE.tsv
else
  echo "~~~ File $workdir/featureCounts/MANE.tsv already exists. No action taken ~~~"
fi

#list all the bams SR: restrict to the probands listed in the "proband" column of candidate_genes.csv
if [ ! -f "$candidate_genes" ]; then
   echo "ERROR: candidate genes file $candidate_genes not found" >&2
   exit 1
fi

probands=$(mktemp)
awk -F',' 'NR==1{for(i=1;i<=NF;i++){gsub(/"|\r/,"",$i); if($i=="proband") c=i}; if(!c){print "ERROR: no proband column" > "/dev/stderr"; exit 1}; next}
           {gsub(/"|\r/,"",$c); if($c!="") print $c}' "$candidate_genes" | sort -u >$probands
[ -s "$probands" ] || { echo "ERROR: no probands read from $candidate_genes" >&2; exit 1; }

# keep a bam only if its name (minus directory and "_sorted*.bam"/".bam") is a listed proband
bams_SR=$(ls -1 $bamdir*bam | awk 'NR==FNR{p[$1]; next} {n=$0; sub(/.*\//,"",n); sub(/_sorted.*$/,"",n); sub(/\.bam$/,"",n); if(n in p) print}' $probands - | tr '\n' ' ')
if [ -z "${bams_SR// /}" ]; then echo "ERROR: no bams match the probands" >&2; rm $probands; exit 1; fi

echo "~~~ $(wc -w <<< "$bams_SR") bams kept for $(wc -l < $probands) probands ~~~"

rm $probands

#per gene and exon (or called feature_counts_pergene_SR.txt)
if [ ! -f "$workdir/featureCounts/feature_counts_pergene.txt" ]; then
  featureCounts -a $genome_in -o $workdir/featureCounts/feature_counts_pergene.txt -T $cpu -p -B -C -g gene_id -t exon -O --fraction $bams_SR
  featureCounts -a $genome_out -o $workdir/featureCounts/feature_counts_perexon_pertranscript.txt -T $cpu  -p -B -C -g exon_id -t exon -O $bams_SR
else
  echo "~~~ File $workdir/featureCounts/feature_counts_pergene.txt already exists. No action taken ~~~"
fi 


###LONG READ DATA
###list all the bams LR
#ls -1 $bamdir/long_read/*bam >bamlist
#bams_LR=$(tr '\n' ' ' < "bamlist")
#rm bamlist

#per gene and exon
#if [ ! -f "$workdir/featureCounts/feature_counts_pergene_LR.txt" ]; then
#  featureCounts -a $genome_in -o $workdir/featureCounts/feature_counts_pergene_LR.txt -T $cpu -L -g gene_id -t exon -O --fraction $bams_LR
#  featureCounts -a $genome_out -o $workdir/featureCounts/feature_counts_perexon_pertranscript_LR.txt -T $cpu  -L -g exon_id -t exon -O $bams_LR
#else
#  echo "~~~ File $workdir/featureCounts/feature_counts_pergene_LR.txt already exists. No action taken ~~~"
#fi 


#add the expression values of LR files to SR dataset
#paste $workdir/featureCounts/feature_counts_pergene_SR.txt  <(cut -f7-  $workdir/featureCounts/feature_counts_pergene_LR.txt) >$workdir/featureCounts/feature_counts_pergene.txt
#paste $workdir/featureCounts/feature_counts_perexon_pertranscript_SR.txt  <(cut -f7-  $workdir/featureCounts/feature_counts_perexon_pertranscript_LR.txt) >$workdir/featureCounts/feature_counts_perexon_pertranscript.txt

#Rscript for data clean-up
Rscript featureCounts.R $workdir $ens_gene $masterlog $fc_exons $fc_genes $genome_in

#echo $workdir $ens_gene $masterlog $fc_exons $fc_genes $genome_in
echo 'Done featureCounts.R Rscript'

