#!/bin/bash

#arguments
gene=$1
transcript=$2

ref_file=$3
ref_annot=$4

bam_dir=$5
bam_file=$6
consensus_dir=$7
region=$8

#1.extract protein_coding or lncRNA exons from gtf. 
awk 'NR==1 || $3 == "exon"' $ref_annot | awk '$24 ~ /protein_coding|lncRNA/' >$bam_dir'exons.gtf'
cut -f1,4,5,9 $bam_dir'exons.gtf' >$bam_dir'exons.bed'

#2extract ONLY the region of interest AND only the transcript of interest
grep $transcript $bam_dir'exons.bed' >$bam_dir$transcript'.bed'

#3.produce vcf then a tsv
bcftools mpileup -I -r $region  -O v -f $ref_file $bam_dir$bam_file | \
bcftools call -m -O v | \
bcftools query -HH -f '%CHROM\t%POS\t%DP\t%REF\t%ALT\n' >$bam_dir'variants.tsv'

#4.convert variants to bed file (and remove the 'chr' from the chromosome number).
awk 'BEGIN{OFS="\t"} {print $1, $2, $2, $3, $4, $5}' $bam_dir'variants.tsv' | awk 'BEGIN{FS=OFS="\t"} {gsub(/chr/, "", $1); print}' >$bam_dir'variants.bed'

#5.annotate exons
bedtools intersect -header -loj -wa -wb -a $bam_dir'variants.bed' -b $bam_dir$transcript'.bed' | awk '{print $1,$3,$4,$5,$6,$15,$19,$21}' OFS="\t" >$bam_dir'gene'$gene'variants_annotated.tsv'

#cleanup
rm $bam_dir$transcript'.bed'  $bam_dir'exons.gtf'  $bam_dir'exons.bed'

