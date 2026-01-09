#!/bin/bash

#arguments
chromosome=$1
start=$2
stop=$3
geneID=$4
proband=$5

#data directories
rnasplice_bamdir="$HOME/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/"
fraser_temp_bamdir="$HOME/scratch/raredisease_rnaseq/FRASER/bams_temp_subset/"
fraser_bamdir="$HOME/scratch/raredisease_rnaseq/FRASER/bams_subset/"
fraser_perregion=$fraser_bamdir'gene'$geneID'_chr'$chromosome'_'$start'_'$stop

mkdir -p $fraser_temp_bamdir
mkdir -p $fraser_bamdir
mkdir -p $fraser_perregion

#subset all the big .bam file to keep only the region of interest.
for bam_in in $rnasplice_bamdir/*sorted.bam
  do
    bam_chr_out=${bam_in//'.bam'/'_chrN.bam'}
    bam_chr_out=${bam_chr_out//$rnasplice_bamdir/$fraser_temp_bamdir}
    samtools view -b "$bam_in" $chromosome:$start-$stop > "$bam_chr_out"

    file_out=${bam_chr_out//'temp_'}
    file_out=${file_out//$fraser_bamdir/$fraser_perregion}  

    #2. Convert BAM to SAM
    samtools view -h $bam_chr_out >tempFRASER.sam

    #3. Prefix chromosome names with 'chr' (skip header lines starting with '@')
    awk 'BEGIN{OFS="\t"} 
     /^@SQ/ {
             # Edit SN (sequence name) in the header
             sub(/^SN:/, "SN:chr", $2)
             print
             next
         }
         /^@/ {print; next}  # keep other headers unchanged
         {
             # Edit reference column in alignment lines
             if ($3 != "*" && $3 !~ /^chr/) $3 = "chr"$3
             print
         }' tempFRASER.sam > tempFRASER_chr.sam

    #4. Convert back to BAM and sort
    samtools view -bS tempFRASER_chr.sam | samtools sort -o $file_out

    #5. Index the sorted BAM
    samtools index $file_out

    #6. Clean up temporary files
    rm tempFRASER.sam tempFRASER_chr.sam
  done

#run coverage once all bam files have been generated for a particular gene
samtools depth -H -a $fraser_perregion/*sorted_chrN.bam -r chr$chromosome:$start-$stop | awk 'NR % 5 == 1' >'$fraser_perregion'/gene_'$geneID'_'$proband'_depth5.csv

#clean up
rm -r $fraser_temp_bamdir
echo 'ALL DONE'



