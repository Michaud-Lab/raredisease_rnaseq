#!/bin/bash

#arguments
chromosome=$1
rnasplice_bamdir=$2
fraser_chr_bamdir=$3

#chr 16 is problematic because of hemoglobin gene (HBA), so let's downsample it. 
if [ $chromosome = "16" ]; then
    fraction=" -s 0.25 "
else
    fraction=""
fi


#subset all the big .bam file to keep only the region of interest.
for bam_in in $rnasplice_bamdir/*sorted.bam
  do
    bam_chr_out=${bam_in//'sorted'/$chromosome}
    bam_chr_out=${bam_chr_out//$rnasplice_bamdir/$fraser_chr_bamdir}
    samtools view $fraction --threads 16 -b "$bam_in" $chromosome > "$bam_chr_out"

    #5. Index the sorted BAM
    samtools index $bam_chr_out
  done

echo "Done chromosome:" $chromosome
