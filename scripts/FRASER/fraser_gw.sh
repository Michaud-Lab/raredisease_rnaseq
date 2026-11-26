#!/bin/bash

#arguments
chromosome=$1
rnasplice_bamdir=$2
fraser_chr_bamdir=$3

#chr 16 is problematic because of hemoglobin gene (HBA), so let's downsample it. 
if [ "$chromosome" = "16" ]; then
    fraction=" -s 0.10 "
elif [ "$chromosome" = "21" ]; then
    fraction=" -s 0.10 "
else
    fraction=""
fi


#subset all the big .bam file to keep only the region of interest.
for bam_in in $rnasplice_bamdir/*sorted.bam
  do
    bam_chr_out=${bam_in//'sorted'/$chromosome}
    bam_chr_out=$(basename "$bam_chr_out")
    bam_chr_out=$fraser_chr_bamdir$bam_chr_out
    #bam_chr_out=${bam_chr_out//$rnasplice_bamdir/$fraser_chr_bamdir}
    if [ ! -f "$bam_chr_out" ]; then
      samtools view $fraction --threads 10 -b "$bam_in" $chromosome > "$bam_chr_out"
    else
      echo "~~~ File $bam_chr_out already exists. No action taken ~~~"
    fi

    #5. Index the sorted BAM (only if there is no index yet)
    if [ ! -f "$bam_chr_out.bai" ]; then
      samtools index "$bam_chr_out"
    fi
  done
