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

    ###Down sample region of interest if necessary
    nreads=$(samtools view -c $rnasplice_bamdir'HSJ_001_03_PAX_sorted.bam' $chromosome:$start-$stop)

    if [ "$nreads" -gt 2000000 ]; then
       echo "Very high coverage ($nreads reads) → 0.1%"
       echo "$bam_in"
       frac="0.001"
    elif [ "$nreads" -gt 1000000 ]; then
       echo "High coverage ($nreads reads) → 10%"
       echo "$bam_in"
       frac="0.10"
    else
       echo "Regular coverage ($nreads reads) → no downsampling"
       frac=""
    fi

    if [ -n "$frac" ]; then
       samtools view -s "$frac" -b "$bam_in" $chromosome:$start-$stop > "$bam_chr_out"
    else
       samtools view -b "$bam_in" $chromosome:$start-$stop > "$bam_chr_out"
    fi

    file_out=${bam_chr_out//'temp_'}
    file_out=${file_out//$fraser_bamdir/$fraser_perregion}  

    #2. Convert BAM to SAM
    if [ ! -f "$file_out" ];
      then samtools view -h $bam_chr_out >tempFRASER.sam

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

      #3B. chrMT ->chrM ((skip header lines starting with '@')
      sed 's/\bchrMT\b/chrM/g' tempFRASER_chr.sam >tempFRASER_chrM.sam

      #4. Convert back to BAM and sort
      samtools view -bS tempFRASER_chrM.sam | samtools sort -o $file_out

      #5. Index the sorted BAM
      samtools index $file_out

      #6. Clean up temporary files
      rm tempFRASER.sam tempFRASER_chr.sam tempFRASER_chrM.sam
    fi
  done

#run coverage once all bam files have been generated for a particular gene
[ "$chromosome" = "MT" ] && chromosome="M"
samtools depth -H -a $fraser_perregion/*sorted_chrN.bam -r chr$chromosome:$start-$stop | awk 'NR % 5 == 1' >"$fraser_perregion"/gene_"$geneID"_"$proband"_depth5.csv

#clean up
rm -r $fraser_temp_bamdir
echo 'DONE fraser.sh' $fraser_perregion
