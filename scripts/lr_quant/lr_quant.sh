
lrpath="/home/renaut/scratch/raredisease_rnaseq/"
resultpath="/home/renaut/scratch/raredisease_rnaseq/lr_quant/"
ref_genome="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiensChr.GRCh38.dna.primary_assembly.fa"
gtf="/home/renaut/scratch/reference/Human_hg38_Gencode_v39/gencode.v39.annotation.gtf"
gtf_c="/home/renaut/scratch/reference/Human_hg38_Gencode_v39/gencode.v39.annotation.sorted.gtf"
cpu=$1
rnasplice_bamdir="$HOME/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/long_read/"


mkdir -p $resultpath

#0.generate a subsetted samples (1M reads)
#samtools view -h "$lrpath"1_A01/hifi_reads/m84196_260317_144727_s1.hifi_reads.bcM0001.bam | head -n 600000 | samtools view -b -o "$lrpath"test_data/hifi_reads/test.hifi_reads.bcM0001.bam

#1. Generate segmented reads
#skera split -j $cpu "$lrpath"1_A01/hifi_reads/m84196_260317_144727_s1.hifi_reads.bcM0001.bam "$lrpath"/reference/mas8_primers.fasta "$resultpath"segmented.bam

#2. Primer removal and demux
#lima -j $cpu "$resultpath"segmented.bam "$lrpath"/reference/IsoSeq_v2_primers_12.fasta "$resultpath"movieX.fl.bam --isoseq --peek-guess

#
samples=$(ls "$resultpath"movieX.fl.IsoSeqX_bc*bam)

for s in $samples
  do
    samp_dir=${s//'movieX.fl.IsoSeqX_'} 
    samp_dir="${samp_dir//'_5p--IsoSeqX_3p.bam'}" 
    bc=${samp_dir//'/home/renaut/scratch/longread_rnaseq/results_6M/'}
    mkdir -p $samp_dir
    echo $samp_dir

    #3. Refine (trimming polyA & concatemer)
    isoseq refine -j $cpu "$resultpath"movieX.fl.IsoSeqX_"$bc"_5p--IsoSeqX_3p.bam "$lrpath"/reference/IsoSeq_v2_primers_12.fasta "$samp_dir"/movieX.flnc.bam --require-polya

    #3b. Merge SMART cells (if necessary)
    #`ls movie*.flnc.bam movie*.flnc.bam movie*.flnc.bam > flnc.fofn`

    #4. Cluster isoforms
    isoseq cluster2 -j $cpu "$samp_dir"/movieX.flnc.bam "$samp_dir"/transcripts.bam

    #5. Check sequences (generate a .fastq)
    #samtools fastq "$samp_dir"/transcripts.bam >"$samp_dir"/transcripts.fastq.gz

    #6a. ppmm2: map clustered reads to human genome
    pbmm2 align -j $cpu --preset ISOSEQ --sort "$samp_dir"/transcripts.bam "$ref_genome" "$samp_dir"/mapped.bam

    #6b. ppmm2: map ALL reads to human genome (I need this for FRASER)
    pbmm2 align -j $cpu --preset ISOSEQ --sort "$samp_dir"/movieX.flnc.bam "$ref_genome" "rnasplice_bamdir"/"$bc"_sorted.bam 

    #7. Collapse into single isoforms (gff contains the sequence annotation. collapsed.flnc_count.txt contains the count of each unique sequence) 
    isoseq collapse "$samp_dir"/mapped.bam "$samp_dir"/movieX.flnc.bam "$samp_dir"/collapsed.gff

    #8. Prepare reference files for pigeon (create a sorted gtf)
    pigeon prepare "$gtf" "$ref_genome"
    pigeon prepare "$samp_dir"/collapsed.gff

    #9. Classify isoforms
    pigeon classify -j $cpu "$samp_dir"/collapsed.sorted.gff "$gtf_c" "$ref_genome"  --fl "$samp_dir"/collapsed.flnc_count.txt  -d "$samp_dir" #add count data outputted from isoseq collapse
    pigeon filter -j $cpu --max-distance 500 "$samp_dir"/collapsed_classification.txt --isoforms "$samp_dir"/collapsed.sorted.gff # filter .gff

    # 11. Gene saturation to check if you sequenced enough
    pigeon report --exclude-singletons "$samp_dir"/collapsed_classification.filtered_lite_classification.txt "$samp_dir"/saturation.txt

    echo "DONE sample ~~~ ""$bc"
  done


