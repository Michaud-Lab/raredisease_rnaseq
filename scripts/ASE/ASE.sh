############################
# STEP 0: DEFINE INPUTS
############################
REF="/home/renaut/scratch/reference/Homo_sapiens_autosomesXYMT/Homo_sapiens.GRCh38.dna.autosomesXYMT.fa"
GTF="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
BAM_DIR="/home/renaut/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/"
VCF_DIR="/home/renaut/scratch/raredisease_rnaseq/ASE/vcf/"
ASE="/home/renaut/scratch/raredisease_rnaseq/ASE/"
export REF ASE BAM_DIR VCF_DIR

mkdir -p $ASE/ase
mkdir -p $ASE/bed
mkdir -p $ASE/plots

echo "DONE ~~~ STEP0 ~~~ define inputs ~~~ $(date)"


#create the ref if it does not exist
#samtools faidx Homo_sapiens.GRCh38.dna.toplevel.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT >Homo_sapiens.GRCh38.dna.autosomesXYMT.fa
#samtools faidx Homo_sapiens.GRCh38.dna.autosomesXYMT.fa
#gatk CreateSequenceDictionary -R ${REF} -O ${REF}.dict

############################
# STEP 1: prepare vcfs
############################
for vcf_file in ${VCF_DIR}/*slivar.vcf.gz
  do
    bcftools view -m2 -M2 -v snps --force-samples -e 'GT="mis"' -s ^samples_to_keep.txt -Oz ${vcf_file} | \
    bcftools annotate --rename-chrs ${ASE}/chr_map.txt -Oz | \
    bcftools norm -d all -O b -o ${vcf_file}.temp
    bcftools index ${vcf_file}.temp
  done

#merge and index all vcf's 
bcftools merge ${VCF_DIR}/*temp -Oz -o ${VCF_DIR}/biallelic_sites.vcf.gz
gatk IndexFeatureFile -I ${VCF_DIR}/biallelic_sites.vcf.gz --verbosity ERROR

echo "DONE ~~~ STEP1 ~~~ bcftools prep ~~~ $(date)"

############################
# STEP 2: ASEReadCounter
############################
cat ${ASE}samples.txt | parallel -j 5 'gatk ASEReadCounter \
-R ${REF} \
-I ${BAM_DIR}{}_sorted.bam \
-V ${VCF_DIR}biallelic_sites.vcf.gz \
-O ${ASE}ase/{}.ase.tsv \
--verbosity ERROR \
--java-options '-Xmx60G' \
--max-depth-per-sample 1000 \
--min-depth 20 \
--min-mapping-quality 20 \
--min-base-quality 10 \
'

echo "DONE ~~~ STEP2 ~~~ ASEReadCounter ~~~ $(date)"

############################
# STEP 3: Convert ASE tables → BED
############################
for s in $(cat ${ASE}/samples.txt);
  do
    awk -v sample=$s 'BEGIN{OFS="\t"}
    NR>1{
    print $1,$2-1,$2,sample,$6,$7
    }' ${ASE}/ase/${s}.ase.tsv > ${ASE}/bed/${s}.bed
  done


echo "DONE ~~~ STEP3 ~~~  ASE -> BED conversio ~~~ $(date)"

############################
# STEP 4: Extract genes from GTF
############################
awk '$3=="gene"' $GTF |
awk 'BEGIN{OFS="\t"}{
match($0,/gene_id "([^"]+)"/,a);
print $1,$4-1,$5,a[1]
}' > ${ASE}/bed/genes.BED

sort -k1,1 -k2,2n ${ASE}/bed/genes.BED > ${ASE}/bed/genes.sorted.BED

echo "DONE ~~~ STEP4 ~~~ Extract genes from GTF ~~~ $(date)"

############################
# STEP 5: Merge all SNPs
############################
cat ${ASE}/bed/*.bed >${ASE}/bed/all_snps.BED
sort -k1,1 -k2,2n ${ASE}/bed/all_snps.BED > ${ASE}/bed/all_snps.sorted.BED

echo "DONE ~~~ STEP5 ~~~ Merge all SNPs ~~~ $(date)"

############################
# STEP 6: SNP → gene overlap
############################
bedtools intersect \
-a ${ASE}/bed/all_snps.sorted.BED \
-b ${ASE}/bed/genes.sorted.BED \
-wa -wb > ${ASE}/overlaps.tsv

echo "DONE ~~~ STEP6 ~~~ dedtools intersect ~~~ $(date)"

############################
# STEP 7:  Deseq2
###########################
Rscript ./ASE.R ${ASE}/overlaps.tsv ${VCF_DIR}/biallelic_sites.vcf.gz

echo "DONE ~~~ STEP7 ~~~ DeSeq2 ~~~ $(date)"

echo "Pipeline finished ~~~ $(date)"
