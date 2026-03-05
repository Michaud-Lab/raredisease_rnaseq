############################
# INPUTS
############################

REF="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"
GTF="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
BAM_DIR="/home/renaut/scratch/raredisease_rnaseq/ASE/bam_subset/"
VCF_DIR="/home/renaut/scratch/raredisease_rnaseq/ASE/vcf_subset/"
ASE="/home/renaut/scratch/raredisease_rnaseq/ASE/"

THREADS=4

samples=(HSJ_001_03_PAX HSJ_003_03_PAX HSJ_004_03_PAX HSJ_006_03_PAX HSJ_007_03_PAX HSJ_014_03_PAX HSJ_015_03_PAX HSJ_016_03_PAX HSJ_017_03_PAX HSJ_018_03_PAX HSJ_019_03_PAX)

mkdir -p $ASE/temp/ase
mkdir -p $ASE/temp/bed
#mkdir -p $ASE/temp/matrix

############################
# STEP 0
# prepare files
############################
#bcftools norm -m -any -f $ref.fa input.vcf.gz | \
#bcftools view -v snps -m 2 -M 2 | \
#bcftools view -i 'GT="het"' -O z -o filtered_output.vcf.gz

#gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.toplevel.fa -O Homo_sapiens.GRCh38.dna.toplevel.dict
#gatk IndexFeatureFile -I SNPs_16.vcf.gz 

#merge vcf / drop genotypes / bi-allelic SNPs/ rename 'chr1' to '1'
#bcftools merge p0*_chr16.vcf.gz -Oz | bcftools view -G -Oz | bcftools view -m2 -M2 -v snps -Oz | bcftools annotate --rename-chrs chr_map.txt -Oz -o biallelic_sites.vcf.gz

#index
#gatk IndexFeatureFile -I biallelic_sites.vcf.gz

############################
# STEP 1
# ASEReadCounter
############################

echo "START ~~~ ASEReadCounter"

for s in $(cat ${ASE}/samples.txt);
  do
    gatk ASEReadCounter \
       -R $REF \
       -I ${BAM_DIR}/${s}_sorted_chr16.bam \
       -V ${VCF_DIR}/biallelic_sites.vcf.gz \
       -O ${ASE}/temp/ase/${s}.ase.tsv \
       --min-depth 4 \
       --min-mapping-quality 20 \
       --min-base-quality 10

   echo 'Done ~~~ '$s
 done

echo "DONE ~~~ ASE counting"

############################
# STEP 2
# Convert ASE tables → BED
############################

for s in $(cat ${ASE}/samples.txt);
  do
    awk -v sample=$s 'BEGIN{OFS="\t"}
    NR>1{
    print $1,$2-1,$2,sample,$6,$7
    }' ${ASE}/temp/ase/${s}.ase.tsv > ${ASE}/temp/bed/${s}.bed
  done


echo "DONE ~~~ ASE -> BED conversion"

############################
# STEP 3
# Extract genes from GTF
############################

awk '$3=="gene"' $GTF |
awk 'BEGIN{OFS="\t"}{
match($0,/gene_id "([^"]+)"/,a);
print $1,$4-1,$5,a[1]
}' > ${ASE}/temp/bed/genes.bed

sort -k1,1 -k2,2n ${ASE}/temp/bed/genes.bed > ${ASE}/temp/bed/genes.sorted.bed

echo "DONE ~~~ Extract genes from GTF"

############################
# STEP 4
# Merge all SNPs
############################

cat ${ASE}/temp/bed/*.bed >${ASE}/temp/bed/all_snps.BED
sort -k1,1 -k2,2n ${ASE}/temp/bed/all_snps.BED > ${ASE}/temp/bed/all_snps.sorted.BED

echo "DONE ~~~ Merge all SNPs"

############################
# STEP 5
# SNP → gene overlap
############################

bedtools intersect \
-a ${ASE}/temp/bed/all_snps.sorted.BED \
-b ${ASE}/temp/bed/genes.sorted.BED \
-wa -wb > ${ASE}/temp/bed/overlaps.tsv

echo "DONE ~~~ Bedtools intersect"

############################
# STEP 6
# Aggregate counts per gene per sample
############################


awk 'BEGIN{OFS="\t"}
{
sample=$4
ref=$5
alt=$6
gene=$10

refSum[sample,gene]+=ref
altSum[sample,gene]+=alt
}
END{
for (k in refSum)
{
split(k,a,SUBSEP)
print a[1],a[2],refSum[k],altSum[k]
}
}' ${ASE}/temp/bed/overlaps.tsv > ${ASE}/temp/gene_counts.tsv


echo "DONE ~~~ Aggregate counts per gene per sample"



#####
#Deseq2
#####

Rscript  ./ASE.R

echo "Pipeline finished"
