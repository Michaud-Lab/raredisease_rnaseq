############################
# INPUTS
############################
REF="/home/renaut/scratch/reference/Homo_sapiens_autosomesXYMT/Homo_sapiens.GRCh38.dna.autosomesXYMT.fa"
GTF="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
#BAM_DIR="/home/renaut/scratch/raredisease_rnaseq/ASE/bam_subset/"
BAM_DIR="/home/renaut/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/"

VCF_DIR="/home/renaut/scratch/raredisease_rnaseq/ASE/vcf/"
ASE="/home/renaut/scratch/raredisease_rnaseq/ASE/"

mkdir -p $ASE/temp/ase
mkdir -p $ASE/temp/bed
mkdir -p $ASE/temp/plots

############################
# STEP 0
# prepare files
############################
#bcftools norm -m -any -f $ref.fa input.vcf.gz | \
#bcftools view -v snps -m 2 -M 2 | \
#bcftools view -i 'GT="het"' -O z -o filtered_output.vcf.gz

#gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.toplevel.fa -O Homo_sapiens.GRCh38.dna.toplevel.dict
#gatk IndexFeatureFile -I SNPs_16.vcf.gz 

#*joint.GRCh38.small_variants.phased.norm.slivar.vcf.gz 

for vcf_file in ${VCF_DIR}/*slivar.vcf.gz
  do
    bcftools view -m2 -M2 -v snps --force-samples -e 'GT="mis"' -s ^samples_to_keep.txt -Oz ${vcf_file} | \
    bcftools annotate --rename-chrs ${ASE}/chr_map.txt -Oz | \
    bcftools norm -d all -O b -o ${vcf_file}.temp
  done

#merge vcf /  bi-allelic SNPs/ rename 'chr1' to '1'
#bcftools merge p0*_chr16.vcf.gz -Oz | bcftools view -Oz | bcftools view -m2 -M2 -v snps -Oz | bcftools annotate --rename-chrs ../chr_map.txt -Oz |  bcftools norm -d all -o biallelic_sites.vcf.gz
if [ ! -f "${VCF_DIR}/biallelic_sites.vcf.gz" ]; then
  #bcftools merge ${VCF_DIR}/*joint*gz -Oz | bcftools view -Oz | bcftools view -m2 -M2 -v snps -Oz | bcftools annotate --rename-chrs ${ASE}/chr_map.txt -Oz |  bcftools norm -d all -o ${VCF_DIR}/biallelic_sites.vcf.gz
  bcftools merge ${VCF_DIR}/*temp -Oz -o ${VCF_DIR}/biallelic_sites.vcf.gz
  echo 'File'${VCF_DIR}'/biallelic_sites.vcf.gz does not exists. Lets create it'
else
  echo 'File'${VCF_DIR}'/biallelic_sites.vcf.gz already exists.'
fi


#index
gatk IndexFeatureFile -I ${VCF_DIR}/biallelic_sites.vcf.gz

############################
# STEP 1
# ASEReadCounter
############################
echo "START ~~~ ASEReadCounter"

cat /home/renaut/scratch/raredisease_rnaseq/ASE/samples.txt | parallel -j 5 'gatk ASEReadCounter \
-R /home/renaut/scratch/reference/Homo_sapiens_autosomesXYMT/Homo_sapiens.GRCh38.dna.autosomesXYMT.fa \
-I /home/renaut/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/{}_sorted.bam \
-V /home/renaut/scratch/raredisease_rnaseq/ASE/vcf_subset/biallelic_sites.vcf.gz \
-O /home/renaut/scratch/raredisease_rnaseq/ASE/temp/ase/{}.ase.tsv \
--java-options '-Xmx60G' \
--max-depth-per-sample 1000 \
--min-depth 20 \
--min-mapping-quality 20 \
--min-base-quality 10 \
'


#for s in $(cat ${ASE}/samples_small.txt);
#  do
#    gatk ASEReadCounter \
#       -R $REF \
#       -I ${BAM_DIR}/${s}_sorted.bam \
#       -V ${VCF_DIR}/biallelic_sites.vcf.gz \
#       -O ${ASE}/temp/ase/${s}.ase.tsv \
#       --java-options '-Xmx100G' \
#       --max-depth-per-sample 1000 \
#       --min-depth 20 \
#       --min-mapping-quality 20 \
#       --min-base-quality 10
#
#   echo 'Done ~~~ '$s
# done

echo "DONE ~~~ ASE counting"

############################
# STEP 2
# Convert ASE tables → BED
############################

for s in $(cat ${ASE}/samples_small.txt);
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
}' > ${ASE}/temp/bed/genes.BED

sort -k1,1 -k2,2n ${ASE}/temp/bed/genes.BED > ${ASE}/temp/bed/genes.sorted.BED

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
-wa -wb > ${ASE}/overlaps.tsv

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
}' ${ASE}/overlaps.tsv > ${ASE}/temp/gene_counts.tsv


echo "DONE ~~~ Aggregate counts per gene per sample"

#####
#Deseq2
#####

Rscript ./ASE2.R ${ASE}/overlaps.tsv ${VCF_DIR}/biallelic_sites.vcf.gz

echo "Pipeline finished"
