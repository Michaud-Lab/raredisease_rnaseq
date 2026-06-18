############################
# STEP 0: DEFINE INPUTS
############################
REF="/home/renaut/scratch/reference/Homo_sapiens_autosomesXYMT/Homo_sapiens.GRCh38.dna.autosomesXYMT.fa"
GTF="/home/renaut/scratch/reference/Homo_sapiens/Homo_sapiens.GRCh38.114.gtf"
BAM_DIR="/project/def-rallard/COMMUN/raredisease_rnaseq/results_nextflow_rnasplice_09_05_2026/star_salmon/"
VCF_DIR="/project/def-rallard/COMMUN/raredisease_rnaseq/ASE/vcf/"
ASE_DIR="/project/def-rallard/COMMUN/raredisease_rnaseq/ASE/"
INPUT_DIR="/project/def-rallard/COMMUN/raredisease_rnaseq/data/input/"
export REF ASE_DIR BAM_DIR VCF_DIR INPUT_DIR

mkdir -p $ASE_DIR/ase
mkdir -p $ASE_DIR/bed
mkdir -p $ASE_DIR/plots

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
    bcftools view -m2 -M2 -v snps --force-samples -e 'GT="mis"' -S <(cut -f2 ${INPUT_DIR}vcf_bam_families.txt) -Oz ${vcf_file} | \
    bcftools annotate --rename-chrs ${INPUT_DIR}chr_map.txt -Oz | \
    bcftools norm -d all -O b -o ${vcf_file}.temp
    bcftools index ${vcf_file}.temp
  done

#merge and index all vcf's  (all VCFs calls merged here, irrespective if it is present in the RNAseq .bam data, because this is filtered further.
bcftools merge ${VCF_DIR}/*temp -Oz -o ${VCF_DIR}/biallelic_sites.vcf.gz
gatk IndexFeatureFile -I ${VCF_DIR}/biallelic_sites.vcf.gz --verbosity ERROR

echo "DONE ~~~ STEP1 ~~~ bcftools prep ~~~ $(date)"

############################
# STEP 2: ASEReadCounter
############################
cut -f1 ${INPUT_DIR}vcf_bam_families.txt | parallel -j 5 'gatk ASEReadCounter \
-R ${REF} \
-I ${BAM_DIR}{}_sorted.bam \
-V ${VCF_DIR}biallelic_sites.vcf.gz \
-O ${ASE_DIR}ase/{}.ase.tsv \
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
for s in $(cut -f1 ${INPUT_DIR}vcf_bam_families.txt);
  do
    awk -v sample=$s 'BEGIN{OFS="\t"}
    NR>1{
    print $1,$2-1,$2,sample,$6,$7
    }' ${ASE_DIR}/ase/${s}.ase.tsv > ${ASE_DIR}/bed/${s}.bed
  done


echo "DONE ~~~ STEP3 ~~~ ASE -> BED conversion ~~~ $(date)"

############################
# STEP 4: Extract genes from GTF
############################
awk '$3=="gene"' $GTF |
awk 'BEGIN{OFS="\t"}{
match($0,/gene_id "([^"]+)"/,a);
print $1,$4-1,$5,a[1]
}' > ${ASE_DIR}/bed/genes.BED

sort -k1,1 -k2,2n ${ASE_DIR}/bed/genes.BED > ${ASE_DIR}/bed/genes.sorted.BED

echo "DONE ~~~ STEP4 ~~~ Extract genes from GTF ~~~ $(date)"

############################
# STEP 5: Merge all SNPs
############################
cat ${ASE_DIR}/bed/*.bed >${ASE_DIR}/bed/all_snps.BED
sort -k1,1 -k2,2n ${ASE_DIR}/bed/all_snps.BED > ${ASE_DIR}/bed/all_snps.sorted.BED

echo "DONE ~~~ STEP5 ~~~ Merge all SNPs ~~~ $(date)"

############################
# STEP 6: SNP → gene overlap
############################
bedtools intersect \
-a ${ASE_DIR}/bed/all_snps.sorted.BED \
-b ${ASE_DIR}/bed/genes.sorted.BED \
-wa -wb > ${ASE_DIR}/overlaps.tsv

echo "DONE ~~~ STEP6 ~~~ bedtools intersect ~~~ $(date)"

############################
# STEP 7:  Deseq2
###########################
Rscript ./ASE.R ${ASE_DIR}/overlaps.tsv ${VCF_DIR}/biallelic_sites.vcf.gz ${INPUT_DIR}/high_confidence_imprinted_genes.csv

echo "DONE ~~~ STEP7 ~~~ DeSeq2 ~~~ $(date)"

echo "Pipeline finished ~~~ $(date)"
