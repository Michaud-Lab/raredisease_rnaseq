module load bcftools

ASE="/home/renaut/scratch/slivar_vcf/"
VCF_DIR="/home/renaut/scratch/slivar_vcf/"

#subset the vcf to create a test dataset
for vcf in $(cut -f3 ${ASE}/vcf_bam_families.txt);
  do
    echo $vcf
    cp /home/renaut/projects/ctb-rallard/COMMUN/PacBioData/OutputFamilies/$vcf/_LAST/out/tertiary_small_variant_filtered_vcf/$vcf.joint.GRCh38.small_variants.phased.norm.slivar.vcf.gz .
  done

#prepare vcf
for vcf_file in ${VCF_DIR}/*slivar.vcf.gz
  do
    bcftools view -m2 -M2 -v snps --force-samples -e 'GT="mis"' -S <(cut -f2 ${VCF_DIR}vcf_bam_families.txt) -Oz ${vcf_file} | \
    bcftools annotate --rename-chrs ${VCF_DIR}chr_map.txt -Oz | \
    bcftools norm -d all -O b -o ${vcf_file}.tempVCFASE
    bcftools index ${vcf_file}.tempVCFASE
  done

bcftools merge ${VCF_DIR}/*tempVCFASE -Oz -o ${VCF_DIR}/biallelic_sites.vcf.gz

rm *tempVCFASE*
