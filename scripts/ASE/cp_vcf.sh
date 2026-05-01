

ASE="/home/renaut/scratch/raredisease_rnaseq/ASE/"


#subset the vcf to create a test dataset
for vcf in $(cut -f2 ${ASE}/samples.txt);;
  do
    echo $vcf
    cp /home/renaut/projects/ctb-rallard/COMMUN/PacBioData/OutputFamilies/$vcf/_LAST/out/tertiary_small_variant_filtered_vcf/$vcf.joint.GRCh38.small_variants.phased.norm.slivar.vcf.gz .
#    gzip -d $vcf.joint.GRCh38.small_variants.phased.norm.slivar.vcf.gz
#    grep -w '^#\|^#CHROM\|^chr16' $vcf.joint.GRCh38.small_variants.phased.norm.slivar.vcf > $vcf'_chr16.vcf'
#    bgzip $vcf'_chr16.vcf'
#    bcftools index -t $vcf'_chr16.vcf.gz'
  done

##move from narval to fir with globus
#all_16.vcf
