# =============================================================================
# ASE.R - Allele-specific expression analysis from ASEReadCounter files
# =============================================================================
source("../rnaseq_helper_functions.R")
load_install_library(c('dplyr', 'tidyr', 'DESeq2', 'vcfR', 'ggplot2',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
params = list(overlap = args[1])
params$vcf = args[2]
params$imprinting = args[3]
params$workdir = dirname(params$overlap)
ase = read.table(params$overlap,header=FALSE)

colnames(ase) = c("CHROM","SNPstart","POS","sampleID",'ref','alt','chr','START','STOP','ensemblID')

# get rid of HSJ_045_04 because it should not be there...
ase = ase[ase$sampleID != 'HSJ_042_04_PAX', ]
ase = ase[ase$sampleID != 'HSJ_083_03_PAX', ]
ase$sample_vcf = gsub('_PAX','',ase$sampleID)

ase = ase %>%
  mutate(total = ref + alt, RNA_ratio = signif(alt / total,2),RNA_DP = paste0(ref,',',alt))

#pivot
pase = ase[,c(1,3,4,14)] %>% 
  distinct(CHROM,POS,sampleID, .keep_all = TRUE) %>%
  pivot_wider(names_from = sampleID, values_from = RNA_DP,values_fill = ".")

#to filter
key_target <- paste(pase$CHROM, pase$POS, sep = ":")


# -----------------------------------------------------------------------------
# 2. Read and filter ASE based on .vcf to check the expected count
# -----------------------------------------------------------------------------
vcf_all = read.vcfR(params$vcf,verbose = TRUE)
key_vcf <- paste(vcf_all@fix[, "CHROM"], vcf_all@fix[, "POS"], sep = ":")
vcf <- vcf_all[key_vcf %in% key_target, ]

#vcf = read.vcfR('test.vcf.gz',verbose = TRUE) 

colnames(vcf@gt) = gsub('-','_',colnames(vcf@gt))


gt <- extract.gt(vcf, element = "GT")
ad <- extract.gt(vcf, element = "AD")
ad_ref <- masplit(ad, record = 1) # reference allele depth
ad_alt <- masplit(ad, record = 2) # alt  
dp  <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)
ratio <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)
vcf_extract = list(gt,ad_ref,ad_alt,dp,gq,ratio,ratio)

for(i in 1:7){
  gt_df <- data.frame(
    CHROM = vcf@fix[, "CHROM"],
    POS   = as.integer(vcf@fix[, "POS"]),
    REF   = vcf@fix[, "REF"],
    ALT   = vcf@fix[, "ALT"],
    gt <- vcf_extract[[i]],
    check.names = FALSE,
    stringsAsFactors = FALSE
    )
  
  #
  if(i==1) WGS_GT = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==2) WGS_REF = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==3) WGS_ALT = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==4) WGS_DP = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==5) WGS_GQ = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==6) WGS_ratio = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
  if(i==7) WGS_pvalue = left_join(pase[,c(1,2)], gt_df, by = c("CHROM", "POS"))
}

WGS_pvalue[,-c(1:4)] = NA

for(i in 1:nrow(WGS_GT))
#  for(i in 1:10000)
{
  for(j in 5:ncol(WGS_GT))
  {
    if(!is.na(WGS_DP[i, j])) {
      if(WGS_GT[i,j] == "0|1" | WGS_GT[i,j] == "1|0" | WGS_GT[i,j] == "0/1" | WGS_GT[i,j] == "1/0" ) {
        WGS_ratio[i,j] = 0.5
        WGS_pvalue[i,j] = signif(binom.test(as.numeric(WGS_REF[i,j]), as.numeric(WGS_DP[i,j]), p = 0.5)$p.value, 2) # binomial test for ASE
      }
    }
  }
    if(i %% 1000 == 0) print(paste0('Done ',i, ' of ' ,nrow(WGS_GT),', Time is: ', Sys.time()))
  }

     
#
WGS_ref_pivoted = pivot_longer(WGS_REF,cols = 5:ncol(WGS_pvalue),values_to = 'REF_DP')
WGS_alt_pivoted = pivot_longer(WGS_ALT,cols = 5:ncol(WGS_pvalue),values_to = 'ALT_DP')
WGS_gq_pivoted = pivot_longer(WGS_GQ,cols = 5:ncol(WGS_pvalue),values_to = 'GQ')
WGS_pvalue_pivoted = pivot_longer(WGS_pvalue,cols = 5:ncol(WGS_pvalue),values_to = 'pvalue')
pase_pivoted = pivot_longer(pase,cols = 3:ncol(pase),values_to = 'RNA_DP')

WGS_ref_pivoted = WGS_ref_pivoted[!is.na(WGS_pvalue_pivoted$pvalue),]
WGS_alt_pivoted = WGS_alt_pivoted[!is.na(WGS_pvalue_pivoted$pvalue),]
WGS_gq_pivoted = WGS_gq_pivoted[!is.na(WGS_pvalue_pivoted$pvalue),]
pase_pivoted = pase_pivoted[!is.na(WGS_pvalue_pivoted$pvalue),]
WGS_pvalue_pivoted = WGS_pvalue_pivoted[!is.na(WGS_pvalue_pivoted$pvalue),]

WGS = data.frame(WGS_pvalue_pivoted,REF_DP = WGS_ref_pivoted$REF_DP, ALT_DP = WGS_alt_pivoted$ALT_DP, GQ = WGS_gq_pivoted$GQ, RNA_DP = pase_pivoted$RNA_DP)
WGS_signif = WGS[WGS$pvalue < 0.1,]

# -----------------------------------------------------------------------------
# 4. Add geneID from hg38 annotation file
# -----------------------------------------------------------------------------
# mapping genes
map = select(org.Hs.eg.db, keys=keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
                keytype="ENTREZID", columns=c("ENSEMBL","SYMBOL"))
map = map[!is.na(map$ENSEMBL),]
colnames(map)[3] = 'geneID'

# add the gene length info
gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene,single.strand.genes.only=FALSE)
gene_locations = as.data.frame(gene_locations)
gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames))<6,]
map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
colnames(map)[c(2,3,5)] = c('ensemblID','geneID','CHROM')
map$CHROM = gsub('chr','',map$CHROM)
map = map[,c(2,3,5:9)]

# -----------------------------------------------------------------------------
# 5. Imprinting and X-linked variants
# -----------------------------------------------------------------------------
high_confidence_imprinted_genes = read.csv(params$imprinting, header =TRUE)

WGS_map = inner_join(WGS,map, by = join_by(CHROM, between(POS, start, end)))
WGS_X = WGS_map[WGS_map$CHROM == 'X',]

if(nrow(WGS_X)==0) WGS_X = WGS_map[1:10,]
WGS_X$Type = 'X'

ase_map_imprinted = WGS_map[WGS_map$geneID %in% high_confidence_imprinted_genes$Gene,]
ase_map_imprinted$Type = 'I'
ase_map_imprinted = rbind(ase_map_imprinted,WGS_X)

write.table(ase_map_imprinted, file.path(params$workdir,"gwImprinted.tsv"),sep = '\t',quote = FALSE)


# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------
signif_map = inner_join(WGS_signif,map, by = join_by(CHROM, between(POS, start, end)))
signif_map = signif_map[!is.na(signif_map$geneID),]


write.table(signif_map, file.path(params$workdir,"gwASE2.tsv"),sep = '\t',quote = FALSE)
