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
colnames(ase) = c("SNPchr","SNPstart","pos","sampleID",'ref','alt','chr','START','STOP','ensemblID')

# get rid of HSJ_045_04 because it should not be there...
ase = ase[ase$sampleID != 'HSJ_042_04_PAX', ]
ase$sample_vcf = gsub('_PAX','',ase$sampleID)
ase$pos = as.character(ase$pos)
ase$pvalue = NA
ase$WGS_GT = NA
ase$WGS_DP = NA
ase$WGS_GQ = NA
ase$WGS_ratio = NA

ase = ase %>%
  mutate(total = ref + alt, RNA_ratio = signif(alt / total,2),RNA_DP = paste0(ref,',',alt))

# get the actual chromosome names from the ase file
aes_rle = rle(ase$SNPchr)
aes_rle = data.frame(lengths = aes_rle$lengths, values = aes_rle$values, c = 1:length(aes_rle$lengths))

# -----------------------------------------------------------------------------
# 2. Read and filter ASE based on .vcf to check the expected count
# -----------------------------------------------------------------------------
vcf = read.vcfR(params$vcf,verbose = TRUE)
colnames(vcf@gt) = gsub('-','_',colnames(vcf@gt))

# split things per chromosome for quicker access
vcf_per_chr = list()

for(c in 1:length(aes_rle$values)){
  vcf_per_chr[[c]] = data.frame(position = vcf@fix[vcf@fix[,1] == aes_rle$values[c],2],vcf@gt[vcf@fix[,1] == aes_rle$values[c],])
}
names(vcf_per_chr) = aes_rle$values

# -----------------------------------------------------------------------------
# 3. Calculate expected count and binomial test for 0.5 HET variants
# -----------------------------------------------------------------------------
for(i in 1:nrow(ase)){
  # temp vcf file when using a new chromosome
  if((i == 1) | ((i!=1) && (ase[i-1,1] != ase[i,1]))) {
    temp_vcf_per_chr = vcf_per_chr[names(vcf_per_chr) == ase$SNPchr[i]][[1]]
  }

  # temp vcf file when using a locus
  if((i == 1) | ((i!=1) && (ase[i-1,2] != ase[i,2]))) {
    temp_geno = temp_vcf_per_chr[temp_vcf_per_chr$position == ase$pos[i],]
  }

  temp_geno2 = 'nosnp'
  if(nrow(temp_geno) >0){
   temp_geno2 = temp_geno[,colnames(temp_geno) == ase$sample_vcf[i]]}

  if(nchar(temp_geno2) > 15) {
    gt = strsplit(temp_geno2, ':')[[1]]
    ase$WGS_GT[i] = gt[1]
    ase$WGS_DP[i] = gt[3]
    ase$WGS_GQ[i] = gt[4]
  }

  if(i %% 10000 == 0) print(paste('Done ',i,' of ',nrow(ase), ', Time is: ',Sys.time()))
}

# subset the hetero SNVs.
ase$WGS_ratio[ase$WGS_GT == "0|1" | ase$WGS_GT == "1|0" | ase$WGS_GT == "0/1" | ase$WGS_GT == "1/0"] = 0.5
ase$WGS_ratio[ase$WGS_GT == ".|." | ase$WGS_GT == "./."] = 1
ase$WGS_ratio[ase$WGS_GT == "0|0" | ase$WGS_GT == "0/0"] = 1
ase$WGS_ratio[ase$WGS_GT == "1|1" | ase$WGS_GT == "1/1"] = 0

ase_signif = ase
ase_signif = ase_signif[!is.na(ase_signif$WGS_GT),]
ase_signif = ase_signif[ase_signif$WGS_GT == "0|1" | ase_signif$WGS_GT == "1|0" | ase_signif$WGS_GT == "0/1" | ase_signif$WGS_GT == "1/0",]

for(i in 1:nrow(ase_signif)){
  ase_signif$WGS_ratio[i] = 0.5 # heterozygous
  ase_signif$pvalue[i] = signif(binom.test(ase_signif$ref[i], ase_signif$total[i], p = 0.5)$p.value, 2) # binomial test for ASE

  if(i %% 5000 == 0) print(paste('Done ',i,' of ',nrow(ase_signif), ', Time is: ',Sys.time()))
}

ase_signif = ase_signif[ase_signif$pvalue < 0.5, ]
ase_signif$pvalue[ase_signif$pvalue < 1e-50] = 1e-50

ase_signif_dedup = ase[!duplicated(ase[,c(3,4)]),]
ase_pivoted = pivot_wider(ase_signif_dedup[,c(3,4,19)],names_from = sampleID,values_from = RNA_DP)

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
colnames(map)[c(2,3,5)] = c('ensemblID','geneID','chr')
map$chr = gsub('chr','',map$chr)

# -----------------------------------------------------------------------------
# 5. Imprinting and X-linked variants
# -----------------------------------------------------------------------------
high_confidence_imprinted_genes = read.csv(params$imprinting, header =TRUE)
ase_map = merge(ase,map,by = 'ensemblID',all.x = TRUE)
ase_map = ase_map[!is.na(ase_map$WGS_ratio),]
ase_map = ase_map[!is.na(ase_map$geneID),]
ase_map_X = ase_map[ase_map$chr.x == 'X',]
ase_map_X = ase_map_X[ase_map_X$WGS_ratio == 0.5, ]
ase_map_X$Type = 'X'

ase_map_imprinted = ase_map[ase_map$geneID %in% high_confidence_imprinted_genes$Gene,]
ase_map_imprinted = ase_map_imprinted[ase_map_imprinted$WGS_ratio == 0.5, ]
ase_map_imprinted$Type = 'I'
ase_map_imprinted = rbind(ase_map_imprinted, ase_map_X)
ase_map_imprinted = ase_map_imprinted[,c(5,1,21,2,4,19:18,13:16,12,28)]
colnames(ase_map_imprinted)[4] = 'chr'

write.table(ase_map_imprinted, file.path(params$workdir,"gwImprinted.tsv"),sep = '\t',quote = FALSE)


# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------
signif_map = merge(ase_signif, map, by = 'ensemblID', all.x = TRUE)
signif_map_ASE = signif_map[, c(5, 1, 21, 2, 4, 19:18, 13:16, 12)]
signif_map_ASE = signif_map_ASE[!is.na(signif_map_ASE$geneID), ]
colnames(signif_map_ASE)[4] = 'chr'

write.table(signif_map_ASE, file.path(params$workdir,"gwASE.tsv"),sep = '\t',quote = FALSE)
