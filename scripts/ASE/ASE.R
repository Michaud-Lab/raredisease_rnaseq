########################################
# ASE analysis from ASEReadCounter files
########################################
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(vcfR)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))


########################################
# 1. Read and filter ASE overlap.txt file
########################################
args = commandArgs(trailingOnly=TRUE)
params = list(overlap = args[1])
params$vcf = args[2]
params$imprinting = args[3]
params$workdir = dirname(params$overlap)
ase <- read.table(params$overlap,header=FALSE)
colnames(ase) <- c("SNPchr","SNPstart","pos","sampleID",'ref','alt','chr','START','STOP','ensemblID')

ase$sample_vcf = gsub('_PAX','',ase$sampleID)
ase$pos = as.character(ase$pos)
ase$pvalue = NA
ase$WGS_GT = NA
ase$WGS_DP = NA
ase$WGS_GQ = NA
ase$WGS_ratio = 1

ase <- ase %>%
  mutate(total = ref + alt, RNA_ratio = signif(alt / total,2),RNA_DP = paste0(ref,',',alt))

#get the actual chromosome names from the ase file
aes_rle = rle(ase$SNPchr)
aes_rle = data.frame(lengths = aes_rle$lengths, values = aes_rle$values, c = 1:length(aes_rle$lengths))

########################################
# 2. Read and filter ASE based on .vcf to check the expected count.
########################################
vcf <- read.vcfR(params$vcf,verbose = F)
colnames(vcf@gt) = gsub('-','_',colnames(vcf@gt))

#split things per chromosome for quicker access
vcf_per_chr = list()

for(c in 1:length(aes_rle$values)){
  vcf_per_chr[[c]] = data.frame(position = vcf@fix[vcf@fix[,1] == aes_rle$values[c],2],vcf@gt[vcf@fix[,1] == aes_rle$values[c],])
}

names(vcf_per_chr) = aes_rle$values



########################################
# 3. Calculate expected count and do a binomial test for the 0.5 HET variants.
########################################

for(i in 1:nrow(ase)){
  #temp vcf file per chromosome and genotype
  temp_vcf_per_chr = vcf_per_chr[names(vcf_per_chr) == ase$SNPchr[i]][[1]]
  temp_geno = temp_vcf_per_chr[temp_vcf_per_chr$position == ase$pos[i],colnames(temp_vcf_per_chr) == ase$sample_vcf[i]]

  if(length(temp_geno) > 0) {
    geno = strsplit(temp_geno, ':')[[1]][1]
    DP = strsplit(temp_geno, ':')[[1]][3] 
    GQ = strsplit(temp_geno, ':')[[1]][4]
    if(geno == ".|." | geno == "./.") ase$WGS_ratio[i] = 1 #missing genotype, so we call it reference
    
    if(geno == "0|0" | geno == "0/0") ase$WGS_ratio[i] = 1 #homozygous reference
    
    if(geno == "0|1" | geno == "0/1" | geno == "1|0" | geno == "1/0" ) {
      ase$WGS_ratio[i] = 0.5 #heterozygous
      ase$pvalue[i] = signif(binom.test(ase$ref[i],ase$total[i], p=0.5)$p.value,2)} #binomial test for ASE
      ase$WGS_GT[i] = geno
      ase$WGS_DP[i] = DP
      ase$WGS_GQ[i] = GQ

    if(geno == "1|1" | geno == "1/1") ase$WGS_ratio[i] = 0 #homo alternate
    
    } else ase$WGS_ratio[i] = NA #sample not present in the .vcf
  
  if (i %% 50000 == 0) print(paste0('Done ',i,' of ',nrow(ase), ' ~~~ Time is: ',Sys.time()))
}

#filter results.
ase_signif = ase[!is.na(ase$WGS_ratio),]
ase_signif = ase_signif[ase_signif$WGS_ratio == 0.5,]

print('Dimension of ase_signif: ')
print(dim(ase_signif))

print('Number of unique genes tested: ')
print(ase_signif %>% group_by(sampleID) %>% summarise(unique_genes = length(unique(ensemblID))))

print('Number of unique SNV tested: ')
print(ase_signif %>% group_by(sampleID) %>% summarise(unique_genes = length(ensemblID)))


ase_signif = ase_signif[ase_signif$pvalue<0.5,]
ase_signif$pvalue[ase_signif$pvalue<1e-50] = 1e-50

########################################
# 4. add geneID from hg38 annotation file.
########################################
#mapping genes
map <- select(org.Hs.eg.db, keys=keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
                keytype="ENTREZID", columns=c("ENSEMBL","SYMBOL"))
map = map[!is.na(map$ENSEMBL),]
colnames(map)[3] = 'geneID'

#add the gene length info
gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene,single.strand.genes.only=FALSE)
gene_locations = as.data.frame(gene_locations)
gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames))<6,]
map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
colnames(map)[c(2,3,5)] = c('ensemblID','geneID','chr')
map$chr = gsub('chr','',map$chr)


########################################
# 5. imprinting and X specific.
########################################
high_confidence_imprinted_genes = read.csv(params$imprinting, header =T)
ase_map = merge(ase,map,by = 'ensemblID',all.x = T)
ase_map = ase_map[!is.na(ase_map$WGS_ratio),]
ase_map = ase_map[!is.na(ase_map$geneID),]
ase_map_X = ase_map[ase_map$chr.x == 'X',]
ase_map_X = ase_map_X[ase_map_X$WGS_ratio == 0.5, ]
ase_map_X$Type = 'X'

ase_map_imprinted = ase_map[ase_map$geneID %in% high_confidence_imprinted_genes$Gene,]
ase_map_imprinted = ase_map_imprinted[ase_map_imprinted$WGS_ratio == 0.5, ]
ase_map_imprinted$Type = 'I'
ase_map_imprinted = rbind(ase_map_imprinted,ase_map_X)
#ase_map_imprinted = ase_map_imprinted[!is.na(ase_map_imprinted$WGS_ratio),]
#ase_map_imprinted = ase_map_imprinted[!is.na(ase_map_imprinted$geneID),]
ase_map_imprinted = ase_map_imprinted[,c(5,1,21,2,4,19:18,13:16,12,28)]
colnames(ase_map_imprinted)[4] = 'chr'

write.table(ase_map_imprinted, file.path(params$workdir,"gwImprinted.tsv"),sep = '\t',quote = F)


########################################
# 5. Save
########################################
signif_map = merge(ase_signif,map,by = 'ensemblID',all.x = T)
signif_map_ASE = signif_map[,c(5,1,21,2,4,19:18,13:16,12)]
#signif_map[,c(5,1,21,2,4,18:17, 19, 13:15,12)]
#signif_map_ASE = signif_map[,c(5,18,1,2,4,6,7,15,13,12)]
signif_map_ASE = signif_map_ASE[!is.na(signif_map_ASE$geneID),]
#signif_map_ASE$Imprinted = 0
#signif_map_ASE = rbind(signif_map_ASE,ase_map_imprinted)
colnames(signif_map_ASE)[4] = 'chr'

write.table(signif_map_ASE, file.path(params$workdir,"gwASE.tsv"),sep = '\t',quote = F)

