########################################
# ASE analysis from ASEReadCounter files
########################################


library(dplyr)
#library(tidyr)
library(DESeq2)
library(vcfR)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


########################################
# 1. Read and filter ASE overlap.txt file
########################################
args = commandArgs(trailingOnly=TRUE)
params = list(overlap = args[1])
params$workdir = dirname(params$overlap)
ase <- read.table("../overlaps.tsv",header=FALSE)
colnames(ase) <- c("SNPchr","SNPstart","SNPstop","sample",'ref','alt','chr','START','STOP','ensemblID')

ase$sample_vcf = gsub('_PAX','',ase$sample)
ase$SNPstop = as.character(ase$SNPstop)
ase$binom_pval = NA

ase <- ase %>%
  mutate(total = ref + alt, ref_fraction_observed = ref / total)

#get the actual chromosome names from the ase file
aes_rle = rle(ase$SNPchr)
aes_rle = data.frame(lengths = aes_rle$lengths, values = aes_rle$values, c = 1:length(aes_rle$lengths))


########################################
# 2. Read and filter ASE based on .vcf to check the expected count.
########################################
vcf <- read.vcfR("allchr_biallelic_sites.vcf.gz")
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
ase$ref_fraction_expected = 1

for(i in 1:nrow(ase)){
  #temp vcf file per chromosome and genotype
  temp_vcf_per_chr = vcf_per_chr[names(vcf_per_chr) == ase$SNPchr[i]][[1]]
  temp_geno = temp_vcf_per_chr[temp_vcf_per_chr$position == ase$SNPstop[i],colnames(temp_vcf_per_chr) == ase$sample_vcf[i]]

  if(length(temp_geno) > 0) {
    geno = strsplit(temp_geno, ':')[[1]][1]
    if(geno == ".|." | geno == "./.") ase$ref_fraction_expected[i] = 1 #missing so ref
    
    if(geno == "0|0" | geno == "0/0") ase$ref_fraction_expected[i] = 1 #homo ref
    
    if(geno == "0|1" | geno == "0/1" | geno == "1|0" | geno == "1/0" ) {
      ase$ref_fraction_expected[i] = 0.5 #het
      ase$binom_pval[i] = binom.test(ase$ref[i],ase$total[i], p=0.5)$p.value} #binomial test for ASE

    if(geno == "1|1" | geno == "1/1") ase$ref_fraction_expected[i] = 0 #homo alt
    
    } else ase$ref_fraction_expected[i] = NA #sample not present in the .vcf
  
  if (i %% 10000 == 0) print(paste0('Done ',i,' of ',nrow(ase), ' ~~~ Time is: ',Sys.time()))
}

#filter results.
ase_signif = ase[!is.na(ase$ref_fraction_expected),]
ase_signif = ase_signif[ase_signif$ref_fraction_expected == 0.5,]
ase_signif = ase_signif[ase_signif$binom_pval<0.05,]

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
# 5. Save
########################################
signif_map = merge(ase_signif,map,by = 'ensemblID',all.x = T)
write.table(signif_map, file.path(params$workdir,"ase_sign_table.tsv"))


########################################
# 6. Plot
########################################
samples = unique(signif_map$sample)
gplot_ASE = list()
#g
for(i in 1:length(samples)){
  gplot_ASE[[i]] = ase %>%
  filter(sample == samples[i]) %>%
  filter(binom_pval < 0.01) %>%
  pivot_longer(5:6, names_to = "allele", values_to = "count")%>%
  ggplot(aes(gene,count, fill=allele)) +
  geom_bar(stat="identity") +
  ggtitle(paste0('ASE (p-val<0.01) ~~~ ', samples[i])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
}  

#
for(i in 1:length(samples)){
  pdf(paste0('ase/gplot_ASE_',samples[i],'.pdf'),width = 14,height = 12) 
  print(gplot_ASE[[i]])
  dev.off()
}



