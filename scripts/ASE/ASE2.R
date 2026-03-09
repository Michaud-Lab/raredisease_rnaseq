########################################
# ASE analysis from ASEReadCounter files
########################################

library(dplyr)
library(tidyr)
library(DESeq2)
library(purrr)
library(ggplot2)

########################################
# 1. Read and filter ASE overlap.txt file
########################################
args = commandArgs(trailingOnly=TRUE)
params = list(overlap = args[1])
params$workdir = dirname(params$overlap)

ase <- read.table(params$overlap,header=FALSE)
colnames(ase) <- c("SNPchr","SNPstart","SNPstop","sample",'ref','alt','chr','START','STOP','gene')

ase <- ase %>%
  mutate(total = ref + alt, ref_fraction = ref / total) %>%
  filter(alt >= 2) %>% #alternatively, you could also filter to check if it is truly a SNP from the .vcf file for that sample.
  filter(ref_fraction < 0.99) %>%
  filter(total >= 10) 

########################################
# 2. Binomial test per SNP per sample
########################################

ase$binom_pval = 0

for(i in 1:nrow(ase)){
    ase$binom_pval[i] = binom.test(ase$ref[i],ase$total[i], p=0.5)$p.value
    if(i%%100==0) print(paste0('Done ',i,' Time is: ', Sys.time()))
  }

########################################
# 3. Plotting
########################################
write.table(ase[ase$binom_pval<0.01,], file.path(params$workdir,"ase_sign_table.tsv"),header=FALSE)

samples = unique(ase$sample)
gplot_ASE = list()
#
for(i in 1:length(samples)){
  gplot_ASE[[i]] = ase[,c(4:6,10)] %>%
  filter(sample == samples[i]) %>%
  pivot_longer(2:3, names_to = "allele", values_to = "count") %>%
  ggplot(aes(gene,count, fill=allele)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
}  

#
for(i in 1:length(samples)){
  pdf(file.path(params$workdir,paste0('gplot_ASE_',samples[i],'.pdf')),width = 14,height = 12) 
  #print(i)
  print(gplot_ASE[[i]])
  dev.off()
}
