########################################
# ASE analysis from ASEReadCounter files
########################################

library(dplyr)
library(tidyr)
library(DESeq2)
library(purrr)
library(tibble)

########################################
# 1. Read and filter ASE overlap.txt file
########################################

ase <- read.table("/home/renaut/scratch/raredisease_rnaseq/ASE/overlaps.tsv",header=FALSE)
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
ase[,c(4:6,10)] %>%
  filter(sample == 'HSJ_001_03_PAX') %>%
  pivot_longer(2:3, names_to = "allele", values_to = "count")%>%
  ggplot(aes(gene,count, fill=allele)) +
  geom_bar(stat="identity") +
 # facet_wrap(~gene) +
  theme_bw()
