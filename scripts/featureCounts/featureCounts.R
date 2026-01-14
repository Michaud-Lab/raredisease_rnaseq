#Libraries
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))

#args and params 
args=commandArgs(trailingOnly=TRUE)

params = list(workdir = args[1])
params$candidate_genes = args[2]
params$FCdir = file.path(params$workdir,'/featureCounts/')
params$datadir = file.path(params$workdir,'data/')
params$ens_gene = args[3]
params$masterlog = args[4]
params$fc_exons = args[5]
params$fc_genes = args[6]

#Ensembl - GeneID correspondance file
ensembl_geneid = read.table(params$ens_gene,header = T)

#candidate genes
candidates = read.csv(params$candidate_genes)

#Clinical data (from Maude)
clinical = readxl::read_xlsx(params$ensembl_geneid, sheet = 'Suivi - RNAseq',skip = 1)
clinical$type =  "Parent"; clinical$type[!is.na(clinical$Mutation)] = 'Proband' #Parent versus Proband
clinical = clinical[order(clinical$`Patient ID`),] #Same order as the transript expression data.
clinical$age = as.numeric(clinical$`Âge (années)`); clinical$age[clinical$`Âge (années)` == '0 (3 mois)'] = 0.25; clinical$age[clinical$`Âge (années)` == '0 (9 mois)'] = 0.75 # Prettify
clinical$PatientID = gsub('_PAX','',clinical$`Patient ID`)

#featureCounts (per gene and per exons)
fc_exons = read.table(params$fc_exons,sep = '\t',header = T,comment.char = '#',check.names = F)
fc_genes = read.table(params$fc_genes,sep = '\t',header = T,comment.char = '#',check.names = F)

#featureCounts average per sample
colsum_genes_counts = colSums(fc_genes[,7:ncol(fc_genes)])
colmean_genes_counts = colsum_genes_counts / mean(colsum_genes_counts)
names(colmean_genes_counts) = gsub("^.*/", "",names(colmean_genes_counts))
names(colmean_genes_counts) = gsub('_sorted.bam','',names(colmean_genes_counts)) 
write.table(colmean_genes_counts,file.path(params$FCdir,'colmean_genes_counts.tsv'))

#EXONS rawcounts and TPM both for genes and exons
fc_exons$ensemblID = sapply(strsplit(fc_exons$Geneid,'_'),'[',1)
fc_exons$transcriptID = sapply(strsplit(fc_exons$Geneid,'_'),'[',2)
fc_exons$exonID = sapply(strsplit(fc_exons$Geneid,'_'),'[',3)

#sum it up per gene in order to find the transcript with the greatest number of reads aligned.
transcripts_summed = fc_exons %>% group_by(ensemblID,transcriptID) %>% summarise(across(c(7:(ncol(fc_exons)-4)), sum))

#keep the MANE instead of the longest gene. 
MANE = read.table(paste0(params$FCdir,"/MANE.tsv"))
fc_exons = fc_exons[fc_exons$transcriptID %in% MANE[,2],]

#clean up
fc_exons = merge(fc_exons,ensembl_geneid,by.x = 'ensemblID',by.y = 'gene_id')
ncol = ncol(fc_exons)
fc_exons_raw_ALL = fc_exons[,c(ncol,1,(ncol-2),(ncol-1),7:(ncol-3))]
colnames(fc_exons_raw_ALL)[1] = 'geneID'
colnames(fc_exons_raw_ALL)[-c(1:5)] = sapply(lapply(strsplit(colnames(fc_exons_raw_ALL),'/'),'['),tail,1)[-c(1:5)]
colnames(fc_exons_raw_ALL)[-c(1:5)] = gsub('_sorted.bam','',colnames(fc_exons_raw_ALL)[-c(1:5)])

fc_exons_tpm = fc_exons_raw_ALL
fc_exons_tpm[,-c(1:5)]  = fc_exons_tpm[,-c(1:5)]  / fc_exons_tpm[,5] * 1000
fc_exons_tpm[,-c(1:5)]  <- lapply(fc_exons_tpm[,-c(1:5)] , function(x) x/sum(x)* 1000000)
fc_exons_tpm[,-c(1:5)] =  round(fc_exons_tpm[,-c(1:5)],2)

fc_exons_raw_ALL = fc_exons_raw_ALL[,colnames(fc_exons_raw_ALL) %in% c('geneID','ensemblID','exonID','transcriptID','Length',clinical$`Patient ID`)]  
fc_exons_tpm = fc_exons_tpm[,colnames(fc_exons_tpm) %in% c('geneID','ensemblID','exonID','transcriptID','Length',clinical$`Patient ID`)] 

colnames(fc_exons_raw_ALL) = gsub('_PAX','',colnames(fc_exons_raw_ALL))
colnames(fc_exons_tpm) = gsub('_PAX','',colnames(fc_exons_tpm))

#ggplot data formatting
fc_exons_tpm_ggplot = merge(fc_exons_tpm,candidates[,c(1,3)],sort = F)
fc_exons_tpm_ggplot$proband = gsub('_PAX','',fc_exons_tpm_ggplot$proband)
fc_exons_tpm_ggplot = fc_exons_tpm_ggplot %>% pivot_longer(cols = c(6:(ncol(fc_exons_tpm_ggplot)-1)), names_to = 'PatientID',values_to = 'expression')
fc_exons_tpm_ggplot = merge(fc_exons_tpm_ggplot,clinical[,colnames(clinical) %in% c('PatientID','Sexe','type','age')])  

#filter out parents and keep only probands
fc_exons_tpm = fc_exons_tpm[,colnames(fc_exons_tpm) %in% c('geneID','ensemblID','transcriptID','exonID','Length',gsub('_PAX','',clinical$`Patient ID`[clinical$type =='Proband']))] #filter patients
fc_exons_raw_ALL = fc_exons_raw_ALL[,colnames(fc_exons_raw_ALL) %in% c('geneID','ensemblID','transcriptID','exonID','Length',gsub('_PAX','',clinical$`Patient ID`[clinical$type =='Proband']))] #filter patients
fc_exons_tpm = fc_exons_tpm[fc_exons_tpm$geneID %in% candidates$geneID, ] #filter ONLY genes of interest
fc_exons_raw = fc_exons_raw_ALL[fc_exons_raw_ALL$geneID %in% candidates$geneID, ] #filter ONLY genes of interest

#GENES
fc_genes = merge(fc_genes,ensembl_geneid,by.x = 'Geneid',by.y = 'gene_id')
ncol = ncol(fc_genes)
fc_genes_raw = fc_genes[,c(ncol,1,6:(ncol-1))]
colnames(fc_genes_raw)[1:2] = c('geneID','ensemblID')
colnames(fc_genes_raw) =  gsub("^.*/", "",colnames(fc_genes_raw))
colnames(fc_genes_raw) = gsub('_sorted.bam','',colnames(fc_genes_raw))

fc_genes_tpm = fc_genes_raw
fc_genes_tpm[,-c(1:3)]  = fc_genes_tpm[,-c(1:3)]  / fc_genes_tpm[,3] * 1000
fc_genes_tpm[,-c(1:3)]  <- lapply(fc_genes_tpm[,-c(1:3)] , function(x) x/sum(x) * 1000000)
fc_genes_tpm[,-c(1:3)] =  round(fc_genes_tpm[,-c(1:3)],2)

fc_genes_tpm = fc_genes_tpm[,colnames(fc_genes_tpm) %in% c('geneID','ensemblID','exonID','Length',clinical$`Patient ID`[clinical$type =='Proband'])] #filter ONLY patients of interest
fc_genes_raw = fc_genes_raw[,colnames(fc_genes_raw) %in% c('geneID','ensemblID','exonID','Length',clinical$`Patient ID`[clinical$type =='Proband'])] #filter ONLY patients of interest

fc_genes_tpm = fc_genes_tpm[fc_genes_tpm$geneID %in% candidates$geneID, ] #filter ONLY genes of interest
fc_genes_raw = fc_genes_raw[fc_genes_raw$geneID %in% candidates$geneID, ] #filter ONLY genes of interest

colnames(fc_genes_raw) = gsub('_PAX','',colnames(fc_genes_raw))
colnames(fc_genes_tpm) = gsub('_PAX','',colnames(fc_genes_tpm))


#gene annotation
source(file.path(params$workdir,'scripts/featureCounts/rnaseq_helper_functions.R'))
gene_annotations = gene_annotation(unique_transcript_id = unique(fc_exons_raw$transcriptID),candidates = candidates)
save(gene_annotations, file= file.path(params$FCdir,"gene_annotations.rda"))

write.table(fc_genes_raw,file.path(params$FCdir,'fc_genes_raw.tsv'),sep = '\t',quote = F)
write.table(fc_genes_tpm,file.path(params$FCdir,'fc_genes_tpm.tsv'),sep = '\t',quote = F)

write.table(fc_exons_raw_ALL,file.path(params$FCdir,'fc_exons_raw.tsv_ALL'),sep = '\t',quote = F)
write.table(fc_exons_raw,file.path(params$FCdir,'fc_exons_raw.tsv'),sep = '\t',quote = F)
write.table(fc_exons_tpm,file.path(params$FCdir,'fc_exons_tpm.tsv'),sep = '\t',quote = F)
write.table(fc_exons_tpm_ggplot,file.path(params$FCdir,'fc_exons_tpm_ggplot.tsv'), sep = '\t',quote = F)
write.table(clinical,file.path(params$FCdir,'clinical.tsv'), sep = '\t',quote = T)

print(paste0('Done write table --- ', Sys.time()))


