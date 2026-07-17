# =============================================================================
# featureCounts.R - Compute raw counts and TPM from featureCounts output
# =============================================================================
# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

params = list(workdir = args[1])
params$candidate_genes = args[2]
params$FCdir = file.path(args[1],'featureCounts/')
params$datadir = file.path(args[1],'data/')
params$ens_gene = args[3]
params$masterlog = args[4]
params$fc_exons = args[5]
params$fc_genes = args[6]
params$candidate_genes_LR = args[7]
params$candidate_genes_extra = args[8]
 
# -----------------------------------------------------------------------------
# 2. Load reference data
# -----------------------------------------------------------------------------
# Libraries
source(file.path(params$workdir,'scripts/rnaseq_helper_functions.R'))
load_install_library(c('readxl', 'tidyr', 'dplyr'))

# Ensembl - GeneID correspondence file
ensembl_geneid = read.table(params$ens_gene,header = TRUE)

# candidate genes
candidates_original = read.csv(params$candidate_genes)
candidates_extra = read.table(params$candidate_genes_extra,comment.char = "#",header = T ,sep = ',');candidates_extra[is.na(candidates_extra)] = ''
candidates_extra = candidates_extra[, colnames(candidates_original)]
candidates_extra$origin = 'extra gene added by data analyst'

candidate_genes_automated_list = list(NULL,NULL,NULL)
gwfiles = c(paste0(params$datadir,c('/gwFRASER.csv','/gw_genes_OUTRIDER.tsv','/gwASE.tsv')))

for(i in 1:3)
{
  if(file.exists(gwfiles[i])){candidate_genes_automated_list[[i]] = candidate_genes_automated(gwfile = gwfiles[i])}
  }

candidates_original$origin = 'original gene'
candidates_extra$origin = 'extra gene'

candidates = rbind(candidates_original,candidates_extra,candidate_genes_automated_list[[1]],candidate_genes_automated_list[[2]],candidate_genes_automated_list[[3]]) %>%
  distinct(geneID,ensembl, proband, .keep_all = TRUE)

write.csv(candidates,file.path(params$datadir, 'input/candidate_genes_ALL.csv'),quote = F, row.names = F)


# Clinical data
clinical = readxl::read_xlsx(params$masterlog, sheet = 'Suivi - RNAseq', skip = 1)
clinical$type = "Parent"
clinical$type[grepl(pattern = '_0[34]_', x = clinical$`Patient ID`)] = 'Proband' # Parent versus Proband
clinical = clinical[order(clinical$`Patient ID`), ] # Same order as the transcript expression data
clinical$age = as.numeric(clinical$`Âge (années)`)
clinical$age[clinical$`Âge (années)` == '0 (3 mois)'] = 0.25
clinical$age[clinical$`Âge (années)` == '0 (9 mois)'] = 0.75
clinical$PatientID = gsub('_PAX', '', clinical$`Patient ID`)
clinical$Notes[is.na(clinical$Notes)] = ''


# -----------------------------------------------------------------------------
# 3. Load featureCounts output
# -----------------------------------------------------------------------------
# featureCounts (per gene and per exon)
fc_exons = read.table(params$fc_exons,sep = '\t',header = TRUE,comment.char = '#',check.names = FALSE)
fc_genes = read.table(params$fc_genes,sep = '\t',header = TRUE,comment.char = '#',check.names = FALSE)

# featureCounts average per sample
colsum_genes_counts = colSums(fc_genes[,7:ncol(fc_genes)])
colmean_genes_counts = colsum_genes_counts / mean(colsum_genes_counts)
names(colmean_genes_counts) = gsub("^.*/", "",names(colmean_genes_counts))
names(colmean_genes_counts) = gsub('_sorted.bam','',names(colmean_genes_counts))
colmean_genes_counts = colmean_genes_counts[order(names(colmean_genes_counts))]
write.table(colmean_genes_counts,file.path(params$FCdir,'colmean_genes_counts.tsv'))

# -----------------------------------------------------------------------------
# 4. Exon-level counts and TPM
# -----------------------------------------------------------------------------
fc_exons$ensemblID = sapply(strsplit(fc_exons$Geneid,'_'),'[',1)
fc_exons$transcriptID = sapply(strsplit(fc_exons$Geneid,'_'),'[',2)
fc_exons$exonID = sapply(strsplit(fc_exons$Geneid,'_'),'[',3)

# sum it up per gene in order to find the transcript with the greatest number of reads aligned.
transcripts_summed = fc_exons %>% group_by(ensemblID,transcriptID) %>% summarise(across(c(7:(ncol(fc_exons)-4)), sum))

# keep the MANE instead of the longest gene.
MANE = read.table(paste0(params$FCdir,"/MANE.tsv"))
MANE = MANE %>% group_by(V2) %>% slice_max(V1, n = 1, with_ties = FALSE)
MANE = as.data.frame(MANE[,c(1,3)])

MANE = rbind(MANE, c(7, "ENST00000374555")) # MDS2 not in MANE (it is a lncRNA) but exon data exists ALTERNATIVELY IF ONE OF THE GENE IS NOT IN THE MANE, WE CAN POTENTIALLY FIND IT USING THE FUNCTION GENE_ANNOTATION...       
MT = fc_exons$transcriptID[fc_exons$Chr == 'MT']
MANE = rbind(MANE, data.frame(V1 = rep(1, length(MT)), V3 = MT)) # add mitochondrial genes

fc_exons = fc_exons[fc_exons$transcriptID %in% MANE[,2],]

# clean up
fc_exons = merge(fc_exons,ensembl_geneid,by.x = 'ensemblID',by.y = 'gene_id')
ncol = ncol(fc_exons)
fc_exons_raw_ALL = fc_exons[,c(ncol,1,(ncol-2),(ncol-1),7:(ncol-3))]
colnames(fc_exons_raw_ALL)[1] = 'geneID'
colnames(fc_exons_raw_ALL)[-c(1:5)] = sapply(lapply(strsplit(colnames(fc_exons_raw_ALL),'/'),'['),tail,1)[-c(1:5)]
colnames(fc_exons_raw_ALL)[-c(1:5)] = gsub('_sorted.bam','',colnames(fc_exons_raw_ALL)[-c(1:5)])

#remove hemoglobing gens
hemo_ensembl = ensembl_geneid$gene_id[ensembl_geneid$gene_name %in% c('HBB','HBA1','HBA2','HBD')]

fc_exons_tpm = fc_exons_raw_ALL[!fc_exons_raw_ALL$ensemblID %in% hemo_ensembl, ]
fc_exons_tpm[,-c(1:5)]  = fc_exons_tpm[,-c(1:5)]  / fc_exons_tpm[,5] * 1000
fc_exons_tpm[,-c(1:5)] = lapply(fc_exons_tpm[,-c(1:5)] , function(x) x/sum(x)* 1000000)
fc_exons_tpm[,-c(1:5)] =  round(fc_exons_tpm[,-c(1:5)],2)

fc_exons_raw_ALL = fc_exons_raw_ALL[,colnames(fc_exons_raw_ALL) %in% c('geneID','ensemblID','exonID','transcriptID','Length',clinical$`Patient ID`)]
fc_exons_tpm = fc_exons_tpm[,colnames(fc_exons_tpm) %in% c('geneID','ensemblID','exonID','transcriptID','Length',clinical$`Patient ID`)]

colnames(fc_exons_raw_ALL) = gsub('_PAX','',colnames(fc_exons_raw_ALL))
colnames(fc_exons_tpm) = gsub('_PAX','',colnames(fc_exons_tpm))

# ggplot data formatting
fc_exons_tpm_ggplot = merge(fc_exons_tpm,candidates[,c(1,3)],sort = FALSE)
fc_exons_tpm_ggplot$proband = gsub('_PAX','',fc_exons_tpm_ggplot$proband)
fc_exons_tpm_ggplot = fc_exons_tpm_ggplot %>% pivot_longer(cols = c(6:(ncol(fc_exons_tpm_ggplot)-1)), names_to = 'PatientID',values_to = 'expression')
fc_exons_tpm_ggplot = merge(fc_exons_tpm_ggplot,clinical[,colnames(clinical) %in% c('PatientID','Sexe','type','age')])

# Filter: probands only, then candidate genes only
proband_ids = gsub('_PAX', '', clinical$`Patient ID`[clinical$type == 'Proband'])
id_cols = c('geneID', 'ensemblID', 'transcriptID', 'exonID', 'Length')
fc_exons_tpm = fc_exons_tpm[, colnames(fc_exons_tpm) %in% c(id_cols, proband_ids)]
fc_exons_raw_ALL = fc_exons_raw_ALL[, colnames(fc_exons_raw_ALL) %in% c(id_cols, proband_ids)]
fc_exons_tpm = fc_exons_tpm[fc_exons_tpm$geneID %in% candidates$geneID, ]
fc_exons_raw = fc_exons_raw_ALL[fc_exons_raw_ALL$geneID %in% candidates$geneID, ]

# -----------------------------------------------------------------------------
# 5. Gene-level counts and TPM
# -----------------------------------------------------------------------------
fc_genes = merge(fc_genes,ensembl_geneid,by.x = 'Geneid',by.y = 'gene_id')
ncol = ncol(fc_genes)
fc_genes_raw = fc_genes[,c(ncol,1,6:(ncol-1))]
colnames(fc_genes_raw)[1:2] = c('geneID','ensemblID')
colnames(fc_genes_raw) =  gsub("^.*/", "",colnames(fc_genes_raw))
colnames(fc_genes_raw) = gsub('_sorted.bam','',colnames(fc_genes_raw))

fc_genes_tpm = fc_genes_raw
fc_genes_tpm[,-c(1:3)]  = fc_genes_tpm[,-c(1:3)]  / fc_genes_tpm[,3] * 1000
fc_genes_tpm[,-c(1:3)] = lapply(fc_genes_tpm[,-c(1:3)] , function(x) x/sum(x) * 1000000)
fc_genes_tpm[,-c(1:3)] =  round(fc_genes_tpm[,-c(1:3)],2)

proband_patient_ids = clinical$`Patient ID`[clinical$type == 'Proband']
fc_genes_tpm = fc_genes_tpm[, colnames(fc_genes_tpm) %in% c('geneID', 'ensemblID', 'exonID', 'Length', proband_patient_ids)]
fc_genes_raw = fc_genes_raw[, colnames(fc_genes_raw) %in% c('geneID', 'ensemblID', 'exonID', 'Length', proband_patient_ids)]

colnames(fc_genes_raw) = gsub('_PAX','',colnames(fc_genes_raw))
colnames(fc_genes_tpm) = gsub('_PAX','',colnames(fc_genes_tpm))

write.table(fc_genes_raw,file.path(params$FCdir,'fc_genes_raw_ALL.tsv'),sep = '\t',quote = FALSE)

fc_genes_tpm = fc_genes_tpm[fc_genes_tpm$geneID %in% candidates$geneID, ]
fc_genes_raw = fc_genes_raw[fc_genes_raw$geneID %in% candidates$geneID, ]

# -----------------------------------------------------------------------------
# 6. Gene model annotation and save outputs
# -----------------------------------------------------------------------------
gene_annotations = gene_annotation(unique_transcript_id = unique(fc_exons_raw$transcriptID),candidates = candidates)
save(gene_annotations, file= file.path(params$FCdir,"gene_annotations.rda"))

write.table(fc_genes_raw,file.path(params$FCdir,'fc_genes_raw.tsv'),sep = '\t',quote = FALSE)
write.table(fc_genes_tpm,file.path(params$FCdir,'fc_genes_tpm.tsv'),sep = '\t',quote = FALSE)

write.table(fc_exons_raw_ALL,file.path(params$FCdir,'fc_exons_raw.tsv_ALL'),sep = '\t',quote = FALSE)
write.table(fc_exons_raw,file.path(params$FCdir,'fc_exons_raw.tsv'),sep = '\t',quote = FALSE)
write.table(fc_exons_tpm,file.path(params$FCdir,'fc_exons_tpm.tsv'),sep = '\t',quote = FALSE)
write.table(fc_exons_tpm_ggplot,file.path(params$FCdir,'fc_exons_tpm_ggplot.tsv'), sep = '\t',quote = FALSE)
write.table(clinical,file.path(params$FCdir,'clinical.tsv'), sep = '\t',quote = TRUE)

print(paste0('Done write table --- ', Sys.time()))
