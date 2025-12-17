params = list(OUTRIDER = file.path(commandArgs(trailingOnly=TRUE)[1],"OUTRIDER/"))

#library
suppressMessages(suppressWarnings(library(OUTRIDER)))
suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))

#parameters
ncores = 8
register(MulticoreParam(ncores, ncores*2, progressbar = TRUE))

#candidate genes
candidates = read.csv(file.path(params$OUTRIDER,'../data/candidate_genes_3.txt'))

###from featurecount gene expression data
fc_exons_raw_ALL = read.table(file.path(params$OUTRIDER,'../featureCounts/fc_exons_raw.tsvALL'),sep = '\t',check.names = F)#;fc_exons_raw[,-c(1:5)] = round(fc_exons_raw[,-c(1:5)])
rownames(fc_exons_raw_ALL) = paste0(fc_exons_raw_ALL$geneID,"_",fc_exons_raw_ALL$ensemblID,"_",fc_exons_raw_ALL$transcriptID,"_",fc_exons_raw_ALL$exonID)
genes_counts = fc_exons_raw_ALL[,-c(1:5)]

#filter dataset to remove very low expression genes
rowmeans = rowMeans(genes_counts)
genes_counts = genes_counts[rowmeans>1,]

#OUTRIDER
ods <- OutriderDataSet(countData=genes_counts)
ods <- OUTRIDER(ods)
dim(ods)

#Result table
table = results(ods,padjCutoff=1)
table$combinedID = table$geneID
table$geneID = sapply(strsplit(table$combinedID,'_'),'[',1)
table$ensemblID = sapply(strsplit(table$combinedID,'_'),'[',2)
table$transcriptID = sapply(strsplit(table$combinedID,'_'),'[',3)
table$exonID = sapply(strsplit(table$combinedID,'_'),'[',4)
table$ensemblID_sampleID =  paste0(table$ensemblID,'_',table$sampleID)
table$pValue = signif(table$pValue,4)

#candidates only
#candidates$ensembl = sapply(strsplit(candidates$ensembl,'.',fixed = T),"[[",1)
candidates$ensembl_proband2 = apply(candidates[,colnames(candidates) %in% c('ensembl','proband2')],1,paste0,collapse ='_')
candidate_table_exon  = table[table$ensemblID_sampleID %in% candidates$ensembl_proband2,]

#write the results
write.table(candidate_table_exon,file.path(params$OUTRIDER,'candidates_perexons_OUTRIDER.tsv'),sep = '\t',quote = F)

#message
print(paste0('Done OUTRIDER per exon--- Time is:',Sys.time()) )
