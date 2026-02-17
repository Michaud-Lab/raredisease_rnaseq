#param arguments
args = commandArgs(trailingOnly=TRUE)
params = list(OUTRIDER = file.path(args[1],"OUTRIDER/"))
params$candidate_genes = args[2]
params$fc_pergene = file.path(args[1],args[3])
params$fc_perexon = file.path(args[1],args[4])

#library
suppressMessages(suppressWarnings(library(OUTRIDER)))
suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))

#parameters
ncores = 8
register(MulticoreParam(ncores, ncores*2, progressbar = TRUE))

#candidate genes
candidates = read.csv(params$candidate_genes)
#candidates$ensembl = sapply(strsplit(candidates$ensembl,'.',fixed = T),"[[",1)
candidates$ensembl_proband = apply(candidates[,colnames(candidates) %in% c('ensembl','proband')],1,paste0,collapse ='_')
#candidates$ensembl_proband2 = sapply(strsplit(candidates$ensembl_proband2,'@'),'[[',1)

##mapping genes
map <- select(org.Hs.eg.db, keys=keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
keytype="ENTREZID", columns=c("ENSEMBL","SYMBOL"))
map = map[!is.na(map$ENSEMBL),]
colnames(map)[3] = 'geneID'

#add the gene length info
gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene,single.strand.genes.only=FALSE)
gene_locations = as.data.table(gene_locations)
gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames))<6,]
map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
colnames(map)[c(2,3,5)] = c('ensemblID','geneID','chr')
map$pos = (map$start+map$end)/2
map$chr = gsub('chr','',map$chr)

###from featurecount gene expression data
genes_counts = read.table(params$fc_pergene,sep = '\t',header = T,comment.char = '#',check.names = F)
rownames(genes_counts) = genes_counts[,1]
genes_counts = genes_counts[,-c(1:6)]
genes_counts = round(genes_counts)
colnames(genes_counts) = sapply(lapply(strsplit(colnames(genes_counts),'/'),'['),tail,1)
colnames(genes_counts) = gsub('_sorted.bam','',colnames(genes_counts))

#probands only (and the LC and F0 from Philippe Campeau, and 04 because that is a twin of a 03)
probands = colnames(genes_counts)[grepl('_0[34]_',colnames(genes_counts)) | grepl('LC_',colnames(genes_counts)) | grepl('F0',colnames(genes_counts)) ]
genes_counts = genes_counts[,colnames(genes_counts) %in% c('gene_id',probands)]

#OUTRIDER
ods <- OutriderDataSet(countData=genes_counts)
ods <- filterExpression(ods, TxDb.Hsapiens.UCSC.hg38.knownGene, mapping=map[,1:2],filterGenes=T, savefpkm=TRUE, fpkmCutoff = 1)
ods <- filterExpression(ods)
ods <- OUTRIDER(ods)
dim(ods)

#Result table
table = as.data.frame(results(ods,padjCutoff=1))
table$pValue = signif(table$pValue,4)
#table$sampleID = gsub('_PAX','',table$sampleID)

significant_table = table[table$pValue<0.05,]
significant_table$padjust = signif(significant_table$padjust,4)
colnames(significant_table)[1] = 'ensemblID' 
significant_table = merge(significant_table,map,by = 'ensemblID')
table$ensemblID_sampleID =  paste0(table$geneID,'_',table$sampleID)

#candidate only
candidate_table  = table[table$ensemblID_sampleID %in% candidates$ensembl_proband,]
colnames(candidate_table)[1] = 'ensemblID'
candidate_table = merge(candidate_table,map,by = 'ensemblID')

#write the results
colnames_ALL = c("sampleID","geneID","ensemblID","chr","start","end","pos","width","pValue","padjust","zScore","l2fc","rawcounts","meanRawcounts","normcounts","meanCorrected")
colnames_candidate_genes = c("sampleID","geneID","ensemblID","pValue","zScore","l2fc","rawcounts","meanRawcounts","normcounts","meanCorrected")
write.table(significant_table[,colnames_ALL],file.path(params$OUTRIDER,'results_OUTRIDER.tsv'),sep = '\t',quote = F)
write.table(candidate_table[,colnames_candidate_genes],file.path(params$OUTRIDER,'candidates_OUTRIDER.tsv'),sep = '\t',quote = F)

#message
print(paste0('Done OUTRIDER --- Time is: ',Sys.time()) )

####PER EXON
#from featurecount gene expression data
fc_exons_raw_ALL = read.table(params$fc_perexon,sep = '\t',check.names = F)
colnames(fc_exons_raw_ALL)[grep('HSJ',colnames(fc_exons_raw_ALL))] = paste0(colnames(fc_exons_raw_ALL)[grep('HSJ',colnames(fc_exons_raw_ALL))],'_PAX')#this is ugly and I should fix this. It's there because HSJ samples have two names, one with '_PAX', one without.
rownames(fc_exons_raw_ALL) = paste0(fc_exons_raw_ALL$geneID,"_",fc_exons_raw_ALL$ensemblID,"_",fc_exons_raw_ALL$transcriptID,"_",fc_exons_raw_ALL$exonID)

#filter dataset to remove very low expression genes
genes_counts = fc_exons_raw_ALL[,-c(1:5)]
genes_counts = genes_counts[rowMeans(genes_counts) > 1,]
#genes_counts = genes_counts[1:10000,]

#OUTRIDER
ods <- OutriderDataSet(countData=genes_counts)
ods <- OUTRIDER(ods)
dim(ods)

#Result table
table = as.data.frame(results(ods,padjCutoff=1))
table$combinedID = table$geneID
table$geneID = sapply(strsplit(table$combinedID,'_'),'[',1)
table$ensemblID = sapply(strsplit(table$combinedID,'_'),'[',2)
table$transcriptID = sapply(strsplit(table$combinedID,'_'),'[',3)
table$exonID = sapply(strsplit(table$combinedID,'_'),'[',4)
table$ensemblID_sampleID =  paste0(table$ensemblID,'_',table$sampleID)
table$pValue = signif(table$pValue,4)

#candidates exons only
colnames_candidate_exon = c("sampleID","geneID","ensemblID","transcriptID","exonID","pValue","zScore","l2fc","rawcounts","meanRawcounts","normcounts","meanCorrected")
candidate_table_exon = table[table$ensemblID_sampleID %in% candidates$ensembl_proband,]
candidate_table_exon = candidate_table_exon[,colnames_candidate_exon]

#all significant exons
colnames_significant_exon = c("sampleID","geneID","ensemblID","transcriptID","exonID","chr","start","end","pos","width","pValue","padjust","zScore","l2fc","rawcounts","meanRawcounts","normcounts","meanCorrected")
table = table[table$pValue < 0.01,]
significant_table_exon = merge(table,map[,-3],by = 'ensemblID')
significant_table_exon = significant_table_exon[,colnames_significant_exon]

#write the results
write.table(significant_table_exon,file.path(params$OUTRIDER,'exons_OUTRIDER.tsv'),sep = '\t',quote = F)
write.table(candidate_table_exon,file.path(params$OUTRIDER,'candidates_perexons_OUTRIDER.tsv'),sep = '\t',quote = F)

#message
print(paste0('Done OUTRIDER per exon --- Time is: ',Sys.time()) )




