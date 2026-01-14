i = as.numeric(commandArgs(trailingOnly=TRUE)[1])
params = list(FRASER = file.path(commandArgs(trailingOnly=TRUE)[2],"FRASER/"))

options(scipen = 999)

#Load required packages
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggtranscript)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(FRASER)))
suppressMessages(suppressWarnings(library(tidyr)))


#fc_exons and genes
fc_exons = read.table(file.path(params$FRASER,'../featureCounts/fc_exons_raw.tsv'),check.names = F)
unique_transcript_id = unique(fc_exons$transcriptID)

#genes_counts = read.table(file.path(params$FRASER,'../featureCounts/feature_counts_pergene.txt'),sep = '\t',header = T,comment.char = '#',check.names = F)
#colnames(genes_counts) = gsub("^.*/", "",colnames(genes_counts))
#colnames(genes_counts) = gsub('_sorted.bam','',colnames(genes_counts))
#colsum_genes_counts = colSums(genes_counts[,7:ncol(genes_counts)])
#colmean_genes_counts = colsum_genes_counts / mean(colsum_genes_counts)
colmean_genes_counts = read.table(file.path(params$FRASER,'../featureCounts/colmean_genes_counts.tsv'))

#candidates
candidates = read.csv(file.path(params$FRASER,'../data/candidate_genes_3.txt'))

#prepare the bam subsetting BASH script
chr = candidates$chromosome[i]
start = candidates$start[i] - 100000; if(start <= 0) start = 1
stop = candidates$stop[i] + 100000
geneID = candidates$geneID[i]
proband = candidates$proband[i]

if(geneID != "") {

command = paste('./fraser.sh', chr, start, stop, geneID,proband)
out_dir = paste0(params$FRASER,'bams_subset/gene',candidates$geneID[i],'_chr',chr,'_',start,'_',stop)
FRASER_test_output = file.path(params$FRASER,paste0('results/output_FRASER_',candidates$proband[i]))
system(command)

#generate a coverage table with samtools for the position of interest
#samtools_depth_cmd = paste0('samtools depth -H -a ',out_dir,'/*sorted_chrN.bam -r chr',chr,':',start,'-',stop,' >',out_dir,'/gene.depth')
#samtools_depth_cmd = paste0("samtools depth -H -a ",out_dir,"/*sorted_chrN.bam -r chr",chr,":",start,"-",stop," | awk 'NR % 5 == 1' >",out_dir,"/gene_",candidates$geneID[i],"_",candidates$proband[i],"_depth5.csv")

#print(samtools_depth_cmd)
#system(samtools_depth_cmd)
    
#FRASER analysis
dir.create(params$FRASER,showWarnings = T)
dir.create(file.path(params$FRASER,'results'),showWarnings = T)
dir.create(FRASER_test_output,showWarnings = T)

sampleTable = data.table(data.frame(sampleID =  gsub('_sorted_chrN.bam','',list.files(paste0(params$FRASER,'bams_subset/gene',geneID,'_chr',chr,'_',start,'_',stop),pattern = '*bam$')),
                         bamFile = list.files(paste0(params$FRASER,'bams_subset/gene',geneID,'_chr',chr,'_',start,'_',stop),pattern = '*bam$',full.names = F),
                         group = 1,
                         pairedEnd = TRUE))

sampleTable$group = 1:nrow(sampleTable)
sampleTable$bamFile = paste0(params$FRASER,'bams_subset/gene',geneID,'_chr',chr,'_',start,'_',stop,'/',sampleTable$bamFile) 

#probands only
probands = c(1:nrow(sampleTable))[grepl('_0[34]_',sampleTable$bamFile) | grepl('LC_',sampleTable$bamFile)| grepl('F0',sampleTable$bamFile)]
settings <- FraserDataSet(colData=sampleTable[probands,], workingDir=FRASER_test_output)

fds <- countRNAData(settings)
fds <- calculatePSIValues(fds)
fds <- annotateRanges(fds,GRCh=38)
fds <- FRASER(fds, q=c(jaccard=2))

res <- results(fds, all=TRUE, padjCutoff=NA, deltaPsiCutoff=NA)
sort(unique(res$hgncSymbol))

res_dt <- as.data.table(res)
res_dt <- res_dt[sampleID == candidates$proband[i],]

#load more annotation packages
suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
txdb_chr <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene,paste0('chr',c(1:22,'X','Y','M')), pruning.mode = "coarse")

temp = c(1:nrow(as.data.frame(res_dt)))[res_dt$hgncSymbol == candidates$geneID[i]]; candidate_gene_localisation = temp[!is.na(temp)][1]
control_samples = sampleTable$sampleID[grepl('_03_',sampleTable$bamFile) | grepl('LC_',sampleTable$bamFile)]
control_samples = control_samples[control_samples != candidates$proband[i]][1:5]

#save it
res_dt_candidate_gene = res_dt[temp[!is.na(temp)],]
res_dt_candidate_gene$mean = (res_dt_candidate_gene$start + res_dt_candidate_gene$end) / 2
res_dt_candidate_gene$minuslogpval = -log(res_dt_candidate_gene$pValue,10)
write.csv(res_dt_candidate_gene,paste0(out_dir,"/gene_",candidates$geneID[i],"_",candidates$proband[i],"_res_dt_candidate_gene.csv"))


#fail-safe in case the plotting does not work.
png(file.path(FRASER_test_output,paste0('gene_',candidates$geneID[i],'_',candidates$proband[i],'_sashimi.png')),width = 1000,height = 600)
plot(0, main = 'failed test')
dev.off()

png(file.path(FRASER_test_output,paste0('gene_',candidates$geneID[i],'_',candidates$proband[i],'_sashimi.png')),width = 1000 + length(temp[!is.na(temp)])*5 ,height = 1200)
try(plotBamCoverageFromResultTable(fds, result=res_dt[candidate_gene_localisation,], show_full_gene=TRUE,txdb=txdb_chr,
orgDb=org.Hs.eg.db,
control_samples = control_samples,
splicegraph_labels="id",
mar = c(1,10,0.1,1))
,silent = T)
dev.off()

}

print(paste0('Done sample ~~~ ', i))

