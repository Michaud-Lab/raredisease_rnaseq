
args = commandArgs(trailingOnly=TRUE)

params = list(FRASER = args[1])
params$workdir = args[2]
i = as.numeric(args[3])
options(scipen = 999)


#Load required packages
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(FRASER)))
suppressMessages(suppressWarnings(library(tidyr)))

#CPUs
ncores = 8
register(MulticoreParam(ncores, ncores*2, progressbar = TRUE))


if(file.exists(file.path(params$workdir,'../raredisease_rnaseq/data/candidate_genes_LR.txt')) !=T) {
  #the samples we processed with LR
  samples = data.frame(long_read_code = c("bc01","bc02","bc03",'bc04'), 
                       proband = c('HSJ_003_03_PAX','HSJ_017_03_PAX','HSJ_018_03_PAX','HSJ_036_03_PAX'),
                       candidate_gene = c('KDM6B','SAV1','SAV1','COMMD9'))
  
  #candidates
  candidates_SR = read.csv(file.path(params$workdir,'../data/candidate_genes_3.txt'))
  candidates_LR = candidates_SR[candidates_SR$proband %in% samples$proband,]
  candidates_LR = merge(candidates_LR,samples[,1:2])
  
  #save the candidate genes with the long read data for the Shiny dashboard
  write.csv(candidates_LR,file.path(params$workdir,'../data/candidate_genes_LR.txt'),row.names = F, quote = T)
} else candidates_LR = read.csv(file.path(params$workdir,'../data/candidate_genes_LR.txt'), header = T)

params$output = paste0(args[1],'/results/output_FRASER_',candidates_LR$long_read_code[i])
dir.create(params$output,showWarnings = T)

#prepare the bam subsetting script
chr = candidates_LR$chromosome[i]
start = candidates_LR$start[i] - 5000; if(start <= 0) start = 1
stop = candidates_LR$stop[i] + 5000
geneID = candidates_LR$geneID[i]
proband = candidates_LR$long_read_code[i]
  
#bam_in = paste0("$HOME/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/long_read/",proband,"_sorted.bam")
bams_in =  data.frame(bams_in = list.files("/home/renaut/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/long_read/", pattern = 'bam$',recursive = F, full.names = T))
bams_in$long_read_code = gsub('_sorted.bam','',basename(bams_in$bams_in))
  
for(b in 1:nrow(bams_in)){
  out_dir = paste0(params$FRASER,'bams_subset/gene',geneID,'_chr',chr,'_',start,'_',stop)
  command_samtools = paste0("samtools view -h ",bams_in$bams_in[b]," chr",chr,":",start,"-",stop, " | sed 's/\bchrMT\b/chrM/g' | samtools view -bS | samtools sort -o ",out_dir,"/",bams_in$long_read_code[b],"_sorted_chrN.bam") 
  command_index = paste0("samtools index ", out_dir,"/",bams_in$long_read_code[b],"_sorted_chrN.bam")
  
  #
  print(command_samtools)
  print(command_index)
  
  system(command_samtools)
  system(command_index)
}

print('DONE subsetting')

  ###coverage updated
  command_coverage = paste0("samtools depth -H -a ",out_dir,"/*sorted_chrN.bam -r chr",chr,":",start,"-",stop," | awk 'NR % 5 == 1' >",out_dir,"/gene_",geneID,"_",proband,"_depth5.csv")
  print(command_coverage)
  system(command_coverage)

  ###FRASER
  sampleTable = data.table(data.frame(sampleID = gsub('_sorted_chrN.bam','',list.files(out_dir,pattern = '*bam$')),
                                      bamFile = list.files(out_dir,pattern = '*bam$',full.names = F),
                                      group = 1))
  
  sampleTable$pairedEnd = TRUE
  sampleTable$pairedEnd[grep('bc',sampleTable$sampleID)] = FALSE
  
  sampleTable$group = 1:nrow(sampleTable)
  sampleTable$bamFile = paste0(out_dir,'/',sampleTable$bamFile) 
  
  #probands only
  probands = c(1:nrow(sampleTable))[grepl('_0[34]_',sampleTable$bamFile) | grepl('LC_',sampleTable$bamFile)| grepl('F0',sampleTable$bamFile) | grepl('bc',sampleTable$bamFile)]
  settings <- FraserDataSet(colData=sampleTable[probands,], workingDir=params$output)
  
  fds <- try(countRNAData(settings,recount=TRUE))
  
  if(class(fds) != "try-error") {
    
    fds <- calculatePSIValues(fds)
    fds <- annotateRanges(fds,GRCh=38)
    fds <- fit(fds, q=c(jaccard=2))
    fds = calculatePvalues(fds)
    fds = calculatePadjValues(fds,method = 'none',geneLevel=F)
    fds = calculateZscore(fds)  
    
    res <- results(fds, all=TRUE, padjCutoff=NA, deltaPsiCutoff=NA)
    sort(unique(res$hgncSymbol))
    
    res_dt <- as.data.table(res)
    res_dt <- res_dt[sampleID == proband,]
    
    #load more annotation packages
    suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
    suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
    
    txdb_chr <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene,paste0('chr',c(1:22,'X','Y','M')), pruning.mode = "coarse")
    
    temp = c(1:nrow(as.data.frame(res_dt)))[res_dt$hgncSymbol == geneID]; candidate_gene_localisation = temp[!is.na(temp)][1]
    
    #save results
    res_dt_candidate_gene = res_dt[temp[!is.na(temp)],]
    res_dt_candidate_gene$mean = (res_dt_candidate_gene$start + res_dt_candidate_gene$end) / 2
    res_dt_candidate_gene$minuslogpval = -log(res_dt_candidate_gene$pValue,10)
    write.csv(res_dt_candidate_gene,paste0(out_dir,"/gene_",geneID,"_",proband,"_res_dt_candidate_gene.csv"))
    
    #
    print(paste('dimensions de res_dt', dim(res_dt))) 
    print(paste0('Done sample ~~~ ',i,' ~~~ Time is: ',Sys.time()))
    print(candidates_LR[i,])
  }
