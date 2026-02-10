args = commandArgs(trailingOnly=TRUE)
params = list(bams_subset = paste0(args[1],'/',args[2]))
params$chromosome = args[2]
params$ncores = args[3]
options(scipen = 999)

#subset chromosome 
dir.create(params$bams_subset,showWarnings = F)
command = paste0('./fraser_gw.sh ', params$chromosome)
if(is.na(list.files(params$bams_subset)[1])) system(command) else{print(paste0('Chromosome ',params$chromosome,' already subsetted'))}

if(file.exists(paste0(params$bams_subset,'res_dt.csv'))==T) print(paste0('ALLREADY DONE CHR ',params$chromosome,': Sys.time is: ', Sys.time())) else{ 

#
options(scipen = 999)

#Load required packages
suppressMessages(suppressWarnings(library(FRASER)))

#parameters
ncores = params$ncores
register(MulticoreParam(ncores, ncores*2, progressbar = TRUE))

sampleTable = data.table(data.frame(sampleID =  gsub('.bam','',list.files(params$bams_subset,pattern = '*bam$')),
                           bamFile = list.files(params$bams_subset,pattern = '*bam$',full.names = F),
                           group = 1,
                           pairedEnd = TRUE))

sampleTable$group = 1:nrow(sampleTable)
sampleTable$bamFile = paste0(params$bams_subset,'/',sampleTable$bamFile) 

#probands only
probands = c(1:nrow(sampleTable))[grepl('_0[34]_',sampleTable$bamFile) | grepl('LC_',sampleTable$bamFile)| grepl('F0',sampleTable$bamFile)]
settings <- FraserDataSet(colData=sampleTable[probands,], workingDir=params$bams_subset)

#FRASER
fds <- countRNAData(settings,recount=FALSE,keepNonStandardChromosomes=F,minExpressionInOneSample = 50,filter = TRUE)
fds <- calculatePSIValues(fds)
fds <- annotateRanges(fds,GRCh=38)
fds <- FRASER(fds, q=c(jaccard=2))

#results  
res <- results(fds, all=TRUE, padjCutoff=NA, deltaPsiCutoff=NA)
res_dt <- as.data.table(res)
res_dt = res_dt[res_dt$pValue < 0.01,]
res_dt$mean = (res_dt$start + res_dt$end) / 2
res_dt$minuslogpval = -log(res_dt$pValue,10)
print(paste0('Dimension de res_dt: ',nrow(res_dt)))
write.csv(res_dt,file.path(params$bams_subset,'res_dt.csv'))
  
print(paste0('DONE: Sys.time is: ', Sys.time()))
}
