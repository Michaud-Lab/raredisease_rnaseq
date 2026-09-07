# =============================================================================
# fraser_gw.R - Genome-wide splicing analysis with FRASER (per chromosome)
source("../rnaseq_helper_functions.R")

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
params = list(bams_subset = paste0(args[1], '/', args[2]))
params$chromosome = args[2]
params$ncores = as.numeric(args[3])
params$rnasplice_bamdir = args[4]
params$fraser_chr_bamdir = args[1]
options(scipen = 999)


# -----------------------------------------------------------------------------
# 2. Subset BAMs to chromosome (skip if already done)
# -----------------------------------------------------------------------------
dir.create(params$bams_subset, showWarnings = FALSE, recursive = TRUE)
command = paste0('./fraser_gw.sh ', params$chromosome, ' ', params$rnasplice_bamdir,
                 ' ', params$fraser_chr_bamdir, '/', params$chromosome)
if (is.na(list.files(params$bams_subset)[1])) {
  system(command)
} else {
  print(paste0('Chromosome ', params$chromosome, ' already subsetted'))
}

# -----------------------------------------------------------------------------
# 3. Run FRASER (skip if results already exist)
# -----------------------------------------------------------------------------
if (file.exists(file.path(params$bams_subset, 'res_dt_ALL.csv')) == TRUE) {
  print(paste0('Already done chr ', params$chromosome, ': Sys.time is: ', Sys.time()))
} else {
  load_install_library(c('FRASER','data.table','TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db','tidyr','dplyr'))
 
     txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
     GenomeInfoDb::seqlevelsStyle(txdb) = "NCBI"
     orgDb = org.Hs.eg.db

  register(MulticoreParam(params$ncores))

  sampleTable = data.table(data.frame(
    sampleID = gsub('.bam', '', list.files(params$bams_subset, pattern = '*bam$')),
    bamFile = list.files(params$bams_subset, pattern = '*bam$', full.names = FALSE),
    group = 1,
    pairedEnd = TRUE
  ))

  sampleTable$group = 1:nrow(sampleTable)
  sampleTable$bamFile = paste0(params$bams_subset, '/', sampleTable$bamFile)

  # Probands only
  probands = which(
    grepl('_0[34]_', sampleTable$bamFile) |
    grepl('LC_', sampleTable$bamFile) |
    grepl('F0', sampleTable$bamFile)
  )
  settings = FraserDataSet(colData = sampleTable[probands, ], workingDir = params$bams_subset)

  fds =                                          countRNAData(settings, recount = FALSE, keepNonStandardChromosomes = FALSE,
                                                      minExpressionInOneSample = 50, filter = TRUE)

#                         ))
  
  fds = calculatePSIValues(fds)
  fds = annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
  fds = FRASER(fds, q = c(jaccard = 2))

  # Filter to significant results and save
  res = results(fds, all = TRUE, padjCutoff = NA, deltaPsiCutoff = NA)
  res_dt = as.data.table(res)
  res_dt_001 = res_dt[res_dt$pValue < 0.01, ]
  res_dt_001$mean = (res_dt_001$start + res_dt_001$end) / 2
  res_dt_001$minuslogpval = -log(res_dt_001$pValue, 10)
  write.csv(res_dt_001, file.path(params$bams_subset, 'res_dt.csv'))

  res_dt$mean = (res_dt$start + res_dt$end) / 2
  res_dt$minuslogpval = -log(res_dt$pValue, 10)
  res_dt = res_dt[!is.na(res_dt$hgncSymbol),]
  res_dt_min = res_dt %>% group_by(hgncSymbol,sampleID) %>% slice_min(pValue) %>% arrange(hgncSymbol,sampleID)
  res_dt_min$sampleID = sub("_[^_]*$","",res_dt_min$sampleID)
  res_dt_min = res_dt_min[!grepl(';',res_dt_min$hgncSymbol),]

  #res_dt_ALL = res_dt[,c("seqnames","start","end","sampleID","hgncSymbol","pValue","padjust")]
  write.csv(res_dt_min, file.path(params$bams_subset, 'res_dt_min.csv'))

  print(paste0('Done: Sys.time is: ', Sys.time()))
}
