# =============================================================================
# fraser_gw.R - Genome-wide splicing analysis with FRASER (per chromosome)
# =============================================================================

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
if (file.exists(paste0(params$bams_subset, 'res_dt.csv')) == TRUE) {
  print(paste0('Already done chr ', params$chromosome, ': Sys.time is: ', Sys.time()))
} else {
  load_install_library('FRASER')

  register(MulticoreParam(params$ncores, params$ncores * 2, progressbar = TRUE))

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

  fds = countRNAData(settings, recount = FALSE, keepNonStandardChromosomes = FALSE,
                      minExpressionInOneSample = 50, filter = TRUE)
  fds = calculatePSIValues(fds)
  fds = annotateRanges(fds, GRCh = 38)
  fds = FRASER(fds, q = c(jaccard = 2))

  # Filter to significant results and save
  res = results(fds, all = TRUE, padjCutoff = NA, deltaPsiCutoff = NA)
  res_dt = as.data.table(res)
  res_dt = res_dt[res_dt$pValue < 0.01, ]
  res_dt$mean = (res_dt$start + res_dt$end) / 2
  res_dt$minuslogpval = -log(res_dt$pValue, 10)
  print(paste0('Dimensions of res_dt: ', nrow(res_dt)))
  write.csv(res_dt, file.path(params$bams_subset, 'res_dt.csv'))

  print(paste0('Done: Sys.time is: ', Sys.time()))
}
