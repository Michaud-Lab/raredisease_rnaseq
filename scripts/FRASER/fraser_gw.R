# =============================================================================
# fraser_gw.R - Genome-wide splicing analysis with FRASER (per chromosome)
#
# Usage: Rscript fraser_gw.R <fraser_chr_bamdir> <ncores> <rnasplice_bamdir> [reference_cohorte]
#   - fraser_gw_count():    subsets ALL .bam to the chromosome and runs countRNAData() on all of them (once)
#   - fraser_gw_pipeline(): runs FRASER on the reference cohort + one proband (every sampleID of make_sample_table()), re-using the counts above
#   [reference_cohorte]  file with one sample ID per line (default: data/input/reference_cohorte.txt)
# =============================================================================
source("../rnaseq_helper_functions.R")

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
params = list(fraser_chr_bamdir = args[1])
params$ncores = as.numeric(args[2])
params$rnasplice_bamdir = args[3]
params$reference_cohorte = if (!is.na(args[4])) args[4] else '../../data/input/reference_cohorte.txt'
# chromosomes to run (each is subsetted, counted and analysed in turn)
params$chromosomes = c(1:22, 'X', 'Y', 'MT')
options(scipen = 999)

load_install_library(c('FRASER', 'data.table', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'tidyr', 'dplyr'))
register(MulticoreParam(params$ncores))

# -----------------------------------------------------------------------------
# 2. Helper: sample table from the chromosome-subsetted .bam files
#    (sample_ids = NULL -> all .bam files; otherwise only .bam whose name starts with '<id>_')
# -----------------------------------------------------------------------------
make_sample_table = function(chromosome, sample_ids = NULL) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)
  bams = list.files(bams_subset, pattern = 'bam$')
  if (!is.null(sample_ids)) {
    keep = Reduce(`|`, lapply(sample_ids, function(id) startsWith(bams, paste0(id, '_'))))
    bams = bams[keep]
  }
  sampleTable = data.table(
    sampleID = gsub('.bam', '', bams),
    bamFile = file.path(bams_subset, bams),
    group = seq_along(bams),
    pairedEnd = TRUE
  )
  sampleTable
}

# -----------------------------------------------------------------------------
# 3. fraser_gw_count(): subset ALL .bam to the chromosome and count them (once)
# -----------------------------------------------------------------------------
fraser_gw_count = function(chromosome) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)
  fds_name = paste0('chr', chromosome, '_ALL')

  # 3a. Subset BAMs to chromosome (skip if already done)
  dir.create(bams_subset, showWarnings = FALSE, recursive = TRUE)
  command = paste("./fraser_gw.sh", chromosome,
                  paste0("'", params$rnasplice_bamdir, "' ",
                         params$fraser_chr_bamdir, '/', chromosome, '/'))

  if (is.na(list.files(bams_subset)[1])) {
    system(command)
  } else {
    print(paste0('Chromosome ', chromosome, ' already subsetted'))
  }

  # 3b. countRNAData() on all .bam (skip if the counts already exist)
  counts_dir = file.path(bams_subset, 'savedObjects', paste0('raw-', fds_name))
  if (dir.exists(counts_dir)) {
    print(paste0('Chromosome ', chromosome, ' already counted'))
    return(invisible(NULL))
  }

  settings = FraserDataSet(colData = make_sample_table(chromosome), workingDir = bams_subset,
                           name = fds_name)
  fds = countRNAData(settings, recount = FALSE, keepNonStandardChromosomes = FALSE,
                     minExpressionInOneSample = 50, filter = TRUE)
  print(paste0('Done counting chr ', chromosome, ' (', ncol(fds), ' samples): Sys.time is: ', Sys.time()))
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# 4. fraser_gw_pipeline(): FRASER on reference cohort + one proband
#    countRNAData() is called on that subset only; per-sample counts cached by
#    fraser_gw_count() (same workingDir) are re-used since recount = FALSE.
# -----------------------------------------------------------------------------
fraser_gw_pipeline = function(chromosome, proband) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)

  res_min_file = file.path(bams_subset, paste0('res_dt_min_', proband, '.csv'))
  if (file.exists(res_min_file)) {
    print(paste0('Already done chr ', chromosome, ', proband ', proband, ': Sys.time is: ', Sys.time()))
    return(invisible(NULL))
  }

  reference_cohorte = read.table(params$reference_cohorte)[, 1]
  # proband = a full sampleID of make_sample_table(); reference cohort IDs are matched on '<id>_'
  sampleTable = unique(rbind(make_sample_table(chromosome, reference_cohorte),
                             make_sample_table(chromosome)[sampleID == proband, ]),
                             by = 'sampleID')
  sampleTable$group = seq_len(nrow(sampleTable))

  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  GenomeInfoDb::seqlevelsStyle(txdb) = "NCBI"
  orgDb = org.Hs.eg.db

  settings = FraserDataSet(colData = sampleTable, workingDir = bams_subset,
                           name = paste0('chr', chromosome, '_', proband))

  fds = countRNAData(settings, recount = FALSE, keepNonStandardChromosomes = FALSE,
                     minExpressionInOneSample = 50, filter = TRUE)

  fds = calculatePSIValues(fds)
  fds = annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgDb)
  fds = FRASER(fds, q = c(jaccard = 2))

  # Keep the proband's results only, filter to significant results and save
  res = suppressMessages(results(fds, all = TRUE, padjCutoff = NA, deltaPsiCutoff = NA))
  res_dt = as.data.table(res)
  res_dt = res_dt[sampleID == proband, ]

  res_dt_001 = res_dt[res_dt$pValue < 0.01, ]
  res_dt_001$mean = (res_dt_001$start + res_dt_001$end) / 2
  res_dt_001$minuslogpval = -log(res_dt_001$pValue, 10)
  write.csv(res_dt_001, file.path(bams_subset, paste0('res_dt_', proband, '.csv')))

  res_dt$mean = (res_dt$start + res_dt$end) / 2
  res_dt$minuslogpval = -log(res_dt$pValue, 10)
  res_dt = res_dt[!is.na(res_dt$hgncSymbol), ]
  res_dt_min = res_dt %>% group_by(hgncSymbol, sampleID) %>% slice_min(pValue) %>% arrange(hgncSymbol, sampleID)
  res_dt_min$sampleID = sub("_[^_]*$", "", res_dt_min$sampleID)
  res_dt_min = res_dt_min[!grepl(';', res_dt_min$hgncSymbol), ]
  write.csv(res_dt_min, res_min_file)

  print(paste0('Done chr ', chromosome, ', proband ', proband, ': Sys.time is: ', Sys.time()))
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# 5. Run
# -----------------------------------------------------------------------------
for (chromosome in params$chromosomes) {
  print(paste0('RUNNING chromosome ', chromosome, ': Sys.time is: ', Sys.time()))
  fraser_gw_count(chromosome)
  probands = make_sample_table(chromosome)$sampleID
  for (proband in probands) fraser_gw_pipeline(chromosome, proband)

  # merge the per-proband results into one res_dt.csv and one res_dt_min.csv for the chromosome
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)
  for (prefix in c('res_dt', 'res_dt_min')) {
    files = file.path(bams_subset, paste0(prefix, '_', probands, '.csv'))
    files = files[file.exists(files)]
    if (length(files) == 0) next
    merged = rbindlist(lapply(files, function(f) fread(f, drop = 1)), fill = TRUE)  # drop = 1: write.csv row-name column
    fwrite(merged, file.path(bams_subset, paste0(prefix, '.csv')))
    print(paste0('Merged ', length(files), ' ', prefix, ' files for chr ', chromosome, ' (', nrow(merged), ' rows)'))
  }
}
