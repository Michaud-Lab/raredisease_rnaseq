# =============================================================================
# fraser_gw.R - Genome-wide splicing analysis with FRASER (per chromosome)
#
# Usage: Rscript fraser_gw.R <fraser_chr_bamdir> <ncores> <rnasplice_bamdir> [reference_cohorte]
#   - fraser_gw_subset_bams(): runs fraser_gw.sh to subset the .bam to the chromosome (new .bam only)
#   - fraser_gw_count():    runs countRNAData() on all of them (again if new .bam)
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
# 3. fraser_gw_subset_bams(): subset ALL .bam to the chromosome (new .bam only)
# -----------------------------------------------------------------------------
fraser_gw_subset_bams = function(chromosome) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)
  dir.create(bams_subset, showWarnings = FALSE, recursive = TRUE)
  command = paste("./fraser_gw.sh", chromosome,
                  paste0("'", params$rnasplice_bamdir, "' ",
                         params$fraser_chr_bamdir, '/', chromosome, '/'))

  # fraser_gw.sh skips the .bam that are already subsetted, so it is always safe to call (new .bam get subsetted)
  system(command)
  print(paste0('Done subsetting chr ', chromosome, ': Sys.time is: ', Sys.time()))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 3. fraser_gw_count(): count them (again only if there are new .bam, otherwise because recount = FALSE, you should be able to load the cached samples. But calculatePSIValues will run again)
# -----------------------------------------------------------------------------
fraser_gw_count = function(chromosome) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)
  fds_name = paste0('chr', chromosome, '_ALL')

  # 3b. countRNAData(), PSI/theta/jaccard and annotation on all .bam, saved once
  # (skip if the saved dataset exists; countRNAData() itself re-uses cached counts if any)
  fds_file = file.path(bams_subset, 'savedObjects', fds_name, 'fds-object.RDS')
  if (file.exists(fds_file)) {
    # up to date only if every subsetted .bam is already in the saved dataset
    saved_ids = samples(loadFraserDataSet(dir = bams_subset, name = fds_name))
    new_ids = setdiff(make_sample_table(chromosome)$sampleID, saved_ids)
    if (length(new_ids) == 0) {
      print(paste0('Chromosome ', chromosome, ' already counted'))
      return(invisible(NULL))
    }
    print(paste0('Chromosome ', chromosome, ': ', length(new_ids), ' new sample(s) (', paste(new_ids, collapse = ', '), '), recounting'))
    unlink(file.path(bams_subset, 'savedObjects', fds_name), recursive = TRUE)  # per-sample counts stay cached in the 'cache' folder
  }

  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  GenomeInfoDb::seqlevelsStyle(txdb) = "NCBI"
  orgDb = org.Hs.eg.db

  settings = FraserDataSet(colData = make_sample_table(chromosome), workingDir = bams_subset,
                           name = fds_name)
  fds = suppressMessages(countRNAData(settings, recount = FALSE, keepNonStandardChromosomes = FALSE,
                                      minExpressionInOneSample = 50, filter = TRUE))
  fds = calculatePSIValues(fds)   # psi5, psi3, theta, jaccard: per sample, independent of the cohort
  fds = annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgDb)
  saveFraserDataSet(fds)
  print(paste0('Done counting chr ', chromosome, ' (', ncol(fds), ' samples): Sys.time is: ', Sys.time()))
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# 4. fraser_gw_pipeline(): FRASER on reference cohort + one proband
#    Loads the dataset saved by fraser_gw_count() (counts, PSI/theta/jaccard, annotation),
#    keeps the reference cohort + the sample, then runs the cohort-dependent part (FRASER fit,
#    p-values, Z-scores).
# -----------------------------------------------------------------------------
fraser_gw_pipeline = function(chromosome, proband) {
  bams_subset = file.path(params$fraser_chr_bamdir, chromosome)

  res_min_file = file.path(bams_subset, paste0('res_dt_min_', proband, '.csv'))
  if (file.exists(res_min_file)) {
    print(paste0('Already done chr ', chromosome, ', proband ', proband, ': Sys.time is: ', Sys.time()))
    return(invisible(NULL))
  }

  reference_cohorte = read.table(params$reference_cohorte)[, 1]
  fds = loadFraserDataSet(dir = bams_subset, name = paste0('chr', chromosome, '_ALL'))

  # proband = a full sampleID of make_sample_table(); reference cohort IDs are matched on '<id>_'
  all_ids = samples(fds)
  keep = (all_ids == proband) |
    Reduce(`|`, lapply(reference_cohorte, function(id) startsWith(all_ids, paste0(id, '_'))))
  fds = fds[, all_ids[keep]]
  # own name, so FRASER's fit results are not written over the shared 'chr<N>_ALL' files
  name(fds) = paste0('chr', chromosome, '_', proband)

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
# 5a. Subset the BAMs and count, for all chromosomes
for (chromosome in params$chromosomes) {
  print(paste0('SUBSETTING AND COUNTING chromosome ', chromosome, ': Sys.time is: ', Sys.time()))
  fraser_gw_subset_bams(chromosome)
  fraser_gw_count(chromosome)
}

# 5b. Rest of the pipeline (per chromosome x proband), in parallel
#     Each (chromosome, proband) job is independent (own saved dataset name, own result files), so they are
#     spread over params$ncores forked workers. Inside a worker FRASER runs serially (no nested parallelism).
#     Memory scales with the number of workers: lower ncores if jobs get killed (out of memory).
jobs = rbindlist(lapply(params$chromosomes, function(chromosome)
  data.table(chromosome = chromosome, proband = make_sample_table(chromosome)$sampleID)))

run_job = function(i) {
  register(SerialParam())
  fraser_gw_pipeline(jobs$chromosome[i], jobs$proband[i])
  TRUE
}
print(paste0('RUNNING ', nrow(jobs), ' chromosome/proband jobs on ', params$ncores, ' cores: Sys.time is: ', Sys.time()))
job_res = parallel::mclapply(seq_len(nrow(jobs)), function(i) try(run_job(i)),
                             mc.cores = params$ncores, mc.preschedule = FALSE)
failed = which(vapply(job_res, function(r) !isTRUE(r), logical(1)))  # error, or worker killed (NULL)
if (length(failed) > 0) {
  print(paste0('FAILED jobs (chr, proband): ', paste0(jobs$chromosome[failed], ', ', jobs$proband[failed], collapse = ' | ')))
}

# 5c. Merge the per-proband results into one res_dt.csv and one res_dt_min.csv per chromosome
for (chromosome in params$chromosomes) {
  probands = make_sample_table(chromosome)$sampleID
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
