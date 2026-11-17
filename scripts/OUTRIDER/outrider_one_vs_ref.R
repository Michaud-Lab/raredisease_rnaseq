# =============================================================================
# outrider_one_vs_ref.R - Aberrant expression (per gene), one sample vs a reference cohort
# =============================================================================
# For every target sample, OUTRIDER is fitted on: reference cohort (data/input/reference_cohorte.txt,
# minus the target if it is part of it) + the target. Only the target's results are kept, so every
# target is compared against the same reference and targets never help fit each other.
# Targets = probands of candidate_genes_ALL.csv present in the count matrix (or the single sample given in args[6]).
# Everything above outrider_pipeline() is shared setup; outrider_pipeline() iterates through all the targets.
# Per-target results are cached in OUTRIDER/one_vs_ref/ and only recomputed when force_outrider is TRUE.

source("../rnaseq_helper_functions.R")
load_install_library(c('dplyr', 'OUTRIDER', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'data.table'))

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
params = list(workdir = args[1])
params$OUTRIDER = file.path(params$workdir, 'OUTRIDER')
params$outdir = file.path(params$OUTRIDER, 'one_vs_ref')
params$fc_pergene = args[2]
params$cpu = as.numeric(args[4])
params$gtf = args[5]
params$force_outrider = ifelse(is.na(args[6]), FALSE, as.logical(args[6]))
params$target = args[7]                      # optional: run a single target sample (e.g. for a SLURM array)
params$min_reference = 20                    # OUTRIDER needs a reasonably large cohort to fit the autoencoder
params$table_genes_file = file.path(params$OUTRIDER, 'table_genes_OUTRIDER_one_vs_ref.rds')

dir.create(params$OUTRIDER, showWarnings = FALSE)
dir.create(params$outdir, showWarnings = FALSE)
register(MulticoreParam(params$cpu, params$cpu * 2, progressbar = FALSE))

# -----------------------------------------------------------------------------
# 2. Counts, reference cohort and targets
# -----------------------------------------------------------------------------
genes_counts = read.table(params$fc_pergene, sep = '\t', header = TRUE, comment.char = '#', check.names = FALSE)
rownames(genes_counts) = genes_counts[, 1]
genes_counts = genes_counts[, -c(1:6)]
genes_counts = round(genes_counts)
colnames(genes_counts) = sapply(lapply(strsplit(colnames(genes_counts), '/'), '['), tail, 1)
colnames(genes_counts) = gsub('_sorted.bam', '', colnames(genes_counts))


# Haemoglobin genes are excluded before running OUTRIDER
# Map Entrez IDs to Ensembl IDs and gene symbols, plus gene location (as in outrider.R)
map = select(
  org.Hs.eg.db,
  keys = keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
  keytype = "ENTREZID",
  columns = c("ENSEMBL", "SYMBOL")
)
map = map[!is.na(map$ENSEMBL), ]
colnames(map)[3] = 'geneID'
gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = FALSE)
gene_locations = as.data.table(gene_locations)
gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames)) < 6, ]
map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
colnames(map)[c(2, 3, 5)] = c('ensemblID', 'geneID', 'chr')
map$pos = (map$start + map$end) / 2
map$chr = gsub('chr', '', map$chr)

hemo = suppressMessages(select(org.Hs.eg.db, keys = c('HBB','HBA1','HBA2','HBD','HLA','HBG','HBQ','HBE','HBZ','HBM','SELPLG'), keytype = 'SYMBOL', columns = 'ENSEMBL'))
genes_counts = genes_counts[!rownames(genes_counts) %in% hemo$ENSEMBL, ]

reference_cohorte = trimws(gsub('["\']', '', readLines(file.path(params$workdir, 'data/input/reference_cohorte.txt'))))
reference_cohorte = reference_cohorte[reference_cohorte != '']
missing_ref = setdiff(reference_cohorte, colnames(genes_counts))
if (length(missing_ref) > 0) print(paste0('Reference samples not in count matrix (ignored): ', paste(missing_ref, collapse = ', ')))
reference_cohorte = intersect(reference_cohorte, colnames(genes_counts))
if (length(reference_cohorte) < params$min_reference) stop('Only ', length(reference_cohorte), ' reference samples found; need at least ', params$min_reference)


candidates = read.csv(file.path(params$workdir, 'data/input/candidate_genes_ALL.csv'))
candidates$ensembl_proband = paste0(candidates$ensembl, '_', candidates$proband)

# -----------------------------------------------------------------------------
# 3. outrider_pipeline function
# -----------------------------------------------------------------------------
# Runs OUTRIDER (per gene) once per proband in candidates$proband - each against the reference cohort - then
# writes the combined table and the significant_table / candidate_table outputs.
# Arguments:
#   candidates - data.frame: candidate_genes_ALL.csv contents (needs proband and ensembl_proband columns)
#   probands   - samples to run; defaults to all unique candidates$proband (those present in the count matrix)
outrider_pipeline = function(candidates = candidates) {
  
  probands = unique(candidates$proband)
  probands = intersect(probands, colnames(genes_counts))
  print(paste0('RUNNING outrider_pipeline() ~~~ ', length(probands), ' proband(s), ', length(reference_cohorte), ' reference samples. ~~~ Time is: ', Sys.time()))

  # One fit per proband
  for (proband in probands) {
    proband_file = file.path(params$outdir, paste0('table_genes_', proband, '.rds'))
    if (file.exists(proband_file) && !params$force_outrider) {
      print(paste0('Cached: ', proband)); next
    }
    cohort = c(setdiff(reference_cohorte, proband), proband)
    ods = OutriderDataSet(countData = genes_counts[, cohort])
    ods = filterExpression(ods, gtfFile = params$gtf, filterGenes = TRUE, savefpkm = TRUE, fpkmCutoff = 0.5)
    ods = OUTRIDER(ods)

    table_proband = as.data.frame(results(ods, padjCutoff = 1))
    table_proband = table_proband[table_proband$sampleID == proband, ]
    table_proband$pValue = signif(table_proband$pValue, 4)
    saveRDS(table_proband, proband_file)
    print(paste0('Done ', proband, ' (', length(cohort), ' samples in fit). Time is: ', Sys.time()))
  }

  # Combine (same columns as table_genes_OUTRIDER.rds from outrider.R)
  table_genes = rbindlist(lapply(probands, function(t) readRDS(file.path(params$outdir, paste0('table_genes_', t, '.rds')))))
  table_genes = as.data.frame(table_genes)
  saveRDS(table_genes, params$table_genes_file)

  # significant_table and candidate_table (same processing as outrider.R)
  significant_table = table_genes[table_genes$pValue < 0.05, ]
  significant_table$padjust = signif(significant_table$padjust, 4)
  colnames(significant_table)[1] = 'ensemblID'
  significant_table = merge(significant_table, map, by = 'ensemblID')
  table_genes$ensemblID_sampleID = paste0(table_genes$geneID, '_', table_genes$sampleID)

  candidate_table = table_genes[table_genes$ensemblID_sampleID %in% candidates$ensembl_proband, ]
  colnames(candidate_table)[1] = 'ensemblID'
  candidate_table = merge(candidate_table, map, by = 'ensemblID')

  print(table(significant_table$sampleID))

  # Per-exon results are not computed here: keep the exon columns (NA) so the output format matches outrider.R
  significant_table$exon_pValue = NA
  significant_table$exon_zScore = NA

  colnames_ALL = c("sampleID", "geneID", "ensemblID", "chr", "start", "end", "pos", "width",
                   "pValue", "padjust", "zScore", "l2fc", "rawcounts", "meanRawcounts",
                   "normcounts", "meanCorrected", "exon_pValue", "exon_zScore")
  colnames_candidate_genes = c("sampleID", "geneID", "ensemblID", "pValue", "zScore", "l2fc",
                               "rawcounts", "meanRawcounts", "normcounts", "meanCorrected")

  write.table(significant_table[, colnames_ALL], file.path(params$OUTRIDER, 'gw_genes_OUTRIDER_one_vs_ref.tsv'), sep = '\t', quote = FALSE)
  write.table(candidate_table[, colnames_candidate_genes], file.path(params$OUTRIDER, 'candidates_OUTRIDER_one_vs_ref.tsv'), sep = '\t', quote = FALSE)

  print(paste0('Done OUTRIDER one-vs-reference --- Time is: ', Sys.time()))
  invisible(list(significant_table = significant_table, candidate_table = candidate_table))
}

# -----------------------------------------------------------------------------
# 4. Run the outrider_pipeline() function
# -----------------------------------------------------------------------------
outrider_pipeline(candidates, probands = if (!is.na(params$target)) params$target else unique(candidates$proband))
