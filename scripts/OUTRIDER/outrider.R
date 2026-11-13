# =============================================================================
# outrider.R - Aberrant expression analysis (per gene and per exon)
# =============================================================================

# Libraries
source("../rnaseq_helper_functions.R")
load_install_library(c('dplyr', 'OUTRIDER', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'org.Hs.eg.db', 'data.table'))

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
params = list(OUTRIDER = file.path(args[1], "OUTRIDER/"))
params$fc_pergene = args[2]
params$fc_perexon = args[3]
params$cpu = as.numeric(args[4])
params$force_outrider = ifelse(is.na(args[6]), FALSE, as.logical(args[6]))
params$table_genes_file = file.path(params$OUTRIDER, 'table_genes_OUTRIDER.rds')
params$table_exons_file = file.path(params$OUTRIDER, 'table_exons_OUTRIDER.rds')
params$gtf = args[5]
params$reference_cohorte = args[7]
# Optional fixed encoding dimension q for OUTRIDER per gene (args[8]); NA/absent = search for the best q each run
params$q_pergene = if (is.na(args[8])) NULL else as.integer(args[8])
# Optional fixed encoding dimension q for OUTRIDER per exon (args[9]); NA/absent = search for the best q each run
params$q_perexon = if (is.na(args[9])) NULL else as.integer(args[9])

dir.create(params$OUTRIDER,showWarnings = F)
register(MulticoreParam(params$cpu, params$cpu * 2, progressbar = FALSE))


# -----------------------------------------------------------------------------
# 2. Load candidate genes
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(args[1], 'data/input/candidate_genes_ALL.csv'))
candidates$ensembl_proband = paste0(candidates$ensembl,'_',candidates$proband)
reference_cohorte = read.table(params$reference_cohorte)[,1]
# -----------------------------------------------------------------------------
# 3. gene annotation map
# -----------------------------------------------------------------------------
# Map Entrez IDs to Ensembl IDs and gene symbols, with gene location info
create_map = function() {
  # One Entrez ID can map to several Ensembl IDs (1:many): expected, so the message is silenced
  map = suppressMessages(select(
    org.Hs.eg.db,
    keys = keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
    keytype = "ENTREZID",
    columns = c("ENSEMBL", "SYMBOL")
  ))
  map = map[!is.na(map$ENSEMBL), ]
  # Keep one row per Ensembl ID so merges by ensemblID downstream do not duplicate rows
  map = map[!duplicated(map$ENSEMBL), ]
  colnames(map)[3] = 'geneID'

  # Add gene location info
  gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = FALSE)
  gene_locations = as.data.table(gene_locations)
  gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames)) < 6, ]
  map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
  colnames(map)[c(2, 3, 5)] = c('ensemblID', 'geneID', 'chr')
  map$pos = (map$start + map$end) / 2
  map$chr = gsub('chr', '', map$chr)
  map
}


# -----------------------------------------------------------------------------
# 4. OUTRIDER per gene (expensive — cached to table_genes_file; only recomputed
#    when that cache is absent OR force_outrider is explicitly requested)
# -----------------------------------------------------------------------------
# If reference_cohorte is given, OUTRIDER is fitted on the reference cohort + the single proband
# candidates$proband[i] only, and only that proband's results are kept. The cache file is then
# specific to that proband. With reference_cohorte = NULL, all probands are fitted together (original behaviour).
# If q is given, OUTRIDER uses that fixed encoding dimension and skips the (slow) findEncodingDim() search.
# With q = NULL, the best q is searched for and printed, so it can be reused on later runs.
outrider_per_gene = function(params, map, hemo_ensembl, candidates, reference_cohorte = NULL, i = 1, q = NULL) {
use_ref = !is.null(reference_cohorte)
if (use_ref) {
  proband = candidates$proband[i]
  ref_dir = file.path(params$OUTRIDER, 'reference_cohorte')
  dir.create(ref_dir, showWarnings = FALSE)
  params$table_genes_file = file.path(ref_dir, paste0('table_genes_OUTRIDER_', proband, '.rds'))
}
if (!file.exists(params$table_genes_file) | params$force_outrider) {
  genes_counts = read.table(params$fc_pergene, sep = '\t', header = TRUE, comment.char = '#', check.names = FALSE)
  rownames(genes_counts) = genes_counts[, 1]
  genes_counts = genes_counts[, -c(1:6)]
  genes_counts = round(genes_counts)
  colnames(genes_counts) = sapply(lapply(strsplit(colnames(genes_counts), '/'), '['), tail, 1)
  colnames(genes_counts) = gsub('_sorted.bam', '', colnames(genes_counts))

  if (use_ref) {
    # Reference cohorte + this proband only
    if (!proband %in% colnames(genes_counts)) stop(proband, ' is not in the gene count matrix.')
    genes_counts = genes_counts[, colnames(genes_counts) %in% unique(c(reference_cohorte, proband))]
  } else {
    # Probands only (includes LC and F0 samples, and _04 twins of _03 probands)
    probands = colnames(genes_counts)[
      grepl('_0[34]_', colnames(genes_counts)) |
      grepl('LC_', colnames(genes_counts)) |
      grepl('F0', colnames(genes_counts))
    ]
    genes_counts = genes_counts[, colnames(genes_counts) %in% c('gene_id', probands)]
  }

  # Remove haemoglobin genes before running OUTRIDER
  genes_counts = genes_counts[!rownames(genes_counts) %in% hemo_ensembl, ]
  ods = OutriderDataSet(countData = genes_counts)
  ods = filterExpression(ods, gtfFile = params$gtf, filterGenes = TRUE, savefpkm = TRUE, fpkmCutoff = 0.5)
  #ods = filterExpression(ods)
  if (is.null(q)) {
    ods = OUTRIDER(ods, iterations = 5)
    } else {
    ods = OUTRIDER(ods, q = q, iterations = 5)
  }

  # Result tables
  table_genes = as.data.frame(results(ods, padjCutoff = 1))
  table_genes$pValue = signif(table_genes$pValue, 4)

  # With a reference cohorte, only keep (and cache) the proband's results, not the reference samples
  if (use_ref) table_genes = table_genes[table_genes$sampleID == proband, ]

  saveRDS(table_genes, params$table_genes_file)
  } else {
    if (!file.exists(params$table_genes_file)) {
    stop(paste0(params$table_genes_file, ' does not exist yet — rerun with force_outrider = TRUE to compute it.'))
  }
    table_genes = readRDS(params$table_genes_file)
    }

  # Keep only the proband's results (reference samples are only there to fit the model)
  
  significant_table = table_genes[table_genes$pValue < 0.05, ]
  significant_table$padjust = signif(significant_table$padjust, 4)
  colnames(significant_table)[1] = 'ensemblID'
  significant_table = merge(significant_table, map, by = 'ensemblID')
  table_genes$ensemblID_sampleID = paste0(table_genes$geneID, '_', table_genes$sampleID)

  candidate_table = table_genes[table_genes$ensemblID_sampleID %in% candidates$ensembl_proband, ]
  colnames(candidate_table)[1] = 'ensemblID'
  candidate_table = merge(candidate_table, map, by = 'ensemblID')
  colnames_candidate_genes = c("sampleID", "geneID", "ensemblID", "pValue", "zScore", "l2fc",
                               "rawcounts", "meanRawcounts", "normcounts", "meanCorrected")
  candidate_table = candidate_table[, colnames_candidate_genes]

  if (file.exists(params$table_genes_file)) {
    print(paste0('Sample ',proband,' already exists. ~~~  Loaded cached table_genes ~~~  Time is: ',Sys.time()))} else {
    print(paste0('Done OUTRIDER against reference cohort for ',proband,'~~~ Time is: ', Sys.time()))
  }

  list(significant_table = significant_table, candidate_table = candidate_table, table_genes = table_genes)
  }

# -----------------------------------------------------------------------------
# 5. OUTRIDER per exon (expensive — cached to table_exons_file; only recomputed
#    when that cache is absent OR force_outrider is explicitly requested)
# -----------------------------------------------------------------------------
# Same logic as outrider_per_gene(): if reference_cohorte is given, OUTRIDER is fitted on the reference cohort +
# the single proband candidates$proband[i] only, and only that proband's rows are kept and cached (one cache file
# per proband). With reference_cohorte = NULL, all samples are fitted together (original behaviour).
# If q is given, a fixed encoding dimension is used and findEncodingDim() is skipped.
# Returns the (per-proband) table_exons; the summary tables are built by outrider_exon_summary() on the combined table.
outrider_per_exon = function(params, map, hemo_ensembl, candidates, reference_cohorte = NULL, i = 1, q = NULL) {
use_ref = !is.null(reference_cohorte)
if (use_ref) {
  proband = candidates$proband[i]
  params$table_exons_file = file.path(params$OUTRIDER, 'reference_cohorte','table_exons_OUTRIDER_', proband, '.rds')
}
if (!file.exists(params$table_exons_file) | params$force_outrider) {
  fc_exons_raw_ALL = read.table(params$fc_perexon, sep = '\t', check.names = FALSE)
  fc_exons_raw_ALL = fc_exons_raw_ALL[, grepl('^bc', colnames(fc_exons_raw_ALL)) == FALSE]

  rownames(fc_exons_raw_ALL) = paste0(
    fc_exons_raw_ALL$geneID, "_",
    fc_exons_raw_ALL$ensemblID, "_",
    fc_exons_raw_ALL$transcriptID, "_",
    fc_exons_raw_ALL$exonID
  )

  # Filter out very low-expression exons and hemo
  genes_counts = fc_exons_raw_ALL[!fc_exons_raw_ALL$ensemblID %in% hemo_ensembl,]
  genes_counts = genes_counts[, -c(1:5)]

  if (use_ref) {
    # Reference cohorte + this proband only
    if (!proband %in% colnames(genes_counts)) stop(proband, ' is not in the exon count matrix.')
    genes_counts = genes_counts[, colnames(genes_counts) %in% unique(c(reference_cohorte, proband))]
  }
  genes_counts = genes_counts[rowMeans(genes_counts) > 1, ]

  ods = OutriderDataSet(countData = genes_counts)
  if (is.null(q)) {
    ods = OUTRIDER(ods)
  } else {
    ods = OUTRIDER(ods, q = q,iterations = 5)
  }

  # Result table
  table_exons = as.data.frame(results(ods, padjCutoff = 1))
  table_exons$combinedID = table_exons$geneID
  table_exons$geneID = sapply(strsplit(table_exons$combinedID, '_'), '[', 1)
  table_exons$ensemblID = sapply(strsplit(table_exons$combinedID, '_'), '[', 2)
  table_exons$transcriptID = sapply(strsplit(table_exons$combinedID, '_'), '[', 3)
  table_exons$exonID = sapply(strsplit(table_exons$combinedID, '_'), '[', 4)
  table_exons$ensemblID_sampleID = paste0(table_exons$ensemblID, '_', table_exons$sampleID)
  table_exons$pValue = signif(table_exons$pValue, 4)

  # With a reference cohorte, only keep (and cache) the proband's results, not the reference samples
  if (use_ref) table_exons = table_exons[table_exons$sampleID == proband, ]

  saveRDS(table_exons, params$table_exons_file)
  
  if (use_ref) print(paste0('Done outrider_per_exon() against reference cohorte for ', proband, ' ~~~ Time is: ', Sys.time()))
  } else {
    table_exons = readRDS(params$table_exons_file)
    print(paste0('Sample ',proband,' already exists. ~~~  Loaded cached table_exons ~~~  Time is: ',Sys.time()))
  }
table_exons
}

# Summary tables / output files from the (combined) per-exon table
outrider_exon_summary = function(params, map, candidates, table_exons, significant_table) {
  # Candidate exons
  colnames_candidate_exon = c("sampleID", "geneID", "ensemblID", "transcriptID", "exonID",
                             "pValue", "zScore", "l2fc", "rawcounts", "meanRawcounts",
                             "normcounts", "meanCorrected")

  candidate_table_exon = table_exons[table_exons$ensemblID_sampleID %in% candidates$ensembl_proband, ]
  candidate_table_exon = candidate_table_exon[, colnames_candidate_exon]

  # All significant exons
  colnames_significant_exon = c("sampleID", "geneID", "ensemblID", "transcriptID", "exonID",
                               "chr", "start", "end", "pos", "width", "pValue", "padjust",
                               "zScore", "l2fc", "rawcounts", "meanRawcounts", "normcounts",
                               "meanCorrected")
  
  table_exons = table_exons[table_exons$pValue < 0.01, ]
  significant_table_exon = merge(table_exons, map[, -3], by = 'ensemblID')
  significant_table_exon = significant_table_exon[, colnames_significant_exon]

  write.table(significant_table_exon, file.path(params$OUTRIDER, 'gw_exons_OUTRIDER.tsv'), sep = '\t', quote = FALSE)
  write.table(candidate_table_exon, file.path(params$OUTRIDER, 'candidates_perexons_OUTRIDER.tsv'), sep = '\t', quote = FALSE)

  # Aggregate min p-value and max z-score per gene per sample
  table_minmax_exons = table_exons %>%
    group_by(geneID, sampleID) %>%
    filter(pValue != Inf) %>%
    summarise(min_pValue = min(pValue), min_zscore = min(zScore), max_zscore = max(zScore))

  table_minmax_exons$exon_zScore = '0'
  colnames(table_minmax_exons)[3] = 'exon_pValue'

  for (i in 1:nrow(table_minmax_exons)) {
    table_minmax_exons$exon_zScore[i] =
      table_minmax_exons[i, 4:5][abs(table_minmax_exons[i, 4:5]) == max(abs(table_minmax_exons[i, 4:5]))][1]
  }

  table_minmax_exons = table_minmax_exons[, -c(4:5)]
  significant_table = merge(
    significant_table, table_minmax_exons,
    by.x = c('geneID', 'sampleID'), by.y = c('geneID', 'sampleID'),
    all.x = TRUE
  )

  colnames_ALL = c("sampleID", "geneID", "ensemblID", "chr", "start", "end", "pos", "width",
                   "pValue", "padjust", "zScore", "l2fc", "rawcounts", "meanRawcounts",
                   "normcounts", "meanCorrected", "exon_pValue", "exon_zScore")
  significant_table = significant_table[, colnames_ALL]

  print(paste0('Done OUTRIDER per exon --- Time is: ', Sys.time()))
  list(significant_table = significant_table, candidate_table_exon = candidate_table_exon)
}

# -----------------------------------------------------------------------------
# 6. Run and save results
# -----------------------------------------------------------------------------
#mapping genes
map = create_map()

# Haemoglobin genes are excluded before running OUTRIDER; also needed by the per-exon section below
hemo_ensembl = map$ensemblID[map$geneID %in% c('HBB','HBA1','HBA2','HBD','HLA','HBG','HBQ','HBE','HBZ','HBM','SELPLG')]

#outrider per gene: one run per unique proband, each against the reference cohorte only
probands_idx = which(!duplicated(candidates$proband))

res_list = list()
table_genes_all = NULL
for (i in probands_idx) {
  res_list[[candidates$proband[i]]] = outrider_per_gene(params, map, hemo_ensembl, candidates,
                                                        reference_cohorte = reference_cohorte, i = i,
                                                        q = params$q_pergene)
  table_genes_all = rbind(table_genes_all, res_list[[candidates$proband[i]]]$table_genes)
}
rownames(table_genes_all) = NULL

res_genes = list(
  significant_table = do.call(rbind, lapply(res_list, `[[`, 'significant_table')),
  candidate_table   = do.call(rbind, lapply(res_list, `[[`, 'candidate_table'))
)

# Save the combined per-proband table (built in the loop above) as the single file read by the Shiny app (and cp_and_cleanup.R)
saveRDS(table_genes_all, params$table_genes_file)

# outrider per exon: same, one run per unique proband against the reference cohorte only
table_exons_list = list()
for (i in probands_idx) {
  table_exons_list[[candidates$proband[i]]] = outrider_per_exon(params, map, hemo_ensembl, candidates,
                                                                reference_cohorte = reference_cohorte, i = i,
                                                                q = params$q_perexon)
}
table_exons_all = do.call(rbind, table_exons_list)
rownames(table_exons_all) = NULL
saveRDS(table_exons_all, params$table_exons_file)

res_exons = outrider_exon_summary(params, map, candidates, table_exons_all, res_genes$significant_table)

write.table(res_exons$significant_table, file.path(params$OUTRIDER, 'gw_genes_OUTRIDER.tsv'), sep = '\t', quote = FALSE)
write.table(res_genes$candidate_table,file.path(params$OUTRIDER, 'candidates_OUTRIDER.tsv'), sep = '\t', quote = FALSE)
