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
params$force_outrider = ifelse(is.na(args[5]), FALSE, as.logical(args[5]))
params$table_genes_file = file.path(params$OUTRIDER, 'table_genes.rds')
params$table_exons_file = file.path(params$OUTRIDER, 'table_exons.rds')

dir.create(params$OUTRIDER)
register(MulticoreParam(params$cpu, params$cpu * 2, progressbar = FALSE))


# -----------------------------------------------------------------------------
# 2. Load candidate genes and gene annotation map
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(args[1], 'data/input/candidate_genes_ALL.csv'))
candidates$ensembl_proband = paste0(candidates$ensembl,'_',candidates$proband)

# Map Entrez IDs to Ensembl IDs and gene symbols
map = select(
  org.Hs.eg.db,
  keys = keys(TxDb.Hsapiens.UCSC.hg38.knownGene, keytype = "GENEID"),
  keytype = "ENTREZID",
  columns = c("ENSEMBL", "SYMBOL")
)
map = map[!is.na(map$ENSEMBL), ]
colnames(map)[3] = 'geneID'

# Add gene location info
gene_locations = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = FALSE)
gene_locations = as.data.table(gene_locations)
gene_locations = gene_locations[nchar(as.character(gene_locations$seqnames)) < 6, ]
map = merge(map, gene_locations, by.y = 'group_name', by.x = 'ENTREZID')
colnames(map)[c(2, 3, 5)] = c('ensemblID', 'geneID', 'chr')
map$pos = (map$start + map$end) / 2
map$chr = gsub('chr', '', map$chr)

# Haemoglobin genes are excluded before running OUTRIDER; also needed by the per-exon section below
hemo_ensembl = map$ensemblID[map$geneID %in% c('HBB','HBA1','HBA2','HBD')]

# -----------------------------------------------------------------------------
# 3. OUTRIDER per gene (expensive — cached to table_genes_file; only recomputed
#    when that cache is absent AND force_outrider is explicitly requested)
# -----------------------------------------------------------------------------
if (!file.exists(params$table_genes_file) && params$force_outrider) {
  genes_counts = read.table(params$fc_pergene, sep = '\t', header = TRUE, comment.char = '#', check.names = FALSE)
  rownames(genes_counts) = genes_counts[, 1]
  genes_counts = genes_counts[, -c(1:6)]
  genes_counts = round(genes_counts)
  colnames(genes_counts) = sapply(lapply(strsplit(colnames(genes_counts), '/'), '['), tail, 1)
  colnames(genes_counts) = gsub('_sorted.bam', '', colnames(genes_counts))

  # Probands only (includes LC and F0 samples, and _04 twins of _03 probands)
  probands = colnames(genes_counts)[
    grepl('_0[34]_', colnames(genes_counts)) |
    grepl('LC_', colnames(genes_counts)) |
    grepl('F0', colnames(genes_counts))
  ]
  genes_counts = genes_counts[, colnames(genes_counts) %in% c('gene_id', probands)]

  # Remove haemoglobin genes before running OUTRIDER
  genes_counts = genes_counts[!rownames(genes_counts) %in% hemo_ensembl, ]
  ods = OutriderDataSet(countData = genes_counts)
  ods = filterExpression(ods, TxDb.Hsapiens.UCSC.hg38.knownGene, mapping = map[, 1:2], filterGenes = TRUE, savefpkm = TRUE, fpkmCutoff = 0.5)
  #ods = filterExpression(ods)
  ods = OUTRIDER(ods)

  # Result tables
  table_genes = as.data.frame(results(ods, padjCutoff = 1))
  table_genes$pValue = signif(table_genes$pValue, 4)

  saveRDS(table_genes, params$table_genes_file)
} else {
  if (!file.exists(params$table_genes_file)) {
    stop(paste0(params$table_genes_file, ' does not exist yet — rerun with force_outrider = TRUE to compute it.'))
  }
  table_genes = readRDS(params$table_genes_file)
  print(paste0('Loaded cached table_genes from ', params$table_genes_file ))
}

significant_table = table_genes[table_genes$pValue < 0.05, ]
significant_table$padjust = signif(significant_table$padjust, 4)
colnames(significant_table)[1] = 'ensemblID'
significant_table = merge(significant_table, map, by = 'ensemblID')
table_genes$ensemblID_sampleID = paste0(table_genes$geneID, '_', table_genes$sampleID)

candidate_table = table_genes[table_genes$ensemblID_sampleID %in% candidates$ensembl_proband, ]
colnames(candidate_table)[1] = 'ensemblID'
candidate_table = merge(candidate_table, map, by = 'ensemblID')


print(table(significant_table$sampleID))

# Output column sets
colnames_ALL = c("sampleID", "geneID", "ensemblID", "chr", "start", "end", "pos", "width",
                 "pValue", "padjust", "zScore", "l2fc", "rawcounts", "meanRawcounts",
                 "normcounts", "meanCorrected", "exon_pValue", "exon_zScore")
colnames_candidate_genes = c("sampleID", "geneID", "ensemblID", "pValue", "zScore", "l2fc",
                              "rawcounts", "meanRawcounts", "normcounts", "meanCorrected")

print(paste0('Done OUTRIDER per gene --- Time is: ', Sys.time()))

# -----------------------------------------------------------------------------
# 4. OUTRIDER per exon (expensive — cached to table_exons_file; only recomputed
#    when that cache is absent AND force_outrider is explicitly requested)
# -----------------------------------------------------------------------------
if (!file.exists(params$table_exons_file) && params$force_outrider) {
  fc_exons_raw_ALL = read.table(params$fc_perexon, sep = '\t', check.names = FALSE)
  fc_exons_raw_ALL = fc_exons_raw_ALL[, grepl('^bc', colnames(fc_exons_raw_ALL)) == FALSE]

  # HSJ samples appear with two naming conventions; add '_PAX' suffix to align them
  colnames(fc_exons_raw_ALL)[grep('HSJ', colnames(fc_exons_raw_ALL))] =
    paste0(colnames(fc_exons_raw_ALL)[grep('HSJ', colnames(fc_exons_raw_ALL))], '_PAX')

  rownames(fc_exons_raw_ALL) = paste0(
    fc_exons_raw_ALL$geneID, "_",
    fc_exons_raw_ALL$ensemblID, "_",
    fc_exons_raw_ALL$transcriptID, "_",
    fc_exons_raw_ALL$exonID
  )

  # Filter out very low-expression exons and hemo
  genes_counts = fc_exons_raw_ALL[!fc_exons_raw_ALL$ensemblID %in% hemo_ensembl,]
  genes_counts = fc_exons_raw_ALL[, -c(1:5)]
  genes_counts = genes_counts[rowMeans(genes_counts) > 1, ]

  ods = OutriderDataSet(countData = genes_counts)
  ods = OUTRIDER(ods)

  # Result table
  table_exons = as.data.frame(results(ods, padjCutoff = 1))
  table_exons$combinedID = table_exons$geneID
  table_exons$geneID = sapply(strsplit(table_exons$combinedID, '_'), '[', 1)
  table_exons$ensemblID = sapply(strsplit(table_exons$combinedID, '_'), '[', 2)
  table_exons$transcriptID = sapply(strsplit(table_exons$combinedID, '_'), '[', 3)
  table_exons$exonID = sapply(strsplit(table_exons$combinedID, '_'), '[', 4)
  table_exons$ensemblID_sampleID = paste0(table_exons$ensemblID, '_', table_exons$sampleID)
  table_exons$pValue = signif(table_exons$pValue, 4)

  saveRDS(table_exons, params$table_exons_file)
} else {
  if (!file.exists(params$table_exons_file)) {
    stop(paste0(params$table_exons_file, ' does not exist yet — rerun with force_outrider = TRUE to compute it.'))
  }
  table_exons = readRDS(params$table_exons_file)
  print(paste0('Loaded cached table_exons from ', params$table_exons_file))
}

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
print(table(significant_table_exon$sampleID))


write.table(significant_table_exon, file.path(params$OUTRIDER, 'gw_exons_OUTRIDER.tsv'), sep = '\t', quote = FALSE)
write.table(candidate_table_exon, file.path(params$OUTRIDER, 'candidates_perexons_OUTRIDER.tsv'), sep = '\t', quote = FALSE)

# Aggregate min p-value and max z-score per gene per sample
table_minmax_exons = table_exons %>%
  group_by(geneID, sampleID) %>%
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

# -----------------------------------------------------------------------------
# 5. Save results
# -----------------------------------------------------------------------------
write.table(significant_table[, colnames_ALL], file.path(params$OUTRIDER, 'gw_genes_OUTRIDER.tsv'), sep = '\t', quote = FALSE)
write.table(candidate_table[, colnames_candidate_genes], file.path(params$OUTRIDER, 'candidates_OUTRIDER.tsv'), sep = '\t', quote = FALSE)

print(paste0('Done OUTRIDER per exon --- Time is: ', Sys.time()))
