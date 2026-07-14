# =============================================================================
# global.R - Load packages and datasets for the RNAseq Shiny app
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load or Install packages
# -----------------------------------------------------------------------------
source('../rnaseq_helper_functions.R')

load_install_library(c('remotes','BiocManager','GenomeInfoDb','DT', 'plotly', 'tidyr', 'shiny', 'shinyjs', 'jsonlite', 'igvShiny', 'shinymanager',
             'GenomicAlignments', 'dplyr', 'ggtranscript', 'patchwork', 'Hmisc',
             'bslib', 'RColorBrewer', 'ggrepel', 'R.utils', 'logger', 'rtracklayer'))


# -----------------------------------------------------------------------------
# 2. Parameters and theme
# -----------------------------------------------------------------------------
data_dir = if (exists("use_data_minimal") && use_data_minimal) "data_minimal" else "data"
params = list(
  datadir = file.path(getwd(), data_dir)
)

theme = bs_theme(
  version = 5,
  bootswatch = "cosmo",
  primary = "#0d6efd",
  base_font = font_google("Roboto")
)

load(file = file.path(params$datadir, "gene_annotations.rda"))

# -----------------------------------------------------------------------------
# 3. Load datasets
# -----------------------------------------------------------------------------
logger::log_info("Loading datasets")

fc_exons_raw = read.table(file.path(params$datadir, 'fc_exons_raw.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
fc_exons_raw[, -c(1:5)] = round(fc_exons_raw[, -c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir, 'fc_exons_tpm.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
fc_genes_tpm = read.table(file.path(params$datadir, 'fc_genes_tpm.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
fc_genes_raw_ALL = read.table(file.path(params$datadir, 'fc_genes_raw_ALL.tsv'), sep = '\t', check.names = FALSE, header = TRUE)

gwOUTRIDER = read.csv(file.path(params$datadir, 'gw_genes_OUTRIDER.tsv'), sep = '\t', check.names = FALSE)
gwOUTRIDER$chr = factor(gwOUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

if (file.exists(file.path(params$datadir, 'gwASE.tsv'))) {
  gwASE = read.csv(file.path(params$datadir, 'gwASE.tsv'), row.names = 1, sep = '\t')
  gwASE$chr = factor(gwASE$chr, levels = c(1:22, 'X', 'Y', 'MT'))
} else {
  gwASE = NULL
  logger::log_warn("gwASE.tsv not found — ASE results will not be available.")
}

if (file.exists(file.path(params$datadir, 'gwImprinted.tsv'))) {
  gwASE_IMX = read.csv(file.path(params$datadir, 'gwImprinted.tsv'), row.names = 1, sep = '\t')
} else {
  gwASE_IMX = NULL
  logger::log_warn("gwImprinted.tsv not found — imprinted/X-linked ASE results will not be available.")
}

significant_perexons_OUTRIDER = read.csv(file.path(params$datadir, 'gw_exons_OUTRIDER.tsv'), sep = '\t', check.names = FALSE)
significant_perexons_OUTRIDER$chr = factor(significant_perexons_OUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

candidates_OUTRIDER = read.csv(file.path(params$datadir, 'candidates_OUTRIDER.tsv'), sep = '\t', check.names = FALSE)
candidates_perexons_OUTRIDER = read.csv(file.path(params$datadir, 'candidates_perexons_OUTRIDER.tsv'), sep = '\t', check.names = FALSE)

transcripts_named_filtered = read.csv(file.path(params$datadir, 'transcripts_named_filtered.tsv'), sep = '\t', check.names = FALSE)
transcripts_named_filtered_ggplot = read.csv(file.path(params$datadir, 'transcripts_named_filtered_ggplot.tsv'), sep = '\t', check.names = FALSE)

fc_exons_tpm_ggplot = read.csv(file.path(params$datadir, 'fc_exons_tpm_ggplot.tsv'), sep = '\t', check.names = FALSE)
candidates = read.csv(file.path(params$datadir, 'input/candidate_genes.csv'))
candidates_extra = read.table(file.path(params$datadir, 'input/candidate_genes_extra.csv'),comment.char = "#",header = T ,sep = ',');candidates_extra[is.na(candidates_extra)] = ''
candidates_extra = candidates_extra[, colnames(candidates)]
candidates = candidates[!grepl('bc',candidates$proband),]

if(file.exists(file.path(params$datadir, 'input/candidate_genes_automated.csv'))) {
  candidate_genes_automated = read.csv(file.path(params$datadir, 'input/candidate_genes_automated.csv'));candidate_genes_automated[is.na(candidate_genes_automated)] = ''
} else {candidate_genes_automated = NULL}

candidates = rbind(candidates,candidates_extra,candidate_genes_automated) %>%
  distinct(geneID,ensembl, proband, .keep_all = TRUE)



clinical = read.csv(file.path(params$datadir, 'clinical.tsv'), sep = '\t', check.names = FALSE)
html_files = list.files(params$datadir,pattern = 'multiqc_report',full.names = T) 
gwFRASER = read.csv(file.path(params$datadir, 'gwFRASER.csv'), row.names = 1)

report_version = read_json(file.path(params$datadir, '/VERSION.json'))
report_version$data = params$zipfile

logger::log_info("Global resources loaded")
