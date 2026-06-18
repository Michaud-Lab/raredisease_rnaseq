# =============================================================================
# global.R - Load packages and datasets for the RNAseq Shiny app
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------
packages = c('DT', 'plotly', 'tidyr', 'shiny', 'shinyjs', 'jsonlite', 'igvShiny',
             'GenomicAlignments', 'dplyr', 'ggtranscript', 'patchwork', 'Hmisc',
             'bslib', 'RColorBrewer', 'ggrepel', 'R.utils', 'logger', 'rtracklayer')
logger::log_info('Loading required packages')

for (p in 1:length(packages)) {
  if (packages[p] %in% installed.packages()) {
    suppressMessages(suppressWarnings(library(packages[p], character.only = TRUE)))
  } else if (packages[p] == 'ggtranscript') {
    stop('Install ggtranscript with: remotes::install_github("dzhang32/ggtranscript")')
  } else {
    stop(paste0('Error in library(): there is no package called ', packages[p]))
  }
}

# -----------------------------------------------------------------------------
# 2. Parameters and theme
# -----------------------------------------------------------------------------
params = list(
  datadir = file.path(getwd(), "data/")
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

fc_exons_raw = read.table(file.path(params$datadir, 'fc_exons_raw.tsv'), sep = '\t', check.names = FALSE)
fc_exons_raw[, -c(1:5)] = round(fc_exons_raw[, -c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir, 'fc_exons_tpm.tsv'), sep = '\t', check.names = FALSE)
fc_genes_tpm = read.table(file.path(params$datadir, 'fc_genes_tpm.tsv'), sep = '\t', check.names = FALSE)
fc_genes_raw_ALL = read.table(file.path(params$datadir, 'fc_genes_raw_ALL.tsv'), sep = '\t', check.names = FALSE)

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

if (file.exists(file.path(params$datadir, 'input/candidate_genes_LR.csv'))) {
  candidates_LR = read.csv(file.path(params$datadir, 'input/candidate_genes_LR.csv'))
  candidates_LR = candidates_LR[, c(2, 3, 10, 4, 5, 6, 7, 8, 9)]
  colnames(candidates_LR) = colnames(candidates)
  candidates = rbind(candidates, candidates_LR)
}

clinical = read.csv(file.path(params$datadir, 'clinical.tsv'), sep = '\t', check.names = FALSE)
html_file = file.path(params$datadir, 'multiqc_report.html')
gwFRASER = read.csv(file.path(params$datadir, 'gwFRASER.csv'), row.names = 1)

report_version = read_json(file.path(params$datadir, '/VERSION.json'))
report_version$data = params$zipfile

logger::log_info("Global resources loaded")
