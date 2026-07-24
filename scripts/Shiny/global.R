# =============================================================================
# global.R - Load packages and datasets for the RNAseq Shiny app
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load or Install packages
# -----------------------------------------------------------------------------
source(file.path(params$scriptsdir, 'rnaseq_helper_functions.R'))

load_install_library(c('remotes','BiocManager','GenomeInfoDb','DT', 'reactable', 'plotly', 'tidyr', 'shiny', 'shinyjs', 'jsonlite', 'igvShiny', 'shinymanager',
             'GenomicAlignments', 'dplyr', 'ggtranscript', 'patchwork', 'Hmisc',
             'bslib', 'RColorBrewer', 'ggrepel', 'R.utils', 'logger', 'rtracklayer'))


# -----------------------------------------------------------------------------
# 2. Parameters and theme
# -----------------------------------------------------------------------------
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

fc_exons_raw = read.table(file.path(params$datadir, 'fc_exons_raw.tsv'), check.names = FALSE, header = TRUE)
fc_exons_raw[, -c(1:5)] = round(fc_exons_raw[, -c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir, 'fc_exons_tpm.tsv'), check.names = FALSE, header = TRUE)
fc_genes_tpm = read.table(file.path(params$datadir, 'fc_genes_tpm.tsv'), check.names = FALSE, header = TRUE)
fc_genes_raw_ALL = read.table(file.path(params$datadir, 'fc_genes_raw_ALL.tsv'), check.names = FALSE, header = TRUE)

gwOUTRIDER = read.table(file.path(params$datadir, 'gw_genes_OUTRIDER.tsv'), check.names = FALSE)
gwOUTRIDER$chr = factor(gwOUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

if (file.exists(file.path(params$datadir, 'gwASE.tsv'))) {
  gwASE = read.table(file.path(params$datadir, 'gwASE.tsv'), row.names = 1)
  gwASE$chr = factor(gwASE$chr, levels = c(1:22, 'X', 'Y', 'MT'))
} else {
  gwASE = NULL
  logger::log_warn("gwASE.tsv not found — ASE results will not be available.")
}

if (file.exists(file.path(params$datadir, 'gwImprinted.tsv'))) {
  gwASE_IMX = read.table(file.path(params$datadir, 'gwImprinted.tsv'), row.names = 1)
} else {
  gwASE_IMX = NULL
  logger::log_warn("gwImprinted.tsv not found — imprinted/X-linked ASE results will not be available.")
}

significant_perexons_OUTRIDER = read.table(file.path(params$datadir, 'gw_exons_OUTRIDER.tsv'), check.names = FALSE)
significant_perexons_OUTRIDER$chr = factor(significant_perexons_OUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

candidates_OUTRIDER = read.table(file.path(params$datadir, 'candidates_OUTRIDER.tsv'), check.names = FALSE)
candidates_perexons_OUTRIDER = read.table(file.path(params$datadir, 'candidates_perexons_OUTRIDER.tsv'), check.names = FALSE)

transcripts_named_filtered = read.table(file.path(params$datadir, 'transcripts_named_filtered.tsv'), check.names = FALSE)
transcripts_named_filtered_ggplot = read.table(file.path(params$datadir, 'transcripts_named_filtered_ggplot.tsv'), check.names = FALSE)

fc_exons_tpm_ggplot = read.table(file.path(params$datadir, 'fc_exons_tpm_ggplot.tsv'), check.names = FALSE)
candidates = read.csv(file.path(params$datadir, 'input/candidate_genes_ALL.csv'), check.names = FALSE)

clinical = read.table(file.path(params$datadir, 'clinical.tsv'), check.names = FALSE)
html_files = list.files(params$datadir,pattern = 'multiqc_report',full.names = T) 
gwFRASER = read.table(file.path(params$datadir, 'gwFRASER.tsv'), row.names = 1)

report_version = read_json(file.path(params$datadir, '/VERSION.json'))
report_version$data = params$zipfile

