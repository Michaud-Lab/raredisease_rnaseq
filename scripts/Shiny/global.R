
# R libraries that you need
packages = c('DT','plotly','tidyr','shiny','shinyjs','jsonlite','igvShiny','GenomicAlignments','dplyr','ggtranscript','patchwork','Hmisc','bslib','RColorBrewer','ggrepel','R.utils','logger')
logger::log_info('Loading required packages')

for(p in 1:length(packages)) {
  if(packages[p] %in% installed.packages()) {
    suppressMessages(suppressWarnings(library(packages[p],character.only = TRUE)))
  } else if (packages[p] == 'ggtranscript'){
    stop('Install ggtranscript with : remotes::install_github("dzhang32/ggtranscript")')
  } else {
    stop(paste0('Error in library() : there is no package called ', packages[p]))
  }
}



params <- list(
  datadir = file.path(getwd(), "data/")
)

theme <- bs_theme(
  version = 5,
  bootswatch = "cosmo",
  primary = "#0d6efd",
  base_font = font_google("Roboto")
)
#source("R/reactive_data.R")
#source("R/plotly_expression.R")
load(file = file.path(params$datadir,"gene_annotations.rda"))

# -----------------------------------------------------------------------------
# Load datasets here
# -----------------------------------------------------------------------------

logger::log_info("Global resources loaded")
fc_exons_raw = read.table(file.path(params$datadir,'fc_exons_raw.tsv'),sep = '\t',check.names = F);fc_exons_raw[,-c(1:5)] = round(fc_exons_raw[,-c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir,'fc_exons_tpm.tsv'),sep = '\t',check.names = F)
fc_genes_tpm = read.table(file.path(params$datadir,'fc_genes_tpm.tsv'),sep = '\t',check.names = F)

gwOUTRIDER = read.csv(file.path(params$datadir,'gw_OUTRIDER.tsv'),sep = '\t',check.names = F)
gwOUTRIDER$chr = factor(gwOUTRIDER$chr,levels = c(1:22,'X','Y','MT'))

gwASE = read.csv(file.path(params$datadir,'gwASE.tsv'),row.names = 1,sep = '\t')
gwASE$chr = factor(gwASE$chr,levels = c(1:22,'X','Y','MT'))
gwASE_IMX = read.csv(file.path(params$datadir,'gwImprinted.tsv'),row.names = 1,sep = '\t')

significant_perexons_OUTRIDER = read.csv(file.path(params$datadir,'exons_OUTRIDER.tsv'),sep = '\t',check.names = F)
significant_perexons_OUTRIDER$chr = factor(significant_perexons_OUTRIDER$chr,levels = c(1:22,'X','Y','MT'))

candidates_OUTRIDER = read.csv(file.path(params$datadir,'candidates_OUTRIDER.tsv'),sep = '\t',check.names = F)
candidates_perexons_OUTRIDER = read.csv(file.path(params$datadir,'candidates_perexons_OUTRIDER.tsv'),sep = '\t',check.names = F)

transcripts_named_filtered = read.csv(file.path(params$datadir,'transcripts_named_filtered.tsv'),sep = '\t',check.names = F)
transcripts_named_filtered_ggplot = read.csv(file.path(params$datadir,'transcripts_named_filtered_ggplot.tsv'),sep = '\t',check.names = F)

fc_exons_tpm_ggplot = read.csv(file.path(params$datadir,'fc_exons_tpm_ggplot.tsv'),sep = '\t',check.names = F)
candidates = read.csv(file.path(params$datadir,'input/candidate_genes.csv'))

if(file.exists(file.path(params$datadir,'input/candidate_genes_LR.csv'))){
  candidates_LR = read.csv(file.path(params$datadir,'input/candidate_genes_LR.csv'))
  candidates_LR = candidates_LR[,c(2,3,10,4,5,6,7,8,9)]
  colnames(candidates_LR) = colnames(candidates)
  candidates = rbind(candidates,candidates_LR) 
}

clinical = read.csv(file.path(params$datadir,'clinical.tsv'),sep = '\t',check.names = F)
html_file = file.path(params$datadir,'multiqc_report.html')

gwFRASER = read.csv(file.path(params$datadir,'gwFRASER.csv'),row.names = 1)

report_version = read_json(file.path(params$datadir,'/VERSION.json'))
report_version$data = params$zipfile

theme <- bs_theme(
  version = 5,
  bootswatch = "cosmo",
  primary = "#0d6efd",
  base_font = font_google("Roboto")
)
