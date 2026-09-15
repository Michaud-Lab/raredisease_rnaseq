# =============================================================================
# global.R - Load packages and datasets for the RNAseq Shiny app
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load or Install packages
# -----------------------------------------------------------------------------
source(file.path(params$scriptsdir, 'rnaseq_helper_functions.R'))
logger::log_info("Loading libraries")
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

# -----------------------------------------------------------------------------
# 3. Load datasets
# -----------------------------------------------------------------------------
# load_rnaseq_dataset(): reads one dataset directory (e.g. "data") and returns
# everything the app needs for it as a single named list.
load_rnaseq_dataset = function(datadir) {

  logger::log_info(paste0("Loading datasets ~ ",datadir))

  ds = list(datadir = datadir)

  ds$gene_annotations = local({
    e = new.env()
    load(file = file.path(datadir, "gene_annotations.rda"), envir = e)
    e$gene_annotations
  })

  ds$fc_exons_raw = read.table(file.path(datadir, 'fc_exons_raw.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$fc_exons_raw[, -c(1:5)] = round(ds$fc_exons_raw[, -c(1:5)])
  ds$fc_exons_tpm = read.table(file.path(datadir, 'fc_exons_tpm.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$fc_genes_tpm = read.table(file.path(datadir, 'fc_genes_tpm.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$fc_genes_raw_ALL = read.table(file.path(datadir, 'fc_genes_raw_ALL.tsv'), sep = '\t', check.names = FALSE, header = TRUE)

  ds$gwOUTRIDER = read.table(file.path(datadir, 'gw_genes_OUTRIDER.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$gwOUTRIDER$chr = factor(ds$gwOUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

  if (file.exists(file.path(datadir, 'gwASE.tsv'))) {
    ds$gwASE = read.table(file.path(datadir, 'gwASE.tsv'), row.names = 1, sep = '\t', header = TRUE)
    ds$gwASE$chr = factor(ds$gwASE$chr, levels = c(1:22, 'X', 'Y', 'MT'))
  } else {
    ds$gwASE = NULL
    logger::log_warn(paste0("gwASE.tsv not found in ", datadir, " — ASE results will not be available."))
  }

  if (file.exists(file.path(datadir, 'gwImprinted.tsv'))) {
    ds$gwASE_IMX = read.table(file.path(datadir, 'gwImprinted.tsv'), row.names = 1, sep = '\t', header = TRUE)
  } else {
    ds$gwASE_IMX = NULL
    logger::log_warn(paste0("gwImprinted.tsv not found in ", datadir, " — imprinted/X-linked ASE results will not be available."))
  }

  ds$significant_perexons_OUTRIDER = read.table(file.path(datadir, 'gw_exons_OUTRIDER.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$significant_perexons_OUTRIDER$chr = factor(ds$significant_perexons_OUTRIDER$chr, levels = c(1:22, 'X', 'Y', 'MT'))

  # Per-gene, per-sample OUTRIDER statistics (genome-wide, all cohort samples).
  # The .rds stores the Ensembl gene ID under "geneID" — relabel it "ensemblID" and
  # join in the gene symbol (from fc_genes_raw_ALL) as "geneID", matching the
  # geneID=symbol / ensemblID=Ensembl convention used by the other tables above.
  if (file.exists(file.path(datadir, 'table_genes_OUTRIDER.rds'))) {
    ds$table_genes_OUTRIDER = readRDS(file.path(datadir, 'table_genes_OUTRIDER.rds'))
    colnames(ds$table_genes_OUTRIDER)[colnames(ds$table_genes_OUTRIDER) == 'geneID'] = 'ensemblID'
    gene_symbol_map = unique(ds$fc_genes_raw_ALL[, c('geneID', 'ensemblID')])
    ds$table_genes_OUTRIDER = merge(gene_symbol_map, ds$table_genes_OUTRIDER, by = 'ensemblID')
  } else {
    ds$table_genes_OUTRIDER = NULL
    logger::log_warn(paste0("table_genes_OUTRIDER.rds not found in ", datadir, " — per-gene OUTRIDER statistics will not be available."))
  }

  # Per-gene, per-sample minimum FRASER splicing p-values (genome-wide, all cohort samples).
  # Already flattened across chromosomes; "hgncSymbol" can list several overlapping genes
  # separated by ';' per splicing event -- relabel it "geneID" for consistency with the
  # other tables, exact-matched (not token-parsed) against the selected gene in the Shiny app.
  if (file.exists(file.path(datadir, 'gwFRASER_min.rds'))) {
    ds$gwFRASER_min = readRDS(file.path(datadir, 'gwFRASER_min.rds'))
    colnames(ds$gwFRASER_min)[colnames(ds$gwFRASER_min) == 'hgncSymbol'] = 'geneID'
  } else {
    ds$gwFRASER_min = NULL
    logger::log_warn(paste0("gwFRASER_min.rds not found in ", datadir, " — per-gene FRASER statistics will not be available."))
  }

  ds$candidates_OUTRIDER = read.table(file.path(datadir, 'candidates_OUTRIDER.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$candidates_perexons_OUTRIDER = read.table(file.path(datadir, 'candidates_perexons_OUTRIDER.tsv'), sep = '\t', check.names = FALSE, header = TRUE)

  ds$fc_exons_tpm_ggplot = read.table(file.path(datadir, 'fc_exons_tpm_ggplot.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$candidates = read.csv(file.path(datadir, 'input/candidate_genes_ALL.csv'), check.names = FALSE)

  ds$clinical = read.table(file.path(datadir, 'clinical.tsv'), sep = '\t', check.names = FALSE, header = TRUE)
  ds$html_files = list.files(datadir, pattern = 'multiqc_report', full.names = TRUE)
  ds$gwFRASER = read.table(file.path(datadir, 'gwFRASER.tsv'), row.names = 1, sep = '\t', header = TRUE)

  ds$report_version = read_json(file.path(datadir, '/VERSION.json'))
  ds$report_version$data = params$zipfile

  ds
}

