# =============================================================================
# Load or Install packages 
# =============================================================================

load_install_library = function(packages,silent = T) {
  for (p in 1:length(packages)) {
    if (packages[p] %in% installed.packages()) {
      if(silent) suppressMessages(suppressWarnings(library(packages[p], character.only = TRUE))) else library(packages[p], character.only = TRUE)
    } else if (packages[p] == 'ggtranscript') {
      remotes::install_github("dzhang32/ggtranscript")
      if(silent) suppressMessages(suppressWarnings(library(packages[p], character.only = TRUE))) else library(packages[p], character.only = TRUE)
    } else {
      bioc_pkgs = BiocManager::available()
      if(packages[p] %in% bioc_pkgs) {
        BiocManager::install(packages[p])
      } else {
        install.packages(packages[p])
      }
      if(silent) suppressMessages(suppressWarnings(library(packages[p], character.only = TRUE))) else library(packages[p], character.only = TRUE)
    }
    if(p == length(packages)) print('Finished installing required packages')
  }
}

# =============================================================================
# rnaseq_helper_functions.R - Helper functions for featureCounts processing
# =============================================================================

# gene_annotation: fetch exon-level gene models from Ensembl via biomaRt
# Returns a list: [[1]] GRanges of exons, [[2]] GRanges of gene bodies
gene_annotation = function(unique_transcript_id = unique(fc_exons_raw$transcriptID),
                           candidates = candidates){

  load_install_library(c('biomaRt', 'ggbio', 'GenomicAlignments'))

  # Connect to Ensembl and select the human dataset for hg38 (GRCh38)
  ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 38, mirror = 'useast')

  # Get gene annotations
  genes = biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position",
                                "strand", "ensembl_gene_id", "hgnc_symbol", "gene_biotype"),mart = ensembl)

  # Subset to standard chromosomes (exclude scaffolds and patches)
  genes = genes[genes$chromosome_name %in% c(1:22,'X','Y','MT'),]

  # Create GRanges object
  gr = GenomicRanges::GRanges(seqnames = paste0("chr", genes$chromosome_name),
                ranges = IRanges::IRanges(start = genes$start_position, end = genes$end_position),
                strand = ifelse(genes$strand == 1, "+", "-"),
                gene_id = genes$ensembl_gene_id,
                symbol = genes$hgnc_symbol,
                ensembl_id = genes$ensembl_gene_id)

  # Preview result
  wh = gr[gr$symbol %in% candidates$geneID[candidates$geneID !=""]]
  wh = wh[wh@seqnames %in% paste0('chr',candidates$chromosome),]

  # genemodel
  exons_df = biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_exon_id",
      "chromosome_name",
      "exon_chrom_start",
      "exon_chrom_end",
      "strand",
      "rank"
    ),
    filters = "hgnc_symbol",
    values = candidates$geneID,
    mart = ensembl
  )

  # Rename columns for clarity
  colnames(exons_df) = c(
    "gene_id", "transcript_id", "exon_id",
    "chromosome", "start", "end",
    "strand", "exon_rank"
  )

  # Convert strand from numeric (1/-1) to "+" or "-"
  exons_df$strand = ifelse(exons_df$strand == 1, "+", "-")

  # Filter to MANE-selected transcripts
  exons_df = exons_df[exons_df$transcript_id %in% unique_transcript_id, ]

  # Create GRanges object
  gr_exons = GenomicRanges::GRanges(
    seqnames = exons_df$chromosome,
    ranges = IRanges::IRanges(start = exons_df$start, end = exons_df$end),
    strand = exons_df$strand,
    gene_id = exons_df$gene_id,
    transcript_id = exons_df$transcript_id,
    exon_id = exons_df$exon_id,
    exon_rank = exons_df$exon_rank
  )

  # return output
  list(gr_exons,wh)
}
