####
gene_annotation = function(unique_transcript_id = unique(fc_exons_raw$transcriptID),
                           candidates = candidates){

  suppressMessages(suppressWarnings(library(biomaRt)))
  suppressMessages(suppressWarnings(library(ggbio)))
  suppressMessages(suppressWarnings(library(GenomicAlignments)))  

  # Connect to Ensembl and select the human dataset for hg38 (GRCh38)
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 38)
  
  # Get gene annotations
  genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position",
                                "strand", "ensembl_gene_id", "hgnc_symbol", "gene_biotype"),mart = ensembl)
  
  ##subset only  genes on an actual chromosome (i.e. not a scaffold or something weird)
  genes = genes[genes$chromosome_name %in% c(1:22,'X','Y','MT'),]
  
  # Create GRanges object
  gr <- GRanges(seqnames = paste0("chr", genes$chromosome_name),
                ranges = IRanges(start = genes$start_position, end = genes$end_position),
                strand = ifelse(genes$strand == 1, "+", "-"),
                gene_id = genes$ensembl_gene_id,
                symbol = genes$hgnc_symbol,
                ensembl_id = genes$ensembl_gene_id)
  
  # Preview result
  wh = gr[gr$symbol %in% candidates$geneID]
  wh = wh[wh@seqnames@values %in% paste0('chr',candidates$chromosome),]

  #genemodel
  exons_df <- getBM(
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
  
  #Rename columns for clarity
  colnames(exons_df) <- c(
    "gene_id", "transcript_id", "exon_id",
    "chromosome", "start", "end",
    "strand", "exon_rank"
  )
  
  # Convert strand from numeric (1/-1) to "+" or "-"
  exons_df$strand <- ifelse(exons_df$strand == 1, "+", "-")
  
  #filter proper Chr and transcript
  exons_df = exons_df[exons_df$transcript_id %in% unique_transcript_id,]
  
  # Create GRanges object
  gr_exons <- GRanges(
    seqnames = exons_df$chromosome,
    ranges = IRanges(start = exons_df$start, end = exons_df$end),
    strand = exons_df$strand,
    gene_id = exons_df$gene_id,
    transcript_id = exons_df$transcript_id,
    exon_id = exons_df$exon_id,
    exon_rank = exons_df$exon_rank
  )
  
  #return output
  list(gr_exons,wh)
}
