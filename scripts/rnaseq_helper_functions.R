# =============================================================================
# Generate automatically new candidate genes based on FRASER, OUTRIDER, ASE.
# =============================================================================
candidate_genes_automated = function(gwfile = file.path(params$datadir, 'gwFRASER.csv')){
  #Generate the gene annotation  
  gene_annotations = gene_annotation(full =T)
  
  #Load files
  if(grepl('gwFRASER',gwfile)) {sep = ','} else {sep = '\t'}
  
  gw = read.csv(gwfile, row.names = 1,sep = sep)

  #get the candidate genes
  if(grepl('gwFRASER',gwfile)) {
    gw_top = gw %>%
      filter(!grepl("^HBA|^HBB|^HLA|^HBG|^HBD|^HBB|^HBQ|^HBE|^HBZ|^HBM", hgncSymbol), !is.na(hgncSymbol)) %>%
      group_by(sampleID) %>%
      filter(padjust < 0.0001) %>%
      slice_min(padjust, n = 5) %>%
      distinct(hgncSymbol, sampleID,.keep_all = T)
  
      gw_top$hgncSymbol = sapply(strsplit(gw_top$hgncSymbol,';'),'[[',1)
      gw_top$origin = paste('gw FRASER',gw_top$padjust,gw_top$deltaPsi,sep = ', ')
      gw_top = gw_top[,c('hgncSymbol','sampleID','origin')]
  }
  
  #get the candidate genes
  if(grepl('OUTRIDER',gwfile)) {
    gw_top = gw %>%
      filter(!grepl("^HBA|^HBB|^HLA|^HBG|^HBD|^HBB|^HBQ|^HBE|^HBZ|^HBM", geneID), !is.na(geneID)) %>%
      group_by(sampleID) %>%
      filter(pValue < 0.000001) %>%
      slice_min(pValue, n = 2) %>%
      distinct(geneID, sampleID,.keep_all = T)
    
    gw_top$geneID = sapply(strsplit(gw_top$geneID,';'),'[[',1)
    gw_top$origin = paste('gw OUTRIDER',gw_top$pValue,gw_top$l2fc,sep = ', ')
    gw_top = gw_top[,c('geneID','sampleID','origin')]
  }
  
  #get the candidate genes
  if(grepl('gwASE',gwfile)) {
    gw_top = gw %>%
      filter(!grepl("^HBA|^HBB|^HLA|^HBG|^HBD|^HBB|^HBQ|^HBE|^HBZ|^HBM", geneID), !is.na(geneID)) %>%
      group_by(sampleID,geneID) %>%
      filter(pvalue < 1e-49) %>%
      dplyr::filter(dplyr::n() >= 2) %>%
      distinct(geneID, sampleID,.keep_all = T)
    
    gw_top$geneID = sapply(strsplit(gw_top$geneID,';'),'[[',1)
    gw_top$origin = paste('gw ASE (at least 2 markers)',gw_top$pvalue,gw_top$RNA_DP,sep = ', ')
    gw_top = gw_top[,c('geneID','sampleID','origin')]
  }
  
  #format them to the candidate format.
  candidates_automated = data.frame(matrix(ncol = 10, nrow = 0))
  colnames(candidates_automated) = c('geneID','ensembl','proband','chromosome','start','stop','proband2','mutation','position','origin')
  candidates_automated[1:nrow(gw_top),c(1,3,10)] = gw_top
  candidates_automated[,4:6] = 1
  candidates_automated[,8:9] = ''
  candidates_automated$proband2 = gsub('_PAX','',candidates_automated$proband)
  
  for(i in 1:nrow(candidates_automated))
  {
    temp = gene_annotations[[2]][gene_annotations[[2]]$symbol == candidates_automated[i,1],]
    if(length(temp@seqnames)==1) {
      candidates_automated$start[i] = temp@ranges@start
      candidates_automated$stop[i] = temp@ranges@start + temp@ranges@width
      candidates_automated$ensembl[i] = temp$gene_id
      candidates_automated$chromosome[i] = as.numeric(gsub('chr','',temp@seqnames@values))}
  }
  candidates_automated = candidates_automated[!is.na(candidates_automated$chromosome),]
  candidates_automated = candidates_automated[candidates_automated$start != 1,]
  print(paste0('Done candidates_automated, found ', nrow(candidates_automated), ' new candidates'))
  return(candidates_automated)
}

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
    if((p == length(packages)) & silent == F) print(paste0('Finished installing/loading required packages: ',paste0(packages,collapse=', ')))
  }
}

# =============================================================================
# rnaseq_helper_functions.R - Helper functions for featureCounts processing
# =============================================================================
# gene_annotation: fetch exon-level gene models from Ensembl via biomaRt
# Returns a list: [[1]] GRanges of exons, [[2]] GRanges of gene bodies
gene_annotation = function(unique_transcript_id = unique(fc_exons_raw$transcriptID),
                           candidates = candidates, full = F){

  load_install_library(c('biomaRt', 'ggbio', 'GenomicAlignments'))

  # Connect to Ensembl and select the human dataset for hg38 (GRCh38)
  ensembl = biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')

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
  if(full == F) {
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
    )} else {wh = gr; gr_exons = NULL}

    # return output
    list(gr_exons,wh)
}
