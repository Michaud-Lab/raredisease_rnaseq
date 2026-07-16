# =============================================================================
# consensus.R - Build consensus/alternate sequence for splicing visualisation
# =============================================================================
source("../rnaseq_helper_functions.R")
load_install_library('dplyr')

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
params = list(workdir = args[1])
params$consensus = file.path(args[1], 'consensus/')
params$FRASER = file.path(args[1], 'FRASER')
params$fc_exons = args[2]
params$ref_file = args[3]
params$ref_annot = args[4]

# -----------------------------------------------------------------------------
# 2. candidate genes to run
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(params$workdir, 'data/input/candidate_genes_ALL.csv'))
#candidates_extra = read.table(file.path(params$workdir, 'data/input/candidate_genes_extra.csv'),comment.char = "#",header = T ,sep = ',');candidates_extra[is.na(candidates_extra)] = ''
#candidates_extra = candidates_extra[, colnames(candidates)]

#if(file.exists(file.path(params$workdir, 'data/input/candidate_genes_automated.csv'))) {
#  candidate_genes_automated = read.csv(file.path(params$workdir, 'data/input/candidate_genes_automated.csv'));candidate_genes_automated[is.na(candidate_genes_automated)] = ''
#} else {candidate_genes_automated = NULL}

#candidates = rbind(candidates,candidates_extra,candidate_genes_automated) %>%
#  distinct(geneID,ensembl, proband, .keep_all = TRUE)


# -----------------------------------------------------------------------------
# 3. consensus_pipeline function
# -----------------------------------------------------------------------------
consensus_pipeline = function(candidates = candidates, i = 1){

  # Identify the MANE transcript for the candidate gene
  dir.create(params$consensus, showWarnings = FALSE)
  fc_exons_tpm = read.table(params$fc_exons, sep = '\t', check.names = FALSE)
  unique_transcript = unique(fc_exons_tpm[, 1:3])
  params$transcript = unique_transcript$transcriptID[unique_transcript$geneID == candidates$geneID[i]]
  
  chr = candidates$chromosome[i]
  start = candidates$start[i] - 5000
  if (start <= 0) start = 1
  stop = candidates$stop[i] + 5000
  params$geneID = candidates$geneID[i]
  params$out_dir = paste0(params$FRASER,'/bams_subset/gene',candidates$geneID[i],'_chr',chr,'_',start,'_',stop,'/')
  params$bam_file = paste0(candidates$proband[i],'_sorted_chrN.bam')
  params$gene_variants_annotated = paste0(params$out_dir,'/gene',candidates$geneID[i],'variants_annotated.tsv')
  params$region = paste0('chr',chr,":",candidates$start[i],"-",candidates$stop[i])
  params$fasta_out = paste0(params$consensus,'gene',params$geneID,'_',candidates$proband[i],".fasta")
  
  # -----------------------------------------------------------------------------
  # 2. Build consensus sequence
  # -----------------------------------------------------------------------------
  if (params$geneID != "") {
    command = paste('./consensus.sh', params$geneID, params$transcript, params$ref_file,
                    params$ref_annot, params$out_dir, params$bam_file, params$consensus, params$region)
    
    if (!file.exists(params$fasta_out)) {
      load_install_library(c('dplyr', 'seqinr'))
      
      # Write placeholder FASTA to avoid re-running if the analysis fails
      write.table(c('>reference', 'ACTG'), params$fasta_out, row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      system(command)
      
      # format the outputs of the bash script
      variant_annotated = read.table(params$gene_variants_annotated,header = TRUE,comment.char ='',check.names=FALSE,na.strings = "",sep = '\t')
      colnames(variant_annotated) = c('chr','position','depth','reference','alternate','transcript','exon','gene')
      variant_annotated$alternate[variant_annotated$alternate =='.'] = variant_annotated$reference[variant_annotated$alternate =='.']
      variant_annotated$transcript = gsub(';','',variant_annotated$transcript)
      variant_annotated$transcript[is.na(variant_annotated$transcript)] = 'intron'
      variant_annotated$exon = as.numeric(gsub(';','',variant_annotated$exon))
      variant_annotated$gene = gsub(';','',variant_annotated$gene)
      
      # annotate the introns
      variant_rle  = rle(variant_annotated$transcript)
      variant_rle = data.frame(values = variant_rle$values, cumsum = cumsum(variant_rle$lengths),lengths = variant_rle$lengths,event = variant_rle$values)
      count = 1
      for (n in 1:nrow(variant_rle)) {
        if (variant_rle$event[n] == "intron") {
          variant_rle$values[n] = paste0('intron', count)
          count = count + 1
        }
      }
      
      variant_rle = variant_rle[c(1,1:nrow(variant_rle)),]
      variant_rle[1,2:3] = c(0,0)
      
      for(n in 2:nrow(variant_rle)){
        if(variant_rle$event[n] == "intron") {variant_annotated$exon[(variant_rle$cumsum[n-1]+1) : variant_rle$cumsum[n]] = variant_rle$values[n]}
      }
      
      # mean expression per intronic position
      mean_intronic_depth = variant_annotated %>% filter(transcript == 'intron') %>% group_by(exon) %>% summarise(mean = mean(depth))
      
      # mean expression per exonic position
      exonic = variant_annotated[variant_annotated$transcript !='intron',]
      exclude_first_and_last = unique(exonic$exon)
      exonic = exonic[exonic$exon %in% exclude_first_and_last[-c(1,length(exclude_first_and_last))],]
      exonic_depth = exonic %>% group_by(exon) %>% summarise(mean = mean(depth))
      mean_exonic_depth = mean(exonic_depth$mean)
      if(is.na(mean_exonic_depth)) mean_exonic_depth = 0
      
      # Flag intron retention: intron mean depth > 30% of mean exon depth
      variant_annotated$retention_event = 0
      if (nrow(mean_intronic_depth) > 0) {
        for (r in 1:nrow(mean_intronic_depth)) {
          if (mean_intronic_depth$mean[r] > (mean_exonic_depth * 0.3)) {
            variant_annotated$retention_event[variant_annotated$exon == mean_intronic_depth$exon[r]] = 1
          }
        }
      }
      
      # Flag exon skipping: exon mean depth < 30% of mean exon depth
      variant_annotated$skipping_event = 0
      if (nrow(exonic_depth) > 0) {
        for (r in 1:nrow(exonic_depth)) {
          if (exonic_depth$mean[r] < (mean_exonic_depth * 0.3)) {
            variant_annotated$skipping_event[variant_annotated$exon == exonic_depth$exon[r]] = 1
          }
        }
      }
      
      # Highlight variant position, retained introns (lowercase), and skipped exons
      if (candidates$position[i] != '') {
        variant_annotated$alternate[variant_annotated$position == candidates$position[i]] =
          paste0("<span style='color: blue;'>",
                 variant_annotated$alternate[variant_annotated$position == candidates$position[i]],
                 "</span>")
      }
      variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate] = paste0("<b>",variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate],"</b>")
      variant_annotated$alternate[variant_annotated$retention_event == 1] = tolower(variant_annotated$alternate[variant_annotated$retention_event == 1])
      variant_annotated$reference[variant_annotated$retention_event == 1] = '-'
      
      # skipping events
      variant_annotated$reference[variant_annotated$skipping_event == 1] = tolower(variant_annotated$reference[variant_annotated$skipping_event == 1])
      variant_annotated$alternate[variant_annotated$skipping_event == 1] = '-'
      
      # sequences
      reference_sequence = variant_annotated$reference[variant_annotated$transcript != 'intron' | variant_annotated$retention_event == 1]
      alternate_sequence = variant_annotated$alternate[variant_annotated$transcript != 'intron' | variant_annotated$retention_event == 1]
      
      print(paste0('Length REF: ', length(reference_sequence[reference_sequence != '-'])))
      print(paste0('Length ALT: ', length(alternate_sequence[alternate_sequence != '-'])))
      
      # Write FASTA output
      my_sequences = list(reference = reference_sequence, alternate = alternate_sequence)
      write.fasta(sequences = my_sequences, names = names(my_sequences),
                  nbchar = 100, file.out = params$fasta_out, open = "w")
      
      print(paste0('Done sample ', i, ' ~~~  ', params$geneID, ' ~~~ Time is: ', Sys.time()))
      
    } else {
      print(paste0('Sample ~~~ ', i, ' already exists ~~~ Time is: ', Sys.time()))
    }
  } else {
    print(paste0('Sample ~~~ ', i, ' did not contain a candidate gene ~~~ Time is: ', Sys.time()))
  }
}


###run the consensus_pipeline()
for(i in 1 : nrow(candidates)){consensus_pipeline(candidates,i=i)}




