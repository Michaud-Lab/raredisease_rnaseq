#plotting
plotting_coverage = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             xlims = c(10,100),
                             gene_annotations='gene_annotations') {

print(paste0('Defined zoom limits: ', xlims[1],' --- ',xlims[2], ' Kb'))
      
#colmean 
colmean_genes_counts = read.table(colmean_genes_counts_file) 
  
#filter proper Chr and transcript
gr_exons = gene_annotations[[1]]
gr_exons = gr_exons[gr_exons$gene_id == candidate$ensembl,]
gr_exons = gr_exons[order(gr_exons$exon_rank),]
gr_exons$exonsnb = paste0('e',1:length(gr_exons@strand))

#wh object for autoplot coverage further down
wh = gene_annotations[[2]]
wh = wh[wh$gene_id == candidate$ensembl,]


  candidate_limit = c(floor(candidate$start/1000),ceiling(candidate$stop/1000))

  #write the gene model
  merged_exons_df <- as.data.frame(gr_exons)[, c("seqnames", "start", "end", "strand","exonsnb")]
  colnames(merged_exons_df) <- c("chromosome", "start", "end", "strand","exonsnb")
  merged_exons_df$start = merged_exons_df$start/1000
  merged_exons_df$end = merged_exons_df$end/1000
  
  if(all.equal(xlims,candidate_limit) ==TRUE) merged_exons_df$exonsnb[!(1:nrow(merged_exons_df) %in% seq(1,nrow(merged_exons_df),by = 5))]= ""
  
  mut_pos = as.numeric(strsplit(candidate$position,'_')[[1]])
  xintercept = mut_pos/1000
  if(length(mut_pos) == 0) {xintercept = (min(merged_exons_df$start)+max(merged_exons_df$end))/2}
  if(length(mut_pos) == 0) mut_pos = xintercept
  mutation = candidate$mutation
 
    #gene model
    candidate_gene_model = merged_exons_df %>%
      ggplot(aes(xstart = start,xend = end,y = '')) +
      geom_range(aes(fill = 'red')) +
      geom_intron(data = to_intron(merged_exons_df),aes(strand = strand),arrow.min.intron.length = 100) +
      geom_text(aes(x = end,vjust = -3, label = exonsnb),size = 6) +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      ylab(candidate$geneID) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
      annotate('text', x = xintercept,y = 0.6, label = ifelse(mutation == '','',paste0(mutation,' (Position: ',mut_pos,')')), col = 'darkblue', vjust = 0, hjust = 0.8, size = 5) +
      xlim(xlims) +
      ggtitle('Gene Model') +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
    #plots
    res_dt_candidate_gene = read.csv(res_dt_candidate_gene_file,row.names = 1)
    res_dt_candidate_gene$mean = res_dt_candidate_gene$mean/1000
    res_dt_candidate_gene$start = res_dt_candidate_gene$start/1000
    res_dt_candidate_gene$end = res_dt_candidate_gene$end/1000
      
    signif = ggplot(res_dt_candidate_gene,aes(x = mean, y = minuslogpval,color = minuslogpval)) +
      geom_point() +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      xlim(xlims) +
      ylab(bquote(-log[10]~(italic(p-value)))) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +    
      scale_color_continuous(palette = c('black','red')) +
      ggtitle('Aberrant splicing') +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
    
    #if the region is less than 5kb lets make the plot a bit prettier
      if(xlims[2]-xlims[1] < 6) {
        res_dt_candidate_gene_subset = res_dt_candidate_gene[res_dt_candidate_gene$start > (xlims[1]) & res_dt_candidate_gene$end < (xlims[2]),]
 
        signif = ggplot(res_dt_candidate_gene_subset,aes(x = start, xend = end, yend = minuslogpval, y = minuslogpval, color = minuslogpval)) +
          geom_segment(linewidth = 1) +
          geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
          xlim(xlims) +
          ylim(c(0,max(res_dt_candidate_gene$minuslogpval))) +
          ylab(bquote(-log[10]~(italic(p-value))))+
          scale_color_continuous(palette = c('black','red'),limits = c(0,max(res_dt_candidate_gene$minuslogpval))) +
          xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
          ggtitle('Aberrant splicing (regions)') +
          theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
      }
    
    #depth
    depth = read.table(depth_file, header = T, check.names = F,comment.char= '')
    depth$POS = depth$POS / 1000
    colmean_genes_counts = colmean_genes_counts
    colnames(depth) = gsub("^.*/", "",colnames(depth))
    colnames(depth) = gsub('_sorted_chrN.bam','',colnames(depth))
    
    for(d in 1:nrow(colmean_genes_counts)){
      depth[,d+2] = depth[,d+2]/colmean_genes_counts[d,1]
    }
    
    depth_filtered = depth %>% filter(POS >= min(xlims),POS < max(xlims))
    if(nrow(depth_filtered) == 0) depth_filtered = head(depth,100)
    
    ylims = c(0,max(c(rowMedians(as.matrix(depth_filtered[,-c(1:2)])),depth_filtered[,colnames(depth_filtered) %in% candidate$proband])))
    
    #pivoted
    depth_pivoted = depth_filtered %>% pivot_longer(cols = c(3:ncol(depth_filtered)), names_to = 'PatientID',values_to = 'Coverage')

    #plot
    plotCov_v2 =
      ggplot(depth_pivoted[depth_pivoted$PatientID != candidate$proband,],aes(x = POS, y = Coverage)) +
      stat_summary(geom="ribbon", fun.data=median_hilow,fun.args = list(conf.int=.5),fill=alpha('darkorange1', alpha =0.7),col= alpha('darkorange1', alpha =0.7)) +
      geom_line(data = depth_pivoted[depth_pivoted$PatientID == candidate$proband,],aes(x = POS, y = Coverage),col = 'black') +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      xlim(xlims) +
      ylab('Normalised coverage') +
      xlab(paste0('Chromosome ',merged_exons_df[1,1], ' Kb')) +
      ggtitle('Coverage (proband in black, 25-75th reference percentiles in orange)') +
      theme(plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
  #list outputs as patchwork ggplot2 object.
  list(list((candidate_gene_model/signif/plotCov_v2) + plot_layout(heights = c(1,1,3))))
  }
  


