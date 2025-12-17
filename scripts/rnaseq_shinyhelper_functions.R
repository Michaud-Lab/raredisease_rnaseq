#plotting
plotting_coverage = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             zoom='full gene',
                             gene_annotations='gene_annotations') {
    
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

  #store the results
  candidate_gene_model = list()
  signif = list()
  plotCov = list()
  plotCov_v2 = list()
  plotNULL = ggplot() +  theme_void() + geom_text(aes(0,0,label='no mutation')) + xlab(NULL);
  
  #for(p in zoom){
  p = c(1,1000,5000)[zoom == c('full gene','mutation zoom (+/-1kb)','mutation zoom (+/-5kb)')]

  #write the gene model
  merged_exons_df <- as.data.frame(gr_exons)[, c("seqnames", "start", "end", "strand","exonsnb")]
  colnames(merged_exons_df) <- c("chromosome", "start", "end", "strand","exonsnb")
      
  if(p == 1) merged_exons_df$exonsnb[!(1:nrow(merged_exons_df) %in% seq(1,nrow(merged_exons_df),by = 5))]= ""

  #check if you have 0,1 or 2 mutations
  mut_pos = as.numeric(strsplit(candidate$position,'_')[[1]])

  if(length(mut_pos) == 0) {xintercept = (min(merged_exons_df$start)+max(merged_exons_df$end))/2}
  if(length(mut_pos) == 0) mut_pos = xintercept
  mutation = candidate$mutation
  
  #for loop for the number of mutations
  for(m in 1:length(mut_pos)) {
    if(p == 1) xintercept = mut_pos
    if(p >1) xintercept = mut_pos[m]
    xlims = c(min(merged_exons_df$start),max(merged_exons_df$end))
    if(p >1) {xlims =  c(max(mut_pos[m]-p,xlims[1]),min(mut_pos[m]+p,xlims[2]))}
    xbreaks = signif(c(xlims[1],mean(xlims),xlims[2]),6)
    
    #gene model
    candidate_gene_model[[m]] = merged_exons_df %>%
      ggplot(aes(xstart = start,xend = end,y = '')) +
      geom_range(aes(fill = 'red')) +
      geom_intron(data = to_intron(merged_exons_df),aes(strand = strand),arrow.min.intron.length = 100) +
      geom_text(aes(x = end,vjust = -3, label = exonsnb),size = 6) +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      ylab(candidate$geneID) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1])) +
      annotate('text', x = xintercept,y = 0.6, label = ifelse(mutation == '','',paste0(mutation,' (Position: ',xintercept,')')), col = 'darkblue', vjust = 0, hjust = 0.8, size = 5) +
      xlim(xlims) +
      ggtitle('Gene Model') +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
    #plots
    #read data
    res_dt_candidate_gene = read.csv(res_dt_candidate_gene_file,row.names = 1)
    
    signif[[m]] = ggplot(res_dt_candidate_gene,aes(x = mean, y = minuslogpval,color = minuslogpval)) +
      geom_point() +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      xlim(xlims) +
      ylab(bquote(-log[10]~(italic(p-value)))) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1])) +    
      scale_color_continuous(palette = c('black','red')) +
      ggtitle('Aberrant splicing') +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
    
    #if p == 2 lets make the plot a bit prettier
      if(p > 1) {
        res_dt_candidate_gene_subset = res_dt_candidate_gene[res_dt_candidate_gene$start > (xintercept-p) & res_dt_candidate_gene$end < (xintercept+p),]
 
        signif[[m]] = ggplot(res_dt_candidate_gene_subset,aes(x = start, xend = end, yend = minuslogpval, y = minuslogpval, color = minuslogpval)) +
          geom_segment(linewidth = 1) +
          geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
          xlim(xlims) +
          ylim(c(0,max(res_dt_candidate_gene$minuslogpval))) +
          ylab(bquote(-log[10]~(italic(p-value))))+
          scale_color_continuous(palette = c('black','red'),limits = c(0,max(res_dt_candidate_gene$minuslogpval))) +
          xlab(paste0('Chromosome ',merged_exons_df[1,1])) +
          ggtitle('Aberrant splicing (regions)') +
          theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
      }
    
    #depth
    depth = read.table(depth_file, header = T, check.names = F,comment.char= '')
    colnames(depth) = gsub("^.*/", "",colnames(depth))
    colnames(depth) = gsub('_sorted_chrN.bam','',colnames(depth))
    
    for(d in 1:nrow(colmean_genes_counts)){
      depth[,d+2] = depth[,d+2]/colmean_genes_counts[d,1]
    }
    
    depth = depth %>% filter(POS >= min(xlims),POS < max(xlims))
    ylims = c(0,max(c(rowMedians(as.matrix(depth[,-c(1:2)])),depth[,colnames(depth) %in% candidate$proband])))
    
    #pivoted
    depth_pivoted = depth %>% pivot_longer(cols = c(3:ncol(depth)), names_to = 'PatientID',values_to = 'Coverage')

    #plot
    plotCov_v2[[m]] =
      ggplot( depth_pivoted[depth_pivoted$PatientID != candidate$proband,],aes(x = POS, y = Coverage)) +
      stat_summary(geom="ribbon", fun.data=median_hilow,fun.args = list(conf.int=.5),fill=alpha('darkorange1', alpha =0.7),col= alpha('darkorange1', alpha =0.7)) +
      geom_line(data = depth_pivoted[depth_pivoted$PatientID == candidate$proband,],aes(x = POS, y = Coverage),col = 'black') +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      ylab('Normalised coverage') +
      xlab(paste0('Chromosome ',merged_exons_df[1,1])) +
      ggtitle('Coverage (proband in black, 25-75th reference percentiles in orange)') +
      scale_x_continuous(breaks = xbreaks) +
      theme(plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
    #fail safe in case there is no mutation
    if(mutation == "" & p>1) {
      candidate_gene_model[[m]] = plotNULL;
      signif[[m]] = plotNULL;
      plotCov[[m]] = plotNULL;
      plotCov_v2[[m]] = plotNULL}
  }
  
  #list outputs
  if(m==1) {output = list(list((candidate_gene_model[[1]]/signif[[1]]/plotCov_v2[[1]]) + plot_layout(heights = c(1,1,3))))}
  if(m == 2 & p>1 ) {output = list(list((candidate_gene_model[[1]]/signif[[1]]/plotCov_v2[[1]] + plot_layout(heights = c(1,1,3))) | (candidate_gene_model[[2]] / signif[[2]] / plotCov_v2[[2]] + plot_layout(heights = c(1,1,3)))))}
  
  if(m == 2 & p==1) {output = list(list((candidate_gene_model[[1]]/signif[[1]]/plotCov_v2[[1]]) + plot_layout(heights = c(1,1,3))))}

  output
  }
  


