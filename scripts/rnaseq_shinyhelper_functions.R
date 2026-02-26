###
###plotting coverage
###
plotting_coverage = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             xlims = c(10,100),
                             gene_annotations='gene_annotations') {

if(candidate$geneID != ""){

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
    if(file.exists(res_dt_candidate_gene_file)) {
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
    } else {
      signif = ggplot() + geom_blank()
      print('No FRASER pvalues available. Likely because expression is too low')
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
      xlab(paste0('Chromosome ',merged_exons_df[1,1], ' (Kb)')) +
      ggtitle('Coverage (proband in black, 25-75th reference percentiles in orange)') +
      theme(plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
    
  #list outputs
  output = list(list((candidate_gene_model/signif/plotCov_v2) + plot_layout(heights = c(1,1,3))))
  output
  } else {
  plot(0, main = 'no gene available');print('No gene available for gene model plot')
  }
}



###
###this is to generate a table of the top p-values of the splicing test.
###
gwFRASER_table = function(res_dt=gwFRASER,sample = 'HSJ_036_03_PAX',pcutoff=0.05,pvalue='padjust', geneID = 'hgncSymbol'){
 
  #factorize
  res_dt$chr = factor(res_dt$chr, levels = c(1:22,'X','Y','MT')) 
  
  #keep sample of interest
  res_dt_sampleID = res_dt[res_dt$sampleID == sample,]
  
  #keep only the top (but remove duplicated splicing events)
  gwFRASER_top = res_dt_sampleID[order(res_dt_sampleID[[pvalue]]),]
  gwFRASER_top = gwFRASER_top[gwFRASER_top[[pvalue]]<0.05,]
  gwFRASER_top[[geneID]][is.na(gwFRASER_top[[geneID]])] = 'na'
  gwFRASER_top = gwFRASER_top[!duplicated(gwFRASER_top[[geneID]], incomparables = 'na'),]
  gwFRASER_top = gwFRASER_top[order(gwFRASER_top$chr),c(1:5,7,9:10,13:14)]
  
  #return
  return(gwFRASER_top)
}


###
###this is to do a manhattan plot of the p-values of the splicing test.
###
manhattan_plot = function(res_dt=gwFRASER,sample = 'HSJ_036_03_PAX',top=25,pcutoff=0.05, pvalue='padjust', geneID = 'hgncSymbol',shape = FALSE){
  
  #shape factor 
  res_dt$shape = 'splicing'
  if(shape) res_dt$shape = ifelse(res_dt$l2fc >0,'over','under')

  #factorize
  res_dt$chr = factor(res_dt$chr, levels = c(1:22,'X','Y','MT')) 

  #keep only samples of interest
  res_dt_sampleID = res_dt[res_dt$sampleID == sample,]

  ##calculate cumulative chromosome sizes
  chr_size <- res_dt %>% 
    group_by(chr) %>% 
    summarise(chr_len=max(end)) %>%
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    select(-chr_len) 
  
  #Add this info to the initial dataset
  results_subset_forggplot = chr_size %>% 
    left_join(res_dt_sampleID, ., by=c("chr"="chr")) %>%
    arrange(chr, pos) %>%
    mutate(BPcum=pos+tot)

  #set colors
  colors = c(brewer.pal(8,'Set2'),brewer.pal(9,'Set1'),brewer.pal(9,'Set3'))
  colors = c(colors[c(1,9,2,11,3,10,4,12,5,6,13,7,15,8,16:18,20:26)],'black')
  
  ###prepare axis labels.
  axisdf = results_subset_forggplot %>%
    group_by(chr) %>%
    dplyr::summarize(center=( max(BPcum) + min(BPcum)+1 ) / 2 )
  
  #keep only the top (but remove duplicated splicing events)
  gw_top = results_subset_forggplot[order(results_subset_forggplot[[pvalue]]),]
  gw_top = gw_top[gw_top[[pvalue]]<pcutoff,]
  gw_top[[geneID]][is.na(gw_top[[geneID]])] = 'na'
  gw_top = gw_top[!duplicated(gw_top[[geneID]], incomparables = 'na'),]  
  gw_top = head(gw_top,top)
  gw_top = gw_top[order(gw_top$chr),]
  
  ###Make the ggplot:
  man_gplot = ggplot(results_subset_forggplot, aes(x = BPcum, y = -log10(.data[[pvalue]]))) +
    
    #Labels
    geom_label_repel(data = gw_top, aes(label = .data[[geneID]], x = BPcum, y = -log10(.data[[pvalue]])), col = 'black',  size = 4,box.padding = 2, max.overlaps = 50) +    
    
    #Show all points
    geom_point( aes(color=chr,shape=shape), alpha=1, size=1.6) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c("over"=17, "under"=6, "splicing"=19)) +
    
    #Custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = c(0.01,0.01)) +
    xlab('Chromosomes') + 
    ylab(paste('-log10 (',pvalue,')')) +
    ylim(0,NA) +
    
    #Customize theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 20)) +
    ggtitle(paste0('Outlier events ~ ',sample))
  
  return(man_gplot)
}


###
###Generate a gene prioritisation data.frame based on gene lists (HPO, outrider, fraser) 
###
gene_prioritization = function(sample = 'HSJ_001_03_PAX',top=100,hpo_sample=clinical,hpo_all=file.path(params$datadir,'genes_to_phenotype.txt'),fraser="",outrider="",custom_genes=""){
  
  #hpo 
  if(!file.exists(hpo_all)) {download.file(url='https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt',dest=file.path(hpo_all))} else {print(paste0('file ', hpo_all,' exists'))}
  if(!exists('hpo')) hpo = read.delim(hpo_all)

  #unique hpo terms
  hpo_un = hpo[!duplicated(hpo$hpo_name),]
  hpo_un$hpo_name_shortened = substring(hpo_un$hpo_name,1,29)
  hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)] = paste0(hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)],'.')
  
  #generate a named list that contains all the genes in all the HPO terms.    
  temp = strsplit(hpo_sample$`HPO terms`[hpo_sample$`Patient ID` == sample],split = '||',fixed = T)[[1]]
  temp = unlist(strsplit(temp,split = ': '))
  temp = temp[grepl('HP:',temp)]

  if(length(temp)>0) {
   temp = gsub(' ','',temp)
   hpo_genes = as.list(temp)
   hpo_terms = hpo_un[hpo_un$hpo_id %in% temp,]
   names(hpo_genes) = paste0(hpo_terms$hpo_id," ",hpo_terms$hpo_name_shortened)
   for(h in 1:length(hpo_genes)){hpo_genes[[h]] = unique(hpo$gene_symbol[hpo$hpo_id == hpo_terms[h,3]])}} else hpo_genes = list(hp_NULL='none')

  #add outrider
  if(!is.null(dim(outrider))) {
   outrider_temp = outrider[outrider$sampleID==sample,colnames(outrider) %in% c('geneID','pValue','zScore','exon_zScore','exon_pValue')]
   outrider_temp = outrider_temp[!is.na(outrider_temp$geneID),]    
   colnames(outrider_temp) = c('geneID','OUTRIDER gene pValue','OUTRIDER gene zScore','OUTRIDER exon pValue','OUTRIDER exon zScore')
  }
  
  #add fraser
  if(!is.null(dim(fraser))) {
    fraser_temp = fraser[fraser$sampleID==sample,colnames(fraser) %in% c('hgncSymbol','pValue')]
    fraser_temp = fraser_temp[!is.na(fraser_temp$hgncSymbol),]    
    colnames(fraser_temp) = c('geneID','FRASER gene pValue')
    fraser_temp = fraser_temp[!duplicated(fraser_temp$geneID),]
  }
  
  #merge Outliers
  outlier_temp = merge(fraser_temp,outrider_temp,by= 'geneID',all=T, sort=F)
   
  #generate a big table stating which gene is listed where.
  all_genes = unique(unlist(hpo_genes))
  table = sapply(hpo_genes, function(x) all_genes %in% x)
  table = data.frame('gene score' = rowSums(table),table*1, check.names = F)
  table$geneID = all_genes
  table = merge(table,outlier_temp,by= 'geneID',all.x= T, sort=F)
  table$`OUTRIDER gene pValue`[is.na(table$`OUTRIDER gene pValue`)] = 'ns'
  table$`OUTRIDER gene zScore`[is.na(table$`OUTRIDER gene zScore`)] = 'ns'
  table$`OUTRIDER exon pValue`[is.na(table$`OUTRIDER exon pValue`)] = 'ns'
  table$`OUTRIDER exon zScore`[is.na(table$`OUTRIDER exon zScore`)] = 'ns'
  table$`FRASER gene pValue`[is.na(table$`FRASER gene pValue`)] = 'ns'
  table = table[,c(1:2,ncol(table):3)]
  table = table[,c(1,2,5,6,7,3,4,8:ncol(table))]
  table = table[order(table$`gene score`,decreasing = T),]
 

  #return top hist
  return(head(table,top))
}
