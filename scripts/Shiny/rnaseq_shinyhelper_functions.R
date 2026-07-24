# plot_hb_fraction: stacked bar chart of haemoglobin gene read fractions per sample.
# Arguments:
#   fc_genes_raw_ALL - data.frame: raw gene counts with columns geneID, ensemblID, then one column per sample
#   hb_genes         - character vector: gene names to include (default: HBA1, HBA2, HBB, HBG1, HBG2)
plot_hb_fraction = function(fc_genes_raw_ALL,
                            hb_genes = c('HBA1', 'HBA2', 'HBB', 'HBG1', 'HBG2','HBD','HBM','HBZ','HBE','HBQ')) {
  non_sample_cols = c('geneID', 'ensemblID', 'Length')
  sample_cols = !colnames(fc_genes_raw_ALL) %in% non_sample_cols
  hb = fc_genes_raw_ALL[fc_genes_raw_ALL$geneID %in% hb_genes, ]
  col_totals = colSums(fc_genes_raw_ALL[, sample_cols], na.rm = TRUE)
  hb_frac = hb
  hb_frac[, sample_cols] = sweep(hb[, sample_cols], 2, col_totals, '/') * 100
  hb_long = pivot_longer(hb_frac, cols = all_of(colnames(hb_frac)[sample_cols]),
                                names_to = 'sample', values_to = 'fraction')
  hb_counts_long = pivot_longer(hb, cols = all_of(colnames(hb)[sample_cols]),
                                       names_to = 'sample', values_to = 'raw_counts')
  plot_ly(hb_long, x = ~sample, y = ~fraction, color = ~geneID, type = 'bar',
          customdata = ~ round(hb_counts_long$raw_counts,0),
          hovertemplate = paste(
            '<b>Sample</b>: %{x}',
            '<br><b>Percentage</b>: %{y:.2f}%',
            '<br><b>Raw counts</b>: %{customdata}',
            '<extra></extra>'
          )) %>%
    layout(
      barmode = 'stack',
      xaxis = list(title = 'Sample', tickangle = -45),
      yaxis = list(title = '% of total reads'),
      legend = list(title = list(text = 'Gene'), orientation = 'h',
                    x = 0, xanchor = 'left', y = 1, yanchor = 'bottom')
    )
}

# plot_total_reads: bar chart of total read counts per sample from fc_genes_raw_ALL.
# Arguments:
#   fc_genes_raw_ALL - data.frame: raw gene counts with columns geneID, ensemblID, then one column per sample
plot_total_reads = function(fc_genes_raw_ALL) {
  non_sample_cols = c('geneID', 'ensemblID', 'Length')
  sample_cols = !colnames(fc_genes_raw_ALL) %in% non_sample_cols
  col_totals = colSums(fc_genes_raw_ALL[, sample_cols], na.rm = TRUE)
  df = data.frame(sample = names(col_totals), total_reads = as.numeric(col_totals))
  plot_ly(df, x = ~sample, y = ~total_reads, type = 'bar',
          hovertemplate = paste(
            '<b>Sample</b>: %{x}',
            '<br><b>Total reads</b>: %{y:,}',
            '<extra></extra>'
          )) %>%
    layout(
      xaxis = list(title = 'Sample', tickangle = -45),
      yaxis = list(title = 'Total reads'),
      showlegend = FALSE
    )
}

# plot_expression_cohort: per-exon TPM box plots for the cohort (children / adults)
#   overlaid with the proband's individual data points.
# Arguments:
#   data_child   - data.frame: cohort rows for patients < 18 yrs
#   data_adults  - data.frame: cohort rows for patients >= 18 yrs
#   data_patient - data.frame: rows for the selected proband
plot_expression_cohort = function(data_child, data_adults, data_patient) {
  n_exons = max(data_patient$exonID, 1)
  plot_ly() %>%
    add_trace(
      data = data_child,
      x = ~exonID,
      y = ~expression,
      type = "box",
      boxpoints = FALSE,
      name = 'Average (<18yrs)',
      color = I('darkblue'),
      hoverinfo = 'none',
      marker = list(size = 12)
    ) %>%
    add_trace(
      data = data_adults,
      x = ~exonID,
      y = ~expression,
      type = "box",
      boxpoints = FALSE,
      name = 'Average (>18yrs)',
      color = I('darkgreen'),
      hoverinfo = 'none',
      marker = list(size = 12)
    ) %>%
    add_trace(
      data = data_patient,
      x = ~exonID,
      y = ~expression,
      type = 'scatter',
      mode = 'markers',
      hoverinfo = 'text',
      marker = list(opacity = 1),
      color = I('darkorange'),
      size = 22,
      name = ~PatientID,
      customdata = ~Sexe,
      text = ~as.factor(age),
      hovertemplate = paste(
        '<b>Exon number</b>: %{x}',
        '<br><b>Proband Age</b>: %{text}',
        '<br><b>Sex</b>: %{customdata}<br>'
      )
    ) %>%
    layout(
      xaxis = list(
        title = "Exon Number",
        ticktext = 1:n_exons,
        tickvals = as.list(1:n_exons),
        tickmode = "array"
      ),
      yaxis = list(title = "Normalised Expression (TPM)"),
      title = paste0('Cohort expression (per exon) for ', data_patient$geneID[1])
    )
}


# plot_expression_family: per-exon TPM scatter plot for all members of the proband's family.
# Arguments:
#   data_family  - data.frame: rows for all family members (proband + parents)
#   data_patient - data.frame: rows for the selected proband (used for gene name in title)
plot_expression_family = function(data_family, data_patient) {
  n_exons = max(data_family$exonID, 1)
  n_members = length(unique(data_family$PatientID))
  plot_ly() %>%
    add_trace(
      data = data_family,
      x = ~exonID,
      y = ~expression,
      type = 'scatter',
      mode = 'markers',
      hoverinfo = 'text',
      marker = list(opacity = 1),
      color = ~PatientID,
      colors = c('darkorange', 'black', 'blue')[n_members:1],
      size = 22,
      name = ~PatientID,
      customdata = ~Sexe,
      text = ~as.factor(age),
      hovertemplate = paste(
        '<b>Exon number</b>: %{x}',
        '<br><b>Age</b>: %{text}',
        '<br><b>Sex</b>: %{customdata}<br>'
      )
    ) %>%
    layout(
      xaxis = list(
        title = "Exon Number",
        ticktext = 1:n_exons,
        tickvals = as.list(1:n_exons),
        tickmode = "array"
      ),
      yaxis = list(title = "Normalised Expression (TPM)"),
      title = paste0('Family expression (per exon) for ', data_patient$geneID[1])
    )
}


# genemodel_plot: builds the gene-model, splicing significance, and coverage plot stack for a candidate gene.
# Arguments:
#   candidate                  - data.frame row: the candidate gene/proband row (geneID, ensembl, start, stop, position, mutation, proband)
#   res_dt_candidate_gene_file - character: path to the FRASER per-candidate-gene results CSV
#   depth_file                 - character: path to the per-base coverage depth file
#   bam_file                   - character: path to the proband's BAM file (kept for interface consistency)
#   colmean_genes_counts_file  - character: path to the column-mean gene counts file (used to normalise coverage)
#   xlims                      - numeric vector length 2: zoom window in Kb (default: c(10, 100))
#   gene_annotations           - list: gene annotation object; [[1]] exons GRanges, [[2]] gene GRanges ("wh")
genemodel_plot = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             xlims = c(10,100),
                             gene_annotations='gene_annotations') {

if(candidate$geneID != "" & file.exists(res_dt_candidate_gene_file)){

print(paste0('Defining zoom limits: ', xlims[1],' --- ',xlims[2], ' Kb'))

# colmean
colmean_genes_counts = read.table(colmean_genes_counts_file)

# filter proper Chr and transcript
gr_exons = gene_annotations[[1]]
gr_exons = gr_exons[gr_exons$gene_id == candidate$ensembl,]
gr_exons = gr_exons[order(gr_exons$exon_rank),]
gr_exons$exonsnb = paste0('e',1:length(gr_exons@strand))

# wh object for autoplot coverage further down
wh = gene_annotations[[2]]
wh = wh[wh$gene_id == candidate$ensembl,]


  candidate_limit = c(floor(candidate$start/1000),ceiling(candidate$stop/1000))

  # write the gene model
  merged_exons_df = as.data.frame(gr_exons)[, c("seqnames", "start", "end", "strand","exonsnb")]
  colnames(merged_exons_df) = c("chromosome", "start", "end", "strand","exonsnb")
  merged_exons_df$start = merged_exons_df$start/1000
  merged_exons_df$end = merged_exons_df$end/1000

  if(all.equal(xlims,candidate_limit) ==TRUE) merged_exons_df$exonsnb[!(1:nrow(merged_exons_df) %in% seq(1,nrow(merged_exons_df),by = 5))]= ""

  mut_pos = as.numeric(strsplit(as.character(candidate$position), '_')[[1]])
  mut_pos = mut_pos[!is.na(mut_pos)]
  xintercept = mut_pos/1000
  if(length(mut_pos) == 0) {xintercept = (min(merged_exons_df$start)+max(merged_exons_df$end))/2}
  if(length(mut_pos) == 0) mut_pos = xintercept
  mutation = candidate$mutation

    # gene model
    candidate_gene_model = merged_exons_df %>%
      ggplot(aes(xstart = start,xend = end,y = '')) +
      geom_range(aes(fill = 'red')) +
      geom_intron(data = to_intron(merged_exons_df),aes(strand = strand),arrow.min.intron.length = 100) +
      geom_text(aes(x = end,vjust = -3, label = exonsnb),size = 6) +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      ylab(candidate$geneID) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
      annotate('text', x = xintercept,y = 0.6, label = ifelse(mutation == '','',paste0(mutation,' (Position: ',mut_pos,')')), col = 'darkblue', vjust = 0, hjust = 0.8, size = 5) +
      coord_cartesian(xlim = xlims) +
      ggtitle('Gene Model') +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))

    # plots
    if(file.exists(res_dt_candidate_gene_file) && R.utils::countLines(res_dt_candidate_gene_file)>1) {
    res_dt_candidate_gene = read.csv(res_dt_candidate_gene_file,row.names = 1)
    res_dt_candidate_gene$mean = res_dt_candidate_gene$mean/1000
    res_dt_candidate_gene$start = res_dt_candidate_gene$start/1000
    res_dt_candidate_gene$end = res_dt_candidate_gene$end/1000

    # min pvalue and max deltaPSI in region of interest.
    temp = res_dt_candidate_gene[res_dt_candidate_gene$mean >= xlims[1] & res_dt_candidate_gene$mean <= xlims[2], ]
    pval = signif(min(temp$pValue),2)
    deltaPSI = temp$deltaPsi[abs(temp$deltaPsi) == max(abs(temp$deltaPsi))]

    signif = ggplot(res_dt_candidate_gene,aes(x = mean, y = minuslogpval,color = minuslogpval)) +
      geom_point() +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      coord_cartesian(xlim = xlims) +
      ylab(bquote(-log[10]~(italic(p-value)))) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
      scale_color_continuous(palette = c('black','red')) +
      ggtitle(paste0('Aberrant splicing (min. pvalue = ',pval,', max deltaPSI = ',deltaPSI,')')) +
      theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))


    # if the region is less than 5kb lets make the plot a bit prettier
      if(xlims[2]-xlims[1] < 6) {
        res_dt_candidate_gene_subset = res_dt_candidate_gene[res_dt_candidate_gene$start > (xlims[1]) & res_dt_candidate_gene$end < (xlims[2]),]

        signif = ggplot(res_dt_candidate_gene_subset,aes(x = start, xend = end, yend = minuslogpval, y = minuslogpval, color = minuslogpval)) +
          geom_segment(linewidth = 1) +
          geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
          coord_cartesian(xlim = xlims, ylim = c(0,max(res_dt_candidate_gene$minuslogpval))) +
          ylab(bquote(-log[10]~(italic(p-value))))+
          scale_color_continuous(palette = c('black','red'),limits = c(0,max(res_dt_candidate_gene$minuslogpval))) +
          xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
          ggtitle(paste0('Aberrant splicing (min. p-value = ',pval,', max. deltaPSI = ',deltaPSI,')')) +
          theme(legend.position = 'none',plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))
      }
    } else {
      signif = ggplot() + geom_blank()
      print('No FRASER pvalues available. Likely because expression is too low')
    }

    # depth
    if(file.exists(depth_file) && R.utils::countLines(depth_file)>1) {
    depth = read.table(depth_file, header = TRUE, check.names = FALSE,comment.char= '')
    depth$POS = depth$POS / 1000
    colnames(depth) = gsub("^.*/", "",colnames(depth))
    colnames(depth) = gsub('_sorted_chrN.bam','',colnames(depth))

    for(d in 3:ncol(depth)){
       depth[,d] = depth[,d] / colmean_genes_counts[rownames(colmean_genes_counts) == colnames(depth)[d],1]
     }

    depth_filtered = depth %>% filter(POS >= min(xlims),POS < max(xlims))
    if(nrow(depth_filtered) == 0) depth_filtered = head(depth,100)

    ylims = c(0,max(c(rowMedians(as.matrix(depth_filtered[,-c(1:2)])),depth_filtered[,colnames(depth_filtered) %in% candidate$proband])))

    # pivoted
    depth_pivoted = depth_filtered %>% pivot_longer(cols = c(3:ncol(depth_filtered)), names_to = 'PatientID',values_to = 'Coverage')
    # plot
    plotCov_v2 =
      ggplot(depth_pivoted[depth_pivoted$PatientID != candidate$proband,],aes(x = POS, y = Coverage)) +
      stat_summary(geom="ribbon", fun.data=median_hilow,fun.args = list(conf.int=.5),fill=alpha('darkorange1', alpha =0.7),col= alpha('darkorange1', alpha =0.7)) +
      geom_line(data = depth_pivoted[depth_pivoted$PatientID == candidate$proband,],aes(x = POS, y = Coverage),col = 'black') +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(mutation=='',0,1)) +
      coord_cartesian(xlim = xlims) +
      ylab('Normalised coverage') +
      xlab(paste0('Chromosome ',merged_exons_df[1,1], ' (Kb)')) +
      ggtitle('Coverage (proband in black, 25-75th reference percentiles in orange)') +
      theme(plot.title = element_text(size = 24),axis.title = element_text(size = 18),axis.text = element_text(size = 14))} else {
        plotCov_v2 = ggplot() + geom_blank()
      }

  # list outputs
  output = list(list((candidate_gene_model/signif/plotCov_v2) + plot_layout(heights = c(1,1,3))))
  output
  } else {
  plot(0, main = 'no gene available');print('No gene available for gene model plot')
  }
}


# gwFRASER_table: table of the top significant splicing events (FRASER) for one sample.
# Arguments:
#   res_dt   - data.frame: genome-wide FRASER results across all samples
#   sample   - character: sampleID to filter to
#   pcutoff  - numeric: p-value cutoff (kept for interface consistency; filtering below is hardcoded to 0.05)
#   pvalue   - character: name of the p-value column to sort/filter on (default: 'padjust')
#   geneID   - character: name of the gene ID column (default: 'hgncSymbol')
gwFRASER_table = function(res_dt=gwFRASER,sample = 'HSJ_036_03_PAX',pcutoff=0.05,pvalue='padjust', geneID = 'hgncSymbol'){

  # factorize
  res_dt$chr = factor(res_dt$chr, levels = c(1:22,'X','Y','MT'))

  # keep sample of interest
  res_dt_sampleID = res_dt[res_dt$sampleID == sample,]

  # keep only the top (but remove duplicated splicing events)
  gwFRASER_top = res_dt_sampleID[order(res_dt_sampleID[[pvalue]]),]
  gwFRASER_top = gwFRASER_top[gwFRASER_top[[pvalue]]<0.05,]
  gwFRASER_top[[geneID]][is.na(gwFRASER_top[[geneID]])] = 'na'
  gwFRASER_top = gwFRASER_top[!duplicated(gwFRASER_top[[geneID]], incomparables = 'na'),]
  gwFRASER_top = gwFRASER_top[order(gwFRASER_top$chr),c(1:5,7,9:10,13:14)]

  # return
  return(gwFRASER_top)
}


# manhattan_plot: genome-wide Manhattan plot of -log10(p-value) for one sample, with top genes labelled.
# Arguments:
#   res_dt   - data.frame: genome-wide results across all samples (must contain chr, pos, sampleID, and the pvalue/geneID columns)
#   sample   - character: sampleID to plot
#   top      - numeric: number of top genes to label (default: 25)
#   pcutoff  - numeric: p-value cutoff used to select labelled/highlighted genes (default: 0.05)
#   pvalue   - character: name of the p-value column (default: 'padjust')
#   geneID   - character: name of the gene ID column (default: 'hgncSymbol')
#   shape    - logical: if TRUE, point shape encodes over/under-expression direction via l2fc (default: FALSE)
#   end      - character: name of the column giving each event's genomic end position, used for chromosome sizes (default: 'end')
manhattan_plot = function(res_dt=gwFRASER,sample = 'HSJ_036_03_PAX',top=25,pcutoff=0.05, pvalue='padjust', geneID = 'hgncSymbol',shape = FALSE, end = 'end'){

  if (is.null(res_dt) || nrow(res_dt) == 0 || !sample %in% res_dt$sampleID) {
    return(plot(0, main = 'no data available'))
  }

  # title definition
  summarise_outliers = res_dt %>%
    group_by(sampleID) %>%
    summarise(signif = sum(.data[[pvalue]] < pcutoff))

  title = paste0(sample,' ~ ',summarise_outliers$signif[summarise_outliers$sampleID==sample], ' outliers events with ',pvalue,' < ',pcutoff,'. Median number is: ', median(summarise_outliers$signif))

  # shape factor
  res_dt$shape = 'splicing'
  if(shape) res_dt$shape = ifelse(res_dt$l2fc >0,'over','under')

  # factorize
  res_dt$chr = factor(res_dt$chr, levels = c(1:22,'X','Y','MT'))

  ## calculate cumulative chromosome sizes
  chr_size = res_dt %>%
    group_by(chr) %>%
    summarise(chr_len=max(.data[[end]])) %>%
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    dplyr::select(-chr_len)

  # Add this info to the initial dataset
  results_subset_forggplot = chr_size %>%
    left_join(res_dt, ., by=c("chr"="chr")) %>%
    arrange(chr, pos) %>%
    mutate(BPcum=pos+tot)

  ### prepare axis labels.
  axisdf = results_subset_forggplot %>%
    group_by(chr) %>%
    dplyr::summarize(center=( max(BPcum) + min(BPcum)+1 ) / 2 )

  # keep only samples of interest
  results_subset_forggplot = results_subset_forggplot[results_subset_forggplot$sampleID == sample,]

  # set colors
  colors = c(brewer.pal(8,'Set2'),brewer.pal(9,'Set1'),brewer.pal(9,'Set3'))
  colors = c(colors[c(1,9,2,11,3,10,4,12,5,6,13,7,15,8,16:18,20:26)],'black')

  # keep only the top (but remove duplicated splicing events)
  gw_top = results_subset_forggplot[order(results_subset_forggplot[[pvalue]]),]
  gw_top = gw_top[gw_top[[pvalue]]<pcutoff,]
  gw_top[[geneID]][is.na(gw_top[[geneID]])] = 'na'
  gw_top = gw_top[!duplicated(gw_top[[geneID]], incomparables = 'na'),]
  gw_top = head(gw_top,top)
  gw_top = gw_top[order(gw_top$chr),]

  ### Make the ggplot:
  man_gplot = ggplot(results_subset_forggplot, aes(x = BPcum, y = -log10(.data[[pvalue]]))) +

    # Labels
    geom_label_repel(data = gw_top, aes(label = .data[[geneID]], x = BPcum, y = -log10(.data[[pvalue]])), col = 'black',  size = 4,box.padding = 2, max.overlaps = 50) +

    # Show all points
    geom_point( aes(color=chr,shape=shape), alpha=1, size=1.6) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c("over"=17, "under"=6, "splicing"=19)) +

    # Custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = c(0.01,0.01)) +
    xlab('Chromosomes') +
    ylab(paste('-log10 (',pvalue,')')) +
    ylim(0,NA) +

    # Customize theme:
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, color = "darkred", face = "bold")) +
      ggtitle(title)

  if(nrow(man_gplot@data)==0) man_gplot = plot(0, main = 'no gene available')

  return(man_gplot)
}


# candidates_summary_reactable: builds the Summary-tab reactable, one row per proband, expandable to that proband's candidate genes.
# Arguments:
#   candidates - data.frame: candidate_genes_ALL.csv contents (one row per proband-gene pair, with clinical/annotation columns)
candidates_summary_reactable = function(candidates) {
  detail_cols = c('geneID', 'Chr', 'position', 'Criteria', 'Hypothèse','HPO proband-gene matches','Mutation','FRASER','OUTRIDER','ASE')

  full = candidates %>%
    select(proband, geneID, Criteria, Age = `Âge (années)`, Sexe, Hypothèse, `HPO terms`, Mutation,`HPO proband-gene matches`,Chr = chromosome, start, stop,FRASER,OUTRIDER,ASE) %>%
    arrange(proband, geneID) %>%
    mutate(across(c(Age, Sexe, Hypothèse, `HPO terms`, Mutation), ~ ifelse(is.na(.x), '', .x))) %>%
    mutate(position = paste0(round((start + stop) / 2000000,2),' Mb'))

  first_non_na = function(x) {
    valid = x[!is.na(x) & x != '']
    if (length(valid) == 0) '' else valid[1]
  }

  summary_tbl = full %>%
    group_by(proband) %>%
    summarise(
      Genes = n(),
      Age = first_non_na(Age),
      Sexe = first_non_na(Sexe),
      `HPO terms` = first_non_na(`HPO terms`),
      .groups = 'drop'
    ) %>%
    arrange(proband)

  reactable(
    summary_tbl,
    columns = list(
      proband = colDef(name = 'Proband',width = 150),
      Genes   = colDef(name = 'Genes', align = 'left', width = 100),
      Age     = colDef(name = 'Age', align = 'left', width = 100),
      Sexe    = colDef(name = 'Sexe', align = 'left', width = 100),
      `HPO terms` = colDef(name = 'HPO terms', align = 'left')
    ),
    details = function(index) {
      proband_genes = full[full$proband == summary_tbl$proband[index], detail_cols]
      htmltools::div(
        style = "padding: 8px 12px 8px 40px",
        reactable(proband_genes, outlined = TRUE, fullWidth = TRUE, rowStyle = list(background = "#fcc95b"),
                  columns = list(
                    Chr = colDef(width = 50),
                    position   = colDef(width = 100),
                    geneID     = colDef(width = 100),
                    Criteria   = colDef(width = 150),
                    Mutation   = colDef(width = 150),
                    `HPO proband-gene matches` = colDef(width = 200),
                    FRASER = colDef(width = 150),
                    OUTRIDER = colDef(width = 150),
                    ASE = colDef(width = 150)
                  ))
      )
    },
    searchable = TRUE,
    striped = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    defaultPageSize = 100
  )
}


# gene_prioritization: ranks candidate genes for a sample by combining HPO term matches with OUTRIDER/FRASER outlier evidence.
# Arguments:
#   sample       - character: sampleID (PatientID) to prioritise genes for
#   top          - numeric: number of top-ranked genes to return (default: 100)
#   hpo_sample   - data.frame: clinical table containing 'Patient ID' and 'HPO terms' columns (default: clinical)
#   hpo_all      - character: path/filename of the HPO gene-to-phenotype annotation file, downloaded if missing (default: 'genes_to_phenotype.txt')
#   fraser       - data.frame: genome-wide FRASER results (splicing outliers)
#   outrider     - data.frame: genome-wide OUTRIDER results (expression outliers)
#   geneprior_rm - character: column name; rows with NA in this column are removed before ranking (default: 'gene score')
gene_prioritization = function(sample = 'HSJ_001_03_PAX',top=100,hpo_sample=clinical,hpo_all='genes_to_phenotype.txt',fraser="",outrider="",geneprior_rm = "gene score"){

  # hpo
  hpo_all = file.path("tmp",hpo_all)
  if(!file.exists(hpo_all)) {
        dir.create("tmp",showWarnings = FALSE)
        download.file(url='https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt',dest=hpo_all)
        } else {print(paste0('file ', hpo_all,' exists'))}
  if(!exists('hpo')) hpo = read.delim(hpo_all)

  # unique hpo terms
  hpo_un = hpo[!duplicated(hpo$hpo_name),]
  hpo_un$hpo_name_shortened = substring(hpo_un$hpo_name,1,29)
  hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)] = paste0(hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)],'.')

  # generate a named list that contains all the genes in all the HPO terms.
  temp = strsplit(hpo_sample$`HPO terms`[hpo_sample$`Patient ID` == sample],split = '||',fixed = TRUE)[[1]]
  temp = unlist(strsplit(temp,split = ': '))
  temp = temp[grepl('HP:',temp)]

  if(length(temp)>0) {
   temp = gsub(' ','',temp)
   hpo_genes = as.list(temp)
   hpo_terms = hpo_un[hpo_un$hpo_id %in% temp,]
   names(hpo_genes) = paste0(hpo_terms$hpo_id," ",hpo_terms$hpo_name_shortened)
   for(h in 1:length(hpo_genes)){hpo_genes[[h]] = unique(hpo$gene_symbol[hpo$hpo_id == hpo_terms[h,3]])}} else {hpo_genes = list(hp_1=c('none','0'));temp = c('','')}

  # add outrider
  if(!is.null(dim(outrider))) {
   outrider_temp = outrider[outrider$sampleID==sample,colnames(outrider) %in% c('geneID','pValue','zScore','exon_zScore','exon_pValue')]
   outrider_temp = outrider_temp[!is.na(outrider_temp$geneID),]
   colnames(outrider_temp) = c('geneID','OUTRIDER gene pValue','OUTRIDER gene zScore','OUTRIDER exon pValue','OUTRIDER exon zScore')
  }

  # add fraser
  if(!is.null(dim(fraser))) {
    fraser_temp = fraser[fraser$sampleID==sample,colnames(fraser) %in% c('hgncSymbol','pValue')]
    fraser_temp = fraser_temp[!is.na(fraser_temp$hgncSymbol),]
    colnames(fraser_temp) = c('geneID','FRASER gene pValue')
    fraser_temp = fraser_temp[!duplicated(fraser_temp$geneID),]
  }

  # merge Outliers
  outlier_temp = merge(fraser_temp,outrider_temp,by= 'geneID',all=TRUE, sort=FALSE)

  # generate a big table stating which gene is listed where.
  all_genes = unique(unlist(hpo_genes))
  table = sapply(hpo_genes, function(x) all_genes %in% x)
  table = data.frame('gene score' = rowSums(table),table*1, check.names = FALSE)
  table$geneID = all_genes
  table = merge(table,outlier_temp,by= 'geneID',all.x= TRUE, sort=FALSE)
  table$`OUTRIDER gene pValue`[is.na(table$`OUTRIDER gene pValue`)] = NA
  table$`OUTRIDER gene zScore`[is.na(table$`OUTRIDER gene zScore`)] = NA
  table$`OUTRIDER exon pValue`[is.na(table$`OUTRIDER exon pValue`)] = NA
  table$`OUTRIDER exon zScore`[is.na(table$`OUTRIDER exon zScore`)] = NA
  table$`FRASER gene pValue`[is.na(table$`FRASER gene pValue`)] = NA
  table = table[,c(1:2,ncol(table):3)]
  table = table[,c(1,2,5,6,7,3,4,8:ncol(table))]
  table = table[!is.na(table[,3]) | !is.na(table[,4]) | !is.na(table[,5]) | !is.na(table[,6]) | !is.na(table[,7]),]
  table = table[order(table$`gene score`,decreasing = TRUE),]
  table = table[!is.na(table[,colnames(table) == geneprior_rm]),]

  # return top hits
  return(head(table,top))
}
