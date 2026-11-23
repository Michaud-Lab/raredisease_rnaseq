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

# plot_expression_cohort: per-exon normcounts (OUTRIDER) box plots for the cohort (children / adults)
#   with a small jittered point for every cohort observation (hover shows PatientID), overlaid with the
#   proband's individual data points.
# Arguments:
#   data_child   - data.frame: cohort rows for patients < 18 yrs
#   data_adults  - data.frame: cohort rows for patients >= 18 yrs
#   data_patient - data.frame: rows for the selected proband
#   unit         - label for the expression unit (default 'normcounts'; 'TPM' when falling back to TPM data)
plot_expression_cohort = function(data_child, data_adults, data_patient, unit = 'normcounts') {
  n_exons = max(data_patient$exonID, 1)
  plot_ly() %>%
    add_trace(
      data = data_child,
      x = ~exonID,
      y = ~expression,
      type = "box",
      boxpoints = 'all',
      jitter = 0.5,
      pointpos = 0,
      name = 'Average (<18yrs)',
      color = I('darkblue'),
      hoveron = 'points',
      hoverinfo = 'text',
      text = ~paste0('<b>Sample</b>: ', PatientID, '<br><b>', unit, '</b>: ', round(expression, 1)),
      marker = list(size = 4, opacity = 0.5)
    ) %>%
    add_trace(
      data = data_adults,
      x = ~exonID,
      y = ~expression,
      type = "box",
      boxpoints = 'all',
      jitter = 0.5,
      pointpos = 0,
      name = 'Average (>18yrs)',
      color = I('darkgreen'),
      hoveron = 'points',
      hoverinfo = 'text',
      text = ~paste0('<b>Sample</b>: ', PatientID, '<br><b>', unit, '</b>: ', round(expression, 1)),
      marker = list(size = 4, opacity = 0.5)
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
      yaxis = list(title = paste0("Normalised Expression (", unit, ")")),
      title = paste0('Cohort expression (per exon, ', unit, ') for ', data_patient$geneID[1])
    )
}


# plotly_message: empty plotly figure showing a single message (used when a plot cannot be drawn).
plotly_message = function(msg) {
  plot_ly() %>% layout(
    xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
    annotations = list(list(text = msg, x = 0.5, y = 0.5, xref = 'paper', yref = 'paper', showarrow = FALSE, font = list(size = 18))))
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
#   conf.int                   - numeric in (0,1): width of the reference coverage interval passed to
#                                 median_hilow() in the coverage ribbon (default: 0.5, i.e. 25th-75th percentiles)
genemodel_plot = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             xlims = c(10,100),
                             gene_annotations='gene_annotations',
                             conf.int = 0.5) {

if(candidate$geneID != "" & file.exists(res_dt_candidate_gene_file)){

print(paste0('Defining zoom limits: ', xlims[1],' --- ',xlims[2], ' Kb'))

# colmean
colmean_genes_counts = read.table(colmean_genes_counts_file)

# filter proper Chr and transcript
gr_exons = gene_annotations[[1]]
gr_exons = gr_exons[gr_exons$gene_id == candidate$ensembl,]

# no exon annotation for this gene (e.g. a manually-added candidate not yet
# included in gene_annotations.rda, or a wrong/outdated ensembl id) -- bail
# out instead of crashing on 1:length(gr_exons@strand) below when it's 0.
if(length(gr_exons) == 0) {
  plot(0, main = 'no gene annotation available')
  print(paste0('No exon annotation found in gene_annotations.rda for ', candidate$geneID, ' (ensembl id: ', candidate$ensembl, '). Regenerate gene_annotations.rda to include this gene.'))
  return(invisible(NULL))
}

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

  mut_pos = as.numeric(strsplit(as.character(candidate$position), '___')[[1]])
  mut_pos = mut_pos[!is.na(mut_pos)]
  xintercept = mut_pos/1000
  if(length(mut_pos) == 0) {xintercept = (min(merged_exons_df$start)+max(merged_exons_df$end))/2}
  if(length(mut_pos) == 0) mut_pos = xintercept
  mutation = strsplit(as.character(candidate$mutation), '___')[[1]]
  if(length(mutation) > 0) {annotate_label = paste0(mutation,' (Position: ',mut_pos,')')}else annotate_label = ''
  
    # gene model
    candidate_gene_model = merged_exons_df %>%
      ggplot(aes(xstart = start,xend = end,y = '')) +
      geom_range(aes(fill = 'red')) +
      geom_intron(data = to_intron(merged_exons_df),aes(strand = strand),arrow.min.intron.length = 100) +
      geom_text(aes(x = end,vjust = -3, label = exonsnb),size = 6) +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(length(mutation)==0,0,1)) +
      ylab(candidate$geneID) +
      xlab(paste0('Chromosome ',merged_exons_df[1,1],' (Kb)')) +
      annotate('text', x = xintercept,y = 0.6, label = annotate_label, col = 'darkblue', vjust = 0, hjust = 0.8, size = 5) +
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
    pval = min(temp$pValue)
    deltaPSI = temp$deltaPsi[temp$pValue == pval]
    pval = signif(pval,2)
    
    signif = ggplot(res_dt_candidate_gene,aes(x = mean, y = minuslogpval,color = minuslogpval)) +
      geom_point() +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(length(mutation)==0,0,1)) +
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
          geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(length(mutation)==0,0,1)) +
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
      stat_summary(geom="ribbon", fun.data=median_hilow,fun.args = list(conf.int=conf.int),fill=alpha('darkorange1', alpha =0.7),col= alpha('darkorange1', alpha =0.7)) +
      geom_line(data = depth_pivoted[depth_pivoted$PatientID == candidate$proband,],aes(x = POS, y = Coverage),col = 'black') +
      geom_vline(xintercept = xintercept,col = 'darkblue',linewidth = 0.5,linetype = "dashed",alpha = ifelse(length(mutation)==0,0,1)) +
      coord_cartesian(xlim = xlims) +
      ylab('Normalised coverage') +
      xlab(paste0('Chromosome ',merged_exons_df[1,1], ' (Kb)')) +
      ggtitle(paste0('Coverage (proband in black, ',100*(1-conf.int)/2,'-',100*(1-(1-conf.int)/2),'th reference percentiles in orange)')) +
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


# genemodel_prepare: exon table and mutation position shared by the splicing_statistics_plotly() panels.
# Arguments:
#   candidate        - data.frame row: the candidate gene/proband row (geneID, ensembl, start, stop, position, mutation)
#   gene_annotations - list: gene annotation object; [[1]] exons GRanges
#   xlims            - numeric(2): zoom window in Kb (exon numbers are thinned out when it is the full gene + buffer)
# Returns a list(merged_exons_df, xintercept, mutation, annotate_label), or NULL if the gene has no exon annotation.
genemodel_prepare = function(candidate, gene_annotations, xlims) {

  # filter proper Chr and transcript
  gr_exons = gene_annotations[[1]]
  gr_exons = gr_exons[gr_exons$gene_id == candidate$ensembl,]

  # no exon annotation for this gene (e.g. a manually-added candidate not yet included in gene_annotations.rda, or a wrong/outdated ensembl id)
  # -- bail out instead of crashing on 1:length(gr_exons@strand) below when it's 0.
  if(length(gr_exons) == 0) {
    print(paste0('No exon annotation found in gene_annotations.rda for ', candidate$geneID, ' (ensembl id: ', candidate$ensembl, '). Regenerate gene_annotations.rda to include this gene.'))
    return(NULL)
  }

  gr_exons = gr_exons[order(gr_exons$exon_rank),]
  gr_exons$exonsnb = paste0('exon ',1:length(gr_exons@strand))

  # full-gene view includes a 5 Kb buffer on each side (the bams_subset data are extracted with the same +/- 5000 bp)
  candidate_limit = c(floor(candidate$start/1000) - 5, ceiling(candidate$stop/1000) + 5)

  merged_exons_df = as.data.frame(gr_exons)[, c("seqnames", "start", "end", "strand","exonsnb")]
  colnames(merged_exons_df) = c("chromosome", "start", "end", "strand","exonsnb")
  merged_exons_df$start = merged_exons_df$start/1000
  merged_exons_df$end = merged_exons_df$end/1000
  merged_exons_df$exon_hover = merged_exons_df$exonsnb   # full exon labels (exonsnb below may be thinned out for the on-plot text)

  if(all.equal(xlims,candidate_limit) ==TRUE) merged_exons_df$exonsnb[!(1:nrow(merged_exons_df) %in% seq(1,nrow(merged_exons_df),by = 5))] = ""

  mut_pos = as.numeric(strsplit(as.character(candidate$position), '___')[[1]])
  mut_pos = mut_pos[!is.na(mut_pos)]
  xintercept = mut_pos/1000
  if(length(mut_pos) == 0) {xintercept = (min(merged_exons_df$start)+max(merged_exons_df$end))/2}
  if(length(mut_pos) == 0) mut_pos = xintercept
  mutation = strsplit(as.character(candidate$mutation), '___')[[1]]
  if(length(mutation) > 0) {annotate_label = paste0(mutation,' (Position: ',mut_pos,')')}else annotate_label = ''

  list(merged_exons_df = merged_exons_df, xintercept = xintercept, mutation = mutation, annotate_label = annotate_label)
}


# panel_title: add a bold title at the top-left corner of a plotly panel's plotting area (paper-referenced annotation).
# plotly::subplot() maps it into the panel's domain, so panels titled this way line up identically when combined.
panel_title = function(p, text) {
  p %>% layout(annotations = list(
    text = paste0('<b>', text, '</b>'), x = 0, y = 1, xref = 'paper', yref = 'paper',
    xanchor = 'left', yanchor = 'bottom', yshift = 4, showarrow = FALSE, font = list(size = 16)))
}


# gene_model_plotly: native plotly gene model (intron line + one filled rectangle per exon, hover = exon label), titled 'Gene Model'.
# Arguments:
#   gene         - list returned by genemodel_prepare()
#   gene_name    - character: y-axis title (gene symbol)
#   xlims        - numeric(2): zoom window in Kb
gene_model_plotly = function(gene, gene_name, xlims) {
  merged_exons_df = gene$merged_exons_df
  xintercept = gene$xintercept; mutation = gene$mutation

  gene_model_plot = plot_ly() %>%
    add_lines(x = c(min(merged_exons_df$start), max(merged_exons_df$end)), y = c(0, 0),
              line = list(color = 'black', width = 1), hoverinfo = 'none', showlegend = FALSE)
  for(r in seq_len(nrow(merged_exons_df))) {
    gene_model_plot = gene_model_plot %>% add_trace(
      x = with(merged_exons_df[r, ], c(start, end, end, start, start)), y = c(-0.3, -0.3, 0.3, 0.3, -0.3),
      type = 'scatter', mode = 'lines', fill = 'toself', fillcolor = 'grey', line = list(color = 'black', width = 1),
      hoveron = 'fills', hoverinfo = I('text'), name = merged_exons_df$exon_hover[r],
      text = I(paste0('<b>', merged_exons_df$exon_hover[r], '</b><br>', round(merged_exons_df$start[r], 2), '-', round(merged_exons_df$end[r], 2), ' Kb')),
      showlegend = FALSE)
  }

  if(length(mutation) > 0) {
    gene_model_plot = gene_model_plot %>%
      add_lines(x = rep(xintercept[1], 2), y = c(-0.8, 0.8), line = list(color = 'darkblue', dash = 'dash', width = 2),
                hoverinfo = 'none', showlegend = FALSE) %>%
      add_text(x = xintercept[1], y = 0.6, text = gene$annotate_label, textposition = 'bottom left',
               textfont = list(size = 14, color = 'darkblue'), hoverinfo = 'none', showlegend = FALSE)
  }
  gene_model_plot %>% layout(
    xaxis = list(title = paste0('Chromosome ',merged_exons_df[1,1],' (Kb)'), range = xlims),
    yaxis = list(title = gene_name, range = c(-0.8, 0.8), showticklabels = FALSE, zeroline = FALSE, showgrid = FALSE),
    hoverlabel = list(bgcolor = "white")) %>%
    panel_title('Gene Model')
}


# fraser_plotly: native plotly FRASER panel (-log10 p-value of each splice junction, black-to-red colour scale).
# Junctions with -log10(p) > 0 are drawn as segments (start to end), zero-valued ones as black markers.
# Arguments:
#   res_dt_candidate_gene_file - character: path to the FRASER per-candidate-gene results CSV
#   xlims                      - numeric(2): zoom window in Kb
#   gene                       - list returned by genemodel_prepare() (chromosome label, mutation line)
# Returns the plotly object, titled with the min. p-value and max. deltaPSI in the window.
fraser_plotly = function(res_dt_candidate_gene_file, xlims, gene) {
  xintercept = gene$xintercept; mutation = gene$mutation

  if(!(file.exists(res_dt_candidate_gene_file) && R.utils::countLines(res_dt_candidate_gene_file)>1)) {
    print('No FRASER pvalues available. Likely because expression is too low')
    return(plotly_message('No FRASER pvalues available (likely because expression is too low)') %>% panel_title('Aberrant splicing'))
  }

  res_dt_candidate_gene = read.csv(res_dt_candidate_gene_file,row.names = 1)
  res_dt_candidate_gene$mean = res_dt_candidate_gene$mean/1000
  res_dt_candidate_gene$start = res_dt_candidate_gene$start/1000
  res_dt_candidate_gene$end = res_dt_candidate_gene$end/1000

  # min pvalue and max deltaPSI in region of interest.
  temp = res_dt_candidate_gene[res_dt_candidate_gene$mean >= xlims[1] & res_dt_candidate_gene$mean <= xlims[2], ]
  pval = min(temp$pValue)
  deltaPSI = temp$deltaPsi[temp$pValue == pval]
  pval = signif(pval,2)

  title = paste0('Aberrant splicing (min. p-value = ',pval,', max. deltaPSI = ',deltaPSI,')')
  ymax = max(res_dt_candidate_gene$minuslogpval, na.rm = TRUE)
  col_pal = colorRampPalette(c('black','red'))(100)
  col_idx = function(x) pmin(100, pmax(1, ceiling(100 * x / ymax)))   # clamp: x/ymax can round just above 1 (index 101 = NA colour)
  hover = function(d) paste0('<b>Junction</b>: ', round(d$start, 2), '-', round(d$end, 2), ' Kb',
                             '<br><b>p-value</b>: ', signif(d$pValue, 3),
                             '<br><b>deltaPSI</b>: ', d$deltaPsi,
                             '<br><b>Reads spanning/flanking junction (Proband)</b>: ', d$counts,' / ',d$totalCounts,
                             '<br><b>Reads spanning/flanking junction (Cohort mean)</b>: ', round(d$meanCounts, 1),' / ',round(d$meanTotalCounts,1))
  res_plot = res_dt_candidate_gene[!is.na(res_dt_candidate_gene$minuslogpval), ]   # NA p-values have no colour (plotly would draw them in its default blue)
  fraser_plot = plot_ly()

  # junctions with -log10(p) > 0: a segment from the junction start to its end (any window size)
  d_pos = res_plot[res_plot$minuslogpval > 0, ]
  for(r in seq_len(nrow(d_pos))) {
    fraser_plot = fraser_plot %>% add_segments(
      x = d_pos$start[r], xend = d_pos$end[r], y = d_pos$minuslogpval[r], yend = d_pos$minuslogpval[r],
      line = list(color = col_pal[col_idx(d_pos$minuslogpval[r])], width = 3),
      hoverinfo = 'text', text = hover(d_pos[r, ]), showlegend = FALSE)
  }
  # junctions with -log10(p) = 0: markers at the junction midpoint
  d_zero = res_plot[res_plot$minuslogpval <= 0, ]
  if(nrow(d_zero) > 0) {
    fraser_plot = fraser_plot %>% add_markers(
      data = d_zero, x = ~mean, y = ~minuslogpval,
      marker = list(color = col_pal[1], size = 7),
      hoverinfo = 'text', text = hover(d_zero), showlegend = FALSE)
  }
  if(length(mutation) > 0) {
    fraser_plot = fraser_plot %>% add_lines(x = rep(xintercept[1], 2), y = c(-0.2, ymax),
      line = list(color = 'darkblue', dash = 'dash', width = 2), hoverinfo = 'none', showlegend = FALSE)
  }
  fraser_plot = fraser_plot %>% layout(
    xaxis = list(title = paste0('Chromosome ',gene$merged_exons_df[1,1],' (Kb)'), range = xlims),
    yaxis = list(title = '-log10(p-value)', range = c(-0.2, ymax * 1.1), zeroline = FALSE),
    hoverlabel = list(bgcolor = "white")) %>%   # no default line at y = 0
    panel_title(title)
}


# coverage_plotly: native plotly coverage panel. The reference band (percentiles of the cohort, proband excluded) is computed
# directly on the wide per-sample matrix, the proband is drawn as a black line.
# Arguments:
#   depth_file                - character: path to the per-base coverage depth file
#   colmean_genes_counts_file - character: path to the column-mean gene counts file (used to normalise coverage)
#   proband                   - character: sampleID of the proband (column of the depth file)
#   xlims                     - numeric(2): zoom window in Kb
#   gene                      - list returned by genemodel_prepare() (mutation line)
#   conf.int                  - numeric in (0,1): width of the reference interval (0.5 = 25th-75th percentiles)
# Returns the plotly object, titled with the reference percentiles shown.
coverage_plotly = function(depth_file, colmean_genes_counts_file, proband, xlims, gene, conf.int = 0.5) {
  xintercept = gene$xintercept; mutation = gene$mutation

  if(!(file.exists(depth_file) && R.utils::countLines(depth_file)>1)) {
    return(plotly_message('No coverage data available') %>% panel_title('Coverage'))
  }

  colmean_genes_counts = read.table(colmean_genes_counts_file)
  depth = read.table(depth_file, header = TRUE, check.names = FALSE,comment.char= '')
  depth$POS = depth$POS / 1000
  colnames(depth) = gsub("^.*/", "",colnames(depth))
  colnames(depth) = gsub('_sorted_chrN.bam','',colnames(depth))

  for(d in 3:ncol(depth)){
    depth[,d] = depth[,d] / colmean_genes_counts[rownames(colmean_genes_counts) == colnames(depth)[d],1]
  }

  depth_filtered = depth %>% filter(POS >= min(xlims),POS < max(xlims))
  if(nrow(depth_filtered) == 0) depth_filtered = head(depth,100)

  # reference band computed directly on the wide matrix (no pivot_longer / ggplot / ggplotly: ~3 million rows for a full gene)
  is_proband = colnames(depth_filtered) == proband
  ref_mat = as.matrix(depth_filtered[, -c(1:2)][, !is_proband[-c(1:2)], drop = FALSE])
  ref_q = matrixStats::rowQuantiles(ref_mat, probs = c(0.5 - conf.int/2, 0.5 + conf.int/2), na.rm = TRUE)
  proband_cov = depth_filtered[, is_proband][1:nrow(depth_filtered)]
  orange = 'rgba(255,140,0,0.7)'   # darkorange1, alpha 0.7
  coverage_plot = plot_ly() %>%
    add_trace(x = depth_filtered$POS, y = ref_q[, 2], type = 'scatter', mode = 'lines', line = list(color = orange, width = 1),
              hoverinfo = 'none', showlegend = FALSE) %>%
    add_trace(x = depth_filtered$POS, y = ref_q[, 1], type = 'scatter', mode = 'lines', line = list(color = orange, width = 1),
              fill = 'tonexty', fillcolor = orange, hoverinfo = 'none', showlegend = FALSE) %>%
    add_trace(x = depth_filtered$POS, y = proband_cov, type = 'scatter', mode = 'lines', line = list(color = 'black', width = 1),
              hoverinfo = 'text', showlegend = FALSE,   # hover on the proband line also reports the reference band at that position
              text = paste0('<b>Position</b>: ', round(depth_filtered$POS, 3), ' Kb',
                            '<br><b>Proband</b>: ', round(proband_cov, 2),
                            '<br><b>Reference ', 100*(0.5 - conf.int/2), '-', 100*(0.5 + conf.int/2), 'th percentiles</b>: ',
                            round(ref_q[, 1], 2), ' - ', round(ref_q[, 2], 2)))
  if(length(mutation) > 0) {
    coverage_plot = coverage_plot %>% add_lines(x = rep(xintercept[1], 2), y = c(0, max(ref_q[, 2], proband_cov, na.rm = TRUE)),
      line = list(color = 'darkblue', dash = 'dash', width = 1), hoverinfo = 'none', showlegend = FALSE)
  }
  title = paste0('Coverage (proband in black, ',100*(1-conf.int)/2,'-',100*(1-(1-conf.int)/2),'th reference percentiles in orange)')
  coverage_plot %>% layout(xaxis = list(title = paste0(sub('chr','Chromosome ',depth[1,1]), ' (Kb)'), range = xlims), yaxis = list(title = 'Normalised coverage')) %>%
    panel_title(title)
}


# splicing_statistics_plotly: plotly version of genemodel_plot() (gene model, FRASER significance and coverage for a
# candidate gene). Combines the panels built by gene_model_plotly(), fraser_plotly(), coverage_plotly() and
# junctions_plotly() into ONE plotly figure (four stacked panels, linked x axes in Kb).
# Arguments:
#   candidate                  - data.frame row: the candidate gene/proband row (geneID, ensembl, start, stop, position, mutation, proband)
#   res_dt_candidate_gene_file - character: path to the FRASER per-candidate-gene results CSV
#   depth_file                 - character: path to the per-base coverage depth file
#   bam_file                   - character: path to the proband's BAM file (kept for interface consistency)
#   colmean_genes_counts_file  - character: path to the column-mean gene counts file (used to normalise coverage)
#   xlims                      - numeric vector length 2: zoom window in Kb (default: c(10, 100))
#   gene_annotations           - list: gene annotation object; [[1]] exons GRanges
#   conf.int                   - numeric in (0,1): width of the reference coverage interval (default: 0.5, i.e. 25th-75th percentiles)
splicing_statistics_plotly = function(candidate = candidates,
                             res_dt_candidate_gene_file = "resdet",
                             depth_file = "depth",
                             bam_file = "bam",
                             colmean_genes_counts_file = 'colmean_genes_counts.tsv',
                             xlims = c(10,100),
                             gene_annotations='gene_annotations',
                             conf.int = 0.5) {

  if(!(candidate$geneID != "" & file.exists(res_dt_candidate_gene_file))) {
    print('No gene available for gene model plot')
    return(plotly_message('No gene available for gene model plot'))
  }

  gene = genemodel_prepare(candidate, gene_annotations, xlims)
  if(is.null(gene)) return(plotly_message('no gene annotation available'))

  # the four panels
  gene_model_plot = gene_model_plotly(gene, candidate$geneID, xlims)
  fraser_plot = fraser_plotly(res_dt_candidate_gene_file, xlims, gene)
  coverage_plot = coverage_plotly(depth_file, colmean_genes_counts_file, candidate$proband, xlims, gene, conf.int)
  junctions_plot = junctions_plotly(res_dt_candidate_gene_file, xlims = xlims, xintercept = if(length(gene$mutation) > 0) gene$xintercept[1] else NULL)

  # Combine the four panels into one interactive figure (own x axis labels on each panel, ranges linked).
  heights = c(2, 3, 5, 3); margin = 0.04
  hs = heights / sum(heights) * (1 - 3 * margin)

  fig = plotly::subplot(gene_model_plot, fraser_plot, coverage_plot, junctions_plot,
                        nrows = 4, shareX = FALSE, titleX = TRUE, titleY = TRUE, margin = margin, heights = hs) %>%
    layout(
      showlegend = TRUE,   # only the junction panel's two traces ('Proband', 'Cohort mean') have a legend entry
      xaxis = list(range = xlims), xaxis2 = list(range = xlims, matches = 'x'), xaxis3 = list(range = xlims, matches = 'x'), xaxis4 = list(range = xlims, matches = 'x')   # zoom/pan still linked
    ) %>% plotly_build()

  # legend of the junction panel: just above that panel, right-aligned
  fig$x$layout$legend = list(orientation = 'h', x = 1, xanchor = 'right', y = fig$x$layout$yaxis4$domain[2] + 0.005, yanchor = 'bottom')

  fig
}



# junctions_plotly: sashimi-style plot of the splice junctions quantified by FRASER (native plotly).
# Each junction is drawn as an arc from its start to its end; the arc height and line width scale with the number of
# split reads supporting it: the proband's `counts` (blue) and the cohort `meanCounts` (orange), both from the FRASER
# results saved by fraser_pipeline() in gene_<gene>_<proband>_res_dt_candidate_gene.csv (all junctions of the gene,
# with p-value and deltaPSI).
# Arguments:
#   fraser_file - character: FRASER per-candidate-gene results CSV
#   xlims       - numeric(2): window to show, in Kb
#   min_count   - numeric: only draw arcs supported by at least this many reads (default 1)
#   xintercept  - numeric (Kb): optional position of a variant, drawn as a dashed vertical line
# Returns a titled plotly object (or a message figure if there are no FRASER results or no junction in the window).
junctions_plotly = function(fraser_file,
                            xlims = c(0, 100),
                            min_count = 1,
                            xintercept = NULL) {

  if(!(file.exists(fraser_file) && R.utils::countLines(fraser_file) > 1)) {
    return(plotly_message('No FRASER junctions available') %>% panel_title('Splice junctions'))
  }
  fr = read.csv(fraser_file, row.names = 1)
  jt = data.frame(start = fr$start / 1000, end = fr$end / 1000, count = fr$counts, meanCount = fr$meanCounts,
                  pValue = fr$pValue, deltaPsi = fr$deltaPsi)   # Kb, like the other panels
  jt = jt[jt$end > xlims[1] & jt$start < xlims[2], ]
  d_prob = jt[!is.na(jt$count) & jt$count >= min_count, ]
  d_mean = jt[!is.na(jt$meanCount) & jt$meanCount >= min_count, ]
  if(nrow(d_prob) + nrow(d_mean) == 0) return(plotly_message('No splice junctions in this window') %>% panel_title('Splice junctions'))

  # arc = half ellipse, peak height = read count; line width scaled to count (1 - 6 px), same scale for both series
  t = seq(0, 1, length.out = 40)
  all_counts = c(d_prob$count, d_mean$meanCount)
  lwd = function(x) 1 + 5 * (log(x) - log(min(all_counts))) / max(log(max(all_counts)) - log(min(all_counts)), 1e-9)
  col_prob = 'steelblue'; col_mean = 'darkorange'
  jhover = function(d) paste0('<b>Junction</b>: ', round(d$start, 3), '-', round(d$end, 3), ' Kb')
  p = plot_ly()

  # cohort mean first (drawn behind the proband)
  for(r in seq_len(nrow(d_mean))) {
    p = p %>% add_trace(
      x = d_mean$start[r] + t * (d_mean$end[r] - d_mean$start[r]), y = d_mean$meanCount[r] * sin(pi * t),
      type = 'scatter', mode = 'lines', line = list(color = col_mean, width = lwd(d_mean$meanCount[r])),
      name = 'Cohort mean', legendgroup = 'mean', showlegend = (r == 1), hoverinfo = 'text',
      text = paste0(jhover(d_mean[r, ]), '<br><b>Cohort mean reads</b>: ', round(d_mean$meanCount[r], 1)))
  }
  for(r in seq_len(nrow(d_prob))) {
    p = p %>% add_trace(
      x = d_prob$start[r] + t * (d_prob$end[r] - d_prob$start[r]), y = d_prob$count[r] * sin(pi * t),
      type = 'scatter', mode = 'lines', line = list(color = col_prob, width = lwd(d_prob$count[r])),
      name = 'Proband', legendgroup = 'proband', showlegend = (r == 1), hoverinfo = 'text',
      text = paste0(jhover(d_prob[r, ]), '<br><b>Proband reads</b>: ', d_prob$count[r],
                    '<br><b>Cohort mean reads</b>: ', round(d_prob$meanCount[r], 1),
                    '<br><b>p-value</b>: ', signif(d_prob$pValue[r], 3),
                    '<br><b>deltaPSI</b>: ', d_prob$deltaPsi[r]))
  }
  if(length(xintercept) > 0) {
    p = p %>% add_lines(x = rep(xintercept[1], 2), y = c(0, max(all_counts)),
      line = list(color = 'darkblue', dash = 'dash', width = 1), hoverinfo = 'none', showlegend = FALSE)
  }
  p %>% layout(
    xaxis = list(title = paste0(sub('chr', 'Chromosome ', fr$seqnames[1]), ' (Kb)'), range = xlims, tickformat = ',.3~f', exponentformat = 'none'),   # same label as the other panels
    yaxis = list(title = 'Supporting split reads (FRASER)', rangemode = 'tozero', zeroline = FALSE),
    legend = list(orientation = 'h', x = 0, y = 1.08),
    hoverlabel = list(bgcolor = 'white')) %>%
    panel_title('Splice junctions quantified by FRASER (arc height = supporting split reads)')
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

  if (is.null(res_dt) || nrow(res_dt) == 0 || length(sample) != 1 || is.na(sample) || !sample %in% res_dt$sampleID) {
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


# omim_phenotypes_lookup: named character vector (names = Approved Gene Symbol) of OMIM Phenotypes, read once from genemap2 and cached.
# Returns NULL (with a warning) if the genemap2 file is not found.
# Arguments:
#   omim_genemap2_name - file name of the OMIM genemap2 export (cached per file name)
# Looks in <datadir>/input (as set in params by the app) first, then in <wd>/data/input (e.g. when datadir is data_minimal).
omim_phenotypes_lookup = local({
  cache = list()
  function(omim_genemap2_name = 'genemap2.20260911.txt') {
    if (is.null(cache[[omim_genemap2_name]])) {
      candidates = file.path(c(if (exists('params')) params$datadir, file.path(getwd(), 'data')), 'input', omim_genemap2_name)
      found = candidates[file.exists(candidates)]
      if (length(found) == 0) {
        warning('OMIM genemap2 file not found; looked in: ', paste(candidates, collapse = ', '))
        return(NULL)
      }
      
      gm = read.delim(found[1], comment.char = '#', header = FALSE, sep = '\t', quote = '', colClasses = 'character', na.strings = '', fill = TRUE)
      # genemap2 columns: 9 = Approved Gene Symbol, 13 = Phenotypes
      gm = gm[!is.na(gm[[9]]) & !is.na(gm[[13]]) & !duplicated(gm[[9]]), ]
      cache[[omim_genemap2_name]] <<- setNames(gm[[13]], gm[[9]])
    }
    cache[[omim_genemap2_name]]
  }
})

# omim_link_geneID: replaces each value in a data.frame's geneID column with an HTML link to its OMIM search entry,
# and adds a 'Phenotypes' column right after geneID (OMIM genemap2 Phenotypes matched on Approved Gene Symbol; NA if none).
# For use with DT::datatable(..., escape = -which(colnames(data) == 'geneID')).
# Arguments:
#   data - data.frame containing a geneID column
omim_link_geneID = function(data) {
  lookup = omim_phenotypes_lookup()
  data$Phenotypes = if (is.null(lookup)) 'no reference data' else unname(lookup[as.character(data$geneID)])
  data$geneID = ifelse(
    is.na(data$geneID) | data$geneID == '',
    data$geneID,
    paste0("<a href='https://www.omim.org/search/?search=", data$geneID, "&type=entry' target='_blank'>", data$geneID, "</a>")
  )
  gene_pos = which(colnames(data) == 'geneID')
  data[, append(setdiff(colnames(data), 'Phenotypes'), 'Phenotypes', after = gene_pos), drop = FALSE]
}

# candidates_summary_reactable: builds the Summary-tab reactable, one row per proband, expandable to that proband's candidate genes.
# Arguments:
#   candidates - data.frame: candidate_genes_ALL.csv contents (one row per proband-gene pair, with clinical/annotation columns)
candidates_summary_reactable = function(candidates) {
  candidates$Mutation[is.na(candidates$Mutation)] = '~'
  candidates$Hypothèse[is.na(candidates$Hypothèse)] = '~'
 
  detail_cols = c('geneID', 'Chr', 'position', 'Criteria', 'Mutation','Hypothèse','HPO proband-gene matches','FRASER','OUTRIDER','ASE')

  full = candidates %>%
    select(Notes, proband, geneID, Criteria, Age = `Âge (années)`, Sexe, Hypothèse, `HPO terms`, Mutation,`HPO proband-gene matches`,Chr = chromosome, start, stop,FRASER,OUTRIDER,ASE) %>%
    arrange(proband, geneID) %>%
    mutate(across(c(Notes, Age, Sexe, Hypothèse, `HPO terms`, Mutation, Criteria, `HPO proband-gene matches`, FRASER, OUTRIDER, ASE), ~ ifelse(is.na(.x), '', as.character(.x)))) %>%
    mutate(position = paste0(round((start + stop) / 2000000,2),' Mb'))

  first_non_na = function(x) {
    valid = x[!is.na(x) & x != '']
    if (length(valid) == 0) '' else valid[1]
  }

  # Concatenate every candidate gene's detail fields for a proband into one blob so the
  # search box can match on them, even though they only appear in the expandable 'details' row.
  summary_tbl = full %>%
    group_by(proband) %>%
    summarise(
      Notes = first_non_na(Notes),
      Genes = n(),
      Age = first_non_na(Age),
      Sexe = first_non_na(Sexe),
      `HPO terms` = first_non_na(`HPO terms`),
      search_detail = paste(geneID, Criteria, Hypothèse, Mutation, `HPO proband-gene matches`, FRASER, OUTRIDER, ASE, collapse = ' | '),
      .groups = 'drop'
    ) %>%
    relocate(Notes) %>%
    arrange(proband)

  reactable(
    summary_tbl,
    columns = list(
      Notes   = colDef(name = 'Notes', align = 'left', width = 200),
      proband = colDef(name = 'Proband',width = 150),
      Genes   = colDef(name = 'Genes', align = 'left', width = 100),
      Age     = colDef(name = 'Age', align = 'left', width = 100),
      Sexe    = colDef(name = 'Sexe', align = 'left', width = 100),
      `HPO terms` = colDef(name = 'HPO terms', align = 'left'),
      # Kept searchable but visually hidden (show = FALSE would drop it from the search index too)
      search_detail = colDef(name = '', width = 0, style = list(display = 'none'), headerStyle = list(display = 'none'))
    ),
    details = function(index) {
      proband_genes = full[full$proband == summary_tbl$proband[index], detail_cols]
      highlight_if_present = function(value) {
        if (grepl('HP:|padj|pval|>',value)) list(background = "#fa4d4d")
      }
      htmltools::div(
        style = "padding: 8px 12px 8px 40px",
        reactable(proband_genes, outlined = TRUE, fullWidth = TRUE, rowStyle = list(background = "#fcc95b"),
                  columns = list(
                    Chr = colDef(width = 50),
                    position   = colDef(width = 100),
                    geneID     = colDef(width = 100, cell = function(value) {
                      if (value == '') return(value)
                      htmltools::tags$a(
                        href = paste0("https://www.omim.org/search/?search=", value, "&type=entry"),
                        target = "_blank", value
                      )
                    }),
                    Criteria   = colDef(width = 150),
                    `HPO proband-gene matches` = colDef(width = 200, style = highlight_if_present),
                    Mutation   = colDef(width = 150, style = highlight_if_present),
                    FRASER = colDef(width = 150, style = highlight_if_present),
                    OUTRIDER = colDef(width = 150, style = highlight_if_present),
                    ASE = colDef(width = 150, style = highlight_if_present)
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
gene_prioritization = function(sample = 'HSJ_001_03_PAX',pcutoff=0.5,hpo_sample=clinical,hpo_all='genes_to_phenotype.txt',fraser="",outrider="",geneprior_rm = "gene score"){

  # hpo
  #hpo_all = file.path("tmp",hpo_all)
  if(!file.exists(hpo_all)) {
        print(paste0('file ', hpo_all,' downloaded'))
        dir.create("tmp",showWarnings = FALSE)
        download.file(url='https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt',dest=hpo_all)
        }
  if(!exists('hpo')) hpo = read.delim(hpo_all)

  # unique hpo terms
  hpo_un = hpo[!duplicated(hpo$hpo_name),]
  #hpo_un = hpo %>% distinct(hpo_id,gene_symbol, .keep_all = TRUE)
  hpo_un$hpo_name_shortened = substring(hpo_un$hpo_name,1,29)
  hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)] = paste0(hpo_un$hpo_name_shortened[nchar(hpo_un$hpo_name_shortened)!=nchar(hpo_un$hpo_name)],'.')

  # generate a named list that contains all the genes in all the HPO terms.
  if(sample %in% hpo_sample$`Patient ID`) {
    temp = strsplit(hpo_sample$`HPO terms`[hpo_sample$`Patient ID` == sample],split = '||',fixed = TRUE)[[1]]
    temp = unlist(strsplit(temp,split = ': '))
    temp = temp[grepl('HP:',temp)]} else {temp = NULL}
  
  if(length(temp)>0) {
   temp = gsub(' ','',temp)
   hpo_genes = as.list(temp)
   hpo_terms = hpo_un[hpo_un$hpo_id %in% temp,]
   names(hpo_genes) = paste0(hpo_terms$hpo_id," ",hpo_terms$hpo_name_shortened)
   for(h in 1:length(hpo_genes)){hpo_genes[[h]] = unique(hpo$gene_symbol[hpo$hpo_id == hpo_terms[h,3]])}} else {return(data.frame(empty = as.character()))}
   hpo_genes = hpo_genes[!is.na(names(hpo_genes))]

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
  
  if(nrow(table)==0) {return(data.frame(empty = as.character()))}
  
  table$pvalue = 1
  table = head(table,50)
  
  #calculate a pvalue for each (top10) gene:
  for(g in 1:nrow(table)){
    total_nb_of_genes = length(unique(hpo$gene_symbol))
    nb_genes = sapply(hpo_genes,length)
    nb_genes = nb_genes[order(names(nb_genes))]
    
    nb_hpo = table[g,colnames(table) %in% names(nb_genes)]
    if(length(nb_hpo)==1) names(nb_hpo) = names(nb_genes) # this is an edge case, with a single HPO term...
    nb_hpo = nb_hpo[order(names(nb_hpo))]
    
    nb_genes_hpo = ifelse(nb_hpo == 1,nb_genes,5000-nb_genes)
    table$pvalue[g] = signif(prod(nb_genes_hpo) / total_nb_of_genes **length(nb_genes_hpo) * (sum(nb_genes_hpo)),2)
  }
  
  table = table[order(table$pvalue),]
  table = cbind(table[,1:2],table[,ncol(table)],table[,3:(ncol(table)-1)])
  colnames(table)[3] = 'gene score adj. pvalue'
  
  # return top hits
  return(table[table[,3] < pcutoff,])
}
