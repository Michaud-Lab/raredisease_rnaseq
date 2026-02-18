
# params
params = list(datadir = file.path(getwd(),'data/')) 

# R libraries that you need
packages = c('DT','plotly','tidyr','shiny','shinyjs','jsonlite','igvShiny','GenomicAlignments','dplyr','ggtranscript','patchwork','Hmisc','bslib','RColorBrewer','ggrepel')

for(p in 1:length(packages)) {
  if(packages[p] %in% installed.packages()) {
    library(packages[p],character.only = TRUE)
    } else if (packages[p] == 'ggtranscript'){
      stop('Install ggtranscript with : remotes::install_github("dzhang32/ggtranscript")')
    } else {
      stop(paste0('Error in library() : there is no package called ', packages[p]))
    }
}

# load files
source(file.path(params$datadir,'scripts/rnaseq_shinyhelper_functions.R'))
fc_exons_raw = read.table(file.path(params$datadir,'fc_exons_raw.tsv'),sep = '\t',check.names = F);fc_exons_raw[,-c(1:5)] = round(fc_exons_raw[,-c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir,'fc_exons_tpm.tsv'),sep = '\t',check.names = F)
fc_genes_tpm = read.table(file.path(params$datadir,'fc_genes_tpm.tsv'),sep = '\t',check.names = F)

gwOUTRIDER = read.csv(file.path(params$datadir,'results_OUTRIDER.tsv'),sep = '\t',check.names = F)
gwOUTRIDER$chr = factor(gwOUTRIDER$chr,levels = c(1:22,'X','Y','MT'))

significant_perexons_OUTRIDER = read.csv(file.path(params$datadir,'exons_OUTRIDER.tsv'),sep = '\t',check.names = F)
significant_perexons_OUTRIDER$chr = factor(significant_perexons_OUTRIDER$chr,levels = c(1:22,'X','Y','MT'))

candidates_OUTRIDER = read.csv(file.path(params$datadir,'candidates_OUTRIDER.tsv'),sep = '\t',check.names = F)
candidates_perexons_OUTRIDER = read.csv(file.path(params$datadir,'candidates_perexons_OUTRIDER.tsv'),sep = '\t',check.names = F)

transcripts_named_filtered = read.csv(file.path(params$datadir,'transcripts_named_filtered.tsv'),sep = '\t',check.names = F)
transcripts_named_filtered_ggplot = read.csv(file.path(params$datadir,'transcripts_named_filtered_ggplot.tsv'),sep = '\t',check.names = F)

fc_exons_tpm_ggplot = read.csv(file.path(params$datadir,'fc_exons_tpm_ggplot.tsv'),sep = '\t',check.names = F)
candidates = read.csv(file.path(params$datadir,'candidate_genes_3.txt'))
clinical = read.csv(file.path(params$datadir,'clinical.tsv'),sep = '\t',check.names = F)
html_file = file.path(params$datadir,'multiqc_report.html')

gwFRASER = read.csv(file.path(params$datadir,'gwFRASER.csv'),row.names = 1)

report_version = read_json(file.path(params$datadir,'/VERSION.json'))
report_version$data = params$zipfile

load(file = file.path(params$datadir,"gene_annotations.rda"))

theme <- bs_theme(
  version = 5,
  bootswatch = "cosmo",
  primary = "#0d6efd",
  base_font = font_google("Roboto")
)

#####
#####UI
#####    
ui <- page_fluid(
    theme = theme,
   
    # Dark title header
    div(
      class = "bg-dark text-white p-3 mb-4",
      h2("RNAseq Dashboard")
    ),

    tabsetPanel(

      ###  Selection & Info
      tabPanel('Proband Selection',
               card(
                 selectInput(
                    inputId = "proband2",
                    label = h2("Select a patient:"),
                    choices = candidates$proband2,
                    selected = 'HSJ_001_03'),
                 div(
                   imageOutput("mainimage"),
                   class = "text-center" 
                 )),
               card(
                 h2("Information"),
                 htmlOutput("description")),
               card(
                 h2("Software Version"),
                 DTOutput("Version")
                 )
               ),  
      
      ### Data table
      tabPanel("Expression",
               card(
                h3('Gene expression per exon'),
                'Reference transcript based on Matched Annotation from NCBI and EMBL-EBI (MANE)',
                selectInput("table_choice", "Choose an expression metric: ",
                  choices = c("Normalised expression", "Raw counts","Isoform-specific expression"),
                  selected = "Normalised expression"),
                DTOutput("exonTPM")
                )
               ),
    
      ### Plot of Expression
      tabPanel("Plot",
               card(plotlyOutput("Expression",width = "1200px")),
               card(plotlyOutput("Expression_perfamily",width = "1200px"))
              ),
    
    ### Data table
      tabPanel("OUTRIDER",
               card(
                 h3('OUTRIDER - OUTlier in RNA-Seq fInDER'),
                 'Identification of aberrant gene expression in RNA-seq data, Outliers are identified as read counts that significantly deviate from the population',
                 '(all probands, including some adults from the LC & F0 cohorts)'),
               card(
                 h4('Genome-wide significance'),
                 plotOutput("gwOUTRIDER", width = '1200px', height = "400px")),
               card(
                  fill=F, 
                  height='200px',
                  h4('Candidate gene'),
                  DTOutput("candidates_OUTRIDER")),
               card(
                  h4('Candidate exons'),
                  DTOutput("candidates_OUTRIDER_exons")),
               card(
                  h4('All significant genes'), 
                  selectInput("pvalue", "Choose a p-value threshold:",
                              choices = c(0.05,0.01,0.005,0.001,0.0005),
                              selected = 0.005),
                  DTOutput("table_OUTRIDER")),
               card(
                  h4('All significant exons (p<0.01)'),
                  DTOutput("significant_OUTRIDER_exons")
                  )
               ),
    
    ### Plot of Expression
    tabPanel("Gene model",
             card(
              h3('Gene model and significance values for FRASER'),
              column(width = 12,align = "center",uiOutput("genemodel_slider"))),
             card(
               plotOutput("Figure_genemodel", width = '1200px', height = "800px")),
             card(
               htmlOutput('Figure_genemodel_legend'))
             ),
    
  
    ### Plot of Structural variation
    tabPanel("IGV",
             card(
              h3('Integrative Genome Viewer'), 
              useShinyjs(), # Initialize shinyjs
              actionButton("addBamLocalFileButton",
                            "Show gene alignment",
                            style = "background-color: red; color: white; border-color: darkorchid;")),
            card(
              igvShinyOutput('igvShiny',width = "1200px")
             )
            ),
  
    ### Plot of Structural variation
    tabPanel(
             "FRASER",
             card(
               h3('FRASER: Find RAre Splicing Events in RNAseq Data'),
               'Splice site map on top, with new exons and junctions in proband in yellow or red line.',
               'Splice site map generated in an annotation-free fashion based on RNA-seq coverage of split reads.',
               'As such, introns/exons may differ from actual genome annotation.',
               'In red, proband of interest. In blue, five representative samples of the population.',
               imageOutput("Sashimi",width = "1200px")
               )
             ),
     
    ##genome-wide FRASER
    tabPanel(
            "gwFRASER",
            card(
              h4('Genome-wide significance values for FRASER'),
              plotOutput("gwFRASER", width = '1200px', height = "400px")),
            card(
              h4('Significant splicing events (below adj. p-value < 0.05)'),
              DTOutput("gwFRASER_table")
              )
            ),

    ### Plot of Structural variation
    tabPanel("fasta",
             card(
               h3('Transcribed sequences: reference / alternate'),
               em('dashes (---) correspond to intron retention events, variants in ',strong('bold')),
               htmlOutput("fasta")
               )
             ),
    
    ### multiQC
    tabPanel(
      "MultiQC",
      h3('Quality Control metrics'),
      htmlOutput("htmlViewer")
    )
)
)



#####
#####server
#####
server <- function(input, output, session) {
  
    ### Selected variable i
    output$selected_var <- renderText({input$proband2})
    i <- reactive({c(1:nrow(candidates))[candidates$proband2 == input$proband2]})
  
    ### Reactive list
    reactive_inputs <- reactive({
        reactive_i  = strsplit(input$proband2,'@')[[1]][1]
        column = c(1:ncol(transcripts_named_filtered))[colnames(transcripts_named_filtered) == reactive_i]
        column_not = c(1:ncol(transcripts_named_filtered))[colnames(transcripts_named_filtered) != reactive_i]
        column_not = column_not[-c(1,2,length(column_not))] 
        transcripts_reactive = transcripts_named_filtered[transcripts_named_filtered$proband == reactive_i,c(1,2,column,column_not)]
        transcripts_reactive = transcripts_reactive[order(transcripts_reactive[,3],decreasing = T),]
        
        # Exon counts RAW
        column = c(1:ncol(fc_exons_raw))[colnames(fc_exons_raw) == reactive_i]
        column_not = c(1:ncol(fc_exons_raw))[colnames(fc_exons_raw) != reactive_i]
        column_not = column_not[-c(1,2,3,4,5,length(column_not))] 
        fc_exons_raw_reactive = fc_exons_raw[fc_exons_raw$geneID == candidates$geneID[i()], c(1,2,3,4,5,column,column_not)]
        fc_exons_raw_reactive = fc_exons_raw_reactive[order(fc_exons_raw_reactive[,4]),]

        # Exon counts TPM
        fc_exons_tpm_reactive = fc_exons_tpm[fc_exons_tpm$geneID == candidates$geneID[i()], c(1,2,3,4,5,column,column_not)]
        fc_exons_tpm_reactive = fc_exons_tpm_reactive[order(fc_exons_tpm_reactive[,4]),]
        
        # OUTRIDER
        table_OUTRIDER = gwOUTRIDER[gwOUTRIDER$sampleID == candidates$proband[i()],]
        table_OUTRIDER = table_OUTRIDER[order(table_OUTRIDER$chr,table_OUTRIDER$pos),]
        table_OUTRIDER_candidate = candidates_OUTRIDER[candidates_OUTRIDER$sampleID == candidates$proband[i()],]    
        table_OUTRIDER_candidate = table_OUTRIDER_candidate[table_OUTRIDER_candidate$geneID == candidates$geneID[i()],]

        table_perexons_OUTRIDER_candidate = candidates_perexons_OUTRIDER[candidates_perexons_OUTRIDER$sampleID == candidates$proband[i()],]    
        table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[table_perexons_OUTRIDER_candidate$geneID == candidates$geneID[i()],]
        table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[order(table_perexons_OUTRIDER_candidate$exonID),]

        table_perexons_OUTRIDER_significant = significant_perexons_OUTRIDER[significant_perexons_OUTRIDER$sampleID == candidates$proband[i()],]    
        table_perexons_OUTRIDER_significant = table_perexons_OUTRIDER_significant[order(table_perexons_OUTRIDER_significant$chr,table_perexons_OUTRIDER_significant$pos,table_perexons_OUTRIDER_significant$exonID),]
        
        # Plotly family of proband
        fc_exons_ggplot_reactive = fc_exons_tpm_ggplot[fc_exons_tpm_ggplot$proband == reactive_i, ]
        fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$geneID == candidates$geneID[i()], ]
        fc_exons_ggplot_reactive_family = fc_exons_ggplot_reactive[gsub('_0[123]$','',fc_exons_ggplot_reactive$PatientID) == gsub('_03$','',reactive_i),]

        # Plotly proband
        fc_exons_ggplot_reactive_patient = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$PatientID == reactive_i,]
        fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$PatientID != reactive_i,] 
        
            # By age
            fc_exons_ggplot_reactive_child = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$age <  18,]
            fc_exons_ggplot_reactive_adults = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$age >=18,]
        
            # X-axis stagger
            fc_exons_ggplot_reactive_child$exonID = fc_exons_ggplot_reactive_child$exonID - 0.1
            fc_exons_ggplot_reactive_adults$exonID = fc_exons_ggplot_reactive_adults$exonID + 0.1
          
        output = list(transcripts_reactive,table_OUTRIDER,table_OUTRIDER_candidate,fc_exons_ggplot_reactive_child,fc_exons_ggplot_reactive_adults,fc_exons_tpm_reactive,fc_exons_raw_reactive,fc_exons_ggplot_reactive,fc_exons_ggplot_reactive_patient,fc_exons_ggplot_reactive_family,table_perexons_OUTRIDER_candidate,table_perexons_OUTRIDER_significant)
        names(output)  = c('transcripts_reactive','table_OUTRIDER','table_OUTRIDER_candidate','fc_exons_ggplot_reactive_child','fc_exons_ggplot_reactive_adults','fc_exons_tpm_reactive','fc_exons_raw_reactive','fc_exons_ggplot_reactive','fc_exons_ggplot_reactive_patient','fc_exons_ggplot_reactive_family','table_perexons_OUTRIDER_candidate','table_perexons_OUTRIDER_significant')
        output
    })
  
    ### Data table
    output$exonTPM <- renderDT({
      if(input$table_choice == 'Normalised expression') {
        datatable(
          reactive_inputs()$fc_exons_tpm_reactive,
          rownames = FALSE,
          options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
        } else if (input$table_choice == "Raw counts") {
        datatable(
          reactive_inputs()$fc_exons_raw_reactive,
          rownames = FALSE,
          options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
        } else if (input$table_choice == "Isoform-specific expression") {
        datatable(
          reactive_inputs()$transcripts_reactive,
          rownames = FALSE,
          options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
      }
    })
    
    
    ### Data table
    output$Table <- renderDT({
      datatable(
        reactive_inputs()$transcripts_reactive,
        rownames = FALSE,
        options = list(pageLength = 50)
        )
    })
    
    ### OUTRIDER table
    output$table_OUTRIDER <- renderDT({
      datatable(
        reactive_inputs()$table_OUTRIDER[reactive_inputs()$table_OUTRIDER$pValue < as.numeric(input$pvalue),],
        rownames = FALSE,
        options = list(pageLength = 100)
      )
    })
    
    ### OUTRIDER candidate gene table
    output$candidates_OUTRIDER <- renderDT({
      datatable(
        reactive_inputs()$table_OUTRIDER_candidate,
        rownames = FALSE,
        options = list(dom = 'tir',info = FALSE)
        )  
    })
    
    ### OUTRIDER candidate exon table 
    output$candidates_OUTRIDER_exons <- renderDT({
      datatable(
        reactive_inputs()$table_perexons_OUTRIDER_candidate,
        rownames = FALSE,
        options = list(pageLength = 50)
      )  
    })   
  
    ### OUTRIDER significant exon table 
    output$significant_OUTRIDER_exons <- renderDT({
      datatable(
        reactive_inputs()$table_perexons_OUTRIDER_significant,
        rownames = FALSE,
        options = list(pageLength = 50)
      )  
    })
   
    ### Plotly
    output$Expression <- renderPlotly({
      #entire cohort
      plot_ly() %>% 
        add_trace(
          data = reactive_inputs()$fc_exons_ggplot_reactive_child,
          x = ~exonID,
          y = ~expression,
          type = "box",
          boxpoints = FALSE,
          name = 'Average (<18yrs)',
          color = I('darkblue'),
          hoverinfo = 'none',
          marker = list(size = 12)) %>%
        add_trace(
          data = reactive_inputs()$fc_exons_ggplot_reactive_adults,
          x = ~exonID,
          y = ~expression,
          type = "box",
          boxpoints = FALSE,
          name = 'Average (>18yrs)',
          color = I('darkgreen'),
          hoverinfo = 'none',
          marker = list(size = 12))  %>%
        add_trace(
          data = reactive_inputs()$fc_exons_ggplot_reactive_patient,
          x = ~exonID,
          y = ~expression,
          type = 'scatter',
          mode = 'markers',
          hoverinfo = 'text',
          marker = list(opacity = 1),
          color = I('darkorange'),
          size = 22,
          name = ~PatientID,
          customdata= ~Sexe,
          text = ~as.factor(age),
          hovertemplate = paste('<b>Exon number</b>: %{x}',
                                        '<br><b>Proband Age</b>: %{text}',
                                        '<br><b>Sex</b>: %{customdata}<br>')
                  ) %>%
        layout(
          xaxis = list(title = "Exon Number",
                       ticktext = 1:(max(reactive_inputs()$fc_exons_ggplot_reactive_patient$exonID,1)), 
                       tickvals = as.list(1:(max(reactive_inputs()$fc_exons_ggplot_reactive_patient$exonID,1))),
                       tickmode = "array"),
          yaxis = list(title = "Normalised Expression (TPM)"),
          title = paste0('Cohort expression (per exon) for ',reactive_inputs()$fc_exons_ggplot_reactive_patient$geneID[1])
      )
})

      
    ### Plotly per family
    output$Expression_perfamily <- renderPlotly({
      plot_ly() %>% 
        add_trace(
          data = reactive_inputs()$fc_exons_ggplot_reactive_family,
          x = ~exonID,
          y = ~expression,
          type = 'scatter',
          mode = 'markers',
          hoverinfo = 'text',
          marker = list(opacity = 1),
          color = ~PatientID,
          colors = c('darkorange','black','blue')[length(unique(reactive_inputs()$fc_exons_ggplot_reactive_family$PatientID)):1],
          size = 22,
          name = ~PatientID,
          customdata= ~Sexe,
          text = ~as.factor(age),
          hovertemplate = paste('<b>Exon number</b>: %{x}',
                                '<br><b>Age</b>: %{text}',
                                '<br><b>Sex</b>: %{customdata}<br>')
        ) %>%
        layout(
          xaxis = list(title = "Exon Number",
                       ticktext = 1:(max(reactive_inputs()$fc_exons_ggplot_reactive_family$exonID,1)), 
                       tickvals = as.list(1:(max(reactive_inputs()$fc_exons_ggplot_reactive_family$exonID,1))),
                       tickmode = "array"),
          yaxis = list(title = "Normalised Expression (TPM)"),
          title = paste0('Family expression (per exon) for ',reactive_inputs()$fc_exons_ggplot_reactive_patient$geneID[1])
        )
    })

    ### Sashimi plots FRASER
    output$Sashimi = renderImage({
      list(
        src = paste0(params$datadir,"/sashimis/gene_",candidates$geneID[i()],"_",candidates$proband[i()],'_sashimi.png')[1], #path to the file
        contentType = "image/png",
        width = 900,
        alt = "My Figure"
        )}, deleteFile = FALSE)
    
    ### Image
    output$mainimage = renderImage({
      list(
        src = file.path(params$datadir,'CHUSJ_CR_Bioinformatique_V2.png'), #path to the file
        contentType = "image/png",
        width = 400,
        alt = "My Figure"
        )},deleteFile = FALSE) 
    
    ### Description
    output$description <- renderUI({
      print(paste0('Selecting ~~~ ',input$proband2))
      selected_ensembl <- candidates$ensembl[i()]
      selected_geneID <- candidates$geneID[i()]
      #url = paste0("https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", selected_ensembl)
      #url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", selected_geneID)
      url = paste0("https://www.proteinatlas.org/", selected_ensembl)
      HTML(
        paste0("<span><b>Gene description for ",selected_geneID,": </b> <a href='", url, "' target='_self'>",selected_ensembl,"</a>",
               "<br><br><b>Mutations: </b>",clinical$Mutation[clinical$PatientID == input$proband2],
               "<br><br><b>Candidate Gene hypothesis: </b>",clinical$Hypothèse[clinical$PatientID == input$proband2],
               "<br><br><b>HPO terms: </b>",clinical$`HPO terms`[clinical$PatientID == input$proband2],
               "<br><br><b>Bam file location (Fir): </b>",system('echo ${HOME}',intern = T),"/scratch/raredisease_rnaseq/results_06_01_2026/star_salmon/",candidates$proband[i()],".sorted.bam<span>")
        )
      })
    
    ### Figure legend
    output$Figure_genemodel_legend  <- renderUI({
      HTML(
        paste0("<span><b>Figure 2:</b> Visualisation des altérations d’épissage détectées par l'outil FRASER.
        <br><b>A:</b> Carte des introns/exons du gène ",candidates$geneID[i()]," et localisation du/des variant.s (ligne bleue pointillée).
        <br><b>B:</b> Évènements d’épissage aberrant pour chaque région détectée (−log₁₀ p-value).
        <br><b>C:</b>  Couverture de séquençage normalisée pour le probant ainsi que 25-75ième percentile de la population de référence (en orange).</span>")
        )
      })
 
    ### Versioning
    output$Version <- renderDT({
      datatable(
        data.frame(Parameter=names(unlist(report_version)),Value=unlist(report_version)),
        rownames = F,options = list(dom = 'p'))
      })
    
    ### multiQC
    output$htmlViewer <- renderUI({
        HTML(paste(readLines(html_file), collapse = "\n"))
    })
    
    ### FASTA
    output$fasta <- renderUI({
      lines <- readLines(paste0(params$datadir,"/consensus/","gene",candidates$geneID[i()],'_',candidates$proband[i()],'.fasta'))
      HTML(paste0("<span style='color: black; font-family: Courier New; font-size: 16px;'>",paste0(lines,collapse = '<br>'),"</span>",collapse = "\n"))
    })
    
    ### IGV
    observeEvent(input$addBamLocalFileButton, {
      gene_dir = paste0(params$datadir,'/bams_subset/gene',candidates$geneID[i()],'_chr',candidates$chromosome[i()],'_',candidates$start[i()]-5000,'_',candidates$stop[i()]+5000,'/')
      color=
      showGenomicRegion(session, id="igvShiny",paste0("chr",candidates$chromosome[i()],":",candidates$start[i()],"-",candidates$stop[i()]))
      bamFile <- paste0(gene_dir,candidates$proband[i()],"_sorted_chrN.bam")
      bamAlign <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq"))
      loadBamTrackFromLocalData(session, id="igvShiny", trackName=input$proband2, data=bamAlign)
      runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'green';")
      })

    ### Default hg38 view
    output$igvShiny <- renderIgvShiny({
      runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'red';")
      genomeOptions <- parseAndValidateGenomeSpec('hg38',initialLocus=paste0("chr",candidates$chromosome[i()],":",candidates$start[i()],"-",candidates$stop[i()]))
      igvShiny(genomeOptions)
      })
    
    #Reative slider
    output$genemodel_slider <- renderUI({
      cmin <- floor(candidates$start[i()]/1000)
      cmax <- ceiling(candidates$stop[i()]/1000)
      sliderInput(
        inputId = "sliderxlims",
        label = "Select genomic window (Kb)",
        min = cmin,
        max = cmax,
        value = c(cmin, cmax),
        width = '80%'
        )
      })
    
    ### Coverage ggplots
    output$Figure_genemodel = renderPlot({
      if(is.null(input$sliderxlims)) {genemodel = ggplot() +  theme_void() + geom_text(aes(0,0,label='Plotting in ¨Progress')) + xlab(NULL)} else {
        gene_dir = paste0(params$datadir,'/bams_subset/gene',candidates$geneID[i()],'_chr',candidates$chromosome[i()],'_',candidates$start[i()]-5000,'_',candidates$stop[i()]+5000,'/')
        genemodel = plotting_coverage(
          candidate = candidates[i(),],
          depth_file = paste0(gene_dir,"gene_",candidates$geneID[i()],"_",candidates$proband[i()],"_depth5.csv"),
          res_dt_candidate_gene_file = paste0(gene_dir,"gene_",candidates$geneID[i()],"_",candidates$proband[i()],"_res_dt_candidate_gene.csv"),
          bam_file = paste0(gene_dir,candidates$proband[i()],"_sorted_chrN.bam"),
          colmean_genes_counts_file = paste0(params$datadir,'/colmean_genes_counts.tsv'),
          gene_annotations=gene_annotations,
          xlims = input$sliderxlims)}
      
      genemodel
      })

    #genome-wide OUTRIDER
    output$gwOUTRIDER = renderPlot({
      manhattan_plot(res_dt=gwOUTRIDER,sample = candidates$proband[i()],geneID = 'geneID', pvalue = 'pValue',pcutoff = 0.001)
    })
    
    #genome-wide FRASER
    output$gwFRASER = renderPlot({
      manhattan_plot(res_dt=gwFRASER,sample = candidates$proband[i()])
    })
    
    output$gwFRASER_table = renderDT({
      datatable(
        gwFRASER_table(res_dt=gwFRASER,sample = candidates$proband[i()]),
        rownames = F,options = list(dom = 'p'))
    })
  }

######shiny app
app = shinyApp(ui, server, options = list(height = 900))
runApp(app, launch.browser = TRUE)

