#Were is the data.zip located. First check if you specified it as a positional argument. Then, check if it is in the current dir. Then check in my default. Finally error out if you cannnot find it.
params = list(zipfile = commandArgs(trailingOnly=TRUE)[1])

#create a temp directory to store things everytime.
tempdir = 'temp'
tempdata = paste0('temp/temp_',format(Sys.time(), format = "%Y_%m_%d_%H_%M_%S"))
dir.create(tempdir,showWarnings = FALSE)
dir.create(tempdata,showWarnings = FALSE)
params$datadir = file.path(getwd(),tempdata,'data')
#params$datadir= "/Users/sebastienrenaut/Documents/St-Justine/nextflow_rnasplice/temp/temp_2025_12_12_13_30_23/data"


#
if(is.na(params$zipfile)) params$zipfile = file.path(getwd(),system('ls -t data*zip|head -1', intern = T))
if(!file.exists(params$zipfile)) {params$zipdir = file.path("/Users/sebastienrenaut/Documents/St-Justine/nextflow_rnasplice/",system('ls -t data*zip|head -1', intern = T));params$datadir = "/Users/sebastienrenaut/Documents/St-Justine/nextflow_rnasplice/data"}
if(!file.exists(params$zipfile)) stop("Cannnot find the data.zip. Please specify its location")

#R libraries that you need
packages = c('DT','plotly','tidyr','shiny','shinyjs','jsonlite','igvShiny','GenomicAlignments','dplyr','ggtranscript','patchwork','Hmisc')

for(p in 1:length(packages)) {
  if(packages[p] %in% installed.packages()) {
    library(packages[p],character.only = TRUE)
    } else if (packages[p] == 'ggtranscript'){
      stop('Install ggtranscript with : remotes::install_github("dzhang32/ggtranscript")')
    } else {
      stop(paste0('Error in library() : there is no package called ', packages[p]))
    }
}

######
#Unzip and load files
######
print(paste0("Using the latest data.zip: ",params$zipfile,' unzipping in: ',params$datadir))
unzip(params$zipfile,exdir = sub('data','.',params$datadir),overwrite=TRUE)

# Create a virtual path "/pdfs" that maps to that directory
#addResourcePath("pdfs",file.path(params$datadir,'/pdfs'))

#fc_exons
#source(file.path('~/Desktop/rnaseq_shinyhelper_functions.R'))
source(file.path(params$datadir,'scripts/rnaseq_shinyhelper_functions.R'))
fc_exons_raw = read.table(file.path(params$datadir,'fc_exons_raw.tsv'),sep = '\t',check.names = F);fc_exons_raw[,-c(1:5)] = round(fc_exons_raw[,-c(1:5)])
fc_exons_tpm = read.table(file.path(params$datadir,'fc_exons_tpm.tsv'),sep = '\t',check.names = F)
fc_genes_tpm = read.table(file.path(params$datadir,'fc_genes_tpm.tsv'),sep = '\t',check.names = F)

results_OUTRIDER = read.csv(file.path(params$datadir,'results_OUTRIDER.tsv'),sep = '\t',check.names = F)
results_OUTRIDER$Chr = factor(gsub('chr','',results_OUTRIDER$Chr),levels = c(1:22,'X','Y'))

candidates_OUTRIDER = read.csv(file.path(params$datadir,'candidates_OUTRIDER.tsv'),sep = '\t',check.names = F)
candidates_perexons_OUTRIDER = read.csv(file.path(params$datadir,'candidates_perexons_OUTRIDER.tsv'),sep = '\t',check.names = F)

transcripts_named_filtered = read.csv(file.path(params$datadir,'transcripts_named_filtered.tsv'),sep = '\t',check.names = F)
transcripts_named_filtered_ggplot = read.csv(file.path(params$datadir,'transcripts_named_filtered_ggplot.tsv'),sep = '\t',check.names = F)

fc_exons_tpm_ggplot = read.csv(file.path(params$datadir,'fc_exons_tpm_ggplot.tsv'),sep = '\t',check.names = F)
candidates = read.csv(file.path(params$datadir,'candidate_genes_3.txt'))
clinical = read.csv(file.path(params$datadir,'clinical.tsv'),sep = '\t',check.names = F)
html_file = file.path(params$datadir,'multiqc_report.html')

report_version = read_json(file.path(params$datadir,'/VERSION.json'))
report_version$data = params$zipfile

load(file = file.path(params$datadir,"gene_annotations.rda"))

#####
#####UI
#####    
ui <- navbarPage(
    tabsetPanel(
      # Page 1: Selection List
      tabPanel(
        'Proband Selection',
        width = "100%",
        mainPanel(
          br(),
          imageOutput("mainimage"),
        ),
        selectInput(
          inputId = "proband2",
          label   = "Select a patient:",
          choices = candidates$proband2,
          selected = 'HSJ_001_03'
        ),
        mainPanel(
          width = "100%",
          h2("Information"),
          br(),
          htmlOutput("description"),
          br(),
          h2("Software Version"),
          DTOutput("Version")
        ),
      ),  
      
      ### Data table
      tabPanel(
          "Expression",
          h3('Gene expression per exon'),
          'Reference transcript according to: Matched Annotation from NCBI and EMBL-EBI (MANE)',
          br(),
          br(),
          br(),
          width = "100%",
          selectInput("table_choice", "Choose an expression metric: ",
                      choices = c("Normalised expression", "Raw counts","Isoform-specific expression"),
                      selected = "Normalised expression"),
            mainPanel(
          #    width = 9,
          #    heigth = 16, 
                DTOutput("exonTPM"),
            )
      ),
    
      ### Plot of Expression
      tabPanel(width = "100%",
        "Plot",
         h3('Gene expression per exon'),
         br(),
           mainPanel(
                plotlyOutput("Expression",width = "1200px"),
                br(),
                plotlyOutput("Expression_perfamily",width = "1200px"),
        )
    ),
    
    ### Data table
      tabPanel(
        "OUTRIDER",
        h3('OUTRIDER - OUTlier in RNA-Seq fInDER'),
        'Identification of aberrant gene expression in RNA-seq data, Outliers are identified as read counts that significantly deviate from the population',
        br(),
        '(all probands, including adults)',
        br(),
        width = "100%",
        mainPanel(
          width = 9,
          heigth = 16, 
          h4('Candidate gene'),
          DTOutput("candidates_OUTRIDER"),
          br(), 
          h4('Candidate exons'),
          DTOutput("candidates_OUTRIDER_exons"),
          br(), 
          br(),
          h4('All significant genes (p < 0.005)'),
          DTOutput("table_OUTRIDER"),
          br(), 
          
        )
    ),
    
    ### Plot of Expression
    tabPanel(width = "100%",
             "Gene model",
             h3('Gene model and significance values for FRASER'),
             selectInput("zoom_choice", "Gene localisation: ",
                         choices = c('full gene','mutation zoom (+/-1kb)','mutation zoom (+/-5kb)'),
                         selected = 'full gene'),
             br(),
             mainPanel(
               uiOutput("Figure_genemodel_dynamic"),
               br(),
               htmlOutput('Figure_genemodel_legend')
             )
    ),
    
  
    ### Plot of Structural variation
    tabPanel("IGV",
           #  width = 9,
            # heigth = 16,
             useShinyjs(), # Initialize shinyjs
             br(),
             br(),
             actionButton("addBamLocalFileButton",
                          "Show gene alignment",
                          style = "background-color: red; color: white; border-color: darkorchid;"),
             br(),
             br(),
             mainPanel(
              igvShinyOutput('igvShiny',width = "1200px")
             )
    ),
  
    ### Plot of Structural variation
    tabPanel(width = "100%",
             "FRASER",
             h3('FRASER: Find RAre Splicing Events in RNAseq Data'),
             'Splice site map on top, with new exons and junctions in proband in yellow or red line.',
             br(),
             'Splice site map generated in an annotation-free fashion based on RNA-seq coverage of split reads.',
             br(),
             'As such, introns/exons may differ from actual genome annotation.',
             br(),
             'In red, proband of interest. In blue, five representative samples of the population.',
             br(),
             mainPanel(
               imageOutput("Sashimi",width = "1200px")
             )
    ),
    
    ### Plot of Structural variation
    tabPanel(width = "100%",
             "fasta",
             h3('Transcribed sequences: reference / alternate'),
             em('intron retention events: dashes (---), variants in'),
             em(strong('bold')),
             br(),
             br(),
             br(),
             mainPanel(
               width = "100%",
               htmlOutput("fasta")
             )
    ),
    
    
      #multiQC
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
  
    #SELECTED VARIABLE
    output$selected_var <- renderText({input$proband2})
  
    #REACTIVE PLOTS
    reactive_ggplots <- reactive({
    plotting_coverage(
      candidate = candidates[candidates$proband2 == input$proband2,],
      depth_file = paste0(params$datadir,"/gene_statistics/gene_",candidates$geneID[candidates$proband2 == input$proband2],"_",candidates$proband[candidates$proband2 == input$proband2],"_depth5.csv"),
      res_dt_candidate_gene_file = paste0(params$datadir,"/gene_statistics/gene_",candidates$geneID[candidates$proband2 == input$proband2],"_",candidates$proband[candidates$proband2 == input$proband2],"_res_dt_candidate_gene.csv"),
      bam_file = paste0(params$datadir,'/gene_statistics/',candidates$proband[candidates$proband2 == input$proband2],"_sorted_chrN.bam"),
      colmean_genes_counts_file = paste0(params$datadir,'/colmean_genes_counts.tsv'),
      gene_annotations=gene_annotations,
      zoom = input$zoom_choice)
    })
    
    #REACTIVE INPUT (AS A LIST)
    reactive_inputs <- reactive({
        reactive_i  = strsplit(input$proband2,'@')[[1]][1]
        multigene_sample = strsplit(input$proband2,'@')[[1]][2]
        column = c(1:ncol(transcripts_named_filtered))[colnames(transcripts_named_filtered) == reactive_i]
        column_not = c(1:ncol(transcripts_named_filtered))[colnames(transcripts_named_filtered) != reactive_i]
        column_not = column_not[-c(1,2,length(column_not))] 
        transcripts_reactive = transcripts_named_filtered[transcripts_named_filtered$proband == reactive_i,c(1,2,column,column_not)]
        transcripts_reactive = transcripts_reactive[order(transcripts_reactive[,3],decreasing = T),]
        
        #Exon Counts RAW
        column = c(1:ncol(fc_exons_raw))[colnames(fc_exons_raw) == reactive_i]
        column_not = c(1:ncol(fc_exons_raw))[colnames(fc_exons_raw) != reactive_i]
        column_not = column_not[-c(1,2,3,4,5,length(column_not))] 
        fc_exons_raw_reactive = fc_exons_raw[fc_exons_raw$geneID == candidates$geneID[candidates$proband2 == input$proband2], c(1,2,3,4,5,column,column_not)]
        fc_exons_raw_reactive = fc_exons_raw_reactive[order(fc_exons_raw_reactive[,4]),]

        #Exon Counts TPM
        fc_exons_tpm_reactive = fc_exons_tpm[fc_exons_tpm$geneID == candidates$geneID[candidates$proband2 == input$proband2], c(1,2,3,4,5,column,column_not)]
        fc_exons_tpm_reactive = fc_exons_tpm_reactive[order(fc_exons_tpm_reactive[,4]),]
        
        #OUTRIDER
        table_OUTRIDER = results_OUTRIDER[results_OUTRIDER$proband == reactive_i,]
        table_OUTRIDER = table_OUTRIDER[order(table_OUTRIDER$Chr,table_OUTRIDER$start),]
      #  table_OUTRIDER_candidate = candidates_OUTRIDER[(candidates_OUTRIDER$geneID ==  candidates$geneID[candidates$proband2 ==  input$proband2]) & (candidates_OUTRIDER$proband == reactive_i),]    
        table_OUTRIDER_candidate = candidates_OUTRIDER[candidates_OUTRIDER$sampleID == reactive_i,]    
        
        table_perexons_OUTRIDER_candidate = candidates_perexons_OUTRIDER[candidates_perexons_OUTRIDER$sampleID == reactive_i,]    
        table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[order(table_perexons_OUTRIDER_candidate$exonID),]
        table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[,]
        
        #Plotly family of proband
        fc_exons_ggplot_reactive = fc_exons_tpm_ggplot[fc_exons_tpm_ggplot$proband == reactive_i, ]
        if(!is.na(multigene_sample)) fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$geneID == multigene_sample, ]
        fc_exons_ggplot_reactive_family = fc_exons_ggplot_reactive[gsub('_0[123]$','',fc_exons_ggplot_reactive$PatientID) == gsub('_03','',reactive_i),]
        

        #Plotly
        fc_exons_ggplot_reactive_patient = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$PatientID == reactive_i,]
        fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$PatientID != reactive_i,] 
        
            #divide by age
            fc_exons_ggplot_reactive_child = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$age <  18,]
            fc_exons_ggplot_reactive_adults = fc_exons_ggplot_reactive[fc_exons_ggplot_reactive$age >=18,]
        
            #stager the positions
            fc_exons_ggplot_reactive_child$exonID = fc_exons_ggplot_reactive_child$exonID - (1/max(fc_exons_ggplot_reactive_child$exonID)*2)
            fc_exons_ggplot_reactive_adults$exonID = fc_exons_ggplot_reactive_adults$exonID + (1/max(fc_exons_ggplot_reactive_adults$exonID)*2)
        
          
          output = list(transcripts_reactive,table_OUTRIDER,table_OUTRIDER_candidate,fc_exons_ggplot_reactive_child,fc_exons_ggplot_reactive_adults,fc_exons_tpm_reactive,fc_exons_raw_reactive,fc_exons_ggplot_reactive,fc_exons_ggplot_reactive_patient,fc_exons_ggplot_reactive_family,table_perexons_OUTRIDER_candidate)
          names(output)  = c('transcripts_reactive','table_OUTRIDER','table_OUTRIDER_candidate','fc_exons_ggplot_reactive_child','fc_exons_ggplot_reactive_adults','fc_exons_tpm_reactive','fc_exons_raw_reactive','fc_exons_ggplot_reactive','fc_exons_ggplot_reactive_patient','fc_exons_ggplot_reactive_family','table_perexons_OUTRIDER_candidate')
          output
    })
  
    #DATA TABLE
    output$exonTPM <- renderDT({
      if(input$table_choice == 'Normalised expression') {datatable(
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
    
    
    #DATA TABLE
    output$Table <- renderDT({
        datatable(
          reactive_inputs()$transcripts_reactive,
            rownames = FALSE,
            options = list(pageLength = 50)
        )
    })
    
    #OUTRIDER TABLE
    output$table_OUTRIDER <- renderDT({
      datatable(
        reactive_inputs()$table_OUTRIDER,
        rownames = FALSE,
        options = list(pageLength = 100)
      )
    })
    
    #OUTRIDER TABLE CANDIDATE
    output$candidates_OUTRIDER <- renderDT({
      datatable(
        reactive_inputs()$table_OUTRIDER_candidate,
        rownames = FALSE,
        options = list(dom = 'tir')
        )  
    })
    
    #OUTRIDER TABLE CANDIDATE EXONS
    output$candidates_OUTRIDER_exons <- renderDT({
      datatable(
        reactive_inputs()$table_perexons_OUTRIDER_candidate,
        rownames = FALSE,
        options = list(pageLength = 50)
      )  
    })
    
  
  
    #PLOTLY PLOT
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
                       ticktext = 1:(max(reactive_inputs()$fc_exons_ggplot_reactive_patient$exonID)), 
                       tickvals = as.list(1:(max(reactive_inputs()$fc_exons_ggplot_reactive_patient$exonID))),
                       tickmode = "array"),
          yaxis = list(title = "Normalised Expression (TPM)"),
          title = paste0('Cohort expression (per exon) for ',reactive_inputs()$fc_exons_ggplot_reactive_patient$geneID[1])
      )

})

    
      
    #PLOTLY BOXPLOT per family
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
                       ticktext = 1:(max(reactive_inputs()$fc_exons_ggplot_reactive_family$exonID)), 
                       tickvals = as.list(1:(max(reactive_inputs()$fc_exons_ggplot_reactive_family$exonID))),
                       tickmode = "array"),
          yaxis = list(title = "Normalised Expression (TPM)"),
          title = paste0('Family expression (per exon) for ',reactive_inputs()$fc_exons_ggplot_reactive_patient$geneID[1])
        )
    })
  
    
    
    
    
    #SASHIMI
    output$Sashimi <- renderImage({
    # Return a list containing the filename and other details
    list(
      src = paste0(params$datadir,"/sashimis/gene_",candidates$geneID[candidates$proband2 == input$proband2],"_",candidates$proband[candidates$proband2 == input$proband2],'_sashimi.png')[1],   # path to the file
      contentType = "image/png",
      width = 900,
      alt = "My Figure"
    )
      }, deleteFile = FALSE)   # don't delete the original file
    
    #Image
    output$mainimage = 
      renderImage({
        # Return a list containing the filename and other details
        list(
          src = file.path(params$datadir,'CHUSJ_CR_Bioinformatique_V2.png'),   # path to the file
          contentType = "image/png",
          width = 400,
          alt = "My Figure"
        )
      },deleteFile = FALSE) 
    
    #DESCRIPTIONS
    output$description <- renderUI({
      print(paste0('Selecting ~~~ ',input$proband2))
      selected_ensembl <- candidates$ensembl[candidates$proband2 == input$proband2]
      selected_geneID <- candidates$geneID[candidates$proband2 == input$proband2]
      url <- paste0("https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", selected_ensembl)
      HTML(
        paste0("<b>Gene description for ",selected_geneID,"</b>: <a href='", url, "' target='_self'>", selected_ensembl, "</a>",
               "<br><br><b>Mutations: </b>",clinical$Mutation[clinical$PatientID ==  strsplit(input$proband2,'@')[[1]][1]],
               "<br><br><b>Candidate Gene hypothesis: </b>",clinical$Hypothèse[clinical$PatientID ==  strsplit(input$proband2,'@')[[1]][1]],
               "<br><br><b>HPO terms: </b>",clinical$`HPO terms`[clinical$PatientID ==  strsplit(input$proband2,'@')[[1]][1]],
               "<br><br><b>Bam file location: </b>",system('echo ${HOME}',intern = T),"/scratch/nextflow_rnasplice/bams/",candidates$proband[candidates$proband2 == input$proband2],".sorted.bam<br><br>")
      )
    })
    
    #Figure legend
    output$Figure_genemodel_legend  <- renderUI({HTML(
      paste0("<b>Figure 2:</b> Visualisation des altérations d’épissage détectées par l'outil FRASER.
      <br><b>A</b>: Carte des introns/exons du gène ",candidates$geneID[candidates$proband2==input$proband2]," et localisation du/des variant.s (ligne bleue pointillée).
      <br><b>B:</b> Évènements d’épissage aberrant pour chaque région détectée (−log₁₀ p-value).
      <br><b>C:</b>  Couverture de séquençage normalisée pour le probant ainsi que 25-75ième percentile de la population de référence (en orange).
   ")
    )})
 
    #Versioning
    output$Version <- renderDT({
      datatable(
        data.frame(Parameter=names(unlist(report_version)),Value=unlist(report_version)),
        rownames = F,options = list(dom = 'p'))
    })
    
    #MULTI QC reports
    output$htmlViewer <- renderUI({
        HTML(paste(readLines(html_file), collapse = "\n"))
    })
    
    
    #FASTA sequence
    output$fasta <- renderUI({
      # Read the lines of the uploaded text file
      lines <- readLines(paste0(params$datadir,"/consensus/","gene",candidates$geneID[candidates$proband2 == input$proband2],'_',candidates$proband[candidates$proband2 == input$proband2],'.fasta'))
      #Print the lines to the output
      HTML(paste0("<span style='color: black; font-family: Courier New; font-size: 16px;'>",paste0(lines,collapse = '<br>'),"</span>",collapse = "\n"))
    })
    
    #IGV
    observeEvent(input$addBamLocalFileButton, {
      color=
      showGenomicRegion(session, id="igvShiny",paste0("chr",candidates$chromosome[candidates$proband2 == input$proband2],":",candidates$start[candidates$proband2 == input$proband2],"-",candidates$stop[candidates$proband2 == input$proband2]))
      bamFile <- paste0(params$datadir,'/bams_subset/',candidates$proband[candidates$proband2 == input$proband2],"_sorted_chrN.bam")
      bamAlign <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq"))
      loadBamTrackFromLocalData(session, id="igvShiny", trackName=input$proband2, data=bamAlign)
      runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'green';")
      })
    
    #default hg38 view
    output$igvShiny <- renderIgvShiny({
    runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'red';")
     genomeOptions <- parseAndValidateGenomeSpec('hg38',initialLocus=paste0("chr",candidates$chromosome[candidates$proband2 == input$proband2],":",candidates$start[candidates$proband2 == input$proband2],"-",candidates$stop[candidates$proband2 == input$proband2]))
     igvShiny(genomeOptions)
  })
    
    #ggplots for coverage    
    output$Figure_genemodel = renderPlot({reactive_ggplots()})
    
    #ggplots for coverage    
    output$Figure_genemodel_dynamic <- renderUI({
      width = '1200px'
      
      if(input$zoom_choice!='full gene'){
      #Check if you have 1 or 2 mutations.
      mut_pos = candidates$position[candidates$proband2 == input$proband2]
      width = ifelse(length(grep('_',mut_pos))==1,'1600px','1200px')
      }
  
      # Return a plotOutput with the specified width and a fixed height
      plotOutput("Figure_genemodel", width = width,height = "800px")
    })
}

######shiny app
app = shinyApp(ui, server, options = list(height = 900))
runApp(app, launch.browser = TRUE)

