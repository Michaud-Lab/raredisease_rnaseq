# Parse command-line arguments
# Usage: Rscript RNAseq_shiny_v2.5.R --data_minimal --use_password
args = commandArgs(trailingOnly = TRUE)
use_data_minimal = "--data_minimal" %in% args
use_password     = "--use_password"  %in% args

# Resolve directories:
data_dir = if (exists("use_data_minimal") && use_data_minimal) "data_minimal" else "data"
scripts_dir = if (dir.exists(file.path(data_dir,"scripts"))) {
  file.path(data_dir,"scripts")} else {
  "scripts"
}

# source
source(file.path(scripts_dir, "Shiny/global.R"))
source(file.path(scripts_dir, "Shiny/rnaseq_shinyhelper_functions.R"))
source(file.path(scripts_dir, "Shiny/reactive_module.R"))


#####
##### Credentials (SQLite database — run scripts/Shiny/create_credentials_db.R once to create it)
#####
if (use_password) {
  credentials_db = file.path(data_dir, "credentials.sqlite")
  if (!file.exists(credentials_db))
    stop("Credentials database not found. Run: Rscript scripts/Shiny/create_credentials_db.R")
}

#####
##### UI
#####
logger::log_info('Defining UI')
app_ui = page_fluid(
  theme = theme,

  # Dark title header
  title = "RNAseq dashboard",
  div(
    class = "bg-dark text-white p-3 mb-4",
    uiOutput("dynamic_title")
  ),

  tabsetPanel(

    ###  Selection & Info
    tabPanel('Proband Selection',
             card(
               layout_columns(
                 col_widths = c(7,5),
                 div(
                   selectInput(
                     inputId = "proband",
                     label = h4("Select a proband:"),
                     choices = sort(unique(candidates$proband)),
                     selected = 'HSJ_001_03_PAX'),
                   reactive_data_UI("reactive_data"),
                   div(imageOutput("mainimage"), class = "text-center")
                 ),
                 div(
                   h4("Search genes/probands:"),
                   DTOutput("candidates_table")
                 )
               )
             ),
             card(
               card_header(strong("Information")),
               htmlOutput("description")),
             card(
               card_header(strong("Software Version")),
               DTOutput("Version")
             )
    ),

    ### Data table
    tabPanel("Expression",
             card(
               card_header(strong("Reference transcript based on Matched Annotation from NCBI and EMBL-EBI (MANE)")),
               selectInput("table_choice", "Choose an expression metric: ",
                           choices = c("Normalised expression", "Raw counts","Isoform-specific expression"),
                           selected = "Normalised expression"),
               DTOutput("exonTPM")
             )
    ),

    ### Plot of Expression
    tabPanel("Plot",
             card(plotlyOutput("Expression",width = "1500px")),
             card(plotlyOutput("Expression_perfamily",width = "1500px"))
    ),

    ### Data table
    tabPanel("OUTRIDER",
             card(
               card_header(strong('OUTlier in RNA-Seq fInDER')),
               'Identification of aberrant gene expression in RNA-seq data, Outliers are identified as read counts that significantly deviate from the population',
               '(all probands, including some adults from the LC & F0 cohorts)'),
             card(
               card_header(strong('Genome-wide significance')),
               plotOutput("gwOUTRIDER", width = '1500px', height = "600px")),
             card(
               fill=FALSE,
               height='200px',
               card_header(strong('Candidate gene')),
               DTOutput("candidates_OUTRIDER")),
             card(
               card_header(strong('Candidate exons')),
               DTOutput("candidates_OUTRIDER_exons")),
             card(
               card_header(strong('All significant genes')),
               selectInput("pvalue", "Select a p-value threshold:",
                           choices = c(0.05,0.01,0.005,0.001,0.0005),
                           selected = 0.005),
               DTOutput("table_OUTRIDER")),
             card(
               card_header(strong('All significant exons (p<0.01)')),
               DTOutput("significant_OUTRIDER_exons")
             )
    ),

    ### Plot of Expression
    tabPanel("Gene model",
             card(
               card_header(strong('Gene model and splicing significance values from FRASER')),
               uiOutput("genemodel_slider"),
               plotOutput("Figure_genemodel", width = '1500px', height = "800px")
             ),
             card(
               htmlOutput('Figure_genemodel_legend'))
    ),


    ### Plot of Structural variation
    tabPanel("IGV",
             card(
               card_header(strong('Integrative Genome Viewer')),
               useShinyjs(), # Initialize shinyjs
               actionButton("addBamLocalFileButton",
                            "Show gene alignment",
                            style = "background-color: red; color: white; border-color: darkorchid;")),
             card(
               igvShinyOutput('igvShiny',width = "100%")
             )
    ),

    ### Plot of Structural variation
    tabPanel(
      "FRASER",
      card(
        card_header(strong('Find RAre Splicing Events in RNAseq Data')),
        'Identification of aberrant splicing events from RNAseq. Outliers are identified intron retention / exon skipping events that significantly deviate from the population.'),
      card(
        card_header(strong('Genome-wide significance')),
        plotOutput("gwFRASER", width = '1500px', height = "600px")),
      card(
        card_header(strong('Significant splicing events (below adj. p-value < 0.05)')),
        DTOutput("gwFRASER_table")
      ),
      card(
        card_header(strong('Splice site map')),
        'Splice site map on top, with new exons and junctions in proband in yellow or red line.',
        'Splice site map generated in an annotation-free fashion based on RNA-seq coverage of split reads.',
        'As such, introns/exons may differ from actual genome annotation.',
        'In red, proband of interest. In blue, five representative samples of the population.',
        imageOutput("Sashimi",width = "1200px")
      )
    ),

    ### gene prioritisation
    tabPanel(
      "Gene Prioritization",
      card(card_header(strong("Gene prioritization")),
           downloadButton("gp_download", "Download CSV"),
           selectInput("geneprior_rm", "Remove missing data in column:",
                       choices = c('gene score','OUTRIDER gene zScore','OUTRIDER gene pValue','FRASER gene pValue','OUTRIDER exon zScore','OUTRIDER exon pValue'),
                       selected = 'gene score'),
           DTOutput("gp"),height = "600px")
    ),

    ### ASE
    tabPanel(
      "ASE",
      card(
        card_header(strong('Allele Specific Expression')),
        htmlOutput('ase_legend')),
      card(
        card_header(strong('Genome-wide significance (Manhattan plot)')),
        plotOutput("gwASE", width = '1500px', height = "600px")),
      card(card_header(strong("Genome-wide significance (table, p<0.01)")),
           DTOutput("gwASE_table"),height = "600px"),
      card(card_header(strong("Known Imprinted genes (table)")),
           DTOutput("gwImprinted_table"),height = "600px"),
      card(card_header(strong("X-chromosome (table)")),
           DTOutput("gwX_table"),height = "600px")
    ),


    ### Plot of Structural variation
    tabPanel("fasta",
             card(
               card_header(strong('Transcribed sequences: reference / alternate')),
               em('dashes (---) in reference/alternate correspond to intron retention/exon skipping events, respectively, variants in ',strong('bold'), 'candidate variant in blue.'),
               htmlOutput("fasta")
             )
    ),

    ### Search expression
    tabPanel("Search expression",
             card(
               card_header(strong("Gene expression across samples (fraction of total reads)")),
               selectizeInput("gene_search", label = "Select a gene:",
                              choices = NULL,
                              selected = NULL,
                              options = list(placeholder = 'Type a gene name...')),
               plotlyOutput("searchExpression", height = "600px")
             ),
             card(
               card_header(strong("Haemoglobin gene expression (fraction of total reads)")),
               plotlyOutput("hb_barplot", height = "600px")
             ),
             card(
               card_header(strong("Total reads assigned to a genic region per sample")),
               plotlyOutput("total_reads_barplot", height = "600px")
             )
    ),

    ### multiQC
    tabPanel(
      "MultiQC",
      card(
        card_header(h3('Quality Control metrics')),
        selectInput("multiQCs", "Select a multiQC:",
                    choices = c(1,2),
                    selected = 1),
        height = "200px"),
      card(
        htmlOutput("htmlViewer"))
    )
  )
)
ui = if (use_password) secure_app(app_ui) else app_ui



#####
##### server
#####
logger::log_info('Defining back-end server')
server = function(input, output, session) {

  if (use_password)
    auth = secure_server(check_credentials = check_credentials(credentials_db))

  #####
  ##### reactive data (module)
  #####
  rd = reactive_data_server(
    id           = "reactive_data",
    proband      = reactive(input$proband),
    pvalue       = reactive(input$pvalue),
    geneprior_rm = reactive(input$geneprior_rm)
  )

  #####
  ##### outputs
  #####
  ### dynamic title
  output$dynamic_title = renderUI({
    titlePanel(
      paste0("RNAseq dashboard (",candidates$proband[rd$i()], ' ~~~ ',candidates$geneID[rd$i()],')')
    )
  })

  ### Image
  output$mainimage = renderImage({
    list(
      src = file.path(params$datadir,'input/CHUSJ_CR_Bioinformatique_V2.png'), #path to the file
      contentType = "image/png",
      width = 400,
      alt = "My Figure"
    )},deleteFile = FALSE)

  ### Data table
  output$exonTPM = renderDT({
    if(input$table_choice == 'Normalised expression') {
      datatable(
        rd$reactive_inputs()$fc_exons_tpm_reactive,
        rownames = FALSE,
        options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
    } else if (input$table_choice == "Raw counts") {
      datatable(
        rd$reactive_inputs()$fc_exons_raw_reactive,
        rownames = FALSE,
        options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
    } else if (input$table_choice == "Isoform-specific expression") {
      datatable(
        rd$reactive_inputs()$transcripts_reactive,
        rownames = FALSE,
        options = list(pageLength = 100,lengthChange = FALSE,info = FALSE))
    }
  })

  ### Data table
  output$Table = renderDT({
    datatable(
      rd$reactive_inputs()$transcripts_reactive,
      rownames = FALSE,
      options = list(pageLength = 50)
    )
  })
  ### Plotly expression
  output$Expression = renderPlotly({
    plot_expression_cohort(
      data_child   = rd$reactive_inputs()$fc_exons_ggplot_reactive_child,
      data_adults  = rd$reactive_inputs()$fc_exons_ggplot_reactive_adults,
      data_patient = rd$reactive_inputs()$fc_exons_ggplot_reactive_patient
    )
  })

  ### Plotly expression per family (only when n_members > 1)
  output$Expression_perfamily_ui = renderUI({
    data_family = rd$reactive_inputs()$fc_exons_ggplot_reactive_family
    n_members = length(unique(data_family$PatientID))
    req(n_members > 1)
    card(plotlyOutput("Expression_perfamily", width = "1500px"))
  })

  ### OUTRIDER table
  output$table_OUTRIDER = renderDT({
    datatable(
      rd$reactive_inputs()$table_OUTRIDER[rd$reactive_inputs()$table_OUTRIDER$pValue < as.numeric(input$pvalue),],
      rownames = FALSE,
      options = list(pageLength = 100)
    )
  })

  ### OUTRIDER candidate gene table
  output$candidates_OUTRIDER = renderDT({
    datatable(
      rd$reactive_inputs()$table_OUTRIDER_candidate,
      rownames = FALSE,
      options = list(dom = 'tir',info = FALSE)
    )
  })

  ### OUTRIDER candidate exon table
  output$candidates_OUTRIDER_exons = renderDT({
    datatable(
      rd$reactive_inputs()$table_perexons_OUTRIDER_candidate,
      rownames = FALSE,
      options = list(pageLength = 50)
    )
  })

  ### OUTRIDER significant exon table
  output$significant_OUTRIDER_exons = renderDT({
    datatable(
      rd$reactive_inputs()$table_perexons_OUTRIDER_significant,
      rownames = FALSE,
      options = list(pageLength = 50)
    )
  })

  ### Sashimi plots FRASER
  output$Sashimi = renderImage({
    list(
      src = paste0(params$datadir,"/sashimis/gene_",candidates$geneID[rd$i()],"_",candidates$proband[rd$i()],'_sashimi.png')[1], #path to the file
      contentType = "image/png",
      width = 900,
      alt = "My Figure"
    )}, deleteFile = FALSE)

  ### Reactive slider
  output$genemodel_slider = renderUI({
    cmin = floor(candidates$start[rd$i()]/1000)
    cmax = ceiling(candidates$stop[rd$i()]/1000)
    sliderInput(
      inputId = "sliderxlims",
      label = "Select genomic window (Kb)",
      min = cmin,
      max = cmax,
      value = c(cmin, cmax),
      width = '100%'
    )
  })

  ### Coverage ggplots
  output$Figure_genemodel = renderPlot({
    req(candidates$geneID[rd$i()]!="")
    if(is.null(input$sliderxlims)) {genemodel = ggplot() +  theme_void() + geom_text(aes(0,0,label='Plotting in ¨Progress')) + xlab(NULL)} else {
      gene_dir = paste0(params$datadir,'/bams_subset/gene',candidates$geneID[rd$i()],'_chr',candidates$chromosome[rd$i()],'_',candidates$start[rd$i()]-5000,'_',candidates$stop[rd$i()]+5000,'/')
      genemodel = genemodel_plot(
        candidate = candidates[rd$i(),],
        depth_file = paste0(gene_dir,"gene_",candidates$geneID[rd$i()],"_",candidates$proband[rd$i()],"_depth5.csv"),
        res_dt_candidate_gene_file = paste0(gene_dir,"gene_",candidates$geneID[rd$i()],"_",candidates$proband[rd$i()],"_res_dt_candidate_gene.csv"),
        bam_file = paste0(gene_dir,candidates$proband[rd$i()],"_sorted_chrN.bam"),
        colmean_genes_counts_file = paste0(params$datadir,'/colmean_genes_counts.tsv'),
        gene_annotations=gene_annotations,
        xlims = input$sliderxlims)}

    genemodel
  })

  ### gene prioritization Table
  output$gp = renderDT({
    datatable(rd$gene_prioritization_data(),
              rownames = TRUE,options = list(pageLength = 100,columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })

  ### DOWNLOAD gene prioritization Table
  output$gp_download = downloadHandler(
    filename = function() {paste0("gene_prioritization_", candidates$geneID[rd$i()],'_',candidates$proband[rd$i()], ".csv")},
    content = function(file) {write.csv(rd$gene_prioritization_data(),file,row.names = FALSE)}
  )

  ### genome-wide OUTRIDER
  output$gwOUTRIDER = renderPlot({
    manhattan_plot(res_dt=gwOUTRIDER,sample = candidates$proband[rd$i()],geneID = 'geneID',pvalue = 'pValue',pcutoff = 0.01,shape = TRUE)
  })

  ### genome-wide FRASER
  output$gwFRASER = renderPlot({
    manhattan_plot(res_dt=gwFRASER,sample = candidates$proband[rd$i()])
  })

  output$gwFRASER_table = renderDT({
    datatable(
      gwFRASER_table(res_dt=gwFRASER,sample = candidates$proband[rd$i()]),
      rownames = FALSE,options = list(pageLength = 100,columnDefs = list(list(className = 'dt-center', targets = "_all"))))
  })

  ### IGV
  observeEvent(input$addBamLocalFileButton, {
    gene_dir = paste0(params$datadir,'/bams_subset/gene',candidates$geneID[rd$i()],'_chr',candidates$chromosome[rd$i()],'_',candidates$start[rd$i()]-5000,'_',candidates$stop[rd$i()]+5000,'/')
    color=
      showGenomicRegion(session, id="igvShiny",paste0("chr",candidates$chromosome[rd$i()],":",candidates$start[rd$i()],"-",candidates$stop[rd$i()]))
    bamFile = paste0(gene_dir,candidates$proband[rd$i()],"_sorted_chrN.bam")
    bamAlign = readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq"))
    loadBamTrackFromLocalData(session, id="igvShiny", trackName=input$proband, data=bamAlign)
    runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'green';")
  })

  ### Default hg38 view
  output$igvShiny = renderIgvShiny({
    runjs("document.getElementById('addBamLocalFileButton').style.backgroundColor = 'red';")
    genomeOptions = parseAndValidateGenomeSpec('hg38',initialLocus=paste0("chr",candidates$chromosome[rd$i()],":",candidates$start[rd$i()],"-",candidates$stop[rd$i()]))
    igvShiny(genomeOptions)
  })

  ### genome-wide ASE table
  output$gwASE_table = renderDT({
    req(!is.null(gwASE))
    datatable(
      rd$reactive_inputs()$gwASE_table,
      rownames = FALSE, options = list(pageLength = 100,lengthChange = FALSE,info = FALSE,columnDefs = list(list(className = 'dt-left', targets = "_all"))))
  })

  ### genome-wide ASE Imprinted table
  output$gwImprinted_table = renderDT({
    req(!is.null(gwASE_IMX))
    datatable(
      rd$reactive_inputs()$gwIMX_table[rd$reactive_inputs()$gwIMX_table$Type == 'I',1:12],
      rownames = FALSE, options = list(pageLength = 100,lengthChange = FALSE,info = FALSE,columnDefs = list(list(className = 'dt-left', targets = "_all"))))
  })

  ### genome-wide ASE X table
  output$gwX_table = renderDT({
    req(!is.null(gwASE_IMX))
    datatable(
      rd$reactive_inputs()$gwIMX_table[rd$reactive_inputs()$gwIMX_table$Type == 'X',1:12],
      rownames = FALSE, options = list(pageLength = 100,lengthChange = FALSE,info = FALSE,columnDefs = list(list(className = 'dt-left', targets = "_all"))))
  })

  ### genome-wide ASE manhattan
  output$gwASE = renderPlot({
    req(!is.null(gwASE))
    manhattan_plot(res_dt=gwASE,sample = candidates$proband[rd$i()],end= 'pos',pcutoff=0.01, pvalue='pvalue',geneID = 'geneID')
  })

  ### FASTA
  output$fasta = renderUI({
    fasta.file = paste0(params$datadir,"/consensus/","gene",candidates$geneID[rd$i()],'_',candidates$proband[rd$i()],'.fasta')
    if(file.exists(fasta.file)) {lines = readLines(fasta.file)} else {lines = 'No gene specified'}
    HTML(paste0("<span style='color: black; font-family: Courier New; font-size: 16px;'>",paste0(lines,collapse = '<br>'),"</span>",collapse = "\n"))
  })

  ### Description
  output$description = renderUI({
    req(candidates$proband[rd$i()]!="")
    selected_ensembl = candidates$ensembl[rd$i()]
    selected_geneID = candidates$geneID[rd$i()]
    selected_patient = candidates$proband[rd$i()]
    selected_clinical = clinical[clinical$`Patient ID` == selected_patient,]
    selected_bam = paste0(params$datadir,'/bams_subset/gene',selected_geneID,'_chr',candidates$chromosome[rd$i()],'_',candidates$start[rd$i()]-5000,'_',candidates$stop[rd$i()]+5000,'/',selected_patient,"_sorted_chrN.bam")
    if(selected_geneID == "") selected_bam = ''
    log_info(paste0('Selecting ~~~ ',selected_patient,' ~~~ ',selected_geneID,' ~~~ ',rd$i()))
    url = paste0("https://www.proteinatlas.org/", selected_ensembl)
    HTML(
      paste0("<span><b>Notes: </b>",selected_clinical$Notes,
             "<br><br><b>Gene description for ",selected_geneID,": </b> <a href='", url, "' target='_self'>",selected_ensembl,"</a>",
             "<br><br><b>Mutations: </b>",selected_clinical$Mutation,
             "<br><br><b>Candidate Gene hypothesis: </b>",selected_clinical$Hypothèse,
             "<br><br><b>HPO terms: </b>",selected_clinical$`HPO terms`,
             "<br><br><b>Sex: </b>",selected_clinical$Sexe,
             "<br><br><b>Gene .bam location: </b>",selected_bam,
             "<br><br><b>Complete .bam location (Fir): </b>",system('echo ${HOME}',intern = TRUE),"/project/def-rallard/COMMUN/raredisease_rnaseq/results_nextflow_rnasplice_09_05_2026/star_salmon/",selected_patient,".sorted.bam<span>")
    )
  })

  ### Figure legends
  output$Figure_genemodel_legend = renderUI({
    HTML(
      paste0("<span><b>Figure 2:</b> Visualisation des altérations d’épissage détectées par l'outil FRASER.
        <br><b>A:</b> Carte des introns/exons du gène ",candidates$geneID[rd$i()]," et localisation du/des variant.s (ligne bleue pointillée).
        <br><b>B:</b> Évènements d’épissage aberrant pour chaque région détectée (−log₁₀ p-value, min. p-value: pvaleur la plus basse dans la région couverte, max. deltaPSI: proportion d'épissage observé moins attendu).
        <br><b>C:</b>  Couverture de séquençage normalisée pour le probant ainsi que 25-75ième percentile de la population de référence (en orange).</span>")
    )
  })

  output$ase_legend = renderUI({
    HTML(
      paste0("<span>Identification of Allele Specific Expression based on SNV genotypes (.vcf) and short read RNAseq alignement (.bam) files.
                <br>Called with GATK - ASEReadCounter. Significance tested with binomial t-tests. Visualised as Manhattan plots and dynamic table.
                <br><b>ref/alt: </b>RNAseq reference/alternate count.
                <br><b>chr: </b>Chromosome.
		<br><b>pos: </b>Position.
		<br><b>RNA_DP: </b>RNA read depth for reference,alternate allele.
		<br><b>RNA_ratio: </b>RNA alternate allele/ (ref+alt).
		<br><b>WGS_GT: </b>WGS genotype.
		<br><b>WGS_DP: </b>WGS read depth for reference,alternate allele.
		<br><b>WGS_GQ: </b>WGS genotype quality.
                <b>WGSratio: </b>WGS alternate allele/ (ref+alt).</span>")
    )
  })

  ### Versioning
  output$Version = renderDT({
    datatable(
      data.frame(Parameter=names(unlist(report_version)),Value=unlist(report_version)),
      rownames = FALSE,options = list(dom = 'p'))
  })

  ### Candidate genes table
  output$candidates_table = renderDT({
    tbl = candidates %>%
      group_by(proband) %>%
      summarise(geneID = paste(geneID, collapse = ', '), .groups = 'drop') %>%
      arrange(proband)
    datatable(
      tbl,
      rownames = FALSE,
      options = list(pageLength = 6,
                     columnDefs = list(list(className = 'dt-left', targets = '_all'))))
  })

  # Search gene expression — populate choices server-side to avoid sending 20k options to browser
  updateSelectizeInput(session, "gene_search",
                       choices = sort(unique(fc_genes_raw_ALL$geneID)),
                       selected = 'LDHA',
                       server = TRUE)

  output$searchExpression = renderPlotly({
    req(input$gene_search)
    plot_hb_fraction(fc_genes_raw_ALL, hb_genes = input$gene_search)
  })


  ### Haemoglobin barplot
  output$hb_barplot = renderPlotly({
    plot_hb_fraction(fc_genes_raw_ALL)
  })

  ### Total reads barplot
  output$total_reads_barplot = renderPlotly({
    plot_total_reads(fc_genes_raw_ALL)
  })

  ### multiQC
  output$htmlViewer = renderUI({
    HTML(paste(readLines(html_files[as.numeric(input$multiQCs)]), collapse = "\n"))
  })
}
  
###### shiny app
logger::log_info('Deploying app')
app = shinyApp(ui, server, options = list(height = 900))
runApp(app, launch.browser = TRUE)

