###############################################################################
# Shiny module: reactive data for RNAseq dashboard
#
# UI:    reactive_data_UI(id)
# Server: reactive_data_server(id, proband, pvalue, geneprior_rm)
#
# Returns a list of reactive expressions:
#   $geneID              – character vector of gene IDs for the selected proband
#   $i                   – integer index into `candidates` for the selection
#   $reactive_inputs     – filtered / reshaped data list consumed by outputs
#   $gene_prioritization_data – gene prioritization table
###############################################################################

reactive_data_UI = function(id) {
  ns = NS(id)
  uiOutput(ns("geneID_list"))
}

reactive_data_server = function(id, proband, pvalue, geneprior_rm) {

  moduleServer(id, function(input, output, session) {

    ### Gene IDs available for the selected proband
    geneID = reactive({
      candidates$geneID[candidates$proband == proband()]
    })

    ### Dynamic gene-selection dropdown (namespaced inside module)
    output$geneID_list = renderUI({
      selectInput(
        inputId  = session$ns("gene_selection"),
        label    = h4("Select a gene:"),
        selected = 'MCM5',
        choices  = geneID()
      )
    })

    ### Row index in `candidates` for the current proband + gene combination
    i = reactive({
      which(candidates$proband == proband() & candidates$geneID == input$gene_selection)
    })

    ### All filtered / reshaped datasets consumed by the outputs
    reactive_inputs = reactive({

      reactive_i = gsub('_PAX', '', proband())

      # Isoform (transcript-level) expression
      column = which(colnames(transcripts_named_filtered) == reactive_i)
      column_not = which(colnames(transcripts_named_filtered) != reactive_i)
      column_not = column_not[-c(1, 2, length(column_not))]
      transcripts_reactive = transcripts_named_filtered[
        transcripts_named_filtered$proband == reactive_i,
        c(1, 2, column, column_not)
      ]
      transcripts_reactive = transcripts_reactive[
        order(transcripts_reactive[, 3], decreasing = TRUE), ]

      # Exon raw counts
      column = which(colnames(fc_exons_raw) == reactive_i)
      column_not = which(colnames(fc_exons_raw) != reactive_i)
      column_not = column_not[-c(1, 2, 3, 4, 5, length(column_not))]
      fc_exons_raw_reactive = fc_exons_raw[
        fc_exons_raw$geneID == candidates$geneID[i()],
        c(1, 2, 3, 4, 5, column, column_not)
      ]
      fc_exons_raw_reactive = fc_exons_raw_reactive[
        order(fc_exons_raw_reactive[, 4]), ]

      # Exon TPM counts
      fc_exons_tpm_reactive = fc_exons_tpm[
        fc_exons_tpm$geneID == candidates$geneID[i()],
        c(1, 2, 3, 4, 5, column, column_not)
      ]
      fc_exons_tpm_reactive = fc_exons_tpm_reactive[
        order(fc_exons_tpm_reactive[, 4]), ]

      # OUTRIDER – genome-wide table for this proband
      table_OUTRIDER = gwOUTRIDER[gwOUTRIDER$sampleID == candidates$proband[i()], ]
      table_OUTRIDER = table_OUTRIDER[order(table_OUTRIDER$chr, table_OUTRIDER$pos), ]

      # OUTRIDER – candidate gene
      table_OUTRIDER_candidate = candidates_OUTRIDER[
        candidates_OUTRIDER$sampleID == candidates$proband[i()], ]
      table_OUTRIDER_candidate = table_OUTRIDER_candidate[
        table_OUTRIDER_candidate$geneID == candidates$geneID[i()], ]

      # OUTRIDER – candidate exons
      table_perexons_OUTRIDER_candidate = candidates_perexons_OUTRIDER[
        candidates_perexons_OUTRIDER$sampleID == candidates$proband[i()], ]
      table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[
        table_perexons_OUTRIDER_candidate$geneID == candidates$geneID[i()], ]
      table_perexons_OUTRIDER_candidate = table_perexons_OUTRIDER_candidate[
        order(table_perexons_OUTRIDER_candidate$exonID), ]

      # OUTRIDER – significant exons
      table_perexons_OUTRIDER_significant = significant_perexons_OUTRIDER[
        significant_perexons_OUTRIDER$sampleID == candidates$proband[i()], ]
      table_perexons_OUTRIDER_significant = table_perexons_OUTRIDER_significant[
        order(table_perexons_OUTRIDER_significant$chr,
              table_perexons_OUTRIDER_significant$pos,
              table_perexons_OUTRIDER_significant$exonID), ]

      # ASE – genome-wide (p < 0.01)
      gwASE_table = gwASE[gwASE$sampleID == candidates$proband[i()], ]
      gwASE_table = gwASE_table[gwASE_table$pvalue < 0.01, ]
      gwASE_table = gwASE_table[
        order(gwASE_table$chr, gwASE_table$pos, gwASE_table$geneID), ]

      # ASE – imprinted & X-linked
      gwIMX_table = gwASE_IMX[gwASE_IMX$sampleID == candidates$proband[i()], ]
      gwIMX_table = gwIMX_table[
        order(gwIMX_table$chr, gwIMX_table$pos, gwIMX_table$geneID), ]

      # Plotly – family and per-age subsets
      fc_exons_ggplot_reactive = fc_exons_tpm_ggplot[
        fc_exons_tpm_ggplot$proband == reactive_i, ]
      fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[
        fc_exons_ggplot_reactive$geneID == candidates$geneID[i()], ]

      fc_exons_ggplot_reactive_family = fc_exons_ggplot_reactive[
        gsub('_0[123]$', '', fc_exons_ggplot_reactive$PatientID) ==
          gsub('_03$', '', reactive_i), ]

      fc_exons_ggplot_reactive_patient = fc_exons_ggplot_reactive[
        fc_exons_ggplot_reactive$PatientID == reactive_i, ]
      fc_exons_ggplot_reactive = fc_exons_ggplot_reactive[
        fc_exons_ggplot_reactive$PatientID != reactive_i, ]

      fc_exons_ggplot_reactive_child = fc_exons_ggplot_reactive[
        fc_exons_ggplot_reactive$age <  18, ]
      fc_exons_ggplot_reactive_adults = fc_exons_ggplot_reactive[
        fc_exons_ggplot_reactive$age >= 18, ]

      # Stagger x-axis positions slightly for readability
      fc_exons_ggplot_reactive_child$exonID = fc_exons_ggplot_reactive_child$exonID  - 0.1
      fc_exons_ggplot_reactive_adults$exonID = fc_exons_ggplot_reactive_adults$exonID + 0.1

      list(
        transcripts_reactive                = transcripts_reactive,
        table_OUTRIDER                      = table_OUTRIDER,
        table_OUTRIDER_candidate            = table_OUTRIDER_candidate,
        fc_exons_ggplot_reactive_child      = fc_exons_ggplot_reactive_child,
        fc_exons_ggplot_reactive_adults     = fc_exons_ggplot_reactive_adults,
        fc_exons_tpm_reactive               = fc_exons_tpm_reactive,
        fc_exons_raw_reactive               = fc_exons_raw_reactive,
        fc_exons_ggplot_reactive            = fc_exons_ggplot_reactive,
        fc_exons_ggplot_reactive_patient    = fc_exons_ggplot_reactive_patient,
        fc_exons_ggplot_reactive_family     = fc_exons_ggplot_reactive_family,
        table_perexons_OUTRIDER_candidate   = table_perexons_OUTRIDER_candidate,
        table_perexons_OUTRIDER_significant = table_perexons_OUTRIDER_significant,
        gwASE_table                         = gwASE_table,
        gwIMX_table                         = gwIMX_table
      )
    })

    ### Gene prioritization table
    gene_prioritization_data = reactive({
      gene_prioritization(
        sample      = candidates$proband[i()],
        top         = 100,
        hpo_sample  = clinical,
        hpo_all     = 'genes_to_phenotype.txt',
        fraser      = gwFRASER,
        outrider    = gwOUTRIDER,
        geneprior_rm = geneprior_rm()
      )
    })

    ### Return all reactives for use by the main server
    list(
      geneID                 = geneID,
      i                      = i,
      reactive_inputs        = reactive_inputs,
      gene_prioritization_data = gene_prioritization_data
    )
  })
}
