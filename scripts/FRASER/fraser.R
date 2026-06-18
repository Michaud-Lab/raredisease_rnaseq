# =============================================================================
# fraser.R - Per-candidate-gene splicing analysis with FRASER
# =============================================================================

# Libraries (loaded inside the main block to avoid loading on already-done runs)
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(googlesheets4))); gs4_deauth()

# -----------------------------------------------------------------------------
# 1. Arguments and parameters
# -----------------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
params = list(workdir = args[2])
params$FRASER = file.path(params$workdir, 'FRASER')
params$rnasplice_bamdir = args[3]
params$fraser_temp_bamdir = args[4]
params$fraser_bamdir = args[5]
params$candidate_genes_extra = args[6]
options(scipen = 999)

# -----------------------------------------------------------------------------
# 2. Setup output directories
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(params$workdir, 'data/input/candidate_genes.csv'))
candidates_extra = read_sheet(params$candidate_genes_extra, skip = 1)
candidates = rbind(candidates, candidates_extra)
params$output = paste0(params$FRASER, '/results/output_FRASER_', candidates$proband[i])

dir.create(params$FRASER, showWarnings = FALSE)
dir.create(file.path(params$FRASER, 'results'), showWarnings = FALSE)
dir.create(params$output, showWarnings = FALSE)

chr = candidates$chromosome[i]
start = candidates$start[i] - 5000
if (start <= 0) start = 1
stop = candidates$stop[i] + 5000
geneID = candidates$geneID[i]
proband = candidates$proband[i]

# -----------------------------------------------------------------------------
# 3. Run FRASER
# -----------------------------------------------------------------------------
if (geneID != "") {
  out_dir = paste0(params$FRASER, '/bams_subset/gene', geneID, '_chr', chr, '_', start, '_', stop)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  command = paste('./fraser.sh', chr, start, stop, geneID, proband,
                  params$rnasplice_bamdir, params$fraser_temp_bamdir, params$fraser_bamdir)

  res_dt_outfile = paste0(out_dir, "/gene_", geneID, "_", proband, "_res_dt_candidate_gene.csv")

   if (!file.exists(res_dt_outfile) | (file.exists(res_dt_outfile) && length(readLines(res_dt_outfile)) == 0)) {
    # Create an empty placeholder in case the analysis fails, to avoid re-running
    system(paste0('touch ', res_dt_outfile))

    # Subset BAM to gene region
    system(command)

    suppressMessages(suppressWarnings(library(dplyr)))
    suppressMessages(suppressWarnings(library(patchwork)))
    suppressMessages(suppressWarnings(library(FRASER)))
    suppressMessages(suppressWarnings(library(tidyr)))

    sampleTable = data.table(data.frame(
      sampleID = gsub('_sorted_chrN.bam', '', list.files(out_dir, pattern = '*bam$')),
      bamFile = list.files(out_dir, pattern = '*bam$', full.names = FALSE),
      group = 1,
      pairedEnd = TRUE
    ))
    sampleTable$group = 1:nrow(sampleTable)
    sampleTable$bamFile = paste0(out_dir, '/', sampleTable$bamFile)

    # Probands only
    probands = which(
      grepl('_0[34]_', sampleTable$bamFile) |
      grepl('LC_', sampleTable$bamFile) |
      grepl('F0', sampleTable$bamFile)
    )
    settings = FraserDataSet(colData = sampleTable[probands, ], workingDir = params$output)

    fds = tryCatch({
      invisible(capture.output(
        fds <- suppressMessages(suppressWarnings(
          countRNAData(settings, recount = TRUE)
        )),
        file = nullfile()
      ))
      fds = calculatePSIValues(fds)
      fds = annotateRanges(fds, GRCh = 38)
      fds = fit(fds, q = c(jaccard = 2))
      fds = calculatePvalues(fds)
      fds = calculatePadjValues(fds, method = 'none', geneLevel = FALSE)
      fds = calculateZscore(fds)
      fds
    }, error = function(e) {
      print(paste0('FRASER failed for sample ~~~ ', i, ' ~~~ ', geneID,
                   ' ~~~ ','likely because there were no splice junctions ~~~ ', conditionMessage(e)))
      # Write an empty CSV (header only) so this sample is not re-run
      write.csv(data.frame(), res_dt_outfile)
      NULL
    })

    if (is.null(fds)) {
      print(paste0('Skipping sample ~~~ ', i, ' ~~~ Time is: ', Sys.time()))
      stop()
    }

    res = results(fds, all = TRUE, padjCutoff = NA, deltaPsiCutoff = NA)
    res_dt = as.data.table(res)
    res_dt = res_dt[sampleID == proband, ]

    print(paste0('Dimensions of res_dt: ', paste(dim(res_dt), collapse = ' x ')))

    # Annotation packages
    suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
    suppressMessages(suppressWarnings(library(org.Hs.eg.db)))

    txdb_chr = keepSeqlevels(
      TxDb.Hsapiens.UCSC.hg38.knownGene,
      paste0('chr', c(1:22, 'X', 'Y', 'M')),
      pruning.mode = "coarse"
    )

    temp = which(res_dt$hgncSymbol == geneID)
    candidate_gene_localisation = temp[!is.na(temp)][1]
    control_samples = sampleTable$sampleID[grepl('_03_', sampleTable$bamFile) | grepl('LC_', sampleTable$bamFile)]
    control_samples = control_samples[control_samples != proband][1:5]

    # Save FRASER results for this candidate gene
    res_dt_candidate_gene = res_dt[temp[!is.na(temp)], ]
    res_dt_candidate_gene$mean = (res_dt_candidate_gene$start + res_dt_candidate_gene$end) / 2
    res_dt_candidate_gene$minuslogpval = -log(res_dt_candidate_gene$pValue, 10)
    write.csv(res_dt_candidate_gene, res_dt_outfile)

    print(paste0('Dimensions of res_dt_candidate_gene: ', paste(dim(res_dt_candidate_gene), collapse = ' x ')))

    # Sashimi plot (fail-safe: initialise with blank plot in case rendering fails)
    sashimi_file = file.path(params$output, paste0('gene_', geneID, '_', proband, '_sashimi.png'))
    png(sashimi_file, width = 1000, height = 600)
    plot(0, main = 'Placeholder - rendering failed')
    dev.off()

    png(sashimi_file, width = 1000 + length(temp[!is.na(temp)]) * 5, height = 1200)
    try(plotBamCoverageFromResultTable(
      fds,
      result = res_dt[candidate_gene_localisation, ],
      show_full_gene = TRUE,
      txdb = txdb_chr,
      orgDb = org.Hs.eg.db,
      control_samples = control_samples,
      splicegraph_labels = "id",
      mar = c(1, 10, 0.1, 1)
    ), silent = TRUE)
    dev.off()

    print(paste0('Done sample ~~~ ', i, ', gene: ',geneID,' ~~~ Time is: ', Sys.time()))

  } else {
    print(paste0('Sample ~~~ ', i, ', gene: ',geneID,' already exists ~~~ Time is: ', Sys.time()))
  }
} else {
  print(paste0('Sample ~~~ ', i, ' did not contain a candidate gene ~~~ Time is: ', Sys.time()))
}
