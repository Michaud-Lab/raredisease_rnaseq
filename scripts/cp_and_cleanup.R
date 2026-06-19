# =============================================================================
# cp_and_cleanup.R - Consolidate all pipeline outputs into the data/ directory
#                    and package for transfer
# =============================================================================

# Libraries
library(readxl)
library(tidyr)
library(data.table)
library(dplyr)
library(zip)
library(jsonlite)

# -----------------------------------------------------------------------------
# 1. Parameters
# -----------------------------------------------------------------------------
params = list(workdir = "/project/def-rallard/COMMUN/raredisease_rnaseq/")
params$datadir = file.path(params$workdir, '/data/')
params$resultdir = file.path(params$workdir, '/rnasplice/')

dir.create(params$datadir, showWarnings = FALSE)
dir.create(paste0(params$datadir, 'bams_subset'), showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Copy ASE results
# -----------------------------------------------------------------------------
cp_ASE = paste0('cp ', params$workdir, '/ASE/gw* ', params$datadir, '/.')
system(cp_ASE)

# -----------------------------------------------------------------------------
# 3. Compile genome-wide FRASER results
# -----------------------------------------------------------------------------
chrs = list.files(file.path(params$workdir, 'FRASER/bams_chr_subset/'), full.names = TRUE)

for (c in 1:length(chrs)) {
  if (file.exists(file.path(chrs[c], 'res_dt.csv'))) {
    gwFRASER_temp = read.csv(file.path(chrs[c], 'res_dt.csv'), header = TRUE, row.names = 1)
    gwFRASER_temp = gwFRASER_temp[gwFRASER_temp$pValue < 0.001, ]
    gwFRASER_temp$sampleID = gsub(paste0("_", gwFRASER_temp$seqnames[1], "$"), "", gwFRASER_temp$sampleID)
    if (c == 1) gwFRASER = gwFRASER_temp
    if (c != 1) gwFRASER = rbind(gwFRASER, gwFRASER_temp)
  }
}

colnames(gwFRASER)[c(1, 20)] = c('chr', 'pos')
write.csv(gwFRASER, file.path(params$datadir, 'gwFRASER.csv'))

# -----------------------------------------------------------------------------
# 4. Copy scripts and project metadata
# -----------------------------------------------------------------------------
cp_scripts = paste0("rsync -av --exclude='deprecat*' --exclude='.*' --exclude='slurm*' ",
                    params$workdir, "/scripts ", params$datadir, "/.")
cp_info = paste0('cp ', paste0(params$workdir,
                               c('/VERSION.json ', '/CHANGELOG.md ', '/README.md ', '/RNAseq_shiny_v* '),
                               collapse = ''), params$datadir, '/.')
system(cp_scripts)
system(cp_info)

# -----------------------------------------------------------------------------
# 5. Copy featureCounts outputs
# -----------------------------------------------------------------------------
cp_data_fc = paste0('cp ', params$workdir, '/featureCounts/*tsv ', params$datadir, '/.')
cp_rda_fc = paste0('cp ', params$workdir, '/featureCounts/gene_annotations.rda ', params$datadir, '/.')
system(cp_data_fc)
system(cp_rda_fc)

# -----------------------------------------------------------------------------
# 6. Copy sashimi plots, consensus sequences and OUTRIDER
# -----------------------------------------------------------------------------
dir.create(file.path(params$datadir, 'sashimis'), showWarnings = FALSE)
system(paste0('cp ', params$workdir, '/FRASER/results/*/*_sashimi.png ', params$datadir, '/sashimis/.'))
system(paste0('cp -r ', params$workdir, '/consensus ', params$datadir, '/.'))
system(paste0('cp -r ', params$workdir, '/OUTRIDER/*OUTRIDER.tsv ', params$datadir, '/.'))

# -----------------------------------------------------------------------------
# 7. Copy per-candidate BAM subsets, depth files, and FRASER results
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(params$datadir, 'input/candidate_genes.csv'))
candidates_LR = read.csv(file.path(params$datadir, 'input/candidate_genes_LR.csv'))
candidates_LR = candidates_LR[, c(2, 3, 10, 4, 5, 6, 7, 8, 9)]
colnames(candidates_LR) = colnames(candidates)
candidates = rbind(candidates, candidates_LR)
configs = read_json(file.path(params$datadir, 'input/configs.json'))
candidates_extra = read.table(configs$general$candidate_genes_extra,comment.char = "#",header = T ,sep = ',');candidates_extra[is.na(candidates_extra)] = ''
candidates = rbind(candidates, candidates_extra)

for (i in 1:nrow(candidates)) {
  gene_dir = paste0('bams_subset/gene', candidates$geneID[i],
                    '_chr', candidates$chromosome[i],
                    '_', candidates$start[i] - 5000,
                    '_', candidates$stop[i] + 5000, '/')
  in_dir = paste0(params$workdir, '/FRASER/', gene_dir)
  out_dir = paste0(params$datadir, gene_dir)

  dir.create(out_dir, showWarnings = TRUE)

  if (length(Sys.glob(paste0(in_dir, candidates$proband[i], '_sorted_chrN.bam*'))) > 0)
    system(paste0('cp ', in_dir, candidates$proband[i], '_sorted_chrN.bam* ', out_dir, '/.'))

  if (length(Sys.glob(paste0(in_dir, '*depth5.csv'))) > 0)
    system(paste0('cp ', in_dir, '*depth5.csv ', out_dir, '/.'))

  if (length(Sys.glob(paste0(in_dir, '*_res_dt_candidate_gene.csv'))) > 0)
    system(paste0('cp ', in_dir, '*_res_dt_candidate_gene.csv ', out_dir, '/.'))
}

# -----------------------------------------------------------------------------
# 8. Fix and copy MultiQC HTML report (remove malformed JS line)
# -----------------------------------------------------------------------------
html_dirs = list.dirs(params$resultdir,recursive = F)
html_files = paste0(html_dirs, '/multiqc/multiqc_report.html')

for(i in 1:length(html_files)){
  grep_cmd = paste0(
    "grep 'this.renderTo.parentNode.insertBefore(this.dataTableDiv,this.renderTo.nextSibling)),",
    "this.dataTableDiv.innerHTML=this.getTable()},a.getOptions().exporting&&",
    "a.getOptions().exporting.buttons.contextButton.menuItems.push({textKey:' ",
    html_files[i], " -n | cut -d: -f1 >grep_problematic_line"
  )
  system(grep_cmd)
  line = read.table('grep_problematic_line')
  system(paste0("sed '", line, "d' ", html_files[i], " >temp.html"))
  file.copy('temp.html', file.path(params$datadir, paste0('multiqc_report_',i,'.html')), overwrite = TRUE)
  
  system('rm temp.html grep_problematic_line')
}
# -----------------------------------------------------------------------------
# 9. Build transcript expression tables
# -----------------------------------------------------------------------------
ensembl_geneid = read.table(file.path(params$datadir, '/input/ensembl_geneid.tsv'), header = TRUE)
clinical = read.table(file.path(params$datadir, 'clinical.tsv'), check.names = FALSE)

transcripts = read.csv(
  file.path(params$resultdir, '/star_salmon/tximport/salmon.merged.transcript_tpm.tsv'),
  sep = '\t', check.names = FALSE
)
transcripts[, -c(1:2)] = round(transcripts[, -c(1:2)], 2)
transcripts_named = merge(transcripts, ensembl_geneid)
transcripts_named = data.frame(
  transcripts_named[, c(ncol(transcripts_named), 2)],
  transcripts_named[, -c(1, 2, ncol(transcripts_named))],
  check.names = FALSE
)
colnames(transcripts_named)[1:2] = c('geneID', 'isoform/transcript')

proband_ids = clinical$`Patient ID`[clinical$type == 'Proband']
transcripts_named_filtered = transcripts_named[transcripts_named$geneID %in% candidates$geneID, ]
transcripts_named_filtered = transcripts_named_filtered[,
  colnames(transcripts_named_filtered) %in% c('geneID', 'isoform/transcript', proband_ids)
]
transcripts_named_filtered = merge(transcripts_named_filtered, candidates[, c(1, 3)], sort = FALSE)
colnames(transcripts_named_filtered) = gsub('_PAX', '', colnames(transcripts_named_filtered))
transcripts_named_filtered$proband = gsub('_PAX', '', transcripts_named_filtered$proband)

transcripts_named_filtered_ggplot = transcripts_named_filtered %>%
  pivot_longer(cols = c(3:(ncol(transcripts_named_filtered) - 1)), names_to = 'PatientID', values_to = 'expression')
transcripts_named_filtered_ggplot = merge(
  transcripts_named_filtered_ggplot,
  clinical[, colnames(clinical) %in% c('PatientID', 'Sexe', 'type', 'age')]
)

write.table(transcripts_named_filtered,
            file.path(params$datadir, 'transcripts_named_filtered.tsv'), sep = '\t', quote = FALSE)
write.table(transcripts_named_filtered_ggplot,
            file.path(params$datadir, 'transcripts_named_filtered_ggplot.tsv'), sep = '\t', quote = FALSE)
write.table(clinical, file.path(params$datadir, 'clinical.tsv'), sep = '\t', quote = TRUE)

# -----------------------------------------------------------------------------
# 10. Zip everything for transfer
# -----------------------------------------------------------------------------
setwd(params$workdir)
zip(zipfile = paste0('tmp/data_', as.character(format(Sys.time(), format = "%Y_%m_%d_%H_%M")), '.zip'),
    files = 'data')

print(paste0('All done, Time is: ', Sys.time()))
