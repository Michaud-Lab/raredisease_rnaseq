# =============================================================================
# cp_and_cleanup.R - Consolidate all pipeline outputs into the data/ directory
#                    and package for transfer
# =============================================================================
source("rnaseq_helper_functions.R")
load_install_library(c('readxl','tidyr','data.table','dplyr','zip','jsonlite'))

# -----------------------------------------------------------------------------
# 1. Parameters
# -----------------------------------------------------------------------------
params = list(workdir = file.path(getwd(),',,'))
params$datadir = file.path(params$workdir, 'data/')
params$resultdir = file.path(params$workdir, 'rnasplice_results')

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
write.table(gwFRASER, file.path(params$datadir, 'gwFRASER.tsv'),sep = '\t')

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
# 7. Annotate candidates with FRASER/OUTRIDER/ASE/HPO results, then copy
#    per-candidate BAM subsets, depth files, and FRASER results
# -----------------------------------------------------------------------------
candidates = read.csv(file.path(params$datadir, 'input/candidate_genes_ALL.csv'), check.names = FALSE)
configs = read_json(file.path(params$datadir, 'input/configs.json'))

gwfiles = paste0(params$datadir, c('/gwFRASER.tsv', '/gw_genes_OUTRIDER.tsv', '/gwASE.tsv', '/clinical.tsv'))
candidatefiles = paste0(params$datadir, c('/gwFRASER.tsv', '/candidates_OUTRIDER.tsv', '/gwASE.tsv'))

candidates = candidate_genes_gw_annotations(
  candidates, gwfiles = gwfiles, candidatefiles = candidatefiles,
  datadir = file.path(params$workdir, 'FRASER'),
  hpo_all = file.path(params$workdir, 'tmp/genes_to_phenotype.txt')
)
write.csv(candidates, file.path(params$datadir, 'input/candidate_genes_ALL.csv'), quote = TRUE, row.names = FALSE)

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
html_files = paste0(list.dirs(params$resultdir,recursive = F), '/multiqc/multiqc_report.html')

for(i in 1:length(html_files)){
  lines = readLines(html_files[i])
  problematic_line = grep("this.renderTo.parentNode.insertBefore", lines, fixed = TRUE)
  lines = lines[-problematic_line]
  writeLines(lines, file.path(params$datadir, paste0('multiqc_report_', i, '.html')))
}

# -----------------------------------------------------------------------------
# 10. Zip everything for transfer
# -----------------------------------------------------------------------------
setwd(params$workdir)
zip(zipfile = paste0('tmp/data_', as.character(format(Sys.time(), format = "%Y_%m_%d_%H_%M")), '.zip'),
    files = 'data')

print(paste0('All done, Time is: ', Sys.time()))
