  #directory were the files are located
  params = list(workdir = "/home/renaut/scratch/raredisease_rnaseq")
  params$datadir = file.path(params$workdir,'/data/')
  params$resultdir = file.path(params$workdir,'/results_07_10_2025/')

  #r, libraries
  library(readxl)
  library(tidyr)
  library(data.table)
  library(dplyr)
  library(zip)

  #cp the data dir and create a new one with current date
  dir.create(params$datadir,showWarnings = FALSE)
  dir.create(paste0(params$datadir,'gene_statistics'),showWarnings = FALSE)  
  dir.create(paste0(params$datadir,'bams_subset'),showWarnings = FALSE)

  #cp the scripts for backup
  cp_scripts = paste0('cp -r ', params$workdir,'/scripts ',params$datadir,'/.')
  cp_info = paste0('cp ',paste0(params$workdir,c('/VERSION.json ','/CHANGELOG.md ', '/README.md ', '/RNAseq_shiny_v* ')),params$datadir,'/.')
  system(cp_scripts);system(cp_info)

  #featureCounts (per gene and per exons)
  cp_data_fc = paste0('cp ', params$workdir,'/featureCounts/*tsv ',params$datadir,'/.')  
  cp_rda_fc = paste0('cp ', params$workdir,'/featureCounts/gene_annotations.rda ',params$datadir,'/.')

  system(cp_data_fc) 
  system(cp_rda_fc)

  #sashimis
  dir.create(file.path(params$datadir,'sashimis'), showWarnings = FALSE)
  cp_sashimi_cmd = paste0('cp ', params$workdir,'/FRASER/results/*/*_sashimi.png ', params$datadir,'/sashimis/.')
  system(cp_sashimi_cmd)

  #consensus sequences
  cp_consensus = paste0('cp -r ', params$workdir,'/consensus ', params$datadir,'/.') 
  system(cp_consensus)

  #outrider
  cp_outrider = paste0('cp -r ', params$workdir,'/OUTRIDER/*OUTRIDER.tsv ', params$datadir,'/.')
  system(cp_outrider)

  #candidate genes
  candidates = read.csv(file.path(params$datadir,'candidate_genes_3.txt'))

  #cp the bam subset for data viz
  for(i in 1:nrow(candidates)){
       dir = paste0(params$workdir, '/FRASER/bams_subset/gene',candidates$geneID[i],'_chr',candidates$chromosome[i],'_',candidates$start[i]-100000,'_',candidates$stop[i]+100000,'/')
       cp_bamsubset = paste0('cp ',dir,candidates$proband[i],'_sorted_chrN.bam* ',params$datadir,'/bams_subset/.')
       cp_genedepth = paste0('cp ',dir,'*depth5.csv ',params$datadir,'/gene_statistics/.')
       cp_res_dt = paste0('cp ',dir,'*_res_dt_candidate_gene.csv ',params$datadir,'/gene_statistics/.')

       system(cp_bamsubset)
       system(cp_genedepth)
       system(cp_res_dt)
        }


  #html multiqc file (remove problematic line)
  html_file = file.path(params$resultdir,'/multiqc/multiqc_report.html')
  grep_cmd = paste0("grep 'this.renderTo.parentNode.insertBefore(this.dataTableDiv,this.renderTo.nextSibling)),this.dataTableDiv.innerHTML=this.getTable()},a.getOptions().exporting&&a.getOptions().exporting.buttons.contextButton.menuItems.push({textKey:' ",html_file," -n | cut -d: -f1 >grep_problematic_line")
  system(grep_cmd)
  line = read.table('grep_problematic_line')
  cleanup_cmd = paste0("sed '",line,"d' ",html_file," >temp.html")
  system(cleanup_cmd)
  file.copy('temp.html',file.path(params$datadir,'multiqc_report.html'),overwrite = T)  
  system('rm temp.html grep_problematic_line')

  #Ensembl - GeneID correspondance file
  ensembl_geneid = read.table(file.path(params$datadir,'ensembl_geneid.tsv'),header = T)

  #Transcripts expression
  transcripts = read.csv(file.path(params$resultdir,'/star_salmon/tximport/salmon.merged.transcript_tpm.tsv'),sep = '\t',check.names = F)
  transcripts[,-c(1:2)] = round(transcripts[,-c(1:2)],2) #prettify
  transcripts_named = merge(transcripts,ensembl_geneid) #get the gene IDs
  transcripts_named = data.frame(transcripts_named[,c(ncol(transcripts_named),2)],transcripts_named[,-c(1,2,ncol(transcripts_named))], check.names = F) #reorder
  colnames(transcripts_named)[1:2] = c('geneID','isoform/transcript') #rename

  #Clinical data (from Maude)
  clinical = readxl::read_xlsx(file.path(params$workdir,'data/CHUSJ\ Master\ Linking\ Log_modif.xlsx'), sheet = 'Suivi - RNAseq',skip = 1)
  clinical$type =  "Parent"; clinical$type[!is.na(clinical$Mutation)] = 'Proband' #Parent versus Proband
  clinical = clinical[order(clinical$`Patient ID`),] #Same order as the transript expression data.
  clinical$age = as.numeric(clinical$`Âge (années)`); clinical$age[clinical$`Âge (années)` == '0 (3 mois)'] = 0.25 # Prettify
  clinical$PatientID = gsub('_PAX','',clinical$`Patient ID`)

  #all.equal(clinical$`Patient ID`,colnames(transcripts_named_filtered)[-c(1:3)]) #sanity check
  #transcripts dataframe formatting for ggplot / plotly  datatable
  transcripts_named_filtered = transcripts_named[transcripts_named$geneID %in% candidates$geneID, ] #filter ONLY genes of interest
  transcripts_named_filtered = transcripts_named_filtered[,colnames(transcripts_named_filtered) %in% c('geneID','isoform/transcript','geneID',clinical$`Patient ID`[clinical$type =='Proband'])]
  transcripts_named_filtered = merge(transcripts_named_filtered,candidates[,c(1,3)],sort = F)
  colnames(transcripts_named_filtered) = gsub('_PAX','',colnames(transcripts_named_filtered))
  transcripts_named_filtered$proband = gsub('_PAX','',transcripts_named_filtered$proband)

  transcripts_named_filtered_ggplot = transcripts_named_filtered %>% pivot_longer(cols = c(3:(ncol(transcripts_named_filtered)-1)), names_to = 'PatientID',values_to = 'expression')
  transcripts_named_filtered_ggplot = merge(transcripts_named_filtered_ggplot,clinical[,colnames(clinical) %in% c('PatientID','Sexe','type','age')])
 
  #write information
  write.table(transcripts_named_filtered,file.path(params$datadir,'transcripts_named_filtered.tsv'),sep = '\t',quote = F)
  write.table(transcripts_named_filtered_ggplot,file.path(params$datadir,'transcripts_named_filtered_ggplot.tsv'), sep = '\t',quote = F)
  write.table(clinical,file.path(params$datadir,'clinical.tsv'), sep = '\t',quote = T)

  #zip everything into one zip for easier transfer
  setwd(params$workdir)
  zip(zipfile = paste0('data_',as.character(format(Sys.time(), format = "%Y_%m_%d_%H_%M")),'.zip'), files = 'data')
  #zip_cmd = paste0('zip -rj data_',as.character(format(Sys.time(), format = "%Y_%m_%d_%H_%M")),'.zip ',params$workdir,'/data')
  #system(zip_cmd)

  #print
  print(paste0('All done, Time is: ',Sys.time()))


