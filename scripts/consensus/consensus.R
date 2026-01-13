library(seqinr)
library(dplyr)

#Define parameters
args = commandArgs(trailingOnly=TRUE)
i = as.numeric(args[1])
params = list(workdir = args[2])
params$candidates = args[3]
params$consensus = file.path(params$workdir,'consensus/')
params$FRASER = file.path(params$workdir,'FRASER')
params$fc_exons = args[4]
params$ref_file = args[5]
params$ref_annot = args[6]

#candidates
candidates = read.csv(params$candidates)

#transcript analysed
fc_exons_tpm = read.table(params$fc_exons,sep = '\t',check.names = F) 
unique_transcript = unique(fc_exons_tpm[,1:3])
params$transcript = unique_transcript$transcriptID[unique_transcript$geneID == candidates$geneID[i]]

#Define parameters
chr = candidates$chromosome[i]
start = candidates$start[i] - 100000; if(start <= 0) start = 1
stop = candidates$stop[i] + 100000
params$geneID = candidates$geneID[i]
params$out_dir = paste0(params$FRASER,'/bams_subset/gene',candidates$geneID[i],'_chr',chr,'_',start,'_',stop,'/')
params$bam_file = paste0(candidates$proband[i],'_sorted_chrN.bam')
params$gene_variants_annotated = paste0(params$out_dir,'/gene',candidates$geneID[i],'variants_annotated.tsv')
params$region = paste0('chr',chr,":",candidates$start[i],"-",candidates$stop[i])

if(params$geneID != "") {

#Prepare the bam subsetting BASH script
command = paste('./consensus.sh',params$geneID, params$transcript,params$ref_file,params$ref_annot, params$out_dir,params$bam_file,params$consensus, params$region)
system(command)

#format the outputs of the bash script 
variant_annotated = read.table(params$gene_variants_annotated,header = T,comment.char ='',check.names=F,na.strings = "",sep = '\t')
colnames(variant_annotated) = c('chr','position','depth','reference','alternate','transcript','exon','gene')
variant_annotated$alternate[variant_annotated$alternate =='.'] = variant_annotated$reference[variant_annotated$alternate =='.']
variant_annotated$transcript = gsub(';','',variant_annotated$transcript)
variant_annotated$transcript[is.na(variant_annotated$transcript)] = 'intron'
variant_annotated$exon = as.numeric(gsub(';','',variant_annotated$exon))
variant_annotated$gene = gsub(';','',variant_annotated$gene)

#annotate the introns
variant_rle  = rle(variant_annotated$transcript)
variant_rle = data.frame(values = variant_rle$values, cumsum = cumsum(variant_rle$lengths),lengths = variant_rle$lengths,event = variant_rle$values)
count=1

for(n in 1:nrow(variant_rle)){
    if(variant_rle$event[n]=="intron") {variant_rle$values[n] = paste0('intron',count);count=count+1}
    }

variant_rle = variant_rle[c(1,1:nrow(variant_rle)),]
variant_rle[1,2:3] = c(0,0)

for(n in 2:nrow(variant_rle)){
    if(variant_rle$event[n] == "intron") {variant_annotated$exon[(variant_rle$cumsum[n-1]+1) : variant_rle$cumsum[n]] = variant_rle$values[n]}
}

#mean expression per intronic position
mean_intronic_depth = variant_annotated %>% filter(transcript == 'intron') %>% group_by(exon) %>% summarise(mean = mean(depth))

#mean expression per exonic position
exonic = variant_annotated[variant_annotated$transcript !='intron',]
exclude_first_and_last = unique(exonic$exon)
exonic = exonic[exonic$exon %in% exclude_first_and_last[-c(1,length(exclude_first_and_last))],] 
mean_exonic_depth = mean(exonic$depth)



#code the retention event if the whole intron is expressed at more than 30% the exon level..
variant_annotated$retention_event = 0

for(r in 1:nrow(mean_intronic_depth)){
    if(mean_intronic_depth$mean[r] > (mean_exonic_depth*0.3) ) variant_annotated$retention_event[variant_annotated$exon == mean_intronic_depth$exon[r]] = 1
}

#variant_annotated$retention_event[variant_annotated$depth > (mean_exonic_depth*0.3)] = 1
#variant_annotated$retention_event[variant_annotated$transcript != 'intron'] = 0 

#tolower the retention and the SNPs / '---' the retention in the reference
#variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate] = tolower(variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate])
variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate] = paste0("<b>",variant_annotated$alternate[variant_annotated$reference != variant_annotated$alternate],"</b>")
variant_annotated$alternate[variant_annotated$retention_event == 1] = tolower(variant_annotated$alternate[variant_annotated$retention_event == 1])
variant_annotated$reference[variant_annotated$retention_event == 1] = '-'

#sequences
reference_sequence = variant_annotated$reference[variant_annotated$transcript != 'intron' | variant_annotated$retention_event == 1]
#reference_sequence = variant_annotated$reference[variant_annotated$transcript != 'intron']
alternate_sequence = variant_annotated$alternate[variant_annotated$transcript != 'intron' | variant_annotated$retention_event == 1]

#Annotation of IUPAC mutations
#iupac = data.frame(symbol = c('r','y','m','k','s','w'), A = c('A','C','A','G','C','A'), B = c('G','T','C','T','G','T'));iupac$AB = paste0(iupac$A,iupac$B)
#for(m in 1:nrow(variant_annotated)){
#  if(variant_annotated$reference[m] != variant_annotated$alternate[m]) {
#    print(variant_annotated[m,])
#    variable_site = paste0(sort(as.character(variant_annotated[m,4:5])),collapse = '')
#    variable_iupac = iupac$symbol[iupac$AB == variable_site]
#    if(length(variable_iupac)==0) variable_iupac = 'N'
#    variant_annotated$alternate[m] = variable_iupac
#  }
#}

#check the retention events and sequence lengths
rle(variant_annotated$retention_event)
print(paste0('length REF: ',length(reference_sequence)))
print(paste0('length ATL: ',length(alternate_sequence)))

# Write to a FASTA file
my_sequences <- list(reference = reference_sequence, alternate = alternate_sequence)
write.fasta(sequences = my_sequences, names = names(my_sequences),nbchar = 100, file.out = paste0(params$consensus,'gene',params$geneID,'_',candidates$proband[i],".fasta"), open = "w")

}

print(paste0('Done sample ',i, ' ~~~  ', params$geneID,' ~~~ Time is: ',Sys.time()))


