

suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))

data <- read.table("/home/renaut/scratch/raredisease_rnaseq/ASE/overlaps.tsv",header=FALSE)

colnames(data) <- c("sample","gene","ref","alt")

#data$condition <- ifelse(grepl("^CTRL",data$sample),"control","disease")
data$condition = c("ctrl","ctrl","ctrl","ctrl","ctrl","disease","disease","disease","disease","disease")

long <- data %>%
pivot_longer(cols=c(ref,alt),
names_to="allele",
values_to="count")

long$sample_allele <- paste(long$sample,long$allele,sep="_")

count_matrix <- long %>%
select(sample_allele,gene,count) %>%
pivot_wider(names_from=sample_allele,values_from=count)

genes <- count_matrix$gene
count_matrix <- as.data.frame(count_matrix[,-1])
rownames(count_matrix) <- genes

samples <- colnames(count_matrix)

condition <- c("ctrl","ctrl","ctrl","ctrl","ctrl","disease","disease","disease","disease","disease")
allele <- ifelse(grepl("ref$",samples),"ref","alt")

coldata <- data.frame(
condition=factor(condition),
allele=factor(allele)
)

rownames(coldata) <- samples

dds <- DESeqDataSetFromMatrix(
countData=count_matrix,
colData=coldata,
design=~ condition + allele + condition:allele
)

dds <- DESeq(dds)

res <- results(dds,name="conditiondisease.allelealt")

write.csv(as.data.frame(res),"ASE_differential_results.csv")

print("DONE ~~~ ASEpipeline.R")
