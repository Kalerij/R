if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

#to Russian language
Sys.setlocale('LC_ALL', 'russian')

#load expression data
raw_counts <- read.csv("C:/Users/rodic/Downloads/raw_counts,bladder_tcga.csv",row.names = 1)
raw_counts <- na.omit(raw_counts)
sum(is.na(raw_counts))
raw_counts[1:3,1:3]
#load metadata
md <- read.csv("C:/Users/rodic/Downloads/metadata,bladder_tcga.csv")

table(md$primary_site)
table(md$sample_type)

metadata_interest <- data.frame("sample_ID"=gsub("-",".",md$entity_submitter_id),"type" = md$sample_type)

metadata_interest$type = as.factor(metadata_interest$type)                           
dds <- DESeqDataSetFromMatrix(countData =raw_counts, colData = metadata_interest, design = ~ type)

dds <- DESeq(dds)
res <- na.omit(as.data.frame(results(dds)))
head(res)
#how many differential genes were found
print(sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1))



#Deseq2 normalization
# 1st way
#calculate differential denes(see above)
normalized_counts1 <- counts(dds, normalized=TRUE)

#2nd way 
normalized_counts2 <-  t(t(as.matrix(raw_counts)) / estimateSizeFactorsForMatrix(as.matrix(raw_counts)))

#3d way
#it is possible BEFORE differential gene search
dds <- estimateSizeFactors(dds)
normalized_counts3 <-  raw_counts/dds$sizeFactor
