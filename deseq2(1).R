BiocManager::install("BiocGenerics")

# Importing libraries required  
library(BiocGenerics)
library(tidyverse)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(matrixStats)
library(MatrixGenerics)
library(SummarizedExperiment)
library(DESeq2)
library(readxl)

# loading data of the study, removing columns not required 
dux4 = read_excel("GSE227982_all_counts.xlsx")
dux4_new = dux4[, !names(dux4) %in% c("Chr", "Start", "End", "Strand")]

# keeping only the expression data, reordering columns as per convenience 
all_data = dux4_new[, c(3:18)]
new_order = c("EV_1", "EV_2", "EV_3", "EV_4",
              "DUX4_1", "DUX4_2", "DUX4_3", "DUX4_4",
              "IGH_1", "IGH_2", "IGH_3", "IGH_4",
              "delta_50_1", "delta_50_2", "delta_50_3", "delta_50_4")
all_data = all_data[new_order]

# defining conditions as per the columns 
condition = rep(c("EV", "DUX4", "IGH", "delta_50"), each = 4)
rnames = new_order
column_data = cbind(rnames, condition)

rownames(column_data) = column_data[, 1]
all(colnames(all_data) == rownames(column_data))

# creating the DESeqDataSet object
dds = DESeqDataSetFromMatrix(countData = all_data, colData = column_data, design = ~ condition)

## set the factor levels
dds$condition = relevel(dds$condition, ref = "EV")

## run deseq for EV vs all
dds = DESeq(dds)
res = results(dds)
summary(res)

## MA plot for EV vs all
plotMA(res)
write.table(res, file = "normal_vs_all_dux4.csv", row.names=F, sep = ",")  

# 1) To compare DUX4 vs. EV
DUX4_vs_EV = results(dds, contrast=c("condition", "DUX4", "EV"))
summary(DUX4_vs_EV)
### MA plot
plotMA(DUX4_vs_EV)
write.table(res, file = "dux4_vs_normal_dux4.csv", row.names=F, sep = ",") 

# 2) To compare IGH vs. EV
IGH_vs_EV = results(dds, contrast=c("condition", "IGH", "EV"))
summary(IGH_vs_EV)
### MA plot
plotMA(IGH_vs_EV)
write.table(res, file = "IGH_vs_normal_dux4.csv", row.names=F, sep = ",") 

# 3) To compare delta_50 vs. EV
delta_50_vs_EV = results(dds, contrast=c("condition", "delta_50", "EV"))
summary(delta_50_vs_EV)
### MA plot
plotMA(delta_50_vs_EV)
write.table(res, file = "delta_50_vs_normal_dux4.csv", row.names=F, sep = ",") 
