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
new_order = c("EV_1", "EV_2", "EV_3", "EV_4",
              "DUX4_1", "DUX4_2", "DUX4_3", "DUX4_4",
              "IGH_1", "IGH_2", "IGH_3", "IGH_4",
              "delta_50_1", "delta_50_2", "delta_50_3", "delta_50_4")

# removing rows with row sum < 20 to remove noise from the data. 20 is just an arbitrary number.
dux4_filtered = dux4_new[which(rowSums(dux4_new[, new_order]) > 20), ]
write.table(dux4_filtered, file = "dux4_filtered_counts.csv", row.names=F, sep = ",")  

# keeping only the expression data, reordering columns as per convenience 
all_data = dux4_filtered[, c(3:18)]
all_data = all_data[new_order]

# defining conditions as per the columns 
condition = factor(rep(c("EV", "DUX4", "IGH", "delta_50"), each = 4))
rnames = new_order
column_data = data.frame(rnames, condition)

# verifying whether the columns of the expression dataframe are as same as the row of column_data
rownames(column_data) = column_data[, 1]
all(colnames(all_data) == rownames(column_data))

# creating the DESeqDataSet object
dds = DESeqDataSetFromMatrix(countData = all_data, colData = column_data, design = ~ condition)
dds = DESeq(dds)

## variance-stabilizing transformation to stabilize variance across the expression values.
vst_res = vst(dds, blind = FALSE)

## PCA Plot (Ideally each group should represent one unique cluster)
plotPCA(vst_res, intgroup = "condition")

## Plot dispersion estimates (Ideally, fitted < dispersion = 1)
plotDispEsts(dds)

# 1) To compare DUX4 vs. EV
DUX4_vs_EV = results(dds, contrast=c("condition", "DUX4", "EV"))
summary(DUX4_vs_EV)
## MA plot
plotMA(DUX4_vs_EV)
### lfcshrinkage for more stable estimates
DUX4_vs_EV = lfcShrink(dds, contrast=c("condition", "DUX4", "EV"), type="apeglm")
write.table(DUX4_vs_EV, file = "dux4_vs_normal_dux4.csv", row.names=F, sep = ",") 

# 2) To compare IGH vs. EV
IGH_vs_EV = results(dds, contrast=c("condition", "IGH", "EV"))
summary(IGH_vs_EV)
### MA plot
plotMA(IGH_vs_EV)
write.table(IGH_vs_EV, file = "IGH_vs_normal_dux4.csv", row.names=F, sep = ",") 

# 3) To compare delta_50 vs. EV
delta_50_vs_EV = results(dds, contrast=c("condition", "delta_50", "EV"))
summary(delta_50_vs_EV)
### MA plot
plotMA(delta_50_vs_EV)
write.table(delta_50_vs_EV, file = "delta_50_vs_normal_dux4.csv", row.names=F, sep = ",")
