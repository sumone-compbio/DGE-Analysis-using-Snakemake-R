if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocGenerics)
library(stats4)
library(readxl)
library(Biobase)
library(ggplot2)
library(ggrepel)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# 1) DUX4 vs EV

data.tab_dux4 = read.table(file.choose(), header = T, sep = ",")
raw.table = read.csv('dux4_filtered_counts.csv')

## add geneid from original table and remove rows with NA values 

data.tab_dux4 = cbind(raw.table[, 1], data.tab_dux4)
colnames(data.tab_dux4) <- c("GeneID", "baseMean", "log2FoldChange", "IfcSE", "pvalue", "padj")
data.tab_dux4 = data.tab_dux4[complete.cases(data.tab_dux4), ]

# no n/a values observed probably because we filtered the data prior to deseq2 analysis.

upregulated_dux4 = data.tab_dux4[data.tab_dux4$log2FoldChange > 2 & data.tab_dux4$padj < 0.0001, ]
downregulated_dux4 = data.tab_dux4[data.tab_dux4$log2FoldChange < -2 & data.tab_dux4$padj < 0.0001, ]

EnhancedVolcano(data.tab_dux4, x = "log2FoldChange", y = "padj", lab = data.tab_dux4$GeneID)

# 2) IGH vs. EV

data.tab_igh = read.table(file.choose(), header = T, sep = ",")

## add geneid from original table and remove rows with NA values 

data.tab_igh = cbind(raw.table[, 1], data.tab_igh)
colnames(data.tab_igh) <- c("GeneID", "baseMean", "log2FoldChange", "IfcSE", "pvalue", "padj")
data.tab_igh = data.tab_igh[complete.cases(data.tab_igh), ]

upregulated_igh = data.tab_igh[data.tab_igh$log2FoldChange > 2 & data.tab_igh$padj < 0.0001, ]
downregulated_igh = data.tab_igh[data.tab_igh$log2FoldChange < -2 & data.tab_igh$padj < 0.0001, ]

EnhancedVolcano(data.tab_igh, x = "log2FoldChange", y = "padj", lab = data.tab_igh$GeneID)

# 3) delta_50 vs. EV

data.tab_d50 = read.table(file.choose(), header = T, sep = ",")

## add geneid from original table and remove rows with NA values 

data.tab_d50 = cbind(raw.table[, 1], data.tab_d50)
colnames(data.tab_d50) <- c("GeneID", "baseMean", "log2FoldChange", "IfcSE", "pvalue", "padj")
data.tab_d50 = data.tab_d50[complete.cases(data.tab_d50), ]

upregulated_d50 = data.tab_d50[data.tab_d50$log2FoldChange > 2 & data.tab_d50$padj < 0.0001, ]
downregulated_d50 = data.tab_d50[data.tab_d50$log2FoldChange < -2 & data.tab_d50$padj < 0.0001, ]

EnhancedVolcano(data.tab_d50, x = "log2FoldChange", y = "padj", lab = data.tab_d50$GeneID)

