library(clusterProfiler)
library(org.Hs.eg.db)

##### GSEA analysis 

##### IGH vs EV

data.tab_igh = read.table(file.choose(), header = T, sep = ",")

data.tab_igh.order <- data.tab_igh[order(-data.tab_igh$log2FoldChange),]

gene_list_igh <- data.tab_igh.order$log2FoldChange
names(gene_list_igh) <- data.tab_igh.order$GeneID
gene_list_igh

gse_igh <- gseGO(gene_list_igh, ont = "BP", keyType = "SYMBOL", OrgDb = org.Hs.eg.db, eps = 1e-300, nPermSimple = 10000)
gse_igh.df <- as.data.frame(gse_igh)

library(ggplot2)
library(enrichplot) 

# dotplot
plot <- dotplot(gse_igh, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Adjusting font sizes, theme settings, space between y-axis labels, and a title

plot + ggtitle("IGH vs EV") + #title
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1), #size and angle for x-axis labels
    axis.text.y = element_text(size = 10), #size for y-axis labels
    strip.text = element_text(size = 12), #size for facet labels
    legend.text = element_text(size = 10), #size for legend text
    legend.title = element_text(size = 12), #size for legend title
    plot.margin = margin(1, 1, 1, 2, "cm"), #margin on the left side to avoid overlapping of texts
    panel.spacing = unit(1, "cm"), #space between facets
    plot.title = element_text(hjust = 0.5, size = 14) #title font size
  ) + 
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.2))) #space on y-axis
