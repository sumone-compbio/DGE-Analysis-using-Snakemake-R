# Differential-Gene-Expression-DGE-Analysis
A simple-to-follow tutorial on differential gene expression analysis in R

## About the case study: 
### Title: DUX4-r neomorphic activity depending on GTF2I in acute lymphoblastic leukemia [RNA-seq]
##### Description:
The rearranged versions of the transcription factor DUX4 (DUX4-r produced by translocations) are one of the most common causes of B-cell lymphoblastic leukemia (B-ALL). The study discovered that such rearrangements can lead to both a loss and a gain of function in DUX4-r.

Loss:  Loss of CBP/EP300 transcriptional co-activators interaction and inability to bind and activate repressed chromatin. The rearranged DUX4-r can still bind to DNA but has alterations in its C-terminal transcription activation domain, affecting its interaction with the co-activators CBP and EP300.

Gain: Gain of interaction with the transcription factor GTF2I redirecting DUX4-r toward leukemogenic targets. Hence, GTF2I can be a potential target to inhibit DUX4-r from causing leukemia.

##### My Aim:
How different variants of DUX4 e.g. DUX4 (wild type), DUX4-IGH, and DUX4-del50 impact different genes and associated ontologies e.g. Biological Pathways.

##### Link to the study: 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227982

## Plots 
### Let's look at the plots of the analysis:

#### 1) PCA Plot
PCA plot of the study. Ideally, each condition should form a distinct cluster as shown here.
![pca_plot_all_condition](https://github.com/sumenties/Differential-Gene-Expression-DGE-Analysis/assets/43076959/ed36cc2f-978f-4fff-81c6-f9f3ef1b8b59)

#### 2) Dispersion Plot
The fitted (red) line should be lower than the dispersion value of 1. 
![dispersion_plot](https://github.com/sumenties/Differential-Gene-Expression-DGE-Analysis/assets/43076959/fdaed0fb-dc4a-4a70-80a0-71ffae5953ad)

#### 3) MA plots for each condition compared to control, i.e., EV in this study. 

##### DUX4 vs EV
We're interested in outliers here. Outliers at the top are overexpressed genes and the outliers at the bottom are underexpressed genes. 
![dux4_vs_ev](https://github.com/sumenties/Differential-Gene-Expression-DGE-Analysis/assets/43076959/78a1a7dd-b751-4fd8-81fc-564c7ed9805f)

Similarly, MA plots for other conditions can be drawn following the code in the deseq2.R file. 

#### 4) Volcano PLot

##### IGH vs EV
Left-hand side: Downregulated genes (log2FoldChange>2, p-adj value < 0.0001)
Right-hand side: Upregulated genes (log2FoldChange<-2, p-adj value < 0.0001)

Significant genes are labeled based on their log2FoldChange and p-adj value.
![igh_vs_ev_vp](https://github.com/sumone-compbio/Differential-Gene-Expression-DGE-Analysis/assets/43076959/5cadb0ff-54fb-4989-a85b-19ae801a20a1)


Significant genes are labelled based on their log2FoldChange and p-adj values. 

###### Next: Gene Ontology (GO), and GSEA to identify the roles and significance of potential biomarkers.
