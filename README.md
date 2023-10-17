# Differential-Gene-Expression-DGE-Analysis
A simple-easy follow tutorial on differential gene expression analysis in R

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

