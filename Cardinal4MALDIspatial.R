###################################
#### This script makes  https://github.com/raulmejia/DSP-Oszwald
####      
####     Input description:
####        expression matrix: 
####     Output:
####
####    Author of the script: Raúl Mejía
####
#### Example:
####      Rscript /PathA/Differential_expression_hybridized-n-scanned-Data-trough-limma.R /PathB/ExpMat_Demo.tsv /PathC/Annot_table_Demo.tsv /Path/code-this-repository/ /Path/yourResults/ label normal 0 0.05 5 Label_4_your_results
#### 
####    
###################################
#### 0) loading and/or installing required libraries
################################### 
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library("tidyverse")
}
if (!require("Cardinal")) {
  BiocManager::install("Cardinal", ask =FALSE)
  library("Cardinal")
}
########################################
#### data given by the user
#########################################
myargs <- commandArgs(trailingOnly = TRUE)

#bug# path_expression_matrix <- "/media/rmejia/mountme88/Projects/Maja-covid/Results/Normalizations/NK_Geo/Majalog2_OAZ1_HPRT1_ABCF1.tsv"
path_expression_matrix <- myargs[1]

data(pig206, package="CardinalWorkflows")
pig206 <- as(pig206, "MSImagingExperiment")
pig206
image(pig206, mz=885.5, plusminus=0.25)

#Preprocessing
pig206_mean <- summarizeFeatures(pig206, "mean")
plot(pig206_mean)



