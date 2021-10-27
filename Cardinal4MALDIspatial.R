###################################
#### This script makes  https://github.com/raulmejia/..
#### Based on https://www.bioconductor.org/packages/release/data/experiment/html/CardinalWorkflows.html
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
# apt-get install libtiff-dev
if (!require("tiff")) {
  BiocManager::install("tiff", ask =FALSE)
  library("tiff")
}
if (!require("fftwtools")) {
  BiocManager::install("fftwtools", ask =FALSE)
  library("fftwtools")
} ## Error :/ 


# tiff, fftwtools, EBImage , Cardinal


############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-m", "--matrix", type="character", 
                    help="path to your expression matrix")
parser$add_argument("-a", "--annotation", type="character", 
                    help="path to your annotation file")
parser$add_argument("-c", "--code", type="character", 
                    help="path to your code")
parser$add_argument("-l", "--label", type="character", 
                    help="label to your results")
parser$add_argument("-g", "--maingroups", type="character", 
                    help="the name of your column to correct / make intrabatch normalization")
parser$add_argument("-s", "--selectedgenesfromrownamesofadataframe", type="character", 
                    help="A data frame that contains the genes that you want to plot in the row names, only the rownames will be used for that tash")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output folder where you want to store your results")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args( )
# print some progress messages to stderr if "quietly" wasn't requested

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



