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
# apt-get install libtiff-dev
if (!require("tiff")) {
  BiocManager::install("tiff", ask =FALSE)
  library("tiff")
}
#sudo apt-get install libfftw3-dev libfftw3-doc
if (!require("fftwtools")) {
  BiocManager::install("fftwtools", ask =FALSE)
  library("fftwtools")
}
if (!require("EBImage")) {
  BiocManager::install("EBImage", ask =FALSE)
  library("EBImage")
}
if (!require("Cardinal")) {
  BiocManager::install("Cardinal", ask =FALSE)
  library("Cardinal")
}
if (!require("CardinalWorkflows")) {
  BiocManager::install("CardinalWorkflows", ask =FALSE)
  library("CardinalWorkflows")
}



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

#######
data(pig206, package="CardinalWorkflows")
pig206 <- as(pig206, "MSImagingExperiment")
pig206
image(pig206, mz=885.5, plusminus=0.25)
image(pig206, mz=150.5, plusminus=0.25)
str(pig206)
str(pig206@featureData)
pig206@featureData@mz
pig206@featureData¨
pig206@featureData@resolution


#Preprocessing
pig206_mean <- summarizeFeatures(pig206, "mean")
plot(pig206_mean)
str(pig206_mean)
?summarizeFeatures
?summarizePixels

pig206_tic <- summarizePixels( pig206, c(tic="sum"))
image(pig206_tic , mz=160, plusminus=0.25)
dev.off()
image(pig206_tic , mz=885.5, plusminus=0.25)
image(pig206_tic)

pig206_ref <- pig206_mean %>%
  peakPick(SNR=3) %>%
  peakAlign(ref="mean",
            tolerance=0.5,
            units="mz") %>%
  peakFilter() %>%
  process()
image(pig206_ref)

pig206_peaks <- pig206 %>%
  normalize(method="tic") %>%
  peakBin(ref=mz(pig206_ref),
          tolerance=0.5,
          units="mz") %>%
  process()
image(pig206_peaks)

pig206_peaks

image(pig206_peaks, mz=187) # heart
image(pig206_peaks, mz=250)
image(pig206_peaks, mz=350)
image(pig206_peaks, mz=450)
image(pig206_peaks, mz=530)
image(pig206_peaks, mz=535)
image(pig206_peaks, mz=537) # Liver
image(pig206_peaks, mz=539)
image(pig206_peaks, mz=550)
image(pig206_peaks, mz=650)
image(pig206_peaks, mz=660)
image(pig206_peaks, mz=655)
image(pig206_peaks, mz=670)
image(pig206_peaks, mz=680)
image(pig206_peaks, mz=690)
dev.off()
image(pig206_peaks, mz=700)
image(pig206_peaks, mz=710)
image(pig206_peaks, mz=720)
image(pig206_peaks, mz=730)
image(pig206_peaks, mz=790)
image(pig206_peaks, mz=800)
image(pig206_peaks, mz=805)
image(pig206_peaks, mz=807)
image(pig206_peaks, mz=820)
image(pig206_peaks, mz=830)
image(pig206_peaks, mz=840)# spinal cord


pig206_pca <- PCA(pig206_peaks, ncomp=3)
str(pig206_pca)
image(pig206_pca, contrast.enhance="histogram", normalize.image="linear")
?Cardinal::image
plot(pig206_pca, lwd=2)

########
## Shrunken centroids
set.seed(1)
pig206_ssc <- spatialShrunkenCentroids(pig206_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=10)
summary(pig206_ssc)

list(s=c(10,15,20,25))

image(pig206_ssc, model=list(s=c(10,15,20,25)))

image(pig206_ssc, model=list(s=20))
plot(pig206_ssc, model=list(s=20), lwd=2)

cols <- discrete.colors(6)
setup.layout(c(3,1))
plot(pig206_ssc, model=list(s=20), column=1, col=cols[1], lwd=2, layout=NULL)
plot(pig206_ssc, model=list(s=20), column=5, col=cols[5], lwd=2, layout=NULL)
plot(pig206_ssc, model=list(s=20), column=6, col=cols[6], lwd=2, layout=NULL)

plot(pig206_ssc, model=list(s=20), values="statistic", lwd=2)

setup.layout(c(3,1))
plot(pig206_ssc, model=list(s=20), values="statistic",
     column=1, col=cols[1], lwd=2, layout=NULL)
plot(pig206_ssc, model=list(s=20), values="statistic",
     column=5, col=cols[5], lwd=2, layout=NULL)
plot(pig206_ssc, model=list(s=20), values="statistic",
     column=6, col=cols[6], lwd=2, layout=NULL)



topFeatures(pig206_ssc, model=list(s=20), class==1)
topFeatures(pig206_ssc, model=list(s=20), class==5)
topFeatures(pig206_ssc, model=list(s=20), class==6)
