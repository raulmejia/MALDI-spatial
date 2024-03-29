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

#hiPSC_200909 <- readMSIData( "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/what_is_here/maldi/hiPSC_200909/hiPSC_200909_converted_fwhm_0.050000.imzML",
             )

hiPSC_200909 <- readMSIData( "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/what_is_here/maldi/hiPSC_200909/hiPSC_200909.imzML")
plot(hiPSC_200909)
plot(hiPSC_200909,x=150,y=150, plusminus=10, fun=mean)
dev.off()
image()
image(hiPSC_200909, mz=299, plusminus=1)
image(hiPSC_200909, mz=1000, plusminus=1)

hiPSC_200909_tich <- summarizePixels( hiPSC_200909 , c(tic="sum"))
str(hiPSC_200909)
str(hiPSC_200909@featureData)

image(hiPSC_200909_tich, mz=299, plusminus=1)
image(hiPSC_200909_tich, mz=700, plusminus=0.25)
image(hiPSC_200909_tich, mz=800, plusminus=0.25)
image(hiPSC_200909_tich, mz=900, plusminus=0.25)
image(hiPSC_200909_tich, mz=1000, plusminus=1)
saveRDS( hiPSC_200909_tich, file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_tich.RDS")
?saveRDS
warnings()
plot(hiPSC_200909_tich)


str(hiPSC_200909_tich)

range( hiPSC_200909@featureData@mz )

hiPSC_200909_mean <- summarizeFeatures(hiPSC_200909, "mean")

hiPSC_200909_ref <- hiPSC_200909_mean %>%
  peakPick(SNR=3) %>%
  peakAlign(ref="mean",
            tolerance=0.5,
            units="mz") %>%
  peakFilter() %>%
  process()

hiPSC_200909_peaks <- hiPSC_200909 %>%
  Cardinal::normalize(method="tic") %>%
  peakBin(ref=mz(hiPSC_200909_ref),
          tolerance=0.5,
          units="mz") %>%
  process()

saveRDS( hiPSC_200909_peaks , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_peaks.RDS")

#?plot
# Hand selected peaks
#image( hiPSC_200909_peaks , mz=537)
image( hiPSC_200909_peaks , mz=600)
image( hiPSC_200909_peaks , mz=700)
image( hiPSC_200909_peaks , mz=750)
image( hiPSC_200909_peaks , mz=782)
image( hiPSC_200909_peaks , mz=18354)

hiPSC_200909_pca <- PCA( hiPSC_200909_peaks , ncomp=3)
saveRDS( hiPSC_200909_pca , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_pca.RDS")

image(hiPSC_200909_pca, contrast.enhance="histogram", normalize.image="linear")
plot(hiPSC_200909_pca , lwd=2)

hiPSC_200909_ssc <- spatialShrunkenCentroids(hiPSC_200909_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=4)


#pig206_ssc <- spatialShrunkenCentroids(hiPSC_200909_peaks, method="adaptive",
#                                       r=2, s=c(0,5,10,15,20,25), k=8)


image( hiPSC_200909_ssc , model=list(s=c(10,15,20,25)))
image( hiPSC_200909_ssc , model=list(s=c(0,5,10,15)))
image(hiPSC_200909_ssc, model=list(s=5))

plot(pig206_ssc, model=list( s=5), lwd=2 )

cols <- discrete.colors(6)
setup.layout(c(2,1))
plot(hiPSC_200909_ssc , model=list(s=5), column=1, col=cols[1], lwd=2, layout=NULL)
plot(hiPSC_200909_ssc, model=list(s=5), column=2, col=cols[2], lwd=2, layout=NULL)



plot(hiPSC_200909_ssc, model=list(s=5), values="statistic", lwd=2)

setup.layout(c(2,1))
plot( hiPSC_200909_ssc, model=list(s=5), values="statistic",
     column=1, col=cols[1], lwd=2, layout=NULL)
plot( hiPSC_200909_ssc, model=list(s=5), values="statistic",
     column=2, col=cols[2], lwd=2, layout=NULL)


topFeatures( hiPSC_200909_ssc, model=list(s=5), class==1)
topFeatures( hiPSC_200909_ssc, model=list(s=5), class==2)
