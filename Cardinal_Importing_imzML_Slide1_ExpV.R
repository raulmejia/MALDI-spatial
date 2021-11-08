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

#hiPSC_200909 <- readMSIData( "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/what_is_here/maldi/hiPSC_200909/hiPSC_200909_converted_fwhm_0.050000.imzML" )

hiPSC_200909 <- readMSIData( "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/what_is_here/maldi/hiPSC_200909/hiPSC_200909.imzML")
str(hiPSC_200909)
plot(hiPSC_200909@elementMetadata@coord[,"x"])
plot(hiPSC_200909@elementMetadata@coord[,"y"])
plot(hiPSC_200909@elementMetadata@coord[,"x"], hiPSC_200909@elementMetadata@coord[,"y"])

plot(hiPSC_200909)
plot(hiPSC_200909,x=150,y=150, plusminus=10, fun=mean)
dev.off()
image()
image(hiPSC_200909, mz=299, plusminus=1)
image(hiPSC_200909, mz=1000, plusminus=1)

hiPSC_200909_tich <- summarizePixels( hiPSC_200909 , c(tic="sum"))
str(hiPSC_200909)
str(hiPSC_200909@featureData)
range(hiPSC_200909@featureData@mz)

image(hiPSC_200909_tich, mz=299, plusminus=1)
image(hiPSC_200909_tich, mz=700, plusminus=0.25)
image(hiPSC_200909_tich, mz=800, plusminus=0.25)
image(hiPSC_200909_tich, mz=900, plusminus=0.25)
image(hiPSC_200909_tich, mz=1000, plusminus=1)
#509 (2)/ 655 (2)/ 
# 656 (1/2)
# 671 (2)
# 676 (1)
# 782 (1)
# 783 (1)

saveRDS( hiPSC_200909_tich, file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_tich.RDS")
warnings()
plot(hiPSC_200909_tich)

str(hiPSC_200909_tich)

range( hiPSC_200909@featureData@mz )

hiPSC_200909_mean <- summarizeFeatures(hiPSC_200909, "mean")
?summarizeFeatures

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

image( hiPSC_200909 , mz=509)
image( hiPSC_200909_peaks , mz=509)
image( hiPSC_200909 , mz=655.92)
image( hiPSC_200909_peaks , mz=655.92) # One distribution
image( hiPSC_200909 , mz=656.93)
image( hiPSC_200909_peaks , mz=656.93)
image( hiPSC_200909 , mz=671)
image( hiPSC_200909_peaks , mz=671)
image( hiPSC_200909 , mz=676)
image( hiPSC_200909_peaks , mz=676) # Another distribution
image( hiPSC_200909 , mz=782)
image( hiPSC_200909_peaks , mz=782)
image( hiPSC_200909 , mz=783)
image( hiPSC_200909_peaks , mz=783)


hiPSC_200909_pca <- PCA( hiPSC_200909_peaks , ncomp=3)
saveRDS( hiPSC_200909_pca , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_pca.RDS")

image( hiPSC_200909_pca, contrast.enhance="histogram", normalize.image="linear" )
plot( hiPSC_200909_pca , lwd=2 )

hiPSC_200909_ssc <- spatialShrunkenCentroids(hiPSC_200909_peaks, method="adaptive",
                                       r=2, s=c(0,5,10,15,20,25), k=4)
saveRDS( hiPSC_200909_ssc , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_200909_ssc.RDS")


#pig206_ssc <- spatialShrunkenCentroids(hiPSC_200909_peaks, method="adaptive",
#                                       r=2, s=c(0,5,10,15,20,25), k=8)


image( hiPSC_200909_ssc , model=list(s=c(10,15,20,25)))
image( hiPSC_200909_ssc , model=list(s=c(0,5,10,15)))
image(hiPSC_200909_ssc, model=list(s=5))


plot(hiPSC_200909_ssc, model=list( s=25), lwd=2 )
plot(hiPSC_200909_ssc, model=list( s=20), lwd=2 )
plot(hiPSC_200909_ssc, model=list( s=15), lwd=2 )
plot(hiPSC_200909_ssc, model=list( s=10), lwd=2 )
plot(hiPSC_200909_ssc, model=list( s=5), lwd=2 )
plot(hiPSC_200909_ssc, model=list( s=0), lwd=2 )


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



