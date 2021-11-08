dir.create( recursive = TRUE,
            "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_VI")

hiPSC_VI <- readMSIData( "/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/what_is_here/maldi/VI_hiPSC/VI_hiPSC.imzML")
plot(hiPSC_VI)
hiPSC_VI_tic <- summarizePixels(  hiPSC_VI  , c(tic="sum"))
?saveRDS

range( hiPSC_VI@featureData@mz )
hiPSC_VI_mean <- summarizeFeatures( hiPSC_VI, "mean")

hiPSC_VI_ref <- hiPSC_VI_mean %>%
  peakPick(SNR=3) %>%
  peakAlign(ref="mean",
            tolerance=0.5,
            units="mz") %>%
  peakFilter() %>%
  process()

hiPSC_VI_peaks <- hiPSC_VI_mean %>%
  Cardinal::normalize(method="tic") %>%
  peakBin(ref=mz( hiPSC_VI_ref ),
          tolerance=0.5,
          units="mz") %>%
  process()

saveRDS( hiPSC_VI_peaks , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_VI/hiPSC_VI_peaks.RDS")

hiPSC_VI_pca <- PCA(  hiPSC_VI_peaks , ncomp=3)
saveRDS( hiPSC_VI_pca ,
         file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_VI_pca.RDS")
image( hiPSC_VI_pca , contrast.enhance="histogram", normalize.image="linear" )
plot( hiPSC_VI_pca , lwd=2 )

hiPSC_VI_ssc <- spatialShrunkenCentroids( hiPSC_VI_peaks , method="adaptive",
                                             r=2, s=c(0,5,10,15,20,25), k=4)

saveRDS( hiPSC_VI_ssc , file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Results/Cardinal_Spectrum/hiPSC_200909/hiPSC_VI_ssc.RDS")


image( hiPSC_200909_ssc , model=list(s=c(10,15,20,25)))
image( hiPSC_200909_ssc , model=list(s=c(0,5,10,15)))