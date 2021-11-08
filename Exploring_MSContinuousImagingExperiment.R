#Building from scratch
# Exploring MSContinuousImagingExperiment

set.seed(2020)
s <- simulateSpectrum(n=9, peaks=10, from=500, to=600)

str(s)
s[["mz"]]
s[["intensity"]]

coord <- expand.grid(x=1:3, y=1:3)
run <- factor(rep("run0", nrow(coord)))

fdata <- MassDataFrame(mz=s$mz)
pdata <- PositionDataFrame(run=run, coord=coord)

out <- MSImagingExperiment(imageData=s$intensity,
                           featureData=fdata,
                           pixelData=pdata)
out

hiPSC_200909
hiPSC_VI

str(hiPSC_200909)
str(hiPSC_VI)

hiPSC_200909@metadata
hiPSC_VI@metadata

str(hiPSC_200909@imageData)
str(hiPSC_VI@imageData)

str(hiPSC_200909@featureData)
str(hiPSC_VI@featureData)

hiPSC_200909@elementMetadata
hiPSC_VI@elementMetadata
str(hiPSC_200909@elementMetadata)
str(hiPSC_VI@elementMetadata)

levels(hiPSC_VI@elementMetadata@run)
levels(hiPSC_200909@elementMetadata@run)

spec_V <- spectra(hiPSC_200909)*.01
plot(colSums(spec_V))
plot(rowSums(spec_V))
colSums(spectra(hiPSC_200909))
plot(rowSums(spectra(hiPSC_200909)))

spectra(hiPSC_VI)
#pixelData(out)
#mz(out)
