#' Script to obtain pseudo MS/MS spectra from a feature of
#' interest from All-ion fragmentation experiments (e.g. MSe, bbCID, AIF)
#'
#' @author  Goncalo Graca
#'
#' 24 January 2025
#'
#' The script uses XCMS functionality to read and peak-pick the data
#' as well as the PseudoMSMS function used to obtain the pseudo-MS/MS 
#' from a feature of interest

### Read the R packages used to load and plot data ###
library(xcms)
library(ggplot2)
library(gridExtra)

# load the PseudoMSMS.R function that should be located in the working directory
source("PseudoMSMS.R")

### Read and peak-pick data ### 
# repeat this every for every new chromatogram: 

# Read data into XCMS...the example here uses MetaboAnnotatoR demo mzML file
# Separate MS1 and MS2 levels
xcmsF1 <- readMSData("Lipid_Positive_QC.mzML", msLevel. = 1, mode = "onDisk")
xcmsF2 <- readMSData("Lipid_Positive_QC.mzML", msLevel. = 2, mode = "onDisk")

# Only MS1 scans can be used further for EIC extraction
# The 'msLevel' label is then changed from 2 to 1 in xcmsF2 object
xcmsF2@featureData@data$msLevel <- 1

# Peak-picking
cwp <- CentWaveParam(ppm = 5, peakwidth = c(5,30)) # using default parameters...can be adjusted!
peaksF1 <- xcms::findChromPeaks(xcmsF1, param = cwp)
peaksF2 <- xcms::findChromPeaks(xcmsF2, param = cwp)

# Create a tables for MS1 and MS2 peaks (features)
ms1_peaks <- chromPeaks(peaksF1)
ms2_peaks <- chromPeaks(peaksF2)

## Input the m/z and retention time for the feature of interest ##
# replace the values below by those of the feature of interest
fmz <- 585.2691869
frt <- 72.79411156

## Run the PseudoMSMS function ##
PseudoMSMS(xcmsF1, xcmsF2, ms1_peaks, ms2_peaks, fmz, frt, dppm = 5, drt = 30, cthr = 0.8)

# if the feature of interest is found and correlated AIF features are found, then a pdf with the plotted results will be created in the working directory, 
# as well as two files (.csv and .mgf) containing the PseudoMSMS m/z and intensity values

# otherwise error messages will be shown together with the warning: "Please choose another feature or modify search criteria."

# some examples and suggestions to solve the issues are below:

# "MS1 feature not found!" ... this means that the criteria for the feature of interest search must be adjusted... try to increase dppm and/or drt

# "More than one MS1 feature found! Adjust your search criteria!" ... try to decrease dppm or drt so that only the feature of interest is captured

# "No MS2 features found at the MS1 feature RT window!" ... the intensity of the feature of interest might be too low, in this case it might be difficult to find fragment peaks in the AIF scans

# "No MS2 features found matching the correlation criteria!" ... the correlation threshold might be set too high ... try to decrease cthr

