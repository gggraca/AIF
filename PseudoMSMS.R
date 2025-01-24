#' Function to obtain pseudo MS/MS spectra from a feature of
#' interest from All-ion fragmentation experiments (e.g. MSe, bbCID, AIF)
#'
#' @author  Goncalo Graca
#' 24 January 2025
#'
#' @param xcmsF1 XCMS object containing MS1 scans
#' @param xcmsF2 XCMS object containing AIF scans
#' @param peaksF1 LC-MS picked peaks from xcmsF1 dataset using XCMS.
#' @param peaksF2 LC-MS picked peaks from xcmsF2 dataset using XCMS.
#' @param fmz The m/z for the feature of interest.
#' @param frt Retention time in seconds for the feature of interest.
#' @param dppm m/z search tolerance for the feature of interest in ppm.
#' @param drt Retention time tolerance for the feature of interest in seconds.
#' @param cthr Pearson correlation coefficient threshold for EIC profile correlation
#' @param savePseudoMSMS if the generated the generated PseudoMSMS should be saved 
 
PseudoMSMS <- function(xcmsF1, xcmsF2, ms1_peaks, ms2_peaks, fmz, frt, dppm = 5, drt = 30, cthr = 0.8, savePseudoMSMS = TRUE){
  # start with empty pseudoSpec and eicStack objects
  pseudoSpec <- NA
  eicStack <- NA
  # search feature in ms1_peaks
  idx <- which(abs(1e6*(ms1_peaks[,"mz"]-fmz)/fmz) < dppm & 
                 abs((ms1_peaks[,"rt"]-frt)) < drt)
  if(length(idx) == 0){
    stop("MS1 feature not found!")
  } else if(length(idx) > 1) {
    stop("More than one MS1 feature found! Adjust your search criteria!")
  } else if(length(idx) == 1) {
    f <- ms1_peaks[idx,]
  }
  
  # get the MS1 feature EIC
  eicMS1 <- chromatogram(xcmsF1, 
                         mz = c(f["mzmin"], f["mzmax"]),
                         rt = c(f["rtmin"], f["rtmax"]))
  
  # get feature EIC as data frame
  feic <- cbind(intensity(eicMS1[1,1]), rtime(eicMS1[1,1]))
  colnames(feic) <- c("intensity", "rt")
  feic <- as.data.frame(feic)
  
  # get MS2 features within the same RT window of MS1 feature
  idx <- which(ms2_peaks[,"rt"] >= f["rtmin"] &
                 ms2_peaks[,"rt"] <= f["rtmax"])
  
  # save table of selected features
  if(length(idx) == 0) {
    stop("No MS2 features found at the MS1 feature RT window!")
  } else {
    selMS2 <- ms2_peaks[idx,]
  }
  # add filter for features with mz > fmz
  selMS2 <- selMS2[which(selMS2[,"mz"] <= fmz),]
  
  if(exists("selMS2")){
    # create a matrix of intensities matched to the MS1 feature RTs
    # the number of rows is the number of intensity points of MS1 feature
    # the number of columns is the number of MS2 features selected
    # the NA will be replaced by the intensity of the RT-matched MS2 feature
    m <- matrix(NA, nrow = nrow(feic), ncol = nrow(selMS2))
    
    #for(i in 1:nrow(selMS2)){
    
    getEICs <- function(x, xcmsF2, selMS2) {
      eicMS2 <- chromatogram(xcmsF2, 
                             mz = c(selMS2[x,"mzmin"], selMS2[x,"mzmax"]),
                             rt = c(selMS2[x,"rtmin"], selMS2[x,"rtmax"]))
      
      # get feature EIC as data frame
      #frag <- rep(paste(round(selMS2[x,"mz"], 4), "m/z"), 
      #            length(rtime(eicMS2[1,1])))
      #f2eic <- cbind(intensity(eicMS2[1,1]), rtime(eicMS2[1,1]), frag)
      f2eic <- data.frame(intensity = intensity(eicMS2[1,1]), 
                          rtime = rtime(eicMS2[1,1]),
                          mz = rep(paste(round(selMS2[x,"mz"], 4), "m/z"), 
                                   length(rtime(eicMS2[1,1]))))
      colnames(f2eic) <- c("intensity", "rt")
      return(f2eic)
    }
    
    f2eics  <- lapply(1:nrow(selMS2), function(x) getEICs(x, xcmsF2, selMS2))
    
    # Get the intensity closest to the RTs of the MS1 feature from each AIF ion:
    for(i in 1:length(f2eics)){
      for(j in 1:nrow(f2eics[[i]])){
        tmp <- which.min(abs(feic[,2]-f2eics[[i]][j,2]))
        m[tmp,i] <- f2eics[[i]][j,1]
      }
    }
    
    # calculate correlations and select features for pseudo-MSMS
    c <- cor(feic[,1], m, use = "pairwise.complete.obs")
    idx <- which(c > cthr)
    if(length(idx) == 0){
      print("No MS2 features found matching the correlation criteria!")
    } else {
      pseudoSpec <- selMS2[idx,c("mz","into")]
      f2eics <- f2eics[idx]
      # labels for the EICs should be tack from pseudoSpec mz variables
      
    }
    # List of objects to return
    r <- list(fmz, frt, f, pseudoSpec, c, selMS2, feic, f2eics)
    
    # plot EICs and Pseudo-MS/MS spectrum
    df1 <- do.call("rbind", f2eics)
    colnames(df1)[3] <- "mz" 
    df2 <- as.data.frame(pseudoSpec)
    # EICs
    p1 <- ggplot(df1[!is.na(df1$intensity),],
                 aes(x = rt, y = intensity, colour = as.factor(mz))) +
      geom_point() +
      geom_line() +
      labs(x = "RT (s)", y = "Intensity (a.u.)",
           colour = paste("Correlation >", cthr)) +
      ggtitle(paste("Correlated EICs:",
                    round(fmz,4),"m/z",round(frt),"s")) +
      theme(plot.title = element_text(hjust = 0.5))
    # Pseudo MS/MS
    p2 <- ggplot(df2, aes(x = mz, y = into, label = round(mz,4))) +
      geom_segment(aes(xend = mz, yend=0),
                   color="red", lwd=0.5) +
      geom_text(size=3, angle=45, hjust=0, vjust=0) +
      theme_minimal() +
      ggtitle(paste("Pseudo-MS/MS:", round(fmz,4),"m/z",round(frt),"s")) +
      theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ylim(0, max(df2[,2]) + 0.1*max(df2[,2])) +
      xlim(min(df2[,1])-50, max(df2[,1])+50) +
      labs(x = "m/z", y = "Intensity (a.u.)")
    
    pdf(paste("AIF_PseudoMSMS_", round(fmz, 4),
              "mz_",round(frt),"s",".pdf", sep = ""), height = 8, width = 12)
    grid.arrange(p1, p2, nrow = 1)
    dev.off()
    
    # return(r)... don't show any values on screen
    
    # save PseudoMSMS as .csv and .mgf
    if(savePseudoMSMS){
      if(!is.null(pseudoSpec)) {
        # save as .csv
        write.csv(pseudoSpec, paste("pseudoMS_AIF_", 
                                    round(fmz, 3), "mz_",round(frt),"s", 
                                    ".csv", sep = ""), row.names = FALSE)
        # and save as .mgf
        fname <- paste("pseudoMS_AIF_", 
                       round(fmz, 3), "mz_",round(frt),"s", 
                       ".mgf", sep = "")
        if(unique(peaksF1@featureData@data$polarity) == 1) charge = "1+"
        if(unique(peaksF1@featureData@data$polarity) == 2) charge = "1-"
        cat("BEGIN IONS",
            paste("PEPMASS=", fmz, sep = ""),
            paste("CHARGE=", charge, sep = ""),
            paste("TITLE=", "Pseudo-MS/MS from AIF at", frt, "seconds"),
            sep = "\n", file = fname)
        for(i in 1:nrow(pseudoSpec)) {
          cat(paste(pseudoSpec[i,1], " ", pseudoSpec[i,2], sep = ""), "\n", file = fname, append = TRUE)
        }
        cat("END IONS", "\n", file = fname, append = TRUE)
      } else NULL
    }
    
  } else print("Please choose another feature or modify search criteria.")
}
