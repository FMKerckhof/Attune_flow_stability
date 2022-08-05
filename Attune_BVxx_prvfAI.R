### load libraries -------------------------------------------------------------
library(KytoFlow)
library(tidyverse)
library(flowAI)

### custom function ------------------------------------------------------------

flowAIcallR <- function(flowset,param=c("BL1-H","BL2-H","SSC-H","FSC-H"),
                        timesplit=0.1,timechannel="Time",experiment="AttuneYG",
                        minevents=1e3L,decomp=TRUE){
  # based upon original script by Ruben Props
  # Extract sample names
  samplenames <- flowCore::sampleNames(flowset)
  # Only include samples with enough events
  flowset_sel <- flowset[which(fsApply(flowset,function(x)nrow(x))>=minevents)]
  
  # Extract parameters not to base denoising on
  param_f <- BiocGenerics::unique(gsub(param, pattern = "-H|-A|-W",
                                       replacement = ""))
  filter_param <- BiocGenerics::colnames(flowset)
  filter_param <- BiocGenerics::unique(gsub(filter_param,
                                            pattern = "-H|-A|-W",
                                            replacement = ""))
  filter_param <- filter_param[!filter_param %in% param_f &
                                 filter_param != timechannel]
  filter_param <- c(filter_param, "FSC", "SSC")
  
  # Denoise with flowAI
  x <- flowAI::flow_auto_qc(
    flowset_sel,
    alphaFR = 0.01,
    output = 1,
    timeCh = timechannel,
    ChExcludeFM = filter_param,
    ChExcludeFS = filter_param,
    decompFR = decomp,
    second_fractionFR = timesplit,
    fcs_highQ = "HighQ",
    fcs_QC = FALSE,
    folder_results = paste("resultsQC",experiment,sep="_")
  )
  
  # Change characters in parameter description from character back to numeric
  # Otherwise nothing that follows will work
  for(i in 1:length(x)){
    x[[i]]@parameters@data[,3] <- as.numeric(x[[i]]@parameters@data[,3])
    x[[i]]@parameters@data[,4] <- as.numeric(x[[i]]@parameters@data[,4])
    x[[i]]@parameters@data[,5] <- as.numeric(x[[i]]@parameters@data[,5])
    if(x[[i]]@description$FCSversion == "3"){
      x[[i]]@description$ORIGINALGUID <- sampleNames(x)[i]
    } else x[[i]]@description$GUID.original <- sampleNames(x)[i]
  }
  return(x)
}


### load dataset ---------------------------------------------------------------
load("attuneYGBVbeads.trf.prv.beads.Rda")

### flowAI reduced set ---------------------------------------------------------
param_Attune_YGBV <- c("SSC-H","FSC-H","BL1-H")

attuneYGBVbeads.trf.beads.QC.prv <- try(flowAIcallR(attuneYGBVbeads.trf.prv.beads,
                                                param=param_Attune_YGBV,
                                                experiment="AttuneYGBV_preview",
                                                minevents = 1.5e3L))

attuneYGBVbeads.trf.beads.QC.prv_HR <- try(flowAIcallR(attuneYGBVbeads.trf.prv.beads,
                                                   param=param_Attune_YGBV,
                                                   experiment="AttuneYGBV_preview_HR",
                                                   minevents = 1.5e3L,
                                                   timesplit = 0.01))

attuneYGBVbeads.trf.beads.QC.prv_UHR <- try(flowAIcallR(attuneYGBVbeads.trf.prv.beads,
                                                    param=param_Attune_YGBV,
                                                    experiment="AttuneYGBV_preview_UHR",
                                                    minevents = 1.5e3L,
                                                    timesplit = 0.005))

### write data out -------------------------------------------------------------
save(attuneYGBVbeads.trf.beads.QC.prv,attuneYGBVbeads.trf.beads.QC.prv_HR,
     attuneYGBVbeads.trf.beads.QC.prv_UHR,file="attuneYGbeads.prv.Rda")