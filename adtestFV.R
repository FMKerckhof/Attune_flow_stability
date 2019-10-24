### load required libraries ----------------------------------------------------
library(kSamples)

### load data  -----------------------------------------------------------------
load("FACSVerseYGsplitteddf.Rda")

### Run AD test ----------------------------------------------------------------
adresFVList <- lapply(combinedtsVerseyg.md.split,function(x){
  x.split <- split(x,f=x$Gate)
  cat(paste("Processing file",unique(x$Filename),"\n"))
  x.adtestres <- kSamples::ad.test(x.split$ts1$BL1.H,x.split$ts2$BL1.H,
                                   method = "simulated",Nsim=1000)
  return(x.adtestres)
})

### Save results ---------------------------------------------------------------
save(adresFVList,file="AD_test_FACSVerseYG.Rda")