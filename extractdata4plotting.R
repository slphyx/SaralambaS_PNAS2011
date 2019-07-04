## for using with runFixedRing16.R
FittedData <- function(rdsfilename) {
  #Pailin: first position is started at time =0
  ObTimes <- c(1, 2, 4, 6, 8, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78,
               84, 90, 96, 102, 108, 114, 120, 126, 132, 138, 144, 150, 156, 162)
  
    #rdsfilename = datafiles[posIndx];
    cat("Reading your RDS file...")
    
    modfit <- readRDS(rdsfilename)
    
    ex <- rstan::extract(modfit, pars = c("y_pred", "mod_predSplitCT"))

    y_predMed <- apply(ex$y_pred, 2, median)
    y_predCI <- apply(ex$y_pred, 2, quantile, probs = c(0.025, 0.975))


    mod_predSplitCTMed <- apply(ex$mod_predSplitCT, 2, median)
    mod_predSplitCTCI <- apply(ex$mod_predSplitCT, 2, quantile, probs = c(0.025, 0.975))
    obTnp <- bigcounts[[1]]$NPPARA
    obT <- ObTimes[1:obTnp]

    outdataframe <- data.frame(obtime = c(0, obT[2:obTnp]), obdat = bigcounts[[1]]$y, median = y_predMed, 
        below = y_predCI[1,], above = y_predCI[2,], mediansplit = mod_predSplitCTMed[obT],
         belowsplit = mod_predSplitCTCI[1, obT], abovesplit = mod_predSplitCTCI[2, obT])

    write.csv(outdataframe, file = paste0(substr(rdsfilename, 1, 7), ".csv"))
    cat("Output file has been written.\n")
}
