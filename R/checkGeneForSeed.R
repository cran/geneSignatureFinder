checkGeneForSeed <-
function(geneLevels, evaluateBICs = TRUE, cutoff = 1.95)
  {
    aClassify <- classify(geneLevels)
    
    tmp1 <- min(survfit(stData[aClassify$clusters[!aClassify$missing] == 1]~ 1)$surv)
    tmp2 <- min(survfit(stData[aClassify$clusters[!aClassify$missing] == 2]~ 1)$surv)
    if(tmp1 > tmp2) {	
      aClassify$clusters[!aClassify$missing][aClassify$clusters[!aClassify$missing] == 1] <- 0
      aClassify$clusters[!aClassify$missing][aClassify$clusters[!aClassify$missing] == 2] <- 1
    } else aClassify$clusters[!aClassify$missing][aClassify$clusters[!aClassify$missing] == 2] <- 0

    if(evaluateBICs)
      bics <- BICs(geneLevels[!aClassify$missing],
                   aClassify$clusters[!aClassify$missing],
                   cutoff = cutoff)$bics else bics <- c(NA, NA)
    tValue <- survdiff(stData[!aClassify$missing] ~ aClassify$clusters[!aClassify$missing])$chisq
    pValue <- 1 - pchisq(tValue, df = 1)
    ans <- c(tValue, pValue, bics, aClassify$clusters)
    names(ans) <- c("tValue", "pValue", "bic1", "bic2", names(geneLevels))
    return(ans)

  }
