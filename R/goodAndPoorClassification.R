goodAndPoorClassification <-
function(clustering, percentile = .2) {
  missing <- is.na(clustering)
  x <- as.numeric(clustering[!missing])
  levs <- as.numeric(names(table(x)))
  aSurvfit <- survfit(stData[x == levs[1]] ~ 1)
  position <- ceiling(length(aSurvfit$surv) * percentile)
  beta1 <- aSurvfit$surv[length(aSurvfit$surv)-position] - aSurvfit$surv[position]
  beta1 <- beta1/(aSurvfit$time[length(aSurvfit$surv)-position]-aSurvfit$time[position])
  
  aSurvfit <- survfit(stData[x == levs[2]] ~ 1)
  position <- ceiling(length(aSurvfit$surv) * percentile)
  beta2 <- aSurvfit$surv[length(aSurvfit$surv)-position] - aSurvfit$surv[position]
  beta2 <- beta2/(aSurvfit$time[length(aSurvfit$surv)-position]-aSurvfit$time[position])
  
  ans <- rep(NA, length(clustering))
  ans[!missing] <- "good"
  if(beta1 < beta2) 
    ans[!missing][x == levs[1]] <- "poor"  else 
      ans[!missing][x == levs[2]] <- "poor" 
  
  return(as.factor(ans))
}
