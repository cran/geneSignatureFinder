pamUnbiased <-
function(data) {
  n <- length(data)
  missing <- is.na(data)

  wdTmp <- data[!missing]
  nn <- n - sum(missing)

  ansClusters <- rep(NA, n)
  i <- 1
  while(i <= nn) {
    nMissing <- sum(missing)
    clusters <- pam(wdTmp[-i], 2)$clustering
    # qui viene aggiornato l'elengo dei valori mancanti
    # con le osservazioni che fanno fallire PAM (valori anomali che generano cluster con 1 osservazione)
    ckTmp <- rep(FALSE, n)
    if(sum(ckTmp[!missing][-i] <- clusters == 1) < 0.025 * n)
      missing[!missing][-i] <- (missing[!missing][-i] + ckTmp[!missing][-i]) > 0

    ckTmp <- rep(FALSE, n)
    if(sum(ckTmp[!missing][-i] <- clusters == 2) < 0.025 * n)
      missing[!missing][-i] <- (missing[!missing][-i] + ckTmp[!missing][-i]) > 0

    if(nMissing < sum(missing)) {
      nMissing <- sum(missing)
      wdTmp <- data[!missing]
      nn <- n - sum(missing)
      ansClusters <- rep(NA, n)
      i <- 1
      next
    }
    
    mTmp <- m1 <- median(wdTmp[-i][clusters == 1])
    m2 <- median(wdTmp[-i][clusters == 2])
# riallineamento della classificazione dei gruppi (1 = gruppo con media inferiore; 2 = ...)
    if(which.min(c(m1, m2)) != 1) {
      clusters[clusters == 2] <- 0
      clusters <- clusters + 1
      m1 <- m2
      m2 <- mTmp
    }
    
    s1 <- bw.nrd0(wdTmp[-i][clusters == 1])
    s2 <- bw.nrd0(wdTmp[-i][clusters == 2])
# classificazione dell'osservazione in base alla regola di min distanza (standardizzata)
    ansClusters[!missing][i] <- which.min(abs((wdTmp[i] - c(m1, m2))/c(s1, s2)))
    i <- i + 1
  }
  ans <- list()
  ans$clusters <- ansClusters
  ans$missing <- missing
  return(ans)
}
