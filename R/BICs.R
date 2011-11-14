BICs <-
function(data, clusters, cutoff = 2.5) {
  n <- length(data)
  ans0 <- normalApproximation(data, median(data), bw.nrd0(data), cutoff = cutoff)
  if(ans0$cutoff > 0) {
    bic1 <- dnorm(data,  ans0$center, ans0$deviation)
    if(sum(tmp <- bic1 < .Machine$double.eps) > 0)
      bic1[tmp] <- .Machine$double.eps
    bic1 <- -2 * sum(log(bic1)) + 2 * log(n)
  } else bic1 <- NA

  
  center1 <- median(data[clusters == 1]) 
  sd1 <- bw.nrd0(data[clusters == 1])
#  if(sd1 < .Machine$double.eps) sd1 <- .Machine$double.eps
  result1 <- normalApproximation(data[clusters == 1], center1, sd1, cutoff = cutoff)
  
  center2 <- median(data[clusters == 0])
  sd2 <- bw.nrd0(data[clusters == 0])
#  if(sd2 < .Machine$double.eps) sd2 <- .Machine$double.eps
  result2 <- normalApproximation(data[clusters == 0], center2, sd2, cutoff = cutoff)
  
  tau <- sum(clusters == 1)/n
  if(result1$cutoff > 0 && result2$cutoff > 0) {
    bic2 <- dMixNorm(data, result1$center, result1$deviation, result2$center, result2$deviation, tau)
    if(sum(tmp <- bic2 < .Machine$double.eps) > 0)
      bic2[tmp] <- .Machine$double.eps
    bic2 <- -2 * sum(log(bic2)) +  5 * log(n)
  } else bic2 <-  NA

  ret <- list()
  ret$bics <- c(bic1, bic2)
  ret$bic1Parameters <- c(ans0$center, ans0$deviation)
  ret$bic2Parameters <- c(result1$center, result1$deviation, result2$center, result2$deviation, tau)
  return(ret)
}
