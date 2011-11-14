classify <-
function(data) {
  if(is.matrix(data)) {
      ans <- pam(data, 2)
      ans$clusters <- ans$clustering
    } else {ans <- pamUnbiased(data)}
  return(ans)
}
